#include "executor/dft_executor.hpp"

#include "fft/vkfft_vulkan_buffer.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace neo_dfttest {
namespace {

[[nodiscard]] VkDeviceSize checked_bytes(std::size_t count, std::size_t element_size, const char* name) {
  if (count == 0) {
    return 1;
  }
  if (count > std::numeric_limits<VkDeviceSize>::max() / element_size) {
    throw std::runtime_error(std::string("Vulkan workspace size overflow for ") + name);
  }
  return static_cast<VkDeviceSize>(count * element_size);
}

[[nodiscard]] int dft_vkfft_padded_real_stride(const DFTBlockSettings& block, const DFTDerivedGeometry& derived) noexcept {
  const int rows = derived.block_volume / block.spatial_size;
  return rows * (block.spatial_size + 2);
}

[[nodiscard]] VkDeviceSize max_padded_frame_bytes(const DFTKernelContext& context) {
  std::size_t bytes = 0;
  for (int plane = 0; plane < context.format.num_planes; ++plane) {
    bytes = std::max(
      bytes,
      static_cast<std::size_t>(context.planes.pad_block_size[plane]) *
        static_cast<std::size_t>(context.block.temporal_size)
    );
  }
  return checked_bytes(bytes, 1, "padded frame");
}

[[nodiscard]] VkDeviceSize max_accumulation_bytes(const DFTKernelContext& context) {
  std::size_t elements = 0;
  for (int plane = 0; plane < context.format.num_planes; ++plane) {
    elements = std::max(
      elements,
      static_cast<std::size_t>(context.planes.e_stride[plane] / static_cast<int>(sizeof(float))) *
        static_cast<std::size_t>(context.planes.pad_height[plane])
    );
  }
  return checked_bytes(elements, sizeof(float), "accumulation");
}

[[nodiscard]] fft::VulkanDeviceBuffer make_workspace_buffer(fft::VulkanRuntime& runtime, VkDeviceSize size) {
  try {
    return fft::make_vulkan_storage_buffer(runtime, size, VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT);
  } catch (const std::runtime_error&) {
    return fft::make_vulkan_storage_buffer(
      runtime,
      size,
      VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT
    );
  }
}

struct VulkanDftWorkspaceShape {
  int complex_stride {0};
  int fft_storage_stride {0};
  int batch_capacity {0};
  int scratch_slots {0};
  VkDeviceSize padded_frame_bytes {0};
  VkDeviceSize accumulation_bytes {0};

  [[nodiscard]] bool operator==(const VulkanDftWorkspaceShape& other) const noexcept {
    return complex_stride == other.complex_stride &&
      fft_storage_stride == other.fft_storage_stride &&
      batch_capacity == other.batch_capacity &&
      scratch_slots == other.scratch_slots &&
      padded_frame_bytes == other.padded_frame_bytes &&
      accumulation_bytes == other.accumulation_bytes;
  }
};

struct VulkanDftCoefficientShape {
  int block_volume {0};
  int coefficient_count {0};

  [[nodiscard]] bool operator==(const VulkanDftCoefficientShape& other) const noexcept {
    return block_volume == other.block_volume && coefficient_count == other.coefficient_count;
  }
};

[[nodiscard]] VulkanDftWorkspaceShape make_workspace_shape(const DFTKernelContext& context) {
  const int batch_capacity = dft_fft_batch_capacity(context.batch_policy);
  const int scratch_slots = dft_fft_scratch_slots(context.batch_policy);
  return VulkanDftWorkspaceShape{
    dft_scratch_complex_stride(context.derived),
    dft_vkfft_padded_real_stride(context.block, context.derived),
    batch_capacity,
    scratch_slots,
    max_padded_frame_bytes(context),
    max_accumulation_bytes(context)
  };
}

[[nodiscard]] VulkanDftCoefficientShape make_coefficient_shape(const DFTKernelContext& context) noexcept {
  return VulkanDftCoefficientShape{
    context.derived.block_volume,
    context.derived.coefficient_count
  };
}

class VulkanDftWorkspace final {
public:
  explicit VulkanDftWorkspace(fft::VulkanRuntime& runtime) noexcept
    : runtime_(&runtime) {}

  VulkanDftWorkspace(const VulkanDftWorkspace&) = delete;
  VulkanDftWorkspace& operator=(const VulkanDftWorkspace&) = delete;

  void ensure(const DFTKernelContext& context) {
    const VulkanDftWorkspaceShape next = make_workspace_shape(context);
    if (ready_ && next == shape_) {
      return;
    }

    padded_frame_ = make_workspace_buffer(*runtime_, next.padded_frame_bytes);
    fft_storage_ = make_workspace_buffer(
      *runtime_,
      checked_bytes(
        static_cast<std::size_t>(next.fft_storage_stride) * static_cast<std::size_t>(next.scratch_slots),
        sizeof(float),
        "FFT storage"
      )
    );
    removed_mean_ = make_workspace_buffer(
      *runtime_,
      checked_bytes(
        static_cast<std::size_t>(next.complex_stride) * static_cast<std::size_t>(next.scratch_slots),
        sizeof(fft::Complex),
        "removed mean"
      )
    );
    accumulation_ = make_workspace_buffer(*runtime_, next.accumulation_bytes);

    fft::clear_vulkan_buffer(*runtime_, padded_frame_);
    fft::clear_vulkan_buffer(*runtime_, fft_storage_);
    fft::clear_vulkan_buffer(*runtime_, removed_mean_);
    fft::clear_vulkan_buffer(*runtime_, accumulation_);

    shape_ = next;
    ready_ = true;
  }

  [[nodiscard]] DFTMutableRealBatchView fft_real_batch(int slot, int count) const noexcept {
    return dft_device_real_batch_view(fft_storage_view(slot), shape_.fft_storage_stride, count);
  }

  [[nodiscard]] DFTMutableComplexBatchView fft_complex_batch(int slot, int count) const noexcept {
    return dft_device_complex_batch_view(fft_storage_view(slot), shape_.complex_stride, count);
  }

  [[nodiscard]] DFTMutableComplexBatchView removed_mean_batch(int slot, int count) const noexcept {
    return dft_device_complex_batch_view(removed_mean_view(slot), shape_.complex_stride, count);
  }

  [[nodiscard]] fft::DeviceBufferView padded_frame() const noexcept {
    return padded_frame_.view();
  }

  [[nodiscard]] fft::DeviceBufferView accumulation() const noexcept {
    return accumulation_.view();
  }

private:
  [[nodiscard]] std::size_t fft_slot_bytes() const noexcept {
    return static_cast<std::size_t>(shape_.fft_storage_stride) *
      static_cast<std::size_t>(shape_.batch_capacity) *
      sizeof(float);
  }

  [[nodiscard]] std::size_t complex_slot_bytes() const noexcept {
    return static_cast<std::size_t>(shape_.complex_stride) *
      static_cast<std::size_t>(shape_.batch_capacity) *
      sizeof(fft::Complex);
  }

  [[nodiscard]] fft::DeviceBufferView fft_storage_view(int slot) const noexcept {
    auto view = fft_storage_.view();
    view.offset_bytes = fft_slot_bytes() * static_cast<std::size_t>(slot);
    view.size_bytes = fft_slot_bytes();
    return view;
  }

  [[nodiscard]] fft::DeviceBufferView removed_mean_view(int slot) const noexcept {
    auto view = removed_mean_.view();
    view.offset_bytes = complex_slot_bytes() * static_cast<std::size_t>(slot);
    view.size_bytes = complex_slot_bytes();
    return view;
  }

  fft::VulkanRuntime* runtime_ {nullptr};
  VulkanDftWorkspaceShape shape_;
  bool ready_ {false};
  fft::VulkanDeviceBuffer padded_frame_;
  fft::VulkanDeviceBuffer fft_storage_;
  fft::VulkanDeviceBuffer removed_mean_;
  fft::VulkanDeviceBuffer accumulation_;
};

class VulkanDftCoefficientCache final {
public:
  explicit VulkanDftCoefficientCache(fft::VulkanRuntime& runtime) noexcept
    : runtime_(&runtime) {}

  VulkanDftCoefficientCache(const VulkanDftCoefficientCache&) = delete;
  VulkanDftCoefficientCache& operator=(const VulkanDftCoefficientCache&) = delete;

  void ensure(const DFTKernelContext& context) {
    const VulkanDftCoefficientShape next = make_coefficient_shape(context);
    if (ready_ && next == shape_) {
      return;
    }

    window_ = allocate_and_upload(
      context.coefficients.window.data,
      context.coefficients.window.size,
      sizeof(float)
    );
    sigmas_ = allocate_and_upload(
      context.coefficients.sigmas.data,
      context.coefficients.sigmas.size,
      sizeof(float)
    );
    sigmas2_ = allocate_and_upload(
      context.coefficients.sigmas2.data,
      context.coefficients.sigmas2.size,
      sizeof(float)
    );
    pmins_ = allocate_and_upload(
      context.coefficients.pmins.data,
      context.coefficients.pmins.size,
      sizeof(float)
    );
    pmaxs_ = allocate_and_upload(
      context.coefficients.pmaxs.data,
      context.coefficients.pmaxs.size,
      sizeof(float)
    );
    window_dft_ = allocate_and_upload(
      context.coefficients.window_dft.data,
      context.coefficients.window_dft.size,
      sizeof(fft::Complex)
    );

    shape_ = next;
    ready_ = true;
  }

  [[nodiscard]] fft::DeviceBufferView window() const noexcept {
    return window_.view();
  }

  [[nodiscard]] fft::DeviceBufferView sigmas() const noexcept {
    return sigmas_.view();
  }

  [[nodiscard]] fft::DeviceBufferView sigmas2() const noexcept {
    return sigmas2_.view();
  }

  [[nodiscard]] fft::DeviceBufferView pmins() const noexcept {
    return pmins_.view();
  }

  [[nodiscard]] fft::DeviceBufferView pmaxs() const noexcept {
    return pmaxs_.view();
  }

  [[nodiscard]] fft::DeviceBufferView window_dft() const noexcept {
    return window_dft_.view();
  }

private:
  [[nodiscard]] fft::VulkanDeviceBuffer allocate_and_upload(
    const void* source,
    int count,
    std::size_t element_size
  ) const {
    const VkDeviceSize bytes = checked_bytes(static_cast<std::size_t>(std::max(count, 0)), element_size, "coefficient table");
    auto buffer = make_workspace_buffer(*runtime_, bytes);
    if (count > 0) {
      fft::upload_to_vulkan_buffer(*runtime_, buffer, source, bytes);
    } else {
      fft::clear_vulkan_buffer(*runtime_, buffer);
    }
    return buffer;
  }

  fft::VulkanRuntime* runtime_ {nullptr};
  VulkanDftCoefficientShape shape_;
  bool ready_ {false};
  fft::VulkanDeviceBuffer window_;
  fft::VulkanDeviceBuffer sigmas_;
  fft::VulkanDeviceBuffer sigmas2_;
  fft::VulkanDeviceBuffer pmins_;
  fft::VulkanDeviceBuffer pmaxs_;
  fft::VulkanDeviceBuffer window_dft_;
};

} // namespace

class VulkanDftExecutor final : public DftExecutor {
public:
  VulkanDftExecutor(unsigned opt, DFTClipFormat format, fft::VulkanRuntime* runtime)
    : runtime_(runtime),
      fallback_(create_cpu_dft_executor(opt, format)) {
    if (runtime_ == nullptr) {
      throw std::runtime_error("Vulkan DFT executor requires a VkFFT Vulkan runtime");
    }
    coefficients_ = std::make_unique<VulkanDftCoefficientCache>(*runtime_);
  }

  [[nodiscard]] DftExecutorCapabilities capabilities() const noexcept override {
    return DftExecutorCapabilities{
      fft::MemoryDomain::device,
      true,
      true,
      true,
      true
    };
  }

  [[nodiscard]] DftMemoryPlan memory_plan() const noexcept override {
    return dft_vulkan_memory_plan(true);
  }

  [[nodiscard]] DFTBatchPolicy make_batch_policy(
    const DFTBlockSettings& block,
    const fft::BackendCapabilities& fft_capabilities
  ) const noexcept override {
    return fallback_->make_batch_policy(block, fft_capabilities);
  }

  void copy_pad(DftCopyPadRequest request) override {
    fallback_->copy_pad(request);
  }

  void process_frame(DftFrameProcessRequest request) override {
    prepare_workspace(request.workspace, request.context);
    fallback_->process_frame(request);
  }

  void process(DftProcessRequest request) override {
    prepare_workspace(request.workspace, request.context);
    fallback_->process(request);
  }

  void process_spatial(DftProcessSpatialRequest request) override {
    prepare_workspace(request.workspace, request.context);
    fallback_->process_spatial(request);
  }

  void process_temporal(DftProcessTemporalRequest request) override {
    prepare_workspace(request.workspace, request.context);
    fallback_->process_temporal(request);
  }

private:
  void prepare_workspace(DftWorkspaceLease lease, const DFTKernelContext& context) {
    workspace_for_slot(lease.slot_id).ensure(context);
    prepare_coefficients(context);
  }

  void prepare_coefficients(const DFTKernelContext& context) {
    std::lock_guard lock(coefficients_mutex_);
    coefficients_->ensure(context);
  }

  VulkanDftWorkspace& workspace_for_slot(unsigned int slot_id) {
    std::lock_guard lock(workspaces_mutex_);
    if (slot_id >= workspaces_.size()) {
      workspaces_.resize(static_cast<std::size_t>(slot_id) + 1);
    }

    auto& workspace = workspaces_[slot_id];
    if (!workspace) {
      workspace = std::make_unique<VulkanDftWorkspace>(*runtime_);
    }
    return *workspace;
  }

  fft::VulkanRuntime* runtime_ {nullptr};
  std::unique_ptr<VulkanDftCoefficientCache> coefficients_;
  std::unique_ptr<DftExecutor> fallback_;
  std::mutex workspaces_mutex_;
  std::mutex coefficients_mutex_;
  std::vector<std::unique_ptr<VulkanDftWorkspace>> workspaces_;
};

std::unique_ptr<DftExecutor> create_vulkan_dft_executor(
  unsigned opt,
  DFTClipFormat format,
  fft::VulkanRuntime* runtime
) {
  return std::make_unique<VulkanDftExecutor>(opt, format, runtime);
}

} // namespace neo_dfttest
