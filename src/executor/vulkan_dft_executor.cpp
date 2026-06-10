#include "executor/dft_executor.hpp"

#include "executor/dft_fft_operation.hpp"
#include "fft/vkfft_vulkan_buffer.hpp"
#include "vulkan/vulkan_compute.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace neo_dfttest {
namespace {

constexpr std::uint32_t kLoadWindowSpv[] = {
#include "neo_dfttest/load_window_comp_spv.inc"
};

[[nodiscard]] VkDeviceSize checked_bytes(std::size_t count, std::size_t element_size, const char* name) {
  if (count == 0) {
    return 1;
  }
  if (count > std::numeric_limits<VkDeviceSize>::max() / element_size) {
    throw std::runtime_error(std::string("Vulkan workspace size overflow for ") + name);
  }
  return static_cast<VkDeviceSize>(count * element_size);
}

[[nodiscard]] std::size_t align_storage_bytes(std::size_t bytes) noexcept {
  return (bytes + 3u) & ~std::size_t{3u};
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
  return checked_bytes(align_storage_bytes(bytes), 1, "padded frame");
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

[[nodiscard]] std::uint32_t load_window_sample_kind(const DFTClipFormat& format) {
  if (format.bytes_per_sample == 1 || format.bytes_per_sample == 2 || format.bytes_per_sample == 4) {
    return static_cast<std::uint32_t>(format.bytes_per_sample);
  }
  throw std::runtime_error("Vulkan load-window kernel received an unsupported sample size");
}

[[nodiscard]] std::uint32_t ceil_div_u32(int value, int divisor) noexcept {
  return static_cast<std::uint32_t>((value + divisor - 1) / divisor);
}

[[nodiscard]] std::uint32_t checked_u32(std::size_t value, const char* name) {
  if (value > std::numeric_limits<std::uint32_t>::max()) {
    throw std::runtime_error(std::string("Vulkan load-window value is too large for ") + name);
  }
  return static_cast<std::uint32_t>(value);
}

struct LoadWindowPushConstants {
  std::uint32_t source_base_bytes {0};
  std::uint32_t source_stride_bytes {0};
  std::uint32_t source_x {0};
  std::uint32_t source_y {0};
  std::uint32_t block_size {0};
  std::uint32_t sample_kind {0};
  std::uint32_t window_offset {0};
  std::uint32_t output_offset {0};
  std::uint32_t output_stride {0};
  float divisor {1.0f};
};

static_assert(sizeof(LoadWindowPushConstants) == 40);

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
    const VkDeviceSize fft_storage_bytes = checked_bytes(
      static_cast<std::size_t>(next.fft_storage_stride) * static_cast<std::size_t>(next.batch_capacity),
      sizeof(float),
      "FFT storage"
    );
    const VkDeviceSize removed_mean_bytes = checked_bytes(
      static_cast<std::size_t>(next.complex_stride) * static_cast<std::size_t>(next.batch_capacity),
      sizeof(fft::Complex),
      "removed mean"
    );
    for (auto& buffer : fft_storage_) {
      buffer = make_workspace_buffer(*runtime_, fft_storage_bytes);
    }
    for (auto& buffer : removed_mean_) {
      buffer = make_workspace_buffer(*runtime_, removed_mean_bytes);
    }
    accumulation_ = make_workspace_buffer(*runtime_, next.accumulation_bytes);

    fft::clear_vulkan_buffer(*runtime_, padded_frame_);
    for (auto& buffer : fft_storage_) {
      fft::clear_vulkan_buffer(*runtime_, buffer);
    }
    for (auto& buffer : removed_mean_) {
      fft::clear_vulkan_buffer(*runtime_, buffer);
    }
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

  [[nodiscard]] fft::VulkanDeviceBuffer& padded_frame_buffer() noexcept {
    return padded_frame_;
  }

  [[nodiscard]] fft::DeviceBufferView accumulation() const noexcept {
    return accumulation_.view();
  }

  [[nodiscard]] int fft_storage_stride() const noexcept {
    return shape_.fft_storage_stride;
  }

private:
  [[nodiscard]] fft::DeviceBufferView fft_storage_view(int slot) const noexcept {
    return fft_storage_[static_cast<std::size_t>(slot)].view();
  }

  [[nodiscard]] fft::DeviceBufferView removed_mean_view(int slot) const noexcept {
    return removed_mean_[static_cast<std::size_t>(slot)].view();
  }

  fft::VulkanRuntime* runtime_ {nullptr};
  VulkanDftWorkspaceShape shape_;
  bool ready_ {false};
  fft::VulkanDeviceBuffer padded_frame_;
  std::array<fft::VulkanDeviceBuffer, kDftFftPipelineSlots> fft_storage_;
  std::array<fft::VulkanDeviceBuffer, kDftFftPipelineSlots> removed_mean_;
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

class VulkanLoadWindowKernel final {
public:
  explicit VulkanLoadWindowKernel(fft::VulkanRuntime& runtime)
    : pipeline_(
        runtime,
        std::span<const std::uint32_t>{kLoadWindowSpv, sizeof(kLoadWindowSpv) / sizeof(kLoadWindowSpv[0])},
        3,
        sizeof(LoadWindowPushConstants)
      ) {}

  void dispatch(
    fft::DeviceBufferView source,
    fft::DeviceBufferView window,
    fft::DeviceBufferView output,
    const LoadWindowPushConstants& constants
  ) const {
    const std::array<fft::DeviceBufferView, 3> buffers {source, window, output};
    pipeline_.dispatch(
      buffers,
      &constants,
      sizeof(constants),
      ceil_div_u32(static_cast<int>(constants.output_stride), 16),
      ceil_div_u32(static_cast<int>(constants.block_size), 16)
    );
  }

private:
  vulkan::ComputePipeline pipeline_;
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
    load_window_ = std::make_unique<VulkanLoadWindowKernel>(*runtime_);
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
    for (int plane = 0; plane < request.plane_count; ++plane) {
      const auto& plane_request = request.planes[static_cast<std::size_t>(plane)];
      if (request.context.planes.process[plane] == 3) {
        process(DftProcessRequest{
          request.workspace,
          plane,
          request.mode,
          plane_request.sources,
          plane_request.source_count,
          plane_request.destination,
          request.temporal_position,
          request.context
        });
      } else if (request.context.planes.process[plane] == 2) {
        copy_plane_rows(plane_request, plane, request);
      }
    }
  }

  void process(DftProcessRequest request) override {
    if (request.source_count <= 0 || request.source_count > kMaxDftTemporalFrames) {
      throw std::runtime_error("Vulkan DFT executor received an invalid source frame count");
    }

    const int pad_block_size = request.context.planes.pad_block_size[request.plane];
    const int pad_stride = request.context.planes.pad_stride[request.plane];
    AlignedBuffer<unsigned char> padded(
      static_cast<std::size_t>(pad_block_size) * static_cast<std::size_t>(request.source_count)
    );

    for (int index = 0; index < request.source_count; ++index) {
      copy_pad(DftCopyPadRequest{
        request.plane,
        request.sources[static_cast<std::size_t>(index)],
        DFTMutablePlaneBytes{padded.data() + pad_block_size * index, pad_stride},
        request.context
      });
    }

    if (request.mode == DftProcessMode::spatial) {
      process_spatial(DftProcessSpatialRequest{
        request.workspace,
        request.plane,
        DFTPlaneBytes{padded.data(), pad_stride},
        request.destination,
        request.context
      });
      return;
    }

    process_temporal(DftProcessTemporalRequest{
      request.workspace,
      request.plane,
      DFTPlaneBytes{padded.data(), pad_stride},
      request.destination,
      request.temporal_position,
      request.context
    });
  }

  void process_spatial(DftProcessSpatialRequest request) override {
    VulkanDftWorkspace& workspace = prepare_workspace(request.workspace, request.context);
    stage_spatial_forward_batches(workspace, request);
    fallback_->process_spatial(request);
  }

  void process_temporal(DftProcessTemporalRequest request) override {
    VulkanDftWorkspace& workspace = prepare_workspace(request.workspace, request.context);
    stage_temporal_forward_batches(workspace, request);
    fallback_->process_temporal(request);
  }

private:
  VulkanDftWorkspace& prepare_workspace(DftWorkspaceLease lease, const DFTKernelContext& context) {
    VulkanDftWorkspace& workspace = workspace_for_slot(lease.slot_id);
    workspace.ensure(context);
    prepare_coefficients(context);
    return workspace;
  }

  void prepare_coefficients(const DFTKernelContext& context) {
    std::lock_guard lock(coefficients_mutex_);
    coefficients_->ensure(context);
  }

  void upload_padded_source(VulkanDftWorkspace& workspace, DFTPlaneBytes source, VkDeviceSize source_bytes) {
    fft::upload_to_vulkan_buffer(
      *runtime_,
      workspace.padded_frame_buffer(),
      source.data,
      source_bytes
    );
  }

  void stage_spatial_forward_batches(VulkanDftWorkspace& workspace, const DftProcessSpatialRequest& request) {
    upload_padded_source(
      workspace,
      request.source,
      static_cast<VkDeviceSize>(request.context.planes.pad_block_size[request.plane])
    );

    stage_forward_batches(
      workspace,
      request.plane,
      request.context,
      [&](int y, int x, int batch_index, int slot) {
        dispatch_load_window(
          workspace,
          request.context,
          request.source.stride_bytes,
          0,
          x,
          y,
          0,
          batch_output_offset(workspace, batch_index),
          slot
        );
      }
    );
  }

  void stage_temporal_forward_batches(VulkanDftWorkspace& workspace, const DftProcessTemporalRequest& request) {
    upload_padded_source(
      workspace,
      request.source,
      static_cast<VkDeviceSize>(
        request.context.planes.pad_block_size[request.plane] * request.context.block.temporal_size
      )
    );

    stage_forward_batches(
      workspace,
      request.plane,
      request.context,
      [&](int y, int x, int batch_index, int slot) {
        for (int z = 0; z < request.context.block.temporal_size; ++z) {
          const std::size_t source_base =
            static_cast<std::size_t>(request.context.planes.pad_block_size[request.plane]) *
            static_cast<std::size_t>(z);
          const std::size_t window_offset =
            static_cast<std::size_t>(request.context.derived.block_area) *
            static_cast<std::size_t>(z);
          const std::size_t output_offset =
            batch_output_offset(workspace, batch_index) +
            static_cast<std::size_t>(request.context.block.spatial_size + 2) *
            static_cast<std::size_t>(request.context.block.spatial_size) *
            static_cast<std::size_t>(z);

          dispatch_load_window(
            workspace,
            request.context,
            request.source.stride_bytes,
            source_base,
            x,
            y,
            window_offset,
            output_offset,
            slot
          );
        }
      }
    );
  }

  template<typename LoadBatchBlock>
  void stage_forward_batches(
    VulkanDftWorkspace& workspace,
    int plane,
    const DFTKernelContext& context,
    LoadBatchBlock&& load_batch_block
  ) {
    const int width = context.planes.pad_width[plane];
    const int eheight = context.planes.e_height[plane];
    const int batch_capacity = dft_fft_batch_capacity(context.batch_policy);
    int active_slot = 0;

    for (int y = 0; y < eheight; y += context.derived.step) {
      for (int x = 0; x <= width - context.block.spatial_size;) {
        DFTBlockBatch batch;
        for (
          ;
          batch.count < batch_capacity && x <= width - context.block.spatial_size;
          ++batch.count, x += context.derived.step
        ) {
          batch.jobs[static_cast<std::size_t>(batch.count)] = DFTBlockJob{plane, y, x, 0};
        }

        if (batch.count == 0) {
          continue;
        }

        auto real = workspace.fft_real_batch(active_slot, batch.count);
        auto coefficients = workspace.fft_complex_batch(active_slot, batch.count);

        for (int index = 0; index < batch.count; ++index) {
          const DFTBlockJob& job = dft_block_job(batch, index);
          load_batch_block(job.y, job.x, index, active_slot);
        }

        DftFftOperations{context}.submit_forward(real, coefficients).wait();
        active_slot = 1 - active_slot;
      }
    }
  }

  [[nodiscard]] std::size_t batch_output_offset(VulkanDftWorkspace& workspace, int batch_index) const noexcept {
    return static_cast<std::size_t>(workspace.fft_storage_stride()) * static_cast<std::size_t>(batch_index);
  }

  void dispatch_load_window(
    VulkanDftWorkspace& workspace,
    const DFTKernelContext& context,
    int source_stride_bytes,
    std::size_t source_base_bytes,
    int source_x,
    int source_y,
    std::size_t window_offset,
    std::size_t output_offset,
    int slot
  ) {
    const auto real = workspace.fft_real_batch(slot, 1);

    const LoadWindowPushConstants constants {
      checked_u32(source_base_bytes, "source base offset"),
      static_cast<std::uint32_t>(source_stride_bytes),
      static_cast<std::uint32_t>(source_x),
      static_cast<std::uint32_t>(source_y),
      static_cast<std::uint32_t>(context.block.spatial_size),
      load_window_sample_kind(context.format),
      checked_u32(window_offset, "window offset"),
      checked_u32(output_offset, "output offset"),
      static_cast<std::uint32_t>(context.block.spatial_size + 2),
      context.sample.divisor
    };

    load_window_->dispatch(
      workspace.padded_frame(),
      coefficients_->window(),
      real.device,
      constants
    );
  }

  static void copy_plane_rows(
    const DftFramePlaneRequest& plane_request,
    int plane,
    const DftFrameProcessRequest& frame_request
  ) {
    const int source_index = frame_request.mode == DftProcessMode::temporal
      ? frame_request.temporal_position
      : 0;
    if (source_index < 0 || source_index >= plane_request.source_count) {
      throw std::runtime_error("Vulkan DFT executor copy plane source is missing");
    }

    const DFTPlaneBytes source = plane_request.sources[static_cast<std::size_t>(source_index)];
    const DFTMutablePlaneBytes destination = plane_request.destination;
    const int row_bytes =
      frame_request.context.planes.width[plane] * frame_request.context.format.bytes_per_sample;
    const int height = frame_request.context.planes.height[plane];

    for (int y = 0; y < height; ++y) {
      std::memcpy(
        destination.data + static_cast<std::ptrdiff_t>(destination.stride_bytes) * y,
        source.data + static_cast<std::ptrdiff_t>(source.stride_bytes) * y,
        static_cast<std::size_t>(row_bytes)
      );
    }
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
  std::unique_ptr<VulkanLoadWindowKernel> load_window_;
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
