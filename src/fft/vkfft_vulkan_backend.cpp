#include "fft/fft_backend.hpp"
#include "fft/vkfft_vulkan_buffer.hpp"
#include "fft/vkfft_vulkan_runtime.hpp"

#include <algorithm>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

#if !defined(VKFFT_BACKEND)
#define VKFFT_BACKEND 0
#endif

#if VKFFT_BACKEND != 0
#error "vkfft_vulkan_backend.cpp requires VKFFT_BACKEND=0"
#endif

#include <vulkan/vulkan.h>
#include "vkFFT.h"

namespace neo_dfttest::fft {
namespace {

void check_vk(VkResult result, const char* action) {
  if (result != VK_SUCCESS) {
    throw std::runtime_error(std::string(action) + ": Vulkan error " + std::to_string(static_cast<int>(result)));
  }
}

void check_vkfft(VkFFTResult result, const char* action) {
  if (result != VKFFT_SUCCESS) {
    throw std::runtime_error(std::string(action) + ": " + getVkFFTErrorString(result));
  }
}

void validate_batch_capacity(const Plan& plan, int count) {
  if (count < 0 || count > plan.layout().max_batch) {
    throw std::runtime_error("FFT batch count exceeds plan capacity");
  }
}

std::size_t checked_extent(int extent) {
  if (extent <= 0) {
    throw std::runtime_error("invalid VkFFT shape extent");
  }
  return static_cast<std::size_t>(extent);
}

VkBuffer vk_buffer(DeviceBufferView view, VkDeviceSize required_size, const char* name) {
  if (view.handle == nullptr) {
    throw std::runtime_error(std::string("missing Vulkan device buffer for ") + name);
  }
  if (view.offset_bytes != 0) {
    throw std::runtime_error(std::string("VkFFT Vulkan backend does not currently support non-zero device buffer offsets for ") + name);
  }
  if (view.size_bytes != 0 && view.size_bytes < required_size) {
    throw std::runtime_error(std::string("Vulkan device buffer is too small for ") + name);
  }
  return reinterpret_cast<VkBuffer>(view.handle);
}

std::size_t real_transform_elements(TransformShape shape) {
  if (shape.rank < 1 || shape.rank > 3) {
    throw std::runtime_error("unsupported VkFFT rank");
  }

  std::size_t elements = 1;
  for (int index = 0; index < shape.rank; ++index) {
    elements *= checked_extent(shape.extents[static_cast<std::size_t>(index)]);
  }
  return elements;
}

std::size_t complex_transform_elements(TransformShape shape) {
  if (shape.rank < 1 || shape.rank > 3) {
    throw std::runtime_error("unsupported VkFFT rank");
  }

  std::size_t elements = 1;
  for (int index = 0; index < shape.rank; ++index) {
    std::size_t extent = checked_extent(shape.extents[static_cast<std::size_t>(index)]);
    if (index == shape.rank - 1) {
      extent = extent / 2 + 1;
    }
    elements *= extent;
  }
  return elements;
}

std::size_t padded_real_elements(TransformShape shape) {
  const std::size_t real_elements = real_transform_elements(shape);
  const std::size_t x_extent = checked_extent(shape.extents[static_cast<std::size_t>(shape.rank - 1)]);
  return (real_elements / x_extent) * (x_extent + 2);
}

std::size_t batch_storage_elements(int max_batch, std::size_t per_transform_elements) {
  if (max_batch < 1) {
    throw std::runtime_error("VkFFT batch layout must allow at least one transform");
  }
  return static_cast<std::size_t>(max_batch) * per_transform_elements;
}

bool requires_queue_idle_after_submit(VkPhysicalDevice physical_device) noexcept {
  VkPhysicalDeviceProperties properties {};
  vkGetPhysicalDeviceProperties(physical_device, &properties);
  return properties.deviceType == VK_PHYSICAL_DEVICE_TYPE_CPU;
}

class GlslangProcess final {
public:
  GlslangProcess() {
    if (!glslang_initialize_process()) {
      throw std::runtime_error("failed to initialize glslang for VkFFT");
    }
  }

  GlslangProcess(const GlslangProcess&) = delete;
  GlslangProcess& operator=(const GlslangProcess&) = delete;

  ~GlslangProcess() {
    glslang_finalize_process();
  }
};

class VulkanContext final : public VulkanRuntime {
public:
  VulkanContext() {
    create_instance();
    select_device();
    create_device();
    create_command_pool();
    create_fence();
  }

  VulkanContext(const VulkanContext&) = delete;
  VulkanContext& operator=(const VulkanContext&) = delete;

  ~VulkanContext() {
    if (device_ != VK_NULL_HANDLE) {
      vkDeviceWaitIdle(device_);
      if (fence_ != VK_NULL_HANDLE) {
        vkDestroyFence(device_, fence_, nullptr);
      }
      if (command_pool_ != VK_NULL_HANDLE) {
        vkDestroyCommandPool(device_, command_pool_, nullptr);
      }
      vkDestroyDevice(device_, nullptr);
    }
    if (instance_ != VK_NULL_HANDLE) {
      vkDestroyInstance(instance_, nullptr);
    }
  }

  [[nodiscard]] VkInstance instance() const noexcept {
    return instance_;
  }

  [[nodiscard]] VkPhysicalDevice physical_device() const noexcept override {
    return physical_device_;
  }

  [[nodiscard]] VkDevice device() const noexcept override {
    return device_;
  }

  [[nodiscard]] VkQueue queue() const noexcept override {
    return queue_;
  }

  [[nodiscard]] VkCommandPool command_pool() const noexcept override {
    return command_pool_;
  }

  [[nodiscard]] VkFence fence() const noexcept override {
    return fence_;
  }

  [[nodiscard]] std::uint32_t queue_family_index() const noexcept override {
    return queue_family_index_;
  }

  [[nodiscard]] std::uint32_t find_memory_type(std::uint32_t type_bits, VkMemoryPropertyFlags properties, VkDeviceSize size) const override {
    VkPhysicalDeviceMemoryProperties memory_properties {};
    vkGetPhysicalDeviceMemoryProperties(physical_device_, &memory_properties);

    for (std::uint32_t index = 0; index < memory_properties.memoryTypeCount; ++index) {
      const bool type_matches = (type_bits & (1u << index)) != 0;
      const auto& memory_type = memory_properties.memoryTypes[index];
      const auto& heap = memory_properties.memoryHeaps[memory_type.heapIndex];
      if (type_matches && (memory_type.propertyFlags & properties) == properties && heap.size >= size) {
        return index;
      }
    }

    throw std::runtime_error("failed to find a compatible Vulkan memory type");
  }

  [[nodiscard]] std::mutex& submission_mutex() const noexcept override {
    return mutex;
  }

  mutable std::mutex mutex;

private:
  struct DeviceCandidate {
    VkPhysicalDevice device {VK_NULL_HANDLE};
    std::uint32_t queue_family_index {0};
    int score {0};
  };

  void create_instance() {
    VkApplicationInfo application_info {VK_STRUCTURE_TYPE_APPLICATION_INFO};
    application_info.pApplicationName = "neo-dfttest";
    application_info.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    application_info.pEngineName = "neo-dfttest";
    application_info.engineVersion = VK_MAKE_VERSION(1, 0, 0);
    application_info.apiVersion = VK_API_VERSION_1_1;

    VkInstanceCreateInfo create_info {VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO};
    create_info.pApplicationInfo = &application_info;
    check_vk(vkCreateInstance(&create_info, nullptr, &instance_), "create Vulkan instance");
  }

  static int device_score(VkPhysicalDeviceType type) noexcept {
    switch (type) {
      case VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU:
        return 4;
      case VK_PHYSICAL_DEVICE_TYPE_INTEGRATED_GPU:
        return 3;
      case VK_PHYSICAL_DEVICE_TYPE_VIRTUAL_GPU:
        return 2;
      case VK_PHYSICAL_DEVICE_TYPE_CPU:
        return 1;
      default:
        return 0;
    }
  }

  std::vector<VkPhysicalDevice> enumerate_devices() const {
    std::uint32_t device_count = 0;
    check_vk(vkEnumeratePhysicalDevices(instance_, &device_count, nullptr), "enumerate Vulkan devices");
    if (device_count == 0) {
      throw std::runtime_error("VkFFT Vulkan backend requires a Vulkan compute device");
    }

    std::vector<VkPhysicalDevice> devices(device_count);
    check_vk(vkEnumeratePhysicalDevices(instance_, &device_count, devices.data()), "enumerate Vulkan devices");
    return devices;
  }

  static std::vector<VkQueueFamilyProperties> queue_families(VkPhysicalDevice device) {
    std::uint32_t family_count = 0;
    vkGetPhysicalDeviceQueueFamilyProperties(device, &family_count, nullptr);
    std::vector<VkQueueFamilyProperties> families(family_count);
    vkGetPhysicalDeviceQueueFamilyProperties(device, &family_count, families.data());
    return families;
  }

  void select_device() {
    DeviceCandidate best {};
    for (VkPhysicalDevice device : enumerate_devices()) {
      const auto families = queue_families(device);
      for (std::uint32_t index = 0; index < families.size(); ++index) {
        if ((families[index].queueFlags & VK_QUEUE_COMPUTE_BIT) == 0) {
          continue;
        }

        VkPhysicalDeviceProperties properties {};
        vkGetPhysicalDeviceProperties(device, &properties);
        const int score = device_score(properties.deviceType);
        if (best.device == VK_NULL_HANDLE || score > best.score) {
          best = DeviceCandidate{device, index, score};
        }
      }
    }

    if (best.device == VK_NULL_HANDLE) {
      throw std::runtime_error("VkFFT Vulkan backend could not find a compute queue");
    }

    physical_device_ = best.device;
    queue_family_index_ = best.queue_family_index;
  }

  void create_device() {
    const float priority = 1.0f;
    VkDeviceQueueCreateInfo queue_info {VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO};
    queue_info.queueFamilyIndex = queue_family_index_;
    queue_info.queueCount = 1;
    queue_info.pQueuePriorities = &priority;

    VkDeviceCreateInfo create_info {VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO};
    create_info.queueCreateInfoCount = 1;
    create_info.pQueueCreateInfos = &queue_info;

    check_vk(vkCreateDevice(physical_device_, &create_info, nullptr, &device_), "create Vulkan device");
    vkGetDeviceQueue(device_, queue_family_index_, 0, &queue_);
  }

  void create_command_pool() {
    VkCommandPoolCreateInfo create_info {VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO};
    create_info.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    create_info.queueFamilyIndex = queue_family_index_;
    check_vk(vkCreateCommandPool(device_, &create_info, nullptr, &command_pool_), "create Vulkan command pool");
  }

  void create_fence() {
    VkFenceCreateInfo create_info {VK_STRUCTURE_TYPE_FENCE_CREATE_INFO};
    check_vk(vkCreateFence(device_, &create_info, nullptr, &fence_), "create Vulkan fence");
  }

  VkInstance instance_ {VK_NULL_HANDLE};
  VkPhysicalDevice physical_device_ {VK_NULL_HANDLE};
  VkDevice device_ {VK_NULL_HANDLE};
  VkQueue queue_ {VK_NULL_HANDLE};
  VkCommandPool command_pool_ {VK_NULL_HANDLE};
  VkFence fence_ {VK_NULL_HANDLE};
  std::uint32_t queue_family_index_ {0};
};

class VkfftPlan final : public Plan::State {
public:
  VkfftPlan(std::shared_ptr<VulkanContext> context, TransformDirection direction, TransformShape shape, BatchLayout layout)
    : context_(std::move(context)),
      direction_(direction),
      shape_(shape),
      layout_(layout),
      real_elements_(real_transform_elements(shape)),
      complex_elements_(complex_transform_elements(shape)),
      padded_real_elements_(padded_real_elements(shape)),
      x_extent_(checked_extent(shape.extents[static_cast<std::size_t>(shape.rank - 1)])),
      row_count_(real_elements_ / x_extent_),
      buffer_size_(batch_storage_elements(layout.max_batch, padded_real_elements_) * sizeof(Real)),
      buffer_(make_vulkan_storage_buffer(
        *context_,
        buffer_size_,
        VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT
      )) {
    wait_queue_idle_after_submit_ = requires_queue_idle_after_submit(context_->physical_device());
    validate_layout();
    allocate_command_buffer();
    initialize_app();
  }

  VkfftPlan(const VkfftPlan&) = delete;
  VkfftPlan& operator=(const VkfftPlan&) = delete;

  ~VkfftPlan() override {
    if (app_initialized_) {
      deleteVkFFT(&app_);
      app_initialized_ = false;
    }
    if (command_buffer_ != VK_NULL_HANDLE) {
      vkFreeCommandBuffers(context_->device(), context_->command_pool(), 1, &command_buffer_);
      command_buffer_ = VK_NULL_HANDLE;
    }
  }

  [[nodiscard]] TransformDirection direction() const noexcept {
    return direction_;
  }

  void submit_r2c(R2CBatch batch) const {
    std::lock_guard lock(context_->mutex);
    zero_buffer();

    auto* mapped = static_cast<Real*>(buffer_.map());
    copy_real_input_to_padded_buffer(mapped, batch);
    buffer_.unmap();

    run(buffer_.buffer(), false);

    const auto* output = static_cast<const Complex*>(buffer_.map());
    copy_complex_output_from_buffer(output, batch);
    buffer_.unmap();
  }

  void submit_c2r(C2RBatch batch) const {
    std::lock_guard lock(context_->mutex);
    zero_buffer();

    auto* mapped = static_cast<Complex*>(buffer_.map());
    copy_complex_input_to_buffer(mapped, batch);
    buffer_.unmap();

    run(buffer_.buffer(), true);

    const auto* output = static_cast<const Real*>(buffer_.map());
    copy_real_output_from_padded_buffer(output, batch);
    buffer_.unmap();
  }

  void submit_r2c_device(R2CBatch batch) const {
    std::lock_guard lock(context_->mutex);
    const VkBuffer input = vk_buffer(batch.input.device, buffer_size_, "r2c input");
    const VkBuffer output = vk_buffer(batch.output.device, buffer_size_, "r2c output");
    if (input != output) {
      throw std::runtime_error("VkFFT Vulkan device r2c currently requires an in-place buffer");
    }
    run(input, false);
  }

  void submit_c2r_device(C2RBatch batch) const {
    std::lock_guard lock(context_->mutex);
    const VkBuffer input = vk_buffer(batch.input.device, buffer_size_, "c2r input");
    const VkBuffer output = vk_buffer(batch.output.device, buffer_size_, "c2r output");
    if (input != output) {
      throw std::runtime_error("VkFFT Vulkan device c2r currently requires an in-place buffer");
    }
    run(input, true);
  }

private:
  void validate_layout() const {
    if (layout_.max_batch < 1) {
      throw std::runtime_error("VkFFT batch layout must allow at least one transform");
    }
    if (layout_.real_stride_elements < static_cast<std::ptrdiff_t>(real_elements_)) {
      throw std::runtime_error("VkFFT real stride is smaller than one transform");
    }
    if (layout_.complex_stride_elements < static_cast<std::ptrdiff_t>(complex_elements_)) {
      throw std::runtime_error("VkFFT complex stride is smaller than one transform");
    }
  }

  void allocate_command_buffer() {
    VkCommandBufferAllocateInfo allocate_info {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
    allocate_info.commandPool = context_->command_pool();
    allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    allocate_info.commandBufferCount = 1;
    check_vk(vkAllocateCommandBuffers(context_->device(), &allocate_info, &command_buffer_), "allocate Vulkan command buffer");
  }

  void initialize_app() {
    VkFFTConfiguration configuration {};
    configuration.FFTdim = static_cast<pfUINT>(shape_.rank);
    for (int axis = 0; axis < shape_.rank; ++axis) {
      configuration.size[static_cast<std::size_t>(axis)] =
        static_cast<pfUINT>(shape_.extents[static_cast<std::size_t>(shape_.rank - 1 - axis)]);
    }

    configuration.physicalDevice = &physical_device_;
    configuration.device = &device_;
    configuration.queue = &queue_;
    configuration.commandPool = &command_pool_;
    configuration.fence = &fence_;
    configuration.isCompilerInitialized = 1;
    configuration.buffer = &buffer_handle_;
    configuration.bufferSize = &buffer_size_;
    configuration.performR2C = 1;
    configuration.numberBatches = static_cast<pfUINT>(layout_.max_batch);
    configuration.disableSetLocale = 1;
    configuration.disableReorderFourStep = 1;
    configuration.makeForwardPlanOnly = direction_ == TransformDirection::r2c ? 1 : 0;
    configuration.makeInversePlanOnly = direction_ == TransformDirection::c2r ? 1 : 0;

    check_vkfft(initializeVkFFT(&app_, configuration), "initialize VkFFT Vulkan plan");
    app_initialized_ = true;
  }

  void zero_buffer() const {
    void* mapped = buffer_.map();
    std::memset(mapped, 0, static_cast<std::size_t>(buffer_size_));
    buffer_.unmap();
  }

  void copy_real_input_to_padded_buffer(Real* destination, R2CBatch batch) const {
    for (int index = 0; index < batch.count; ++index) {
      const Real* source_block = batch.input.data + batch.input.stride_elements * index;
      Real* destination_block = destination + padded_real_elements_ * static_cast<std::size_t>(index);
      for (std::size_t row = 0; row < row_count_; ++row) {
        std::memcpy(destination_block + row * (x_extent_ + 2), source_block + row * x_extent_, x_extent_ * sizeof(Real));
      }
    }
  }

  void copy_complex_output_from_buffer(const Complex* source, R2CBatch batch) const {
    for (int index = 0; index < batch.count; ++index) {
      const Complex* source_block = source + complex_elements_ * static_cast<std::size_t>(index);
      Complex* destination_block = batch.output.data + batch.output.stride_elements * index;
      std::memcpy(destination_block, source_block, complex_elements_ * sizeof(Complex));
    }
  }

  void copy_complex_input_to_buffer(Complex* destination, C2RBatch batch) const {
    for (int index = 0; index < batch.count; ++index) {
      const Complex* source_block = batch.input.data + batch.input.stride_elements * index;
      Complex* destination_block = destination + complex_elements_ * static_cast<std::size_t>(index);
      std::memcpy(destination_block, source_block, complex_elements_ * sizeof(Complex));
    }
  }

  void copy_real_output_from_padded_buffer(const Real* source, C2RBatch batch) const {
    for (int index = 0; index < batch.count; ++index) {
      const Real* source_block = source + padded_real_elements_ * static_cast<std::size_t>(index);
      Real* destination_block = batch.output.data + batch.output.stride_elements * index;
      for (std::size_t row = 0; row < row_count_; ++row) {
        std::memcpy(destination_block + row * x_extent_, source_block + row * (x_extent_ + 2), x_extent_ * sizeof(Real));
      }
    }
  }

  void run(VkBuffer fft_buffer, bool inverse) const {
    check_vk(vkResetCommandBuffer(command_buffer_, 0), "reset Vulkan command buffer");

    VkCommandBufferBeginInfo begin_info {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    check_vk(vkBeginCommandBuffer(command_buffer_, &begin_info), "begin Vulkan command buffer");

    VkMemoryBarrier before_fft {VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    before_fft.srcAccessMask = VK_ACCESS_HOST_WRITE_BIT | VK_ACCESS_TRANSFER_WRITE_BIT | VK_ACCESS_SHADER_WRITE_BIT;
    before_fft.dstAccessMask = VK_ACCESS_SHADER_READ_BIT | VK_ACCESS_SHADER_WRITE_BIT;
    vkCmdPipelineBarrier(
      command_buffer_,
      VK_PIPELINE_STAGE_HOST_BIT | VK_PIPELINE_STAGE_TRANSFER_BIT | VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
      VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
      0,
      1,
      &before_fft,
      0,
      nullptr,
      0,
      nullptr
    );

    VkFFTLaunchParams launch {};
    VkBuffer buffer_handle = fft_buffer;
    launch.commandBuffer = &command_buffer_;
    launch.buffer = &buffer_handle;
    check_vkfft(VkFFTAppend(&app_, inverse ? 1 : 0, &launch), "execute VkFFT Vulkan plan");

    VkMemoryBarrier after_fft {VK_STRUCTURE_TYPE_MEMORY_BARRIER};
    after_fft.srcAccessMask = VK_ACCESS_SHADER_WRITE_BIT;
    after_fft.dstAccessMask =
      VK_ACCESS_HOST_READ_BIT |
      VK_ACCESS_TRANSFER_READ_BIT |
      VK_ACCESS_SHADER_READ_BIT |
      VK_ACCESS_SHADER_WRITE_BIT;
    vkCmdPipelineBarrier(
      command_buffer_,
      VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
      VK_PIPELINE_STAGE_HOST_BIT | VK_PIPELINE_STAGE_TRANSFER_BIT | VK_PIPELINE_STAGE_COMPUTE_SHADER_BIT,
      0,
      1,
      &after_fft,
      0,
      nullptr,
      0,
      nullptr
    );

    check_vk(vkEndCommandBuffer(command_buffer_), "end Vulkan command buffer");

    VkSubmitInfo submit_info {VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit_info.commandBufferCount = 1;
    submit_info.pCommandBuffers = &command_buffer_;
    check_vk(vkQueueSubmit(context_->queue(), 1, &submit_info, context_->fence()), "submit Vulkan command buffer");
    check_vk(vkWaitForFences(context_->device(), 1, &fence_, VK_TRUE, std::numeric_limits<std::uint64_t>::max()), "wait for Vulkan fence");
    check_vk(vkResetFences(context_->device(), 1, &fence_), "reset Vulkan fence");
    if (wait_queue_idle_after_submit_) {
      check_vk(vkQueueWaitIdle(context_->queue()), "wait for Vulkan CPU queue idle after VkFFT");
    }
  }

  std::shared_ptr<VulkanContext> context_;
  TransformDirection direction_;
  TransformShape shape_;
  BatchLayout layout_;
  std::size_t real_elements_ {0};
  std::size_t complex_elements_ {0};
  std::size_t padded_real_elements_ {0};
  std::size_t x_extent_ {0};
  std::size_t row_count_ {0};
  pfUINT buffer_size_ {0};
  VulkanDeviceBuffer buffer_;
  VkPhysicalDevice physical_device_ {context_->physical_device()};
  VkDevice device_ {context_->device()};
  VkQueue queue_ {context_->queue()};
  VkCommandPool command_pool_ {context_->command_pool()};
  VkFence fence_ {context_->fence()};
  VkBuffer buffer_handle_ {buffer_.buffer()};
  mutable VkCommandBuffer command_buffer_ {VK_NULL_HANDLE};
  mutable VkFFTApplication app_ {};
  mutable bool app_initialized_ {false};
  bool wait_queue_idle_after_submit_ {false};
};

const VkfftPlan& vkfft_plan(const Plan& plan) noexcept {
  return static_cast<const VkfftPlan&>(plan.state());
}

class VkfftVulkanBackend final : public Backend {
public:
  void load() override {
    glslang_ = std::make_unique<GlslangProcess>();
    context_ = std::make_shared<VulkanContext>();
    loaded_ = true;
  }

  void unload() noexcept override {
    context_.reset();
    glslang_.reset();
    loaded_ = false;
  }

  bool loaded() const noexcept override {
    return loaded_;
  }

  bool has_threading() const noexcept override {
    return false;
  }

  void set_thread_count(int) override {}

  BackendCapabilities capabilities() const noexcept override {
    return BackendCapabilities{
      true,
      false,
      true,
      8,
      std::numeric_limits<int>::max(),
    };
  }

  VulkanRuntime* vulkan_runtime() noexcept override {
    return context_.get();
  }

  Plan make_plan(
    TransformDirection direction,
    TransformShape shape,
    BatchLayout layout,
    Real*,
    Complex*,
    PlanOptions
  ) override {
    if (!context_) {
      throw std::runtime_error("VkFFT Vulkan backend is not loaded");
    }
    return Plan(std::make_unique<VkfftPlan>(context_, direction, shape, layout), layout);
  }

  Completion submit_r2c(const Plan& plan, R2CBatch batch, SubmitOptions) const override {
    validate_batch_capacity(plan, batch.count);
    if (batch.count == 0) {
      return Completion::completed();
    }

    const auto& native = vkfft_plan(plan);
    if (native.direction() != TransformDirection::r2c) {
      throw std::runtime_error("VkFFT Vulkan r2c submitted with a non-r2c plan");
    }
    if (batch.input.domain == MemoryDomain::host && batch.output.domain == MemoryDomain::host) {
      native.submit_r2c(batch);
    } else if (batch.input.domain == MemoryDomain::device && batch.output.domain == MemoryDomain::device) {
      native.submit_r2c_device(batch);
    } else {
      throw std::runtime_error("VkFFT Vulkan r2c requires matching host or device memory domains");
    }
    return Completion::completed();
  }

  Completion submit_c2r(const Plan& plan, C2RBatch batch, SubmitOptions) const override {
    validate_batch_capacity(plan, batch.count);
    if (batch.count == 0) {
      return Completion::completed();
    }

    const auto& native = vkfft_plan(plan);
    if (native.direction() != TransformDirection::c2r) {
      throw std::runtime_error("VkFFT Vulkan c2r submitted with a non-c2r plan");
    }
    if (batch.input.domain == MemoryDomain::host && batch.output.domain == MemoryDomain::host) {
      native.submit_c2r(batch);
    } else if (batch.input.domain == MemoryDomain::device && batch.output.domain == MemoryDomain::device) {
      native.submit_c2r_device(batch);
    } else {
      throw std::runtime_error("VkFFT Vulkan c2r requires matching host or device memory domains");
    }
    return Completion::completed();
  }

private:
  bool loaded_ {false};
  std::unique_ptr<GlslangProcess> glslang_;
  std::shared_ptr<VulkanContext> context_;
};

} // namespace

std::unique_ptr<Backend> create_vkfft_vulkan_backend() {
  return std::make_unique<VkfftVulkanBackend>();
}

} // namespace neo_dfttest::fft
