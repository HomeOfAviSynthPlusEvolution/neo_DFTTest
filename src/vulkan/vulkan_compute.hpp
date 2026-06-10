#pragma once

#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)

#include "fft/fft_backend.hpp"
#include "fft/vkfft_vulkan_runtime.hpp"

#include <cstddef>
#include <cstdint>
#include <mutex>
#include <span>

#include <vulkan/vulkan.h>

namespace neo_dfttest::vulkan {

struct ComputeDispatch {
  const void* push_constants {nullptr};
  std::uint32_t push_constant_bytes {0};
  std::uint32_t groups_x {1};
  std::uint32_t groups_y {1};
  std::uint32_t groups_z {1};
  bool barrier_after {false};
};

class ComputeBinding final {
public:
  ComputeBinding() = default;

  ComputeBinding(const ComputeBinding&) = delete;
  ComputeBinding& operator=(const ComputeBinding&) = delete;

  ComputeBinding(ComputeBinding&& other) noexcept;
  ComputeBinding& operator=(ComputeBinding&& other) noexcept;

  ~ComputeBinding();

  [[nodiscard]] explicit operator bool() const noexcept {
    return descriptor_set_ != VK_NULL_HANDLE;
  }

private:
  friend class ComputePipeline;

  ComputeBinding(
    fft::VulkanRuntime& runtime,
    VkDescriptorPool descriptor_pool,
    std::mutex& descriptor_pool_mutex,
    VkDescriptorSet descriptor_set
  ) noexcept;

  void reset() noexcept;

  fft::VulkanRuntime* runtime_ {nullptr};
  VkDescriptorPool descriptor_pool_ {VK_NULL_HANDLE};
  std::mutex* descriptor_pool_mutex_ {nullptr};
  VkDescriptorSet descriptor_set_ {VK_NULL_HANDLE};
};

class ComputePipeline final {
public:
  ComputePipeline(
    fft::VulkanRuntime& runtime,
    std::span<const std::uint32_t> spirv,
    std::uint32_t storage_binding_count,
    std::uint32_t push_constant_bytes
  );

  ComputePipeline(const ComputePipeline&) = delete;
  ComputePipeline& operator=(const ComputePipeline&) = delete;

  ComputePipeline(ComputePipeline&&) = delete;
  ComputePipeline& operator=(ComputePipeline&&) = delete;

  ~ComputePipeline();

  [[nodiscard]] ComputeBinding bind_storage_buffers(
    std::span<const fft::DeviceBufferView> storage_buffers
  ) const;

  void record_dispatch_many(
    VkCommandBuffer command_buffer,
    const ComputeBinding& binding,
    std::span<const ComputeDispatch> dispatches
  ) const;

  void dispatch(
    std::span<const fft::DeviceBufferView> storage_buffers,
    const void* push_constants,
    std::uint32_t push_constant_bytes,
    std::uint32_t groups_x,
    std::uint32_t groups_y = 1,
    std::uint32_t groups_z = 1
  ) const;

  void dispatch_many(
    std::span<const fft::DeviceBufferView> storage_buffers,
    std::span<const ComputeDispatch> dispatches
  ) const;

private:
  void validate_dispatches(std::span<const ComputeDispatch> dispatches) const;

  fft::VulkanRuntime* runtime_ {nullptr};
  std::uint32_t storage_binding_count_ {0};
  std::uint32_t push_constant_bytes_ {0};
  VkShaderModule shader_ {VK_NULL_HANDLE};
  VkDescriptorSetLayout descriptor_set_layout_ {VK_NULL_HANDLE};
  VkPipelineLayout pipeline_layout_ {VK_NULL_HANDLE};
  VkPipeline pipeline_ {VK_NULL_HANDLE};
  VkDescriptorPool descriptor_pool_ {VK_NULL_HANDLE};
  mutable std::mutex descriptor_pool_mutex_;
};

} // namespace neo_dfttest::vulkan

#endif
