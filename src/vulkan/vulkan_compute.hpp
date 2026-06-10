#pragma once

#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)

#include "fft/fft_backend.hpp"
#include "fft/vkfft_vulkan_runtime.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

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
  fft::VulkanRuntime* runtime_ {nullptr};
  std::uint32_t storage_binding_count_ {0};
  std::uint32_t push_constant_bytes_ {0};
  VkShaderModule shader_ {VK_NULL_HANDLE};
  VkDescriptorSetLayout descriptor_set_layout_ {VK_NULL_HANDLE};
  VkPipelineLayout pipeline_layout_ {VK_NULL_HANDLE};
  VkPipeline pipeline_ {VK_NULL_HANDLE};
  VkDescriptorPool descriptor_pool_ {VK_NULL_HANDLE};
  VkDescriptorSet descriptor_set_ {VK_NULL_HANDLE};
};

} // namespace neo_dfttest::vulkan

#endif
