#pragma once

#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)

#include "fft/fft_backend.hpp"
#include "fft/vkfft_vulkan_runtime.hpp"

#include <cstdint>
#include <span>

#include <vulkan/vulkan.h>

namespace neo_dfttest::fft {

struct VulkanBufferUpload {
  DeviceBufferView destination;
  const void* source {nullptr};
  VkDeviceSize bytes {0};
};

class VulkanDeviceBuffer final {
public:
  VulkanDeviceBuffer() = default;
  VulkanDeviceBuffer(
    VulkanRuntime& runtime,
    VkDeviceSize size,
    VkBufferUsageFlags usage,
    VkMemoryPropertyFlags memory_properties
  );

  VulkanDeviceBuffer(const VulkanDeviceBuffer&) = delete;
  VulkanDeviceBuffer& operator=(const VulkanDeviceBuffer&) = delete;

  VulkanDeviceBuffer(VulkanDeviceBuffer&& other) noexcept;
  VulkanDeviceBuffer& operator=(VulkanDeviceBuffer&& other) noexcept;

  ~VulkanDeviceBuffer();

  void reset() noexcept;

  [[nodiscard]] explicit operator bool() const noexcept {
    return buffer_ != VK_NULL_HANDLE;
  }

  [[nodiscard]] VkBuffer buffer() const noexcept {
    return buffer_;
  }

  [[nodiscard]] VkDeviceMemory memory() const noexcept {
    return memory_;
  }

  [[nodiscard]] VkDeviceSize size() const noexcept {
    return size_;
  }

  [[nodiscard]] DeviceBufferView view() const noexcept;

  [[nodiscard]] void* map() const;
  void unmap() const noexcept;

private:
  VulkanRuntime* runtime_ {nullptr};
  VkBuffer buffer_ {VK_NULL_HANDLE};
  VkDeviceMemory memory_ {VK_NULL_HANDLE};
  VkDeviceSize size_ {0};
};

[[nodiscard]] VulkanDeviceBuffer make_vulkan_storage_buffer(
  VulkanRuntime& runtime,
  VkDeviceSize size,
  VkMemoryPropertyFlags memory_properties
);

void submit_vulkan_commands(
  VulkanRuntime& runtime,
  void (*record)(VkCommandBuffer command_buffer, void* user),
  void* user,
  bool wait_queue_idle_after_submit = false
);
void clear_vulkan_buffer(VulkanRuntime& runtime, VulkanDeviceBuffer& buffer, std::uint32_t value = 0);
void clear_vulkan_buffer_view(VulkanRuntime& runtime, DeviceBufferView view, std::uint32_t value = 0);
void upload_to_vulkan_buffer(VulkanRuntime& runtime, VulkanDeviceBuffer& destination, const void* source, VkDeviceSize bytes);
void upload_to_vulkan_buffer_view(VulkanRuntime& runtime, DeviceBufferView destination, const void* source, VkDeviceSize bytes);
void upload_to_vulkan_buffer_views(VulkanRuntime& runtime, std::span<const VulkanBufferUpload> uploads);
void download_from_vulkan_buffer(VulkanRuntime& runtime, VulkanDeviceBuffer& source, void* destination, VkDeviceSize bytes);

} // namespace neo_dfttest::fft

#endif
