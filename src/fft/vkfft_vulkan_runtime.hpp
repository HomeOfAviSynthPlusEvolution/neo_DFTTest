#pragma once

#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)

#include <cstdint>
#include <mutex>

#include <vulkan/vulkan.h>

namespace neo_dfttest::fft {

class VulkanRuntime {
public:
  virtual ~VulkanRuntime() = default;

  [[nodiscard]] virtual VkPhysicalDevice physical_device() const noexcept = 0;
  [[nodiscard]] virtual VkDevice device() const noexcept = 0;
  [[nodiscard]] virtual VkQueue queue() const noexcept = 0;
  [[nodiscard]] virtual VkCommandPool command_pool() const noexcept = 0;
  [[nodiscard]] virtual VkFence fence() const noexcept = 0;
  [[nodiscard]] virtual std::uint32_t queue_family_index() const noexcept = 0;
  [[nodiscard]] virtual std::uint32_t find_memory_type(
    std::uint32_t type_bits,
    VkMemoryPropertyFlags properties,
    VkDeviceSize size
  ) const = 0;
  [[nodiscard]] virtual std::mutex& submission_mutex() const noexcept = 0;
};

} // namespace neo_dfttest::fft

#endif
