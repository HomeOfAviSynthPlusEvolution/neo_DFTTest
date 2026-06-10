#include "fft/vkfft_vulkan_buffer.hpp"

#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)

#include <stdexcept>
#include <string>
#include <utility>

namespace neo_dfttest::fft {
namespace {

void check_vk(VkResult result, const char* action) {
  if (result != VK_SUCCESS) {
    throw std::runtime_error(std::string(action) + ": Vulkan error " + std::to_string(static_cast<int>(result)));
  }
}

} // namespace

VulkanDeviceBuffer::VulkanDeviceBuffer(
  VulkanRuntime& runtime,
  VkDeviceSize size,
  VkBufferUsageFlags usage,
  VkMemoryPropertyFlags memory_properties
)
  : runtime_(&runtime),
    size_(size) {
  if (size_ == 0) {
    throw std::runtime_error("invalid zero-sized Vulkan buffer allocation");
  }

  VkBufferCreateInfo buffer_info {VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO};
  buffer_info.size = size_;
  buffer_info.usage = usage;
  buffer_info.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
  buffer_info.queueFamilyIndexCount = 1;
  const std::uint32_t queue_family = runtime_->queue_family_index();
  buffer_info.pQueueFamilyIndices = &queue_family;
  try {
    check_vk(vkCreateBuffer(runtime_->device(), &buffer_info, nullptr, &buffer_), "create Vulkan buffer");

    VkMemoryRequirements requirements {};
    vkGetBufferMemoryRequirements(runtime_->device(), buffer_, &requirements);

    VkMemoryAllocateInfo allocate_info {VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO};
    allocate_info.allocationSize = requirements.size;
    allocate_info.memoryTypeIndex = runtime_->find_memory_type(
      requirements.memoryTypeBits,
      memory_properties,
      requirements.size
    );
    check_vk(vkAllocateMemory(runtime_->device(), &allocate_info, nullptr, &memory_), "allocate Vulkan buffer memory");
    check_vk(vkBindBufferMemory(runtime_->device(), buffer_, memory_, 0), "bind Vulkan buffer memory");
  } catch (...) {
    if (memory_ != VK_NULL_HANDLE) {
      vkFreeMemory(runtime_->device(), memory_, nullptr);
      memory_ = VK_NULL_HANDLE;
    }
    if (buffer_ != VK_NULL_HANDLE) {
      vkDestroyBuffer(runtime_->device(), buffer_, nullptr);
      buffer_ = VK_NULL_HANDLE;
    }
    runtime_ = nullptr;
    size_ = 0;
    throw;
  }
}

VulkanDeviceBuffer::VulkanDeviceBuffer(VulkanDeviceBuffer&& other) noexcept
  : runtime_(std::exchange(other.runtime_, nullptr)),
    buffer_(std::exchange(other.buffer_, VK_NULL_HANDLE)),
    memory_(std::exchange(other.memory_, VK_NULL_HANDLE)),
    size_(std::exchange(other.size_, 0)) {}

VulkanDeviceBuffer& VulkanDeviceBuffer::operator=(VulkanDeviceBuffer&& other) noexcept {
  if (this != &other) {
    reset();
    runtime_ = std::exchange(other.runtime_, nullptr);
    buffer_ = std::exchange(other.buffer_, VK_NULL_HANDLE);
    memory_ = std::exchange(other.memory_, VK_NULL_HANDLE);
    size_ = std::exchange(other.size_, 0);
  }
  return *this;
}

VulkanDeviceBuffer::~VulkanDeviceBuffer() {
  reset();
}

void VulkanDeviceBuffer::reset() noexcept {
  if (!runtime_) {
    return;
  }
  if (buffer_ != VK_NULL_HANDLE) {
    vkDestroyBuffer(runtime_->device(), buffer_, nullptr);
    buffer_ = VK_NULL_HANDLE;
  }
  if (memory_ != VK_NULL_HANDLE) {
    vkFreeMemory(runtime_->device(), memory_, nullptr);
    memory_ = VK_NULL_HANDLE;
  }
  runtime_ = nullptr;
  size_ = 0;
}

DeviceBufferView VulkanDeviceBuffer::view() const noexcept {
  return DeviceBufferView{
    reinterpret_cast<void*>(buffer_),
    0,
    static_cast<std::size_t>(size_)
  };
}

void* VulkanDeviceBuffer::map() const {
  void* data = nullptr;
  check_vk(vkMapMemory(runtime_->device(), memory_, 0, size_, 0, &data), "map Vulkan buffer memory");
  return data;
}

void VulkanDeviceBuffer::unmap() const noexcept {
  vkUnmapMemory(runtime_->device(), memory_);
}

VulkanDeviceBuffer make_vulkan_storage_buffer(
  VulkanRuntime& runtime,
  VkDeviceSize size,
  VkMemoryPropertyFlags memory_properties
) {
  constexpr VkBufferUsageFlags usage =
    VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
    VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
    VK_BUFFER_USAGE_TRANSFER_DST_BIT;
  return VulkanDeviceBuffer(runtime, size, usage, memory_properties);
}

} // namespace neo_dfttest::fft

#endif
