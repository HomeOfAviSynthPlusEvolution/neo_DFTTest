#include "fft/vkfft_vulkan_buffer.hpp"

#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)

#include <cstring>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace neo_dfttest::fft {
namespace {

void check_vk(VkResult result, const char* action) {
  if (result != VK_SUCCESS) {
    throw std::runtime_error(std::string(action) + ": Vulkan error " + std::to_string(static_cast<int>(result)));
  }
}

void validate_transfer_size(VkDeviceSize size, VkDeviceSize capacity, const char* action) {
  if (size > capacity) {
    throw std::runtime_error(std::string(action) + ": Vulkan transfer exceeds buffer size");
  }
}

VkDeviceSize vk_device_size(std::size_t size, const char* action) {
  if (size > std::numeric_limits<VkDeviceSize>::max()) {
    throw std::runtime_error(std::string(action) + ": Vulkan buffer size is too large");
  }
  return static_cast<VkDeviceSize>(size);
}

VkBuffer vk_view_buffer(DeviceBufferView view, const char* action) {
  if (view.handle == nullptr) {
    throw std::runtime_error(std::string(action) + ": missing Vulkan buffer");
  }
  return reinterpret_cast<VkBuffer>(view.handle);
}

VkDeviceSize vk_view_offset(DeviceBufferView view, const char* action) {
  const VkDeviceSize offset = vk_device_size(view.offset_bytes, action);
  if ((offset % 4u) != 0) {
    throw std::runtime_error(std::string(action) + ": Vulkan fill offset must be 4-byte aligned");
  }
  return offset;
}

VkDeviceSize vk_view_range(DeviceBufferView view, const char* action) {
  if (view.size_bytes == 0) {
    return VK_WHOLE_SIZE;
  }
  const VkDeviceSize range = vk_device_size(view.size_bytes, action);
  if ((range % 4u) != 0) {
    throw std::runtime_error(std::string(action) + ": Vulkan fill range must be a multiple of 4 bytes");
  }
  return range;
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

void submit_vulkan_commands(
  VulkanRuntime& runtime,
  void (*record)(VkCommandBuffer command_buffer, void* user),
  void* user
) {
  if (record == nullptr) {
    throw std::runtime_error("missing Vulkan command recorder");
  }

  std::lock_guard lock(runtime.submission_mutex());

  VkCommandBuffer command_buffer = VK_NULL_HANDLE;
  VkCommandBufferAllocateInfo allocate_info {VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO};
  allocate_info.commandPool = runtime.command_pool();
  allocate_info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
  allocate_info.commandBufferCount = 1;
  check_vk(vkAllocateCommandBuffers(runtime.device(), &allocate_info, &command_buffer), "allocate Vulkan command buffer");

  try {
    VkFence fence = runtime.fence();
    check_vk(vkResetFences(runtime.device(), 1, &fence), "reset Vulkan fence");

    VkCommandBufferBeginInfo begin_info {VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO};
    begin_info.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    check_vk(vkBeginCommandBuffer(command_buffer, &begin_info), "begin Vulkan command buffer");
    record(command_buffer, user);
    check_vk(vkEndCommandBuffer(command_buffer), "end Vulkan command buffer");

    VkSubmitInfo submit_info {VK_STRUCTURE_TYPE_SUBMIT_INFO};
    submit_info.commandBufferCount = 1;
    submit_info.pCommandBuffers = &command_buffer;
    check_vk(vkQueueSubmit(runtime.queue(), 1, &submit_info, fence), "submit Vulkan command buffer");
    check_vk(vkWaitForFences(runtime.device(), 1, &fence, VK_TRUE, std::numeric_limits<std::uint64_t>::max()), "wait for Vulkan fence");
  } catch (...) {
    vkFreeCommandBuffers(runtime.device(), runtime.command_pool(), 1, &command_buffer);
    throw;
  }

  vkFreeCommandBuffers(runtime.device(), runtime.command_pool(), 1, &command_buffer);
}

void clear_vulkan_buffer_view(VulkanRuntime& runtime, DeviceBufferView view, std::uint32_t value) {
  struct ClearRequest {
    DeviceBufferView view;
    std::uint32_t value;
  } request {view, value};

  submit_vulkan_commands(
    runtime,
    [](VkCommandBuffer command_buffer, void* user) {
      const auto& clear = *static_cast<ClearRequest*>(user);
      vkCmdFillBuffer(
        command_buffer,
        vk_view_buffer(clear.view, "clear Vulkan buffer view"),
        vk_view_offset(clear.view, "clear Vulkan buffer view"),
        vk_view_range(clear.view, "clear Vulkan buffer view"),
        clear.value
      );
    },
    &request
  );
}

void clear_vulkan_buffer(VulkanRuntime& runtime, VulkanDeviceBuffer& buffer, std::uint32_t value) {
  clear_vulkan_buffer_view(runtime, buffer.view(), value);
}

void upload_to_vulkan_buffer(
  VulkanRuntime& runtime,
  VulkanDeviceBuffer& destination,
  const void* source,
  VkDeviceSize bytes
) {
  upload_to_vulkan_buffer_view(runtime, destination.view(), source, bytes);
}

void upload_to_vulkan_buffer_view(
  VulkanRuntime& runtime,
  DeviceBufferView destination,
  const void* source,
  VkDeviceSize bytes
) {
  const VulkanBufferUpload upload {destination, source, bytes};
  upload_to_vulkan_buffer_views(runtime, std::span<const VulkanBufferUpload>{&upload, 1});
}

void upload_to_vulkan_buffer_views(VulkanRuntime& runtime, std::span<const VulkanBufferUpload> uploads) {
  if (uploads.empty()) {
    return;
  }

  std::vector<VulkanDeviceBuffer> staging_buffers;
  staging_buffers.reserve(uploads.size());
  for (const VulkanBufferUpload& upload : uploads) {
    const VkDeviceSize destination_range = vk_device_size(upload.destination.size_bytes, "upload to Vulkan buffer");
    validate_transfer_size(upload.bytes, destination_range, "upload to Vulkan buffer");
    if (upload.bytes == 0) {
      staging_buffers.emplace_back();
      continue;
    }
    if (upload.source == nullptr) {
      throw std::runtime_error("upload to Vulkan buffer received null source");
    }

    staging_buffers.emplace_back(
      runtime,
      upload.bytes,
      VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
      VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT
    );
    std::memcpy(staging_buffers.back().map(), upload.source, static_cast<std::size_t>(upload.bytes));
    staging_buffers.back().unmap();
  }

  struct CopyRequest {
    const std::vector<VulkanDeviceBuffer>* staging_buffers;
    const VulkanBufferUpload* uploads;
    std::size_t count;
  } request {&staging_buffers, uploads.data(), uploads.size()};

  submit_vulkan_commands(
    runtime,
    [](VkCommandBuffer command_buffer, void* user) {
      const auto& copy = *static_cast<CopyRequest*>(user);
      for (std::size_t index = 0; index < copy.count; ++index) {
        const VulkanBufferUpload& upload = copy.uploads[index];
        if (upload.bytes == 0) {
          continue;
        }

        VkBufferCopy region {};
        region.dstOffset = vk_device_size(upload.destination.offset_bytes, "upload to Vulkan buffer");
        region.size = upload.bytes;
        vkCmdCopyBuffer(
          command_buffer,
          (*copy.staging_buffers)[index].buffer(),
          vk_view_buffer(upload.destination, "upload to Vulkan buffer"),
          1,
          &region
        );
      }
    },
    &request
  );
}

void download_from_vulkan_buffer(
  VulkanRuntime& runtime,
  VulkanDeviceBuffer& source,
  void* destination,
  VkDeviceSize bytes
) {
  validate_transfer_size(bytes, source.size(), "download from Vulkan buffer");
  if (bytes == 0) {
    return;
  }
  if (destination == nullptr) {
    throw std::runtime_error("download from Vulkan buffer received null destination");
  }

  VulkanDeviceBuffer staging(
    runtime,
    bytes,
    VK_BUFFER_USAGE_TRANSFER_DST_BIT,
    VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT
  );

  struct CopyRequest {
    VulkanDeviceBuffer* source;
    VulkanDeviceBuffer* destination;
    VkDeviceSize bytes;
  } request {&source, &staging, bytes};

  submit_vulkan_commands(
    runtime,
    [](VkCommandBuffer command_buffer, void* user) {
      const auto& copy = *static_cast<CopyRequest*>(user);
      VkBufferCopy region {};
      region.size = copy.bytes;
      vkCmdCopyBuffer(command_buffer, copy.source->buffer(), copy.destination->buffer(), 1, &region);
    },
    &request
  );

  std::memcpy(destination, staging.map(), static_cast<std::size_t>(bytes));
  staging.unmap();
}

} // namespace neo_dfttest::fft

#endif
