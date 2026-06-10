#include "vulkan/vulkan_compute.hpp"

#if defined(NEO_DFTTEST_ENABLE_VKFFT_VULKAN)

#include "fft/vkfft_vulkan_buffer.hpp"

#include <limits>
#include <stdexcept>
#include <string>

namespace neo_dfttest::vulkan {
namespace {

void check_vk(VkResult result, const char* action) {
  if (result != VK_SUCCESS) {
    throw std::runtime_error(std::string(action) + ": Vulkan error " + std::to_string(static_cast<int>(result)));
  }
}

VkBuffer vk_buffer(fft::DeviceBufferView view, const char* name) {
  if (view.handle == nullptr) {
    throw std::runtime_error(std::string("missing Vulkan buffer for ") + name);
  }
  return reinterpret_cast<VkBuffer>(view.handle);
}

VkDeviceSize vk_offset(std::size_t offset_bytes, const char* name) {
  if (offset_bytes > std::numeric_limits<VkDeviceSize>::max()) {
    throw std::runtime_error(std::string("Vulkan buffer offset is too large for ") + name);
  }
  return static_cast<VkDeviceSize>(offset_bytes);
}

VkDeviceSize vk_range(fft::DeviceBufferView view) noexcept {
  return view.size_bytes == 0 ? VK_WHOLE_SIZE : static_cast<VkDeviceSize>(view.size_bytes);
}

} // namespace

ComputePipeline::ComputePipeline(
  fft::VulkanRuntime& runtime,
  std::span<const std::uint32_t> spirv,
  std::uint32_t storage_binding_count,
  std::uint32_t push_constant_bytes
)
  : runtime_(&runtime),
    storage_binding_count_(storage_binding_count),
    push_constant_bytes_(push_constant_bytes) {
  if (spirv.empty()) {
    throw std::runtime_error("cannot create Vulkan compute pipeline from empty SPIR-V");
  }
  if (storage_binding_count_ == 0) {
    throw std::runtime_error("Vulkan compute pipeline requires at least one storage binding");
  }

  VkShaderModuleCreateInfo shader_info {VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO};
  shader_info.codeSize = spirv.size_bytes();
  shader_info.pCode = spirv.data();
  check_vk(vkCreateShaderModule(runtime_->device(), &shader_info, nullptr, &shader_), "create Vulkan shader module");

  std::vector<VkDescriptorSetLayoutBinding> bindings(storage_binding_count_);
  for (std::uint32_t index = 0; index < storage_binding_count_; ++index) {
    bindings[index].binding = index;
    bindings[index].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    bindings[index].descriptorCount = 1;
    bindings[index].stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
  }

  VkDescriptorSetLayoutCreateInfo descriptor_layout_info {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO};
  descriptor_layout_info.bindingCount = static_cast<std::uint32_t>(bindings.size());
  descriptor_layout_info.pBindings = bindings.data();
  check_vk(
    vkCreateDescriptorSetLayout(runtime_->device(), &descriptor_layout_info, nullptr, &descriptor_set_layout_),
    "create Vulkan descriptor set layout"
  );

  VkPushConstantRange push_range {};
  push_range.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
  push_range.offset = 0;
  push_range.size = push_constant_bytes_;

  VkPipelineLayoutCreateInfo pipeline_layout_info {VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO};
  pipeline_layout_info.setLayoutCount = 1;
  pipeline_layout_info.pSetLayouts = &descriptor_set_layout_;
  if (push_constant_bytes_ > 0) {
    pipeline_layout_info.pushConstantRangeCount = 1;
    pipeline_layout_info.pPushConstantRanges = &push_range;
  }
  check_vk(
    vkCreatePipelineLayout(runtime_->device(), &pipeline_layout_info, nullptr, &pipeline_layout_),
    "create Vulkan compute pipeline layout"
  );

  VkPipelineShaderStageCreateInfo stage_info {VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO};
  stage_info.stage = VK_SHADER_STAGE_COMPUTE_BIT;
  stage_info.module = shader_;
  stage_info.pName = "main";

  VkComputePipelineCreateInfo pipeline_info {VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO};
  pipeline_info.stage = stage_info;
  pipeline_info.layout = pipeline_layout_;
  check_vk(
    vkCreateComputePipelines(runtime_->device(), VK_NULL_HANDLE, 1, &pipeline_info, nullptr, &pipeline_),
    "create Vulkan compute pipeline"
  );

  VkDescriptorPoolSize pool_size {};
  pool_size.type = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
  pool_size.descriptorCount = storage_binding_count_;

  VkDescriptorPoolCreateInfo pool_info {VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO};
  pool_info.maxSets = 1;
  pool_info.poolSizeCount = 1;
  pool_info.pPoolSizes = &pool_size;
  check_vk(vkCreateDescriptorPool(runtime_->device(), &pool_info, nullptr, &descriptor_pool_), "create Vulkan descriptor pool");

  VkDescriptorSetAllocateInfo allocate_info {VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO};
  allocate_info.descriptorPool = descriptor_pool_;
  allocate_info.descriptorSetCount = 1;
  allocate_info.pSetLayouts = &descriptor_set_layout_;
  check_vk(vkAllocateDescriptorSets(runtime_->device(), &allocate_info, &descriptor_set_), "allocate Vulkan descriptor set");
}

ComputePipeline::~ComputePipeline() {
  if (!runtime_) {
    return;
  }
  if (descriptor_pool_ != VK_NULL_HANDLE) {
    vkDestroyDescriptorPool(runtime_->device(), descriptor_pool_, nullptr);
  }
  if (pipeline_ != VK_NULL_HANDLE) {
    vkDestroyPipeline(runtime_->device(), pipeline_, nullptr);
  }
  if (pipeline_layout_ != VK_NULL_HANDLE) {
    vkDestroyPipelineLayout(runtime_->device(), pipeline_layout_, nullptr);
  }
  if (descriptor_set_layout_ != VK_NULL_HANDLE) {
    vkDestroyDescriptorSetLayout(runtime_->device(), descriptor_set_layout_, nullptr);
  }
  if (shader_ != VK_NULL_HANDLE) {
    vkDestroyShaderModule(runtime_->device(), shader_, nullptr);
  }
}

void ComputePipeline::dispatch(
  std::span<const fft::DeviceBufferView> storage_buffers,
  const void* push_constants,
  std::uint32_t push_constant_bytes,
  std::uint32_t groups_x,
  std::uint32_t groups_y,
  std::uint32_t groups_z
) const {
  if (storage_buffers.size() != storage_binding_count_) {
    throw std::runtime_error("Vulkan compute dispatch received the wrong number of storage buffers");
  }
  if (push_constant_bytes != push_constant_bytes_) {
    throw std::runtime_error("Vulkan compute dispatch received the wrong push constant size");
  }
  if (push_constant_bytes_ > 0 && push_constants == nullptr) {
    throw std::runtime_error("Vulkan compute dispatch received null push constants");
  }

  std::vector<VkDescriptorBufferInfo> buffer_infos(storage_buffers.size());
  std::vector<VkWriteDescriptorSet> writes(storage_buffers.size());
  for (std::uint32_t index = 0; index < storage_binding_count_; ++index) {
    buffer_infos[index].buffer = vk_buffer(storage_buffers[index], "compute storage binding");
    buffer_infos[index].offset = vk_offset(storage_buffers[index].offset_bytes, "compute storage binding");
    buffer_infos[index].range = vk_range(storage_buffers[index]);

    writes[index].sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
    writes[index].dstSet = descriptor_set_;
    writes[index].dstBinding = index;
    writes[index].descriptorCount = 1;
    writes[index].descriptorType = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    writes[index].pBufferInfo = &buffer_infos[index];
  }

  struct DispatchRequest {
    const ComputePipeline* pipeline;
    const VkWriteDescriptorSet* writes;
    std::uint32_t write_count;
    const void* push_constants;
    std::uint32_t push_constant_bytes;
    std::uint32_t groups_x;
    std::uint32_t groups_y;
    std::uint32_t groups_z;
  } request {
    this,
    writes.data(),
    static_cast<std::uint32_t>(writes.size()),
    push_constants,
    push_constant_bytes,
    groups_x,
    groups_y,
    groups_z
  };

  fft::submit_vulkan_commands(
    *runtime_,
    [](VkCommandBuffer command_buffer, void* user) {
      const auto& request = *static_cast<DispatchRequest*>(user);
      const ComputePipeline& pipeline = *request.pipeline;
      vkUpdateDescriptorSets(
        pipeline.runtime_->device(),
        request.write_count,
        request.writes,
        0,
        nullptr
      );
      vkCmdBindPipeline(command_buffer, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline.pipeline_);
      vkCmdBindDescriptorSets(
        command_buffer,
        VK_PIPELINE_BIND_POINT_COMPUTE,
        pipeline.pipeline_layout_,
        0,
        1,
        &pipeline.descriptor_set_,
        0,
        nullptr
      );
      if (request.push_constant_bytes > 0) {
        vkCmdPushConstants(
          command_buffer,
          pipeline.pipeline_layout_,
          VK_SHADER_STAGE_COMPUTE_BIT,
          0,
          request.push_constant_bytes,
          request.push_constants
        );
      }
      vkCmdDispatch(command_buffer, request.groups_x, request.groups_y, request.groups_z);
    },
    &request
  );
}

} // namespace neo_dfttest::vulkan

#endif
