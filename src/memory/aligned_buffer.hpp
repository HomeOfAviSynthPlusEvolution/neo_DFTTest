#ifndef NEO_DFTTEST_ALIGNED_BUFFER_HPP
#define NEO_DFTTEST_ALIGNED_BUFFER_HPP

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <new>

#ifndef FRAME_ALIGN
#define FRAME_ALIGN 64
#endif

#ifndef _WIN32
inline void* dfttest_aligned_malloc(std::size_t size, std::size_t alignment) {
  void* ptr = nullptr;
  return posix_memalign(&ptr, alignment, size) == 0 ? ptr : nullptr;
}
#define _aligned_malloc(a,b) dfttest_aligned_malloc((a),(b))
#define _aligned_free free
#endif

namespace neo_dfttest {

struct AlignedFree {
  void operator()(void* p) const noexcept {
    _aligned_free(p);
  }
};

template <class T>
class AlignedBuffer {
public:
  AlignedBuffer() = default;

  explicit AlignedBuffer(std::size_t count) {
    reset(count);
  }

  AlignedBuffer(const AlignedBuffer&) = delete;
  AlignedBuffer& operator=(const AlignedBuffer&) = delete;
  AlignedBuffer(AlignedBuffer&&) noexcept = default;
  AlignedBuffer& operator=(AlignedBuffer&&) noexcept = default;

  void reset(std::size_t count = 0) {
    if (count == 0) {
      data_.reset();
      size_ = 0;
      return;
    }

    auto* raw = static_cast<T*>(_aligned_malloc(sizeof(T) * count, FRAME_ALIGN));
    if (!raw) {
      throw std::bad_alloc();
    }

    data_.reset(raw);
    size_ = count;
  }

  T* data() noexcept {
    return data_.get();
  }

  const T* data() const noexcept {
    return data_.get();
  }

  std::size_t size() const noexcept {
    return size_;
  }

  explicit operator bool() const noexcept {
    return data_ != nullptr;
  }

  T& operator[](std::size_t index) noexcept {
    return data_.get()[index];
  }

  const T& operator[](std::size_t index) const noexcept {
    return data_.get()[index];
  }

private:
  std::unique_ptr<T, AlignedFree> data_;
  std::size_t size_ = 0;
};

} // namespace neo_dfttest

#endif
