#ifndef ATOMIC_BUFFER_H
#define ATOMIC_BUFFER_H

#include <memory>
#include <atomic>
#include <utility>

template <typename T>
class atomic_buffer {
public:
  atomic_buffer(int n);
  ~atomic_buffer() = default;

  int size() const noexcept;
  bool empty() const noexcept;
  bool full() const noexcept;

  void put(const T & x, bool last) noexcept;
  std::pair<bool,T> get() noexcept;

private:
  int next_position(int p) const noexcept;

private:
  struct item {
    bool last;
    T value;
  };

  const int size_;
  std::unique_ptr<item[]> buf_;
  alignas(64) std::atomic<int> next_read_ {0};
  alignas(64) std::atomic<int> next_write_ {0};
};

template <typename T>
atomic_buffer<T>::atomic_buffer(int n) :
  size_{n},
  buf_{new item[size_]}
{
}

template <typename T>
int atomic_buffer<T>::size() const noexcept
{
  return size_;
}

template <typename T>
bool atomic_buffer<T>::empty() const noexcept
{
  return next_read_ == next_write_;
}

template <typename T>
bool atomic_buffer<T>::full() const noexcept
{
  const int next = next_position(next_write_.load());
  return next == next_read_.load();
}

template <typename T>
void atomic_buffer<T>::put(const T & x, bool last) noexcept
{
  const int next = next_position(next_write_.load());
  while (next == next_read_.load()) {
    ;
  }
  buf_[next_write_.load()] = item{last,x};
  next_write_.store(next);
}

template <typename T>
std::pair<bool,T> atomic_buffer<T>::get() noexcept
{
  while (empty()) {
    ;
  }
  auto res = buf_[next_read_.load()];
  next_read_.store(next_position(next_read_.load()));
  return std::make_pair(res.last,res.value);
}

template <typename T>
int atomic_buffer<T>::next_position(int p) const noexcept
{
  return p + ((p+1>=size_)?(1-size_):1);
}

#endif
