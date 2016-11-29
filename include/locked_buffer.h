#ifndef LOCKED_BUFFER_H
#define LOCKED_BUFFER_H

#include <memory>
#include <mutex>
#include <condition_variable>
#include <atomic>

template <typename T>
class locked_buffer {
public:
        locked_buffer(int n);
        locked_buffer(const locked_buffer &) = delete;
        ~locked_buffer() = default;

        int size() const noexcept;
        bool empty() const noexcept;
        bool full() const noexcept;

        void put(const T & x, bool last) noexcept;
        std::pair<bool,T> get() noexcept;

private:
        int next_position(int p) const noexcept;
        bool do_empty() const noexcept;
        bool do_full() const noexcept;

private:
        struct item {
                bool last;
                T value;
        };

        const int size_;
        std::unique_ptr<item[]> buf_;
        int next_read_ = 0;
        int next_write_ = 0;

        mutable std::mutex mut_;
        std::condition_variable not_full_;
        std::condition_variable not_empty_;
};

template <typename T>
locked_buffer<T>::locked_buffer(int n) :
        size_ {
        n
},
buf_ {new item[size_]}
{
}

template <typename T>
int locked_buffer<T>::size() const noexcept
{
        return size_;
}

template <typename T>
bool locked_buffer<T>::empty() const noexcept
{
        using namespace std;
        lock_guard<mutex> l {mut_};
        return do_empty();
}

template <typename T>
bool locked_buffer<T>::full() const noexcept
{
        using namespace std;
        lock_guard<mutex> l {mut_};
        return do_full();
}

template <typename T>
void locked_buffer<T>::put(const T & x, bool last) noexcept
{
        using namespace std;
        unique_lock<mutex> l {mut_};
        not_full_.wait(l, [this] { return !do_full(); });
        buf_[next_write_] = item {last,x};
        next_write_ = next_position(next_write_);
        l.unlock();
        not_empty_.notify_one();
}

template <typename T>
std::pair<bool,T> locked_buffer<T>::get() noexcept
{
        using namespace std;
        unique_lock<mutex> l {mut_};
        not_empty_.wait(l, [this] { return !do_empty(); });
        auto res = buf_[next_read_];
        next_read_ = next_position(next_read_);
        l.unlock();
        not_full_.notify_one();
        return make_pair(res.last,res.value);
}

template <typename T>
int locked_buffer<T>::next_position(int p) const noexcept
{
        return p + ((p+1>=size_) ? (1-size_) : 1);
}

template <typename T>
bool locked_buffer<T>::do_empty() const noexcept
{
        return (next_read_ == next_write_);
}

template <typename T>
bool locked_buffer<T>::do_full() const noexcept
{
        const int next = next_position(next_write_);
        return next == next_read_;
}

#endif
