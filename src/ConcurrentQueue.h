//
// Copyright (c) 2013 Juan Palacios juan.palacios.puyana@gmail.com
// Subject to the BSD 2-Clause License
// - see < http://opensource.org/licenses/BSD-2-Clause>
// Based on https://github.com/juanchopanza/cppblog/tree/master/Concurrency/Queue
//

#ifndef CONCURRENT_DEQUE_
#define CONCURRENT_DEQUE_

#include <deque>
#include <thread>
#include <mutex>
#include <condition_variable>

template <typename T>
class ConcurrentDeque
{
 public:
  T pop() 
  {
    std::unique_lock<std::mutex> mlock(mutex_);
    while (deque_.empty())
    {
      cond_.wait(mlock);
    }
    auto val = deque_.front();
    deque_.pop_front();
    return val;
  }
  void push_back(const T item)
  {
    std::unique_lock<std::mutex> mlock(mutex_);
    deque_.push_back(item);
    mlock.unlock();
    cond_.notify_one();
  }
  void push_front(const T item)
  {
    std::unique_lock<std::mutex> mlock(mutex_);
    deque_.push_front(item);
    mlock.unlock();
    cond_.notify_one();
  }
  void push_nolock(const T item)
  {
    deque_.push_back(item);
  }
  ConcurrentDeque()=default;
  ConcurrentDeque(const ConcurrentDeque&) = delete;            // disable copying
  ConcurrentDeque& operator=(const ConcurrentDeque&) = delete; // disable assignment
  
 private:
  std::deque<T> deque_;
  std::mutex mutex_;
  std::condition_variable cond_;
};

#endif
