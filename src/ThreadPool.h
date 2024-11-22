#ifndef NANOFQ_THREADPOOL_H
#define NANOFQ_THREADPOOL_H

#include <queue>
#include <mutex>
#include <thread>
#include <functional>
#include <condition_variable>
#include <future>
#include <iostream>
#include <memory>


class ThreadPool
{
private:
    std::vector<std::thread> m_threads;
    std::queue<std::function<void()>> m_tasks;
    std::mutex m_mtx;
    std::condition_variable m_cv;
    std::atomic<bool> m_running;

public:
    ThreadPool() = delete;
    ThreadPool(const ThreadPool&) = delete;
    ThreadPool(ThreadPool&&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool&&) = delete;

    explicit ThreadPool(int threads_number) : m_running(true)
    {
        if (threads_number < 1) {
            std::cerr << "Thread number must >= 1" << std::endl;
            exit(1);
        }
        m_threads.reserve(threads_number);
        for (int i{0}; i < threads_number; i++) {
            m_threads.emplace_back([this](){
                while (true) {
                    std::function<void()> task;
                    {
                        std::unique_lock lock{m_mtx};
                        if (m_running && !m_tasks.empty()) {
                            task = std::move(m_tasks.front());
                            m_tasks.pop();
                        } else {
                            m_cv.wait(lock, [this](){ return !m_running && !m_tasks.empty(); });
                            if (!m_running) return;
                            task = std::move(m_tasks.front());
                            m_tasks.pop();
                        }
                    }
                    task();
                }
            });
        }
    }

    ~ThreadPool()
    {
        m_running = false;
        m_cv.notify_all();
        for (auto& t : m_threads) {
            if (t.joinable()) t.join();
        }
    }

    template <typename Func, typename... Args>
    auto enqueue(Func&& func, Args&&... args) -> std::future<decltype(func(args...))>
    {
        using return_type = decltype(func(args...));
        auto function = std::make_shared<std::packaged_task<return_type()>>(
            std::bind(std::forward<Func>(func), std::forward<Args>(args)...)
        );
        std::future<return_type> future{function->get_future()};
        {
            std::unique_lock lock{m_mtx};
            m_tasks.emplace([function](){ (*function)(); });
        }
        m_cv.notify_one();
        return future;
    }
};


#endif //NANOFQ_THREADPOOL_H
