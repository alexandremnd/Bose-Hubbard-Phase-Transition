#pragma once

#include <condition_variable>
#include <mutex>
#include <vector>
#include <thread>
#include <queue>
#include <functional>

class ThreadPool {
    public:
        ThreadPool() {}

        void start(int num_threads) {
            if (num_threads == -1) {
                num_threads = std::thread::hardware_concurrency();
            }

            done = false;

            for (int i = 0; i < num_threads; ++i) {
                threads.emplace_back(std::thread(&ThreadPool::thread_loop, this));
            }
        }

        void thread_loop() {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->tasks_mutex);
                    this->condition.wait(lock, [this]() {
                        return this->done || !this->tasks.empty();
                    });

                    if (this->done && this->tasks.empty()) {
                        return;
                    }

                    task = std::move(this->tasks.front());
                    this->tasks.pop();
                }
                task();
            }
        }

        void queue_job(std::function<void()> job) {
            {
                std::unique_lock<std::mutex> lock(this->tasks_mutex);
                this->tasks.push(job);
            }
            this->condition.notify_one();
        }

        bool busy() {
            bool is_busy;
            {
                std::unique_lock<std::mutex> lock(this->tasks_mutex);
                is_busy = !tasks.empty();
            }
            return is_busy;
        }

        void stop() {
            {
                std::unique_lock<std::mutex> lock(this->tasks_mutex);
                done = true;
            }
            condition.notify_all();

            for (std::thread &thread : threads) {
                thread.join();
            }
            threads.clear();
        }

        ~ThreadPool() {
            stop();
        };

    private:
        std::vector<std::thread> threads;
        std::queue<std::function<void()>> tasks;
        std::mutex tasks_mutex;
        std::condition_variable condition;
        bool done;
};