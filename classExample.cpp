#include <iostream>
#include <thread>
#include <mutex>
#include <atomic>

std::mutex mtx;
std::atomic_flag atomFlag = ATOMIC_FLAG_INIT;
std::atomic<int> flag(0);

void print_thread_id(int id){
        int expected =0;
        while(flag.compare_exchange_strong(expected, 1)) expected=0;
        // mtx.lock();
        // while(atomFlag.test_and_set()) ;
        std::cout << "thread @" << id << '\n';
        // atomFlag.clear();
        // mtx.unlock();
}

int main(){
        std::thread threads[10];

        for(int i=0; i<10; i++) threads[i] = std::thread(print_thread_id, i+1);

        for(auto& th : threads) th.join();


}
