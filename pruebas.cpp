#include <iostream>
#include <thread>
#include <mutex>
#include <atomic>

std::mutex mtx;
std::atomic_flag atomFlag = ATOMIC_FLAG_INIT;
std::atomic<int> flag(0);

void one(int i){
        mtx.lock();
        std::cout << "method one-> " << i <<  '\n';
        mtx.unlock();
}
void two(int i){
        mtx.lock();
        std::cout << "method two-> " << i <<  '\n';
        mtx.unlock();
}
void three(int i){
        mtx.lock();
        std::cout << "method three-> " << i <<  '\n';
        mtx.unlock();
}

int main(){
        std::thread threads[3];

        for(int i=0; i<2; i++) {
                threads[0] = std::thread(one, i);
                threads[1] = std::thread(two, i);
                threads[2] = std::thread(three, i);
        }

        for(auto& th : threads) th.join();


}
