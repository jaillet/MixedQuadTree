#include <functional>
#include <chrono>
#include <iostream>

//https://gist.github.com/lizhongz/21e77864aa9d5f66e842
//Modified a little bit

//Only for void functions

/* This file provide a decorator that measure the execution time of a function.

Use it with make_decorator, or count_time. Examples below 
*/

template <class> struct ExeTime;

// Execution time decorator
template <class... Args>
struct ExeTime<void(Args ...)> {
public:
    ExeTime(std::function<void(Args...)> func): f_(func) { } 

    float operator ()(Args ... args) {
        auto start_time = std::chrono::high_resolution_clock::now();
        f_(args...);    
        auto end_time = std::chrono::high_resolution_clock::now();
        
        return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    }   

private:
    std::function<void(Args ...)> f_; 
};

template <class... Args>
ExeTime<void(Args ...)> make_decorator(void (*f)(Args ...)) {
    return ExeTime<void(Args...)>(std::function<void(Args...)>(f));    
}

template <class F, class... Args>
float count_time(F&& f, Args&& ... args) {
    auto et = make_decorator(f);
    return et(std::forward<Args>(args)...);
}


void ex(int a) {
    int b = 0;
    for (int i = 0; i < a; ++i)
    {
        //Do something...
        b += i;
    }
}
void exemple_utilisation() {
    float time = count_time(ex, 5);
}