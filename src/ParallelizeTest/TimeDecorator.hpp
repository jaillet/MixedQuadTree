#include <functional>
#include <chrono>
#include <iostream>
#include <time.h>

//https://gist.github.com/lizhongz/21e77864aa9d5f66e842
//Modified a little bit

//Only for void functions

/* This file provide a decorator that measure the execution time of a function.

Use it with make_decorator, or count_time. Examples below 
*/

//  Windows
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

//  Posix/Linux
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC; // regarder le temps CPU
}
#endif

template <class> struct ExeTime;

// Execution time decorator
template <class... Args>
struct ExeTime<void(Args ...)> {
public:
    ExeTime(std::function<void(Args...)> func): f_(func) { }

    double operator ()(Args ... args) {
        struct timespec start, finish;
        double elapsed;

        clock_gettime(CLOCK_MONOTONIC, &start);
        f_(args...);
        clock_gettime(CLOCK_MONOTONIC, &finish);

        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

        return elapsed * 1000;
    }

    /*
    float operator ()(Args ... args) {
        auto start_time = std::chrono::high_resolution_clock::now();
        f_(args...);    
        auto end_time = std::chrono::high_resolution_clock::now();

        return std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    }
     */

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