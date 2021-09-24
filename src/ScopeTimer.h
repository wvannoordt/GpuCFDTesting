#ifndef SCOPE_TIMER_H
#define SCOPE_TIMER_H
#include <string>
#include <chrono>
#include "print.h"
struct ScopeTimer
{
    std::string task;
    decltype(std::chrono::high_resolution_clock::now()) start;
    ScopeTimer(std::string task_in)
    {
        task = task_in;
        start = std::chrono::high_resolution_clock::now();
    }
    ~ScopeTimer(void)
    {
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> ms_double = end - start;
        print(task, ms_double.count(), "ms");
    }
};

#endif