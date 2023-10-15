#pragma once

#include <chrono>
#include "logger.h"

NAMESPACE_PMT

class Timer
{
public:
    inline Timer() : stopped_(false)
    {
        start();  
    }

    inline Timer(std::string&& msg) :
        msg_(msg), stopped_(false)
    {
        start();
    }

    inline void start()
    {
        start_ = std::chrono::high_resolution_clock::now();
    }

    inline double stop()
    {
        stopped_ = true;

        std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now() - start_;
 
        return diff.count();
    }

    inline ~Timer()
    {
        if (stopped_) return;

        out(std::fixed << stop() << " s: " << msg_);
    }

private:
    std::string msg_;
    std::chrono::time_point<std::chrono::high_resolution_clock> start_;
    bool stopped_;
};

NAMESPACE_PMT_END
