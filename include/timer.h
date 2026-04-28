#ifndef VDPSU_TIMER_H
#define VDPSU_TIMER_H

#include <chrono>
#include <string>

namespace vdpsu {

class Timer {
public:
    void Start() {
        start_ = std::chrono::high_resolution_clock::now();
    }

    double StopMs() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::milli>(end - start_).count();
    }

    double StopUs() {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::micro>(end - start_).count();
    }

private:
    std::chrono::high_resolution_clock::time_point start_;
};

} // namespace vdpsu

#endif // VDPSU_TIMER_H
