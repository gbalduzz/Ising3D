// Simple timing utility. Do not use with multiple threads.
#pragma once

#include <chrono>

class Clock {
public:
  using Timestamp = typeof(std::chrono::high_resolution_clock::now());

  static void start() { t0_ = std::chrono::high_resolution_clock::now(); }

  static double stop() {
    Timestamp tf = std::chrono::high_resolution_clock::now();
    return std::chrono::duration_cast<std::chrono::duration<double>>(tf - t0_)
        .count();
  }

private:
  static Timestamp t0_;
};
