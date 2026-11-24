// Self-test of the solve timer: the same steady_clock and secs() as
// lpsolve.cpp and classic.cpp, checked for resolution, monotonicity
// and agreement with CLOCK_MONOTONIC_RAW, a clock the kernel never
// slews. Exit status 1 on any failure.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <initializer_list>

using Clock = std::chrono::steady_clock;
static double secs(Clock::time_point a) {
    return std::chrono::duration<double>(Clock::now() - a).count();
}

static double raw_now() {
    timespec t;
    clock_gettime(CLOCK_MONOTONIC_RAW, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

static bool check_resolution() {
    double best = 1;
    for (int i = 0; i < 1000000; i++) {
        auto a = Clock::now(), b = Clock::now();
        double d = std::chrono::duration<double>(b - a).count();
        if (d > 0 && d < best) best = d;
    }
    bool ok = best <= 1e-6;
    std::printf("resolution: %.0f ns %s\n", best * 1e9,
                ok ? "ok" : "FAIL");
    return ok;
}

static bool check_monotonic() {
    long back = 0;
    auto prev = Clock::now();
    for (long i = 0; i < 10000000; i++) {
        auto now = Clock::now();
        if (now < prev) back++;
        prev = now;
    }
    std::printf("monotonic: %ld backward steps in 1e7 reads %s\n",
                back, back ? "FAIL" : "ok");
    return back == 0;
}

static bool check_interval(double d) {
    double r0 = raw_now();
    auto t0 = Clock::now();
    while (raw_now() - r0 < d) {}
    double raw = raw_now() - r0, ours = secs(t0);
    double err = ours - raw;
    bool ok = std::fabs(err) <= 1e-5 + 1e-4 * raw;
    std::printf("interval %7.3f s: secs()=%.6f raw=%.6f err=%+.1e s"
                " %s\n", d, ours, raw, err, ok ? "ok" : "FAIL");
    return ok;
}

int main() {
    bool ok = check_resolution() & check_monotonic();
    for (double d : {0.001, 0.01, 0.1, 1.0, 5.0})
        ok &= check_interval(d);
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
