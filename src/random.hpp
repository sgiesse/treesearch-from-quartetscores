#ifndef RANDOM_HPP
#define RANDOM_HPP

namespace {
    std::mt19937 mt;
    bool initialized = false;
}

namespace Random {

    void init() {
        std::random_device rd;
        mt = std::mt19937(rd());
        initialized = true;
    }

    void seed(int s) {
        mt.seed(s);
        initialized = true;
    }

    int get_rand_int(int a, int b) {
        if (!initialized) {
            init();
        }
        std::uniform_int_distribution<int> distr(a,b);
        return distr(mt);
    }

    float get_rand_float(float a, float b) {
        if (!initialized) {
            init();
        }
        std::uniform_real_distribution<> distr(a,b);
        return distr(mt);
    }

    std::mt19937 getMT() {
        if (!initialized) init();
        return mt;
    }
}

#endif
