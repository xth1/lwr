#ifndef RANDOM_GENERATOR
#define RANDOM_GENERATOR 1

#include "HeaderDefault.h"

/* Generator for uniform distributed numbers */
class RandomGenerator {
    /* std::random_device is a uniformly-distributed integer random number 
     * generator that produces non-deterministic random numbers.
     */
    random_device generator;
    uniform_int_distribution<int> distribution_integer;
    uniform_real_distribution<double> distribution_float;

    public:
    RandomGenerator();
    double getFloatPoint01();
    double getFloatPoint(double begin, double end);
    int getInteger(int begin, int end);
};
#endif