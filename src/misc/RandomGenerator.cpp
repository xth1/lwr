#include "RandomGenerator.h"

RandomGenerator::RandomGenerator() {
    double float_upper_bound = 1.0 + numeric_limits<double>::epsilon();
    // build the float point distribution: generates a number in the range
    // [lower_bound, upper_bound)
    this->distribution_float = uniform_real_distribution<double> (
        0.0, //lower_bound
        float_upper_bound
    );

    // build the integer distribution
    this->distribution_integer = uniform_int_distribution<int> (0, INT_MAX);
}

// Produces a float-point number in the range [0,1]
double RandomGenerator::getFloatPoint01() {
    double random_num =  this->distribution_float(this->generator);
    return min(random_num, 1.0);
}

// Produces a float-point number in the range [begin, end]
double RandomGenerator::getFloatPoint(double begin, double end) {
    assert(begin >= 0);
    assert(begin < end);
    double random_num =  this->getFloatPoint01();
    double diff = (end - begin);
    return random_num * diff + begin;
}

// Produces a float-point number in the range [begin, end]
int RandomGenerator::getInteger(int begin, int end) {
    assert(begin >= 0);
    assert(begin < end);
    int random_num = this->distribution_integer(this->generator);
    int diff = end - begin;
    return random_num % (diff + 1) + begin;
}
