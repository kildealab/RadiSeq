#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

/*------------------------------------------------------------------------------------------------------------------------------------------
 This is a header-only implementation of the Mersenne Twister algorithm. You can use this in any code file using the 'rng' namespace.
 Random number generator can be initialized with either a fixed seed or a system generated random seed. If fixed, you are expected to
 pass an unsigned integer seed to the initWithFixedSeed() function. Otherwise, initRandomSeed() will initialize the rng with a seed sequence. 
*------------------------------------------------------------------------------------------------------------------------------------------*/
#include <chrono>
#include <random>


namespace rng{
    
    inline std::mt19937 mt;                                                                     // Initialize a Mersenne Twister rng object

    // Function to initialize mt19937 rng with random seeds
    inline std::mt19937 init(){
        std::random_device rd;                                                                  // Random device will ask the OS to generate random number
        std::seed_seq ss{                                                                       // Create a seed sequence with high-res clock and 7 random numbers from random device
            static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
            return std::mt19937{ss};
    }

    // User function to initialize rng with a fixed seed
    inline void initFixedSeed(unsigned int seed){
        mt = std::mt19937(seed);
    }

    // User function to initialize rng with a random seed
    inline void initRandomSeed(){
        mt = init();
    }

    // Function to generate a random integer number in given range, inclusive of the bounds
    inline int rand_int(int lowerBound,int upperBound){
        std::uniform_int_distribution<int> dist(lowerBound, upperBound);
        return dist(mt);
    }

    // Function to generate a random double number in given range, inclusive of the bounds
    inline double rand_double(double lowerBound,double upperBound){
        std::uniform_real_distribution<double> dist(lowerBound, upperBound);
        return dist(mt);
    }

    // This function calculates the probability of observing k or more successes in n independent Bernoulli trials, 
    // each with a success probability of p. It returns the complement of the cumulative distribution function (CDF), 
    // which represents the probability of observing a value greater than or equal to k.
    inline double binomial_cdf(size_t k, double p, size_t n){
        std::binomial_distribution<size_t> dist(n, p);
        size_t successes = 0;

        for (size_t i = 0; i < n; ++i) {
            if (dist(mt) >= k) successes++;
        }

        double probability = 1.0 - static_cast<double>(successes) / static_cast<double>(n);
        return probability;
    }

    // Function to generate a double value randomly in a gaussian (normal) distribution specified with mean and stdDev
    inline double normal_distribution(double mean, double stdDev){
        std::normal_distribution<double> dist(mean, stdDev);
        return dist(mt);
    }

}


# endif