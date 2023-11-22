#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

/*------------------------------------------------------------------------------------------------------------------------------------------
 This is a header-only implementation of the Mersenne Twister algorithm. You can use this in any code file using the 'rng' namespace.
 Random number generator can be initialized with either a fixed seed or a system generated random seed. If fixed, you are expected to
 pass an unsigned integer seed to the initWithFixedSeed() function. Otherwise, initRandomSeed() will initialize the rng with a seed sequence.
 For thread-local implementation, use initThreadLocalRandomSeed() and initThreadLocalFixedSeed functions respectively.
*------------------------------------------------------------------------------------------------------------------------------------------*/
#include <chrono>
#include <random>
#include <vector>


namespace rng{
    
    inline std::mt19937 mt;                                                                     // Initialize a Mersenne Twister rng object
    inline std::vector<std::mt19937> local_mt;                                                  // Initialize a vector to hold Mersenne Twister rng object for each local thread

    // Function to initialize mt19937 rng with random seeds
    inline std::mt19937 init(){
        std::random_device rd;                                                                  // Random device will ask the OS to generate random number
        std::seed_seq ss{                                                                       // Create a seed sequence with high-res clock and 7 random numbers from random device
            static_cast<unsigned int>(std::chrono::steady_clock::now().time_since_epoch().count()), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
            return std::mt19937{ss};
    }

    // Initialize the thread-local RNG with a fixed seed
    inline void initThreadLocalFixedSeed(unsigned int seed, int nThreads){
        local_mt.resize(nThreads);                                                               // Create a vector to hold nThread instances of the random engine
        for (int i=0; i<nThreads; i++){                                                          // Iteration starts from 1 because 0 correponds to the master thread
            local_mt[i] = std::mt19937(seed);
        }
    }

    // Initialize the thread-local RNG with random seeds
    inline void initThreadLocalRandomSeed(int nThreads){
        local_mt.resize(nThreads);                                                               // Create a vector to hold nThread instances of the random engine
        for (int i=0; i<nThreads; i++){                                                          // Iteration starts from 1 because 0 correponds to the master thread
            local_mt[i] = init();
        }
    }

    // Function to generate a random integer number in given range, inclusive of the bounds
    // Third optional argument can be used to switch to thread-local RNG engine if needed
    inline int rand_int(int lowerBound, int upperBound, int threadID = 0){                       // threadID gets a default value of 0 (for master thread)
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::uniform_int_distribution<int> dist(lowerBound, upperBound);
        return dist(generator);
    }

    // Function to generate a random double number in given range, inclusive of the bounds
    // Third optional argument can be used to switch to thread-local RNG engine if needed
    inline double rand_double(double lowerBound, double upperBound, int threadID = 0){
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::uniform_real_distribution<double> dist(lowerBound, upperBound);
        return dist(generator);
    }

    // This function calculates the probability of observing k or more successes in n independent Bernoulli trials, 
    // each with a success probability of p. It returns the complement of the cumulative distribution function (CDF), 
    // which represents the probability of observing a value greater than or equal to k.
    // Fourth optional argument can be used to switch to thread-local RNG engine if needed
    inline double binomial_cdf(size_t k, double p, size_t n, int threadID = 0){
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::binomial_distribution<size_t> dist(n, p);
        size_t successes = 0;

        for (size_t i = 0; i < n; ++i) {
            if (dist(generator) >= k) successes++;
        }

        double probability = 1.0 - static_cast<double>(successes) / static_cast<double>(n);
        return probability;
    }

    // Function to generate a double value randomly in a gaussian (normal) distribution specified with mean and stdDev
    // Third optional argument can be used to switch to thread-local RNG engine if needed
    inline double normal_distribution(double mean, double stdDev, int threadID = 0){
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::normal_distribution<double> dist(mean, stdDev);
        return dist(generator);
    }

}


# endif