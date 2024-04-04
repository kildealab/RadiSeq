#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

/*------------------------------------------------------------------------------------------------------------------------------------------
 This is a header-only implementation of the Mersenne Twister algorithm. You can use this in any code file using the 'rng' namespace.
 Random number generator can be initialized with either a fixed seed or a system generated random seed. If fixed, you are expected to
 pass an unsigned integer seed to the initWithFixedSeed() function. Otherwise, initRandomSeed() will initialize the rng with a seed sequence.
 For thread-local implementation, use initThreadLocalRandomSeed() and initThreadLocalFixedSeed functions respectively.
-------------------------------------------------------------------------------------------------------------------------------------------*/
#include <chrono>
#include <random>
#include <vector>
#include <type_traits>


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
    // Template is used to switch between int, long, short etc. types
    template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
    inline T rand_int(T lowerBound, T upperBound, int threadID = 0){                            // threadID gets a default value of 0 (for master thread)
        std::mt19937& generator = local_mt[threadID];                                           // Switch to the generator engine instance accordingly based on threadID
        std::uniform_int_distribution<T> dist(lowerBound, upperBound);
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
    /* inline double binomial_cdf(size_t k, double p, size_t n, int threadID = 0){
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::binomial_distribution<size_t> dist(n, p);
        size_t successes = 0;

        for (size_t i = 0; i < n; ++i){                                                          // Perform n independent trials
            if (dist(generator) >= k) successes++;                                               // Find how many were successes
        }

        double probability = 1.0 - static_cast<double>(successes+1) / static_cast<double>(n);    // Probability of observing at most k successes (1 - probability of observing at least k+1 successes)
        return probability;
    } */

    inline size_t binomial_distribution(double p, size_t n, int threadID = 0){
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::binomial_distribution<size_t> dist(n, p);
        return dist(generator);                                                                  // Returns number of successes in n number of trials with a success rate p
    }
    // Function to generate a double value randomly in a gaussian (normal) distribution specified with mean and stdDev
    // Third optional argument can be used to switch to thread-local RNG engine if needed
    /*inline double normal_distribution(double mean, double stdDev, int threadID = 0){
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::normal_distribution<double> dist(mean, stdDev);
        return dist(generator);
    }*/

    // Function for weighted random integer selection. The sum of all weights needs to be 1 and 
    // the returned value will be in the range [1, weights.size()]
    // Third optional argument can be used to switch to thread-local RNG engine if needed
    inline int weighted_rand_int(const std::vector<double>& weights, int threadID = 0){
        std::mt19937& generator = local_mt[threadID];                                            // Switch to the generator engine instance accordingly based on threadID
        std::uniform_real_distribution<double> dist(0.0, 1.0);                                   // Pick a random real number between 0 and 1 
        double randomNumber = dist(generator);

        double cumulativeProbability = 0.0;
        for (size_t i = 0; i < weights.size(); ++i) {
            cumulativeProbability += weights[i];
            if (randomNumber <= cumulativeProbability){
                return (static_cast<int>(i)+1);
            }
        }
        return(0);                                                                               // Properly weighted data should not execute this return
    }

}


# endif