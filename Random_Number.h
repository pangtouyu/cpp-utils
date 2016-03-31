/*
 *  Title: 	randomnumber
 *  Last Mod: 	Fri Mar 18 07:28:57 1988
 *  Author: 	Vincent Broman
 *		<broman@schroeder.nosc.mil>
 */  

void start_random_number(int seed_a, int seed_b);
double next_random_number(void);
    
/*  
 *  This package makes available Marsaglia's highly portable generator 
 *  of uniformly distributed pseudo-random numbers.
 *  
 *  The sequence of 24 bit pseudo-random numbers produced has a period 
 *  of about 2**144, and has passed stringent statistical tests 
 *  for randomness and independence.
 *  
 *  Supplying two seeds to start_random_number is required once
 *  at program startup before requesting any random numbers, like this:
 *      start_random_number(101, 202);
 *      r = next_random_number();
 *  The correspondence between pairs of seeds and generated sequences 
 *  of pseudo-random numbers is many-to-one.
 *  
 *  This package should compile and run identically on any 
 *  machine/compiler which supports >=16 bit integer arithmetic
 *  and >=24 bit floating point arithmetic.
 *  
 *  References:
 *      M G Harmon & T P Baker, ``An Ada Implementation of Marsaglia's
 *      "Universal" Random Number Generator'', Ada Letters, late 1987.
 *      
 *      G Marsaglia, ``Toward a universal random number generator'',
 *      to appear in the Journal of the American Statistical Association.
 *  
 *  George Marsaglia is at the Supercomputer Computations Research Institute
 *  at Florida State University.
 */  
