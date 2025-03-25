/* Header file for uniform random number generator 09/02/06

                $$$$$$$$$$$$$$$$$$$$$$$
                $   URand.h - header  $ 
                $$$$$$$$$$$$$$$$$$$$$$$

   by W.B. Yates
   Copyright (c) W.B. Yates. All rights reserved.
   History:
 
   Uniform Distribution (0,1)
   Wrapper class for STL random
*/

#ifndef __URAND_H__
#define __URAND_H__


#include <iostream>
#include <random>
  
class URand 
{
public:
    
    URand( void );
    explicit URand( unsigned int s );
    
    int
    rndInt( int N ) { ++m_count; return m_rd(m_reng) * N; }
    
    int
    rndInt( int min, int max ) { ++m_count; return min + m_rd(m_reng) * (max - min); }
    
	double 
	rndFloat( void ) { ++m_count; return m_rd(m_reng); }
    
    double 
    rndFloat(double min, double max) { ++m_count; return  ((max - min) * m_rd(m_reng)) + min; }

    // set the seed and reinitialise the sequence
    void
    seed( unsigned int s );

    // what seed is being used
    unsigned int seed( void ) const { return m_seed; }

    // at what position in random sequence seeded by s
    unsigned long int count( void ) const { return m_count; }
    
    // return me to the point p in the sequence seeded by s 
    // reset(s,0) == seed(s)
    void
    reset( unsigned int s, unsigned long p );  
    
    std::default_random_engine&
    engine( void ) { return m_reng;  }
	
private:
    
    unsigned long int m_count;
    unsigned int      m_seed;

    std::default_random_engine m_reng;
    std::uniform_real_distribution<double> m_rd;
};

 
std::ostream&
operator<<(std::ostream& ostr, const URand &r);

std::istream&
operator>>(std::istream& istr, URand &r);

#endif // __URAND_H__
