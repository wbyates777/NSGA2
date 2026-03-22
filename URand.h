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
    
    URand( void ): m_count(0), m_seed(DEFAULT_SEED), m_reng{DEFAULT_SEED}, m_rd{0,1} {}
    explicit URand( unsigned int s ): m_count(0), m_seed(s), m_reng{s}, m_rd{0,1} { seed(s); }
    ~URand( void )=default;
    
    
    int
    rndInt( int N ) { ++m_count; return m_rd(m_reng) * N; } // return [0,...,N]
    
    int
    rndInt( int min, int max ) { ++m_count; return min + m_rd(m_reng) * (max - min); }  // return [min,...,max]
    
	double 
	rndFloat( void ) { ++m_count; return m_rd(m_reng); } // return (0,...,1)
    
    double 
    rndFloat(double min, double max) { ++m_count; return  ((max - min) * m_rd(m_reng)) + min; } // return (min,...,max)

    
    // set the seed and reinitialise the sequence
    void
    seed( unsigned int s )
    {   
        m_seed = s;
        m_count = 0;
        m_reng.seed( m_seed );
    }

    // what seed is being used
    unsigned int seed( void ) const { return m_seed; }

    // at what position in random sequence seeded by s
    unsigned long int count( void ) const { return m_count; }
    
    // return me to the point p in the sequence seeded by s 
    // reset(s,0) == seed(s)
    void
    reset( unsigned int s, unsigned long p )
    {
        seed( s );
        for (long int i = 0; i < p; ++i)
            rndFloat();
    }
    
    std::default_random_engine&
    engine( void ) { return m_reng;  }
	
private:
    
    unsigned long int m_count;
    unsigned int      m_seed;

    std::default_random_engine m_reng;
    std::uniform_real_distribution<double> m_rd;
    
    static constexpr unsigned int DEFAULT_SEED = 13;
};

 

inline std::ostream&
operator<<(std::ostream& ostr, const URand &r)
{
    ostr << r.count() << ' ' << r.seed() << '\n';
    return ostr;
}


inline std::istream&
operator>>(std::istream& istr, URand &r)
{
    unsigned int seed = 0;
    unsigned long count = 0;
    
    istr >> count >> seed;
    
    r.reset( seed, count );
    
    return istr;
}


#endif // __URAND_H__
