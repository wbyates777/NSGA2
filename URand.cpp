/* Code file for uniform distrubution random number generator 09/02/06

                $$$$$$$$$$$$$$$$$$$$$$$
                $   URand.cpp - code  $ 
                $$$$$$$$$$$$$$$$$$$$$$$

   by W.B. Yates
   Copyright (c) W.B. Yates. All rights reserved.
   History:
 
*/

#ifndef __URAND_H__
#include "URand.h"
#endif


const int DEFAULT_SEED = 13;

std::ostream&
operator<<(std::ostream& ostr, const URand &rn)
{
	ostr << rn.count() << ' ' << rn.seed() << '\n';
	return ostr;
}


std::istream&
operator>>(std::istream& istr, URand &rn)
{
    unsigned int seed = 0;
	unsigned long count = 0;

    istr >> count;
	istr >> seed;
	rn.reset( seed, count );
	
	return istr;
}

URand::URand( void ): m_count(0), m_seed(DEFAULT_SEED), m_reng{DEFAULT_SEED}, m_rd{0,1}
{
}

URand::URand( unsigned int s ): m_count(0), m_seed(s), m_reng{s}, m_rd{0,1}
{
    seed(s);
}

void
URand::seed( unsigned int s )
{   
    m_seed = s;
    m_count = 0;
    m_reng.seed( m_seed );
}


void
URand::reset( unsigned int s, unsigned long p )
{
    seed( s );
    for (int i = 0; i < p; ++i)
        rndFloat();
}


