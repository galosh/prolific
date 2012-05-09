/**
 * \file Random.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      Class definition for the Galosh Random class.  This is our local
 *      interface for random number generation, at present using Boost's random
 *      templates.
 * \par Overview:
 *    This file is part of prolific, a library of useful C++ classes for
 *    working with genomic sequence data and Profile HMMs.  Please see the
 *    document CITING, which should have been included with this file.  You may
 *    use at will, subject to the license (Apache v2.0), but *please cite the
 *    relevant papers* in your documentation and publications associated with
 *    uses of this library.  Thank you!
 * \copyright
 *    Copyright &copyl 2007, 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 * \par License:
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    
 *        http://www.apache.org/licenses/LICENSE-2.0
 *    
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 *****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_RANDOM_HPP__
#define __GALOSH_RANDOM_HPP__

#include "Prolific.hpp"

#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/normal_distribution.hpp>

//#include <boost/random/uniform_int.hpp>
//#include <boost/random/uniform_real.hpp>

namespace galosh {

  class Random
  {
    // This is a typedef for a random number generator.
    // Try boost::mt19937 or boost::ecuyer1988 instead of boost::minstd_rand
    //typedef boost::minstd_rand base_generator_type;
    // See http://www.boost.org/libs/random/random-generators.html
    typedef boost::mt19937 base_generator_type;
    /* http://www.boost.org/libs/random/random-generators.html#mt19937 */
    #define DEFAULT_SEED 5489

  public:
    // Define a random number generator and initialize it with a reproducible
    // seed.
    // (The seed is unsigned, otherwise the wrong overload may be selected
    // when using mt19937 as the base_generator_type.)
    uint32_t m_seed;

    base_generator_type m_generator;
    boost::uniform_01<base_generator_type> m_uniform_01;
    boost::uniform_smallint<uint32_t> m_uniform_uint;
    uint32_t m_uniform_uint_max;
    boost::normal_distribution<double> m_normal;

    Random () :
      m_seed( DEFAULT_SEED ),
      m_generator( DEFAULT_SEED ),
      m_uniform_01( m_generator ),
      m_uniform_uint( ( uint32_t )0, ( uint32_t )1 ),
      m_uniform_uint_max( 1 ),
      m_normal( 0, 1 )
    {
      // Do nothing else
    } // <init>()

    Random ( uint32_t const seed ) :
      m_seed( seed ),
      m_generator( seed ),
      m_uniform_01( m_generator ),
      m_uniform_uint( ( uint32_t )0, ( uint32_t )1 ),
      m_uniform_uint_max( ( uint32_t )1 ),
      m_normal( 0, 1 )
    {
      // Do nothing else
    } // <init>( uint32_t )

    double
    nextUniform ()
    {
      return m_uniform_01();
    } // nextUniform()

    uint32_t
    nextUniform ( uint32_t const max )
    {
      if( max != m_uniform_uint_max ) {
        // TODO: REMOVE
        //cout << "about to set the bounds to 0, " << max << endl;
        m_uniform_uint =
          boost::uniform_smallint<uint32_t>( ( uint32_t )0, max );
        m_uniform_uint_max = max;
        // TODO: REMOVE
        //cout << "just set the bounds to 0, " << max << endl;
      }
      return m_uniform_uint( m_generator );
    } // nextUniform( uint32_t )

    double
    nextNormal ()
    {
      return m_normal( m_uniform_01 );
    } // nextNormal()

    /**
     * This method is from:
@article{1979,
  ISSN = {0035-9254},
  abstract = {Gamma variates with index $\alpha > 1$ are produced by combining two adaptations of Kinderman and Monahan's technique for generating random variates by the use of the ratio of uniform variates. Expensive logarithmic/exponential evaluation is frequently avoided. The method is uniformly fast for all $\alpha > 1$, is compact and easy to program.},
  author = {{Cheng, R. C. H.} and {Feast, G. M.}},
  copyright = {Copyright 1979 Royal Statistical Society},
  journal = {Applied Statistics},
  jstor_articletype = {Full Length Article},
  jstor_date = {1979},
  jstor_formatteddate = {1979},
  keywords = {Gamma Distribution, Random Numbers, Rejection Method, Simulation},
  number = {3},
  pages = {290--295},
  publisher = {Royal Statistical Society},
  title = {Some Simple Gamma Variate Generators},
  url = {http://links.jstor.org/sici?sici=0035-9254%281979%2928%3A3%3C290%3ASSGVG%3E2.0.CO%3B2-I},
  volume = {28},
  year = {1979},
}
     */
    template <class RealType>
    RealType
    nextGamma ( RealType const alpha )
    {
      assert( alpha > 0 );
      if( alpha < 1.0 ) {
        // See the note on Page 9 of the above reference.  The method only
        // works for alpha > 1.0, but you can do this trick to work around the
        // problem.
        RealType tmp_gamma = ( RealType )nextGamma( alpha + 1.0 );
        RealType U = ( RealType )nextUniform();
        return ( tmp_gamma * ( pow( U, ( 1.0 / alpha ) ) ) );
      }

      RealType d = alpha - ( 1.0 / 3.0 );
      RealType c = ( 1.0 / sqrt( 9.0 * d ) );
      RealType x,v,U;
      RealType x_squared;
      while( true ) { // Keep trying until we draw a variate.  (Acceptance rate
                      // is high, above .9 for at least alpha < 100)
        // get x and v
        do {
          x = ( RealType )nextNormal();
          v = 1.0 + ( c * x );
        } while( v <= 0.0 );
        v = ( v * v * v );

        // Use them for acceptance-rejection
        x_squared = ( x * x );
        U = ( RealType )nextUniform(); 
        // Cheaper 'squeeze method'
        if( U < ( 1.0 - ( .0331 * ( x_squared * x_squared ) ) ) ) {
          return ( d * v );
        }
        // If that doesn't work, try version with logs.
        if( log( U ) <
            ( ( 0.5 * x_squared ) + ( d * ( 1.0 - v + log( v ) ) ) ) ) {
          return ( d * v );
        }
      } // Until we get a variate...

    } // nextGamma( RealType )

    void
    setSeed ( uint32_t const seed )
    {
      m_seed = seed;
      m_generator.seed( seed );
      // m_uniform_01 stores its own copy of the generator.
      m_uniform_01.base().seed( seed );
    } // setSeed( uint32_t )

    uint32_t
    getSeed ()
    {
      return m_seed;
    } // getSeed()

  }; // End class Random

} // End namespace galosh

#endif // __GALOSH_RANDOM_HPP__
