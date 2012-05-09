/**
 * \file MultinomialDistribution.hpp
 * \author  D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      Class definition for the galosh::MultinomialDistribution class.  It
 *      represents a discrete distribution over a finite set of categories.
 * \par Overview:
 *    This file is part of prolific, a library of useful C++ classes for
 *    working with genomic sequence data and Profile HMMs.  Please see the
 *    document CITING, which should have been included with this file.  You may
 *    use at will, subject to the license (Apache v2.0), but *please cite the
 *    relevant papers* in your documentation and publications associated with
 *    uses of this library.  Thank you!
 *
 * \copyright &copy; 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
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

#ifndef __GALOSH_MULTINOMIALDISTRIBUTION_HPP__
#define __GALOSH_MULTINOMIALDISTRIBUTION_HPP__

#include "Prolific.hpp"

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <limits>
using std::numeric_limits;
#include <math.h>
#include </usr/include/stdint.h>
#include </usr/include/assert.h>

#include "Random.hpp"
#include "Ambiguous.hpp"

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/version.hpp>
#include <boost/lexical_cast.hpp>   ///TAH 2/12

#include <seqan/basic.h>
#include "Algebra.hpp"
using seqan::ordValue;

namespace galosh {

  /**
   * A MultinomialDistribution is a map from an enumerated type to a
   * probability type.  There are requirements on ValueType: it should support:
   * ValueSize<ValueType>::VALUE, and should support ValueType( unsigned ).  See
   * Residue.hpp for examples.
   */
  template <typename ValueType,
            typename ProbabilityType>
  class MultinomialDistribution
  {
    // Boost serialization
  private:
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        ar & BOOST_SERIALIZATION_NVP( m_probs );
      } // serialize( Archive &, const unsigned int )

  public:
    // NOTE These should be protected but have to be public because different
    // template instantiations need to directly access each other's values...
    static uint32_t const m_elementCount = seqan::ValueSize<ValueType>::VALUE;
    ProbabilityType m_probs[ m_elementCount ];
    
    std::string 
    toString () const;
 
    /**
     * Inner class for ambiguous values with callbacks for any change.  Note
     * that the AmbiguousValue will not change if you alter the
     * MultinomialDistribution via a different interface; the callbacks only
     * work the other way: if you alter the AmbiguousValue, the
     * MultinomialDistribution will be altered.  Note also that most methods
     * are not implemented.  At this point, only operator= and operator+= are
     * implemented.
     */
    template <typename AmbiguityCodeType>
    class AmbiguousValue : public ProbabilityType
    {
      AmbiguityCodeType const & m_ambiguityCode;
      MultinomialDistribution & m_callbackDistribution;

    public:
      AmbiguousValue (
        AmbiguityCodeType const & ambiguity_code,
        MultinomialDistribution & callback_distribution
      ) :
        ProbabilityType(),
        m_ambiguityCode( ambiguity_code ),
        m_callbackDistribution( callback_distribution )
      {
        // Do nothing else
      } // <init>( AmbiguityCodeType const &, MultinomialDistribution & )

      template <typename T>
      AmbiguousValue (
        AmbiguityCodeType const & ambiguity_code,
        MultinomialDistribution & callback_distribution,
        T const & v
      ) :
        ProbabilityType( v ),
        m_ambiguityCode( ambiguity_code ),
        m_callbackDistribution( callback_distribution )
      {
        // Do nothing else
      } // <init>( AmbiguityCodeType const &, MultinomialDistribution &, T const & )

      AmbiguousValue (
        AmbiguousValue const & copy_from
      ) :
        ProbabilityType( copy_from.toProbabilityType() ),
        m_ambiguityCode( copy_from.m_ambiguityCode ),
        m_callbackDistribution( copy_from.m_callbackDistribution )
      {
        // Do nothing else
      } // <init>( AmbiguousValue const & )

      template <typename T>
      AmbiguousValue &
      operator= ( T const & v )
      {
        return this->ambiguousAssign( v );
      } // operator= ( T const & )

      template <typename T>
      AmbiguousValue &
      unambiguousAssign ( T const & v )
      {
        ProbabilityType::operator=( v );

        return *this;
      } // unambiguousAssign ( T const & )

      template <typename T>
      AmbiguousValue &
      ambiguousAssign ( T const & v )
      {
        const size_t num_ambiguous_elements =
          ambiguousCount( m_ambiguityCode, ValueType() );

        // Actually divide the value evenly among the ambiguous elements
        ProbabilityType prob_to_assign = v;
        prob_to_assign /= num_ambiguous_elements;

        this->unambiguousAssign( v );
        ValueType possible_value;
        for( size_t i = 0; i < num_ambiguous_elements; i++ ) {
          galosh::ambiguousAssign( possible_value, m_ambiguityCode, i );
          m_callbackDistribution.m_probs[ ordValue( possible_value ) ] =
            prob_to_assign;
        }

        return *this;
      } // ambiguousAssign ( T const & )

      template <typename T>
      AmbiguousValue &
      operator+= ( T const & v )
      {
        return ambiguousIncrement( v );
      }

      template <typename T>
      AmbiguousValue &
      unambiguousIncrement ( T const & v )
      {
        ProbabilityType::operator+=( v );

        return *this;
      } // unambiguousIncrement ( T const & )

      template <typename T>
      AmbiguousValue &
      ambiguousIncrement ( T const & v )
      {
        const size_t num_ambiguous_elements =
          ambiguousCount( m_ambiguityCode, ValueType() );

        // Actually divide the value evenly among the ambiguous elements
        ProbabilityType prob_to_add = v;
        prob_to_add /= num_ambiguous_elements;

        this->unambiguousAssign( 0 );
        ValueType possible_value;
        for( size_t i = 0; i < num_ambiguous_elements; i++ ) {
          galosh::ambiguousAssign( possible_value, m_ambiguityCode, i );
          m_callbackDistribution.m_probs[ ordValue( possible_value ) ] += prob_to_add;
          this->unambiguousIncrement( m_callbackDistribution.m_probs[ ordValue( possible_value ) ] );
        }

        return *this;
      } // ambiguousIncrement ( T const & )

      ProbabilityType
      toProbabilityType () const
      {
        return *this;
      } // toProbabilityType () const

      // NOT (YET) IMPLEMENTED:
      template <typename T>
      AmbiguousValue &
      operator-= ( T const & v )
      {
        return this->ambiguousDecrement( v );
      }
      // NOT (YET) IMPLEMENTED:
      template <typename T>
      AmbiguousValue &
      operator*= ( T const & v )
      {
        return this->ambiguousMultiplyAssign( v );
      }
      // NOT (YET) IMPLEMENTED:
      template <typename T>
      AmbiguousValue &
      operator/= ( T const & v )
      {
        return this->ambiguousDivideAssign( v );
      }
      // NOT (YET) IMPLEMENTED:
      template <typename T>
      AmbiguousValue
      operator+ ( T const & v )
      {
        return this->ambiguousAdd( v );
      }
      // NOT (YET) IMPLEMENTED:
      template <typename T>
      AmbiguousValue
      operator- ( T const & v )
      {
        return this->ambiguousSubtract( v );
      }
      // NOT (YET) IMPLEMENTED:
      template <typename T>
      AmbiguousValue
      operator* ( T const & v )
      {
        return this->ambiguousMultiply( v );
      }
      // NOT (YET) IMPLEMENTED:
      template <typename T>
      AmbiguousValue
      operator/ ( T const & v )
      {
        return this->ambiguousDivide( v );
      }

    }; // Inner class MultinomialDistribution::AmbiguousValue

  public:
    // Starts out with an even distribution
    MultinomialDistribution ();

    // Copy constructor
    template <typename AnyProbabilityType>
    MultinomialDistribution ( MultinomialDistribution<ValueType,AnyProbabilityType> const & copy_from );

    // Start out with an even distribution.
    void
    reinitialize ();

    uint32_t
    size () const;

    /**
     * How many free parameters are there?  If we consider this a
     * MultinomialDistribution, its probabilities should sum to 1, so there are
     * ( size() - 1 ) free parameters.
     */
    uint32_t
    freeParameterCount () const;

    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator= ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist );

    MultinomialDistribution<ValueType,ProbabilityType> &
    operator= ( MultinomialDistribution<ValueType,ProbabilityType> const& other_dist );

    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator= ( AnyProbabilityType const & set_to );

    /**
     * Add to each probability the corresponding probability in the given
     * distribution.  Note that this may violate the rule that the
     * probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator+= ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist );

    /**
     * Add to each probability the given value.  Note that this may violate the
     * rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator+= ( AnyProbabilityType const& value );

    /**
     * Subtract from each probability the corresponding probability in the
     * given distribution.  Note that this may violate the rule that the
     * probabilities are greater than 0.  If it does, you may later run into
     * problems (it is not checked here).
     */
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator-= ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist );

    /**
     * Subtract from each probability the given value.  Note that this may
     * violate the rule that the probabilities are greater than 0.  If it does,
     * you may later run into problems (it is not checked here).
     */
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator-= ( AnyProbabilityType const& value );

    /**
     * Multiply each probability by the corresponding probability in the given
     * distribution.
     */
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator*= (
      MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist
    );

    /**
     * Create a new Multinomial distribution that is the product of this and
     * another.
     */
    // Note that this copies its result to the caller...
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType>
    operator* (
      MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist
    );

    /**
     * Multiply each probability by the given value.
     */
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator*= (
      AnyProbabilityType const& scalar
    );

    /**
     * Divide each probability value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    MultinomialDistribution<ValueType,ProbabilityType> &
    operator/= ( AnyProbabilityType const& denominator );

    template <typename AnyProbabilityType>
    bool
    operator== ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist );

    bool
    operator== ( MultinomialDistribution<ValueType,ProbabilityType> const& other_dist );

    template <typename AnyProbabilityType>
    bool
    operator== ( AnyProbabilityType const & compare_to );

    template <typename AnyType>
    bool
    operator!= ( AnyType const& other_thing );

    ProbabilityType &
    operator[] ( const ValueType & which );

    ProbabilityType const&
    operator[] ( const ValueType & which ) const;

    ProbabilityType &
    operator[] ( const uint32_t which );

    ProbabilityType const&
    operator[] ( const uint32_t which ) const;

    template <typename AnyValueType>
    AmbiguousValue<AnyValueType> const
    operator[] ( const AnyValueType & which ) const
    {
      // This should not compile if attempted with something that is not ambiguous...
      return ambiguousSum( which, typename IsAmbiguous<AnyValueType, ValueType>::Type() );
    } // operator[] ( AnyValueType const & ) const

    template <typename AnyValueType>
    AmbiguousValue<AnyValueType>
    operator[] ( const AnyValueType & which ) 
    {
      // This should not compile if attempted with something that is not ambiguous...
      return ambiguousSum( which, typename IsAmbiguous<AnyValueType, ValueType>::Type() );
    } // operator[] ( AnyValueType const & )

    // Not const, because the return AmbiguousValue can be altered, altering
    // this MultinomialDistribution.
    template <typename AnyValueType>
    AmbiguousValue<AnyValueType>
    ambiguousSum ( const AnyValueType & ambiguity_code, seqan::True const tag = seqan::True() )
    {
      const size_t num_ambiguous_elements =
        ambiguousCount( ambiguity_code, ValueType() );
      AmbiguousValue<AnyValueType> prob( ambiguity_code, *this, 0 );
      ValueType v;
      for( size_t i = 0; i < num_ambiguous_elements; i++ ) {
        galosh::ambiguousAssign( v, ambiguity_code, i );
        prob.unambiguousIncrement( m_probs[ seqan::ordValue( v ) ] );
      }

      return prob;
    } // ambiguousSum ( AnyValueType const & [, seqan::True ] )

    // Not const, because the return AmbiguousValue can be altered, altering
    // this MultinomialDistribution.
    template <typename AnyValueType>
    AmbiguousValue<AnyValueType>
    ambiguousSum ( const AnyValueType & ambiguity_code, seqan::True const tag = seqan::True() ) const
    {
      const size_t num_ambiguous_elements =
        ambiguousCount( ambiguity_code, ValueType() );
      AmbiguousValue<AnyValueType> prob( ambiguity_code, ( *const_cast<MultinomialDistribution *>( this ) ), 0 );
      ValueType v;
      for( size_t i = 0; i < num_ambiguous_elements; i++ ) {
        galosh::ambiguousAssign( v, ambiguity_code, i );
        prob.unambiguousIncrement( m_probs[ ordValue( v ) ] );
      }

      return prob;
    } // ambiguousSum ( AnyValueType const & [, seqan::True ] ) const

    // Evenly distributes value to be added, among all of the values to which
    // the ambiguity code refers, so that the total amount added is the value
    // given.
    template <typename AnyValueType, typename AnyProbabilityType>
    void
    ambiguousIncrement ( const AnyValueType & ambiguity_code, AnyProbabilityType const & value )
    {
      const size_t num_ambiguous_elements =
        ambiguousCount( ambiguity_code, ValueType() );
      // Actually divide the value evenly among the ambiguous elements
      ProbabilityType prob_to_add = value;
      prob_to_add /= num_ambiguous_elements;

      ValueType v;
      for( size_t i = 0; i < num_ambiguous_elements; i++ ) {
        galosh::ambiguousAssign( v, ambiguity_code, i );
        m_probs[ ordValue( v ) ] += prob_to_add;
      }
      return;
    } // ambiguousIncrement ( AnyValueType const &, AnyProbabilityType const & )

    /**
     * Stream reader.
     */
    friend std::istream &
    operator>> (
      std::istream & is,
      MultinomialDistribution<ValueType,ProbabilityType> & md
    )
    {
      md.readMultinomialDistribution( is );

      return is;
    } // friend operator>> ( istream &, MultinomialDistribution<ValueType> & )

    void
    readMultinomialDistribution (
      std::istream& is
    );

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      MultinomialDistribution<ValueType,ProbabilityType> const& md
    )
    {
      md.writeMultinomialDistribution( os );

      return os;
    } // friend operator<< ( basic_ostream &, MultinomialDistribution<ValueType> const& )

    template<class CharT, class Traits>
    void
    writeMultinomialDistribution (
      std::basic_ostream<CharT,Traits>& os
    ) const;

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ();

    /**
     * Set all values to be the same (evenly distributed).
     */
    void
    even ();

    /**
     * Set all values to a random quantity, uniformly distributed over the
     * simplex.
     */
    void
    uniform ( Random & random );

    /**
     * Adjust values such that they sum to one.  This is just like saying *this
     * /= this->total().
     */
    void
    normalize ();

    /**
     * Adjust values such that they sum to one, enforcing the given minimum
     * value.
     */
    template <typename other_type>
    void
    normalize ( other_type const & min );

    /**
     * Adjust values such that they sum to one, enforcing the given minimum
     * value.
     */
    void
    normalize ( ProbabilityType const & min );

    /**
     * We do not enforce that values sum to 1.  This method calculates and
     * returns the sum.
     */
    ProbabilityType
    total () const;

    /**
     * Return the lowest of the probabilities contained herein.
     */
    ProbabilityType
    minimumValue () const;

    /**
     * Return the largest of the probabilities contained herein.
     */
    ProbabilityType
    maximumValue () const;

    /**
     * Return the ValueType with the largest of the probabilities contained
     * herein.
     */
    ValueType
    maximumValueType () const;

    /**
     * Calculate and return the Euclidean distance between this distribution and
     * another distribution (over the same values, with the same probability
     * type).
     *
     * NOTE: this uses only the first (m_elementCount - 1) values, under the
     * assumption that the final value is just 1-(the sum of the rest).  This
     * is arbitrary!  Choosing to use a different subset will result in a
     * different euclideanDistance calculation.  If you know that the
     * distributions aren't normalized, you should either normalize them first,
     * or (if you don't want them to be normalized), add the squared distance
     * between the last two components yourself.
     */
    double
    euclideanDistance (
      MultinomialDistribution<ValueType, ProbabilityType> const& another_dist
    ) const;
  
    /**
     * Calculate and return the square of the Euclidean distance between this
     * distribution and another distribution (over the same values, with the
     * same probability type).
     *
     * NOTE: this uses only the first (m_elementCount - 1) values, under the
     * assumption that the final value is just 1-(the sum of the rest).  This
     * is arbitrary!  Choosing to use a different subset will result in a
     * different euclideanDistance calculation.  If you know that the
     * distributions aren't normalized, you should either normalize them first,
     * or (if you don't want them to be normalized), add the squared distance
     * between the last two components yourself.
     */
    double
    euclideanDistanceSquared (
      MultinomialDistribution<ValueType, ProbabilityType> const& another_dist
    ) const;

    /**
     * Calculate a vector of expected distances, in which expected_distances[ v
     * ] is the expected distance from v to the other values (averaged over
     * this distribution of other values).
     */
    template <typename DistanceMatrixType,
              typename ExpectedDistanceType>
    void
    calculateExpectedDistances (
      MultinomialDistribution<ValueType, MultinomialDistribution<ValueType, DistanceMatrixType> > const & distance_matrix,
      MultinomialDistribution<ValueType, ExpectedDistanceType> & expected_distances
    ) const
    {
      // Note that expected_distances[ i ] will be the expected distance *from*
      // i, averaged over distances *to* another value.
      static MultinomialDistribution<ValueType, ExpectedDistanceType> tmp_vec =
        MultinomialDistribution<ValueType, ExpectedDistanceType>();
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        tmp_vec = *this;
        tmp_vec *=
          distance_matrix[ i ];
        // TODO: REMOVE
        //cout << "calculateExpectedDistances: tmp_vec is " << tmp_vec << endl;
        expected_distances[ i ] = tmp_vec.total();
      }
    } // calculateExpectedDistances( MultinomialDistribution<ValueType, MultinomialDistribution<ValueType, DistanceMatrixType> > const &, MultinomialDistribution<ValueType, ExpectedDistanceType> & ) const

    /**
     * Calculate a vector of expected distances, in which expected_distances[ v
     * ] is the expected distance from v to the other values (averaged over
     * this distribution of other values).
     */
    // overloaded for ExpectedValueType==double as a workaround, since there's
    // no direct casting from Probabilities to doubles.
    template <typename DistanceMatrixType>
    void
    calculateExpectedDistances (
      MultinomialDistribution<ValueType, MultinomialDistribution<ValueType, DistanceMatrixType> > const & distance_matrix,
      MultinomialDistribution<ValueType, double> & expected_distances
    ) const
    {
      // Note that expected_distances[ i ] will be the expected distance *from*
      // i, averaged over distances *to* another value.
      double dist;
      uint32_t j;
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        expected_distances[ i ] = 0.0;
        for( j = 0; j < m_elementCount; j++ ) {
          dist = toDouble( ( *this )[ j ] );
          dist *= toDouble( distance_matrix[ i ][ j ] );
          expected_distances[ i ] += dist;
        }
      }
    } // calculateExpectedDistances( MultinomialDistribution<ValueType, MultinomialDistribution<ValueType, DistanceMatrixType> > const &, MultinomialDistribution<ValueType, double> & ) const

    /**
     * Calculate the expected value, given a vector of values.  This is just
     * the average of the given values, where the average is taken over this
     * distribution.
     */
    template <typename ExpectedValueType>
    ExpectedValueType
    calculateExpectedValue (
      MultinomialDistribution<ValueType, ExpectedValueType> const & values
    ) const
    {
      static MultinomialDistribution<ValueType, ExpectedValueType> tmp_vec =
        MultinomialDistribution<ValueType, ExpectedValueType>();
      tmp_vec = values;
      tmp_vec *= *this;
      // TODO: REMOVE
      //cout << "calculateExpectedValue: tmp_vec is " << tmp_vec << endl;
      return tmp_vec.total();
    } // calculateExpectedValue( MultinomialDistribution<ValueType, ExpectedValuyeType> const & ) const

    /**
     * Calculate the expected value, given a vector of values.  This is just
     * the average of the given values, where the average is taken over this
     * distribution.
     */
    // overloaded for ExpectedValueType==double as a workaround, since there's
    // no direct casting from Probabilities to doubles.
    double
    calculateExpectedValue (
      MultinomialDistribution<ValueType, double> const & values
    ) const
    {
      double total = 0.0;
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        total += ( values[ i ] * toDouble( ( *this )[ i ] ) );
      }
      return total;
    } // calculateExpectedValue( MultinomialDistribution<ValueType, double> const & ) const

    /**
     * Return a value drawn randomly according to this distribution.
     *
     * NOTE: You must ensure that the distribution is normalized first.  For
     * efficiency the test is not performed in this method.  Use (total() == 1),
     * or just call normalize().
     */
    ValueType
    draw ( Random & random ) const;

    /**
     * Change the probabilities of this Multinomial to values drawn from a
     * dirichlet distribution with the given counts.  The counts vector must,
     * of course, be of length m_elementCount.  VectorCountType can be a vector
     * of a Real type (anything coercible to a double using toDouble( count )).
     */
    template <class VectorCountType>
    void
    fromDirichlet ( VectorCountType const & counts, Random & random );

    // Note that this copies its result to the caller.  See the other
    // toBoltzmannGibbs for an alternative.
    template <typename RealType>
    MultinomialDistribution<ValueType, RealType>
    toBoltzmannGibbs ( double const temperature );

    template <typename RealType>
    void
    toBoltzmannGibbs (
      RealType const temperature,
      MultinomialDistribution<ValueType, RealType> & boltzmann_gibbs_to_be_filled
    ) const;

    template <typename RealType>
    void
    fromBoltzmannGibbs (
      RealType const temperature,
      MultinomialDistribution<ValueType, RealType> from_dist
    );


  public:
    class DistanceMatrix :
      public MultinomialDistribution<
               ValueType,
               MultinomialDistribution<ValueType, ProbabilityType> >
    {
    public:
      bool isSymmetric;
      bool isProximityMatrix; // versus a distance matrix
      bool isMultiplicative; // versus additive

      DistanceMatrix () :
        isSymmetric( false ),
        isProximityMatrix( false ), // no, it's a *distance* matrix
        isMultiplicative( false ) // no, it's *additive*
      {
        // Do nothing else
      } // <init>()

      // Why do I have to explicitly delegate this?  I dunno.
      template <typename AnyType>
      MultinomialDistribution<ValueType, ProbabilityType> &
      operator= ( AnyType const& set_to )
      {
        MultinomialDistribution<ValueType, MultinomialDistribution<ValueType, ProbabilityType> >::operator=( set_to );
      } // operator=( AnyType const & )

    }; // End inner class MultinomialDistribution::DistanceMatrix

  }; // End class MultinomialDistribution

  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::MultinomialDistribution ////
  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_INIT
  MultinomialDistribution<ValueType, ProbabilityType>::
    // Starts out with an even distribution
  MultinomialDistribution ()
  {
    // Start out with an even distribution.
    even();
  } // <init>()

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_INIT
  MultinomialDistribution<ValueType, ProbabilityType>::
    // Copy constructor
    MultinomialDistribution ( MultinomialDistribution<ValueType,AnyProbabilityType> const & copy_from )
    {
      // Call operator=
      *this = copy_from;
    } // <init>( MultinomialDistribution<AnyProbabilityType> const & )

  /**
   * \fn std::string MultinomialDistribution::toString() const
   * \brief create a string representation of a MultinomialDistribution object.
   *
   * This has been split off from Paul's original writeMultinomialDistribution, so
   * so that the string for applications other than direct output (e.g. output to
   * somewhere other than cout.
   *
   * TAH 1/12
   */

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_INIT
  std::string
  MultinomialDistribution<ValueType,ProbabilityType>:: 
  toString() const {
    	   if(m_elementCount == 0) return "()";
   	   std::string retVal = "(";
   	   for( uint32_t i = 0; i < m_elementCount; i++ )
   	   {
   	      if( i != 0 ) {
   	         retVal += ",";
   	      }
   	      retVal += ValueType( i );
   	      retVal += "=";
              /// Paul's original note from writeMultnomialDistribution
                /// The cast to double forces the printed output to be how a double
   	        /// prints, not how a ProbabilityType prints.  This is safe because we
   	        /// generally keep values above some minimal value.  See normalize(
   	        /// ProbabilityType ).
   	        /// os << ( double )m_probs[ i ];
                /// 
                /// \todo Here we're using lexical_cast instead of ( double ) - this is after
                /// hours of suffering with "createRandomSequence" which uses "floatrealspace"
                /// for its ProbabilityType and somehow failed to compile.  It would be good
                /// to understand the compilation failure.  
              std::ostringstream oss;
              oss << m_probs[ i ];
       	      retVal += oss.str();

   	   }
   	   retVal += ")";
   	   return retVal;
  }


  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
  reinitialize ()
  {
    // Start out with an even distribution.
    even();
  } // reinitialize()

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
    uint32_t
  MultinomialDistribution<ValueType, ProbabilityType>::
    size () const
    {
      return m_elementCount;
    } // size() const

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
    uint32_t
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * How many free parameters are there?  If we consider this a
     * MultinomialDistribution, its probabilities should sum to 1, so there are
     * ( size() - 1 ) free parameters.
     */
    freeParameterCount () const
    {
      return m_elementCount - 1;
    } // freeParameterCount() const

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  MultinomialDistribution<ValueType, ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator= ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] = other_dist.m_probs[ i ];
      }
      return *this;
    } // operator=( MultinomialDistribution<ValueType,AnyProbabilityType>& )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_COPY
  MultinomialDistribution<ValueType, ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator= ( MultinomialDistribution<ValueType,ProbabilityType> const& other_dist )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] = other_dist.m_probs[ i ];
      }
      return *this;
    } // operator=( MultinomialDistribution<ValueType,ProbabilityType>& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_REINITIALIZE
  MultinomialDistribution<ValueType, ProbabilityType> &
  MultinomialDistribution<ValueType,ProbabilityType>::
    operator= ( AnyProbabilityType const & set_to )
    {
      if( set_to == 0 ) {
        zero();
      } else{
        for( uint32_t i = 0; i < m_elementCount; i++ ) {
          m_probs[ i ] = set_to;
        }
      }
      return *this;
    } // operator=( AnyProbabilityType const & )

  template <typename ValueType,
            typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
  MultinomialDistribution<ValueType,ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Add to each probability the corresponding probability in the given
     * distribution.  Note that this may violate the rule that the
     * probabilities sum to 1.
     */
    operator+= ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] += other_dist.m_probs[ i ];
      }
      return *this;
    } // operator+=( MultinomialDistribution<ValueType,AnyProbabilityType> const& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    MultinomialDistribution<ValueType,ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Add to each probability the given value.  Note that this may violate the
     * rule that the probabilities sum to 1.
     */
    operator+= ( AnyProbabilityType const& value )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] += value;
      }
      return *this;
    } // operator+=( AnyProbabilityType const& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    MultinomialDistribution<ValueType,ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Subtract from each probability the corresponding probability in the
     * given distribution.  Note that this may violate the rule that the
     * probabilities are greater than 0.  If it does, you may later run into
     * problems (it is not checked here).
     */
    operator-= ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] -= other_dist.m_probs[ i ];
      }
      return *this;
    } // operator-=( MultinomialDistribution<ValueType,AnyProbabilityType> const& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    MultinomialDistribution<ValueType,ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Subtract from each probability the given value.  Note that this may
     * violate the rule that the probabilities are greater than 0.  If it does,
     * you may later run into problems (it is not checked here).
     */
    operator-= ( AnyProbabilityType const& value )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] -= value;
      }
      return *this;
    } // operator-=( AnyProbabilityType const& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    MultinomialDistribution<ValueType,ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Multiply each probability by the corresponding probability in the given
     * distribution.
     */
    operator*= (
      MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist
    )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] *= other_dist.m_probs[ i ];
      }
      return *this;
    } // operator*=( MultinomialDistribution const& )


  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    MultinomialDistribution<ValueType,ProbabilityType>
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Create a new Multinomial distribution that is the product of this and
     * another.
     */
    // Note that this copies its result to the caller...
    operator* (
      MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist
    )
    {
      MultinomialDistribution<ValueType,ProbabilityType> md = *this;
      md *= other_dist;
      return md;
    } // operator*( MultinomialDistribution const& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    MultinomialDistribution<ValueType,ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Multiply each probability by the given value.
     */
    operator*= (
      AnyProbabilityType const& scalar
    )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] *= scalar;
      }
      return *this;
    } // operator*=( AnyProbabilityType const& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    MultinomialDistribution<ValueType,ProbabilityType> &
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Divide each probability value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    operator/= ( AnyProbabilityType const& denominator )
    {
      // TODO: REMOVE
      //assert( m_elementCount > 0 );
      //cout << "MultinomialDistribution::operator/=( " << denominator << " )" << endl;
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        //cout << "MultinomialDistribution::operator/=( " << denominator << " ): m_probs[ " << i << " ] is (BEFORE) " << m_probs[ i ] << endl;
        m_probs[ i ] /= denominator;
        //cout << "MultinomialDistribution::operator/=( " << denominator << " ): m_probs[ " << i << " ] is (AFTER) " << m_probs[ i ] << endl;
      }
      return *this;
    } // operator/=( AnyProbabilityType const& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPARE
  bool
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator== ( MultinomialDistribution<ValueType,AnyProbabilityType> const& other_dist )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        if( m_probs[ i ] != other_dist.m_probs[ i ] ) {
          // TODO: REMOVE!
          //cout << "MultinomialDistribution::operator==[1]: diff is at pos i: ours is " << m_probs[ i ] << " != theirs is " << other_dist.m_probs[ i ] << "; diff is " << ( m_probs[ i ] - other_dist.m_probs[ i ] ) << endl;
          return false;
        }
      }
      return true;
    } // operator==( MultinomialDistribution<ValueType,AnyProbabilityType>& )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPARE
  bool
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator== ( MultinomialDistribution<ValueType,ProbabilityType> const& other_dist )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        if( m_probs[ i ] != other_dist.m_probs[ i ] ) {
          // TODO: REMOVE!
          //cout << "MultinomialDistribution::operator==[2]: diff is at pos i: ours is " << m_probs[ i ] << " != theirs is " << other_dist.m_probs[ i ] << "; diff is " << ( m_probs[ i ] - other_dist.m_probs[ i ] ) << endl;
          return false;
        }
      }
      return true;
    } // operator==( MultinomialDistribution<ValueType,ProbabilityType>& )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPARE
  bool
  MultinomialDistribution<ValueType,ProbabilityType>::
    operator== ( AnyProbabilityType const & compare_to )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        if( m_probs[ i ] != compare_to ) {
          // TODO: REMOVE!
          //cout << "MultinomialDistribution::operator==[3]: diff is at pos i: ours is " << m_probs[ i ] << " != " << compare_to << endl; 
          return false;
        }
      }
      return true;
    } // operator==( AnyProbabilityType const & )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename AnyType>
  GALOSH_INLINE_TRIVIAL
  bool
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator!= ( AnyType const& other_thing )
    {
      return !( ( *this ) == other_thing );
    } // operator!=( AnyType const & )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
    ProbabilityType &
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator[] ( const ValueType & which )
    {
      return m_probs[ ( uint32_t )which ];
    } // operator[]( ValueType const & )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
    ProbabilityType const&
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator[] ( const ValueType & which ) const
    {
      return m_probs[ ( uint32_t )which ];
    } // operator[]( ValueType const & ) const

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
    ProbabilityType &
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator[] ( const uint32_t which )
    {
      return m_probs[ which ];
    } // operator[]( uint32_t )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
    ProbabilityType const&
  MultinomialDistribution<ValueType, ProbabilityType>::
    operator[] ( const uint32_t which ) const
    {
      return m_probs[ which ];
    } // operator[]( uint32_t ) const

/**
 * Stream reader.
 */
  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_OSTREAM
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    readMultinomialDistribution (
      std::istream & is
    )
    {
#ifdef NDEBUG
      static const bool do_extra_debugging = false;
#else
      static const bool do_extra_debugging = true;
#endif // NDEBUG .. else ..
      //is >> "(";
      if( do_extra_debugging ) {
        assert( is.get() == '(' );
      } else {
        is.ignore( 1 );
      }
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        if( i != 0 ) {
          //is >> ",";
          if( do_extra_debugging ) {
            assert( is.get() == ',' );
          } else {
            is.ignore( 1 );
          }
        }
        //is >> ValueType( i );
        if( do_extra_debugging ) {
          assert( ( char )is.get() == ( char )ValueType( i ) );
        } else {
          is.ignore( 1 );
        }
        //is >> "=";
        if( do_extra_debugging ) {
          assert( is.get() == '=' );
        } else {
          is.ignore( 1 );
        }
        is >> m_probs[ i ];
        // TODO: REMOVE
        //cout << "reading " << ValueType( i ) << ": " << m_probs[ i ] << endl;
      }
      //is >> ")";
      if( do_extra_debugging ) {
        // TODO: REMOVE
        //cout << "Expecting ')' but getting " << (char)is.peek() << endl;
        assert( is.get() == ')' );
      } else {
        is.ignore( 1 );
      }
    } // readMultinomialDistribution( istream & )

  template <typename ValueType,
            typename ProbabilityType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    writeMultinomialDistribution (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      os << "(";
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        if( i != 0 ) {
          os << ",";
        }
        os << ValueType( i );
        os << "=";
        // The cast to double forces the printed output to be how a double
        // prints, not how a ProbabilityType prints.  This is safe because we
        // generally keep values above some minimal value.  See normalize(
        // ProbabilityType ).
        //os << ( double )m_probs[ i ];
        os << m_probs[ i ];
      }
      os << ")";
    } // writeMultinomialDistribution( basic_ostream & ) const

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
  void
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    zero ()
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] = 0;
      }
    } // zero()

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_INIT
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Set all values to be the same (evenly distributed).
     */
    even ()
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] = ( 1.0 / m_elementCount );
      }
    } // even()

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Set all values to a random quantity, uniformly distributed over the
     * simplex.
     */
    uniform ( Random & random )
    {
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] = random.nextUniform();
      }
      // TODO: REMOVE
      //cout << "uniform() before normalization: " << *this << endl;
      normalize();
      // TODO: REMOVE
      //cout << "uniform() after normalization: " << *this << endl;
    } // uniform( Random & )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Adjust values such that they sum to one.  This is just like saying *this
     * /= this->total().
     */
    normalize ()
    {
      ProbabilityType t = total();
      if( t == 0 ) {
        even();
      } else {
        operator/=( t );
      }
    } // normalize()

  template <typename ValueType,
            typename ProbabilityType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
  void
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Adjust values such that they sum to one, enforcing the given minimum
     * value.
     */
    normalize ( other_type const & min )
    {
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_NORMALIZE_ENSUREMIN
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Adjust values such that they sum to one, enforcing the given minimum
     * value.
     */
    normalize ( ProbabilityType const & min )
    {
      ProbabilityType t = total();
      if( t == 0 ) {
        even();
        return;
      }
      ProbabilityType lowest = minimumValue();
      if( lowest >= ( min * t ) ) {
        operator/=( t );
        return;
      }
      // lowest < (t*min).

      // Add to all values whatever it takes to get lowest >= min.

      // If n is m_elementCount, then we
      // want to_add s.t. min = ( lowest + to_add ) / ( t + n*to_add )
      // t*min + n*to_add*min = lowest + to_add
      // t*min - lowest = to_add - n*to_add*min
      // t*min - lowest = to_add * ( 1 - n*min )
      // (t*min - lowest)/(1-n*min) = to_add

      ProbabilityType vt = ( double )m_elementCount;
      ProbabilityType to_add = ( ( ( t * min ) - lowest ) / ( 1.0 - ( vt * min ) ) );
      
      operator+=( to_add );
      operator/=( t + ( vt * to_add ) );
    } // normalize( ProbabilityType const & )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPLEX_ACCESSOR
    ProbabilityType
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * We do not enforce that values sum to 1.  This method calculates and
     * returns the sum.
     */
    total () const
    {
      ProbabilityType t = m_probs[ 0 ];
      for( uint32_t i = 1; i < m_elementCount; i++ ) {
        t += m_probs[ i ];
      }
      return t;
    } // total()

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPLEX_ACCESSOR
    ProbabilityType
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Return the lowest of the probabilities contained herein.
     */
    minimumValue () const
    {
      ProbabilityType mv = m_probs[ 0 ];
      for( uint32_t i = 1; i < m_elementCount; i++ ) {
        if( mv > m_probs[ i ] ) {
          mv = m_probs[ i ];
        }
      }
      return mv;
    } // minimumValue()

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPLEX_ACCESSOR
    ProbabilityType
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Return the largest of the probabilities contained herein.
     */
    maximumValue () const
    {
      ProbabilityType mv = m_probs[ 0 ];
      for( uint32_t i = 1; i < m_elementCount; i++ ) {
        if( mv < m_probs[ i ] ) {
          mv = m_probs[ i ];
        }
      }
      return mv;
    } // maximumValue()

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPLEX_ACCESSOR
    ValueType
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Return the ValueType with the largest of the probabilities contained
     * herein.
     */
    maximumValueType () const
    {
      // TODO: REMOVE
      //cout << "maximumValueType() called on: " << *this << endl;

      uint32_t mv_which = 0;
      ProbabilityType mv = m_probs[ 0 ];
      for( uint32_t i = 1; i < m_elementCount; i++ ) {
        if( mv < m_probs[ i ] ) {
          mv = m_probs[ i ];
          mv_which = i;
        }
      }
      // TODO: REMOVE
      //cout << "maximum value: " << mv << "(found at " << ValueType( mv_which ) << " == " << mv_which << ")" << endl;

      return ValueType( mv_which );
    } // maximumValueType()

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Calculate and return the Euclidean distance between this distribution and
     * another distribution (over the same values, with the same probability
     * type).
     *
     * NOTE: this uses only the first (m_elementCount - 1) values, under the
     * assumption that the final value is just 1-(the sum of the rest).  This
     * is arbitrary!  Choosing to use a different subset will result in a
     * different euclideanDistance calculation.  If you know that the
     * distributions aren't normalized, you should either normalize them first,
     * or (if you don't want them to be normalized), add the squared distance
     * between the last two components yourself.
     */
    euclideanDistance (
      MultinomialDistribution<ValueType, ProbabilityType> const& another_dist
    ) const
    {
      return sqrt( euclideanDistanceSquared( another_dist ) );
    } // euclideanDistance( MultinomialDistribution<ValueType, ProbabilityType> const& )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_EUCLIDEANDISTSQ
    double
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Calculate and return the square of the Euclidean distance between this
     * distribution and another distribution (over the same values, with the
     * same probability type).
     *
     * NOTE: this uses only the first (m_elementCount - 1) values, under the
     * assumption that the final value is just 1-(the sum of the rest).  This
     * is arbitrary!  Choosing to use a different subset will result in a
     * different euclideanDistance calculation.  If you know that the
     * distributions aren't normalized, you should either normalize them first,
     * or (if you don't want them to be normalized), add the squared distance
     * between the last two components yourself.
     */
    euclideanDistanceSquared (
      MultinomialDistribution<ValueType, ProbabilityType> const& another_dist
    ) const
    {
      double squared_euclidean_distance = 0.0;
      register double tmp_dist;
      for( uint32_t i = 0; i < ( m_elementCount - 1 ); i++ ) {
        if( m_probs[ i ] > another_dist.m_probs[ i ] ) {
          tmp_dist = toDouble( m_probs[ i ] - another_dist.m_probs[ i ] );
        } else {
          tmp_dist = toDouble( another_dist.m_probs[ i ] - m_probs[ i ] );
        }
        squared_euclidean_distance +=
          ( tmp_dist * tmp_dist );
      }
      // TODO: REMOVE!
      //      cout << "Euclidean distance squared between " << *this << " and " << another_dist << " is " << squared_euclidean_distance << endl;

      return squared_euclidean_distance;
    } // euclideanDistanceSquared( MultinomialDistribution<ValueType, ProbabilityType> const& )

  template <typename ValueType,
            typename ProbabilityType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_DRAW
    ValueType
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Return a value drawn randomly according to this distribution.
     *
     * NOTE: You must ensure that the distribution is normalized first.  For
     * efficiency the test is not performed in this method.  Use (total() == 1),
     * or just call normalize().
     */
    draw ( Random & random ) const
    {
      ProbabilityType u = random.nextUniform();
      ProbabilityType next_bin_boundary = 0;
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        next_bin_boundary += m_probs[ i ];
        // TODO: REMOVE
        //        cout << "[MultinomialDistribution::draw(Random&)] u is " << u << ", next_bin_boundary is " << next_bin_boundary << ", i is " << i << ", ValueType( i ) is " << ValueType( i ) << endl;
        if( u < next_bin_boundary ) {
          //cout << "[MultinomialDistribution::draw(Random&)] .. so we are returning " << ValueType( i ) << endl;
          return ValueType( i );
        }
      }
      // We should never reach this point, since the probs should sum to 1.
      assert( false );
      return ValueType( 0 );
    } // draw( Random & ) const

  template <typename ValueType,
            typename ProbabilityType>
    template <class VectorCountType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_FROMDIRICHLET
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    /**
     * Change the probabilities of this Multinomial to values drawn from a
     * dirichlet distribution with the given counts.  The counts vector must,
     * of course, be of length m_elementCount.  VectorCountType can be a vector
     * of a Real type (anything coercible to a double using toDouble( count )).
     */
    fromDirichlet ( VectorCountType const & counts, Random & random )
    {
      assert( counts.size() == m_elementCount );

      double gammas[ m_elementCount ];
      double sum_of_gammas = 0.0;
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        gammas[ i ] = random.nextGamma( toDouble( counts[ i ] ) );
        sum_of_gammas += gammas[ i ];
      }
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        if( sum_of_gammas == 0.0 ) {
          m_probs[ i ] = 0.0;
        } else {
          m_probs[ i ] =
            ( gammas[ i ] / sum_of_gammas );
        }
      }
    } // fromDirichlet( VectorCountType const &, Random & )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename RealType>
  GALOSH_INLINE_TRIVIAL
    MultinomialDistribution<ValueType, RealType>
  MultinomialDistribution<ValueType, ProbabilityType>::
    // Note that this copies its result to the caller.  See the other
    // toBoltzmannGibbs for an alternative.
    toBoltzmannGibbs ( double const temperature )
    {
      MultinomialDistribution<ValueType, RealType> rv =
        MultinomialDistribution<ValueType, RealType>();
      toBoltzmannGibbs( temperature, rv );
      return rv;
    } // toBoltzmannGibbs( double )

  template <typename ValueType,
            typename ProbabilityType>
    template <typename RealType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_TO_BOLTZMANNGIBBS
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    toBoltzmannGibbs (
      RealType const temperature,
      MultinomialDistribution<ValueType, RealType> & boltzmann_gibbs_to_be_filled
    ) const
    {
      register double d;
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        d = toDouble( m_probs[ i ] );
        if( d == 0 ) {
          d = numeric_limits<double>::min();
        }
        boltzmann_gibbs_to_be_filled.m_probs[ i ] =
          static_cast<RealType>( log( d ) );
      }
      boltzmann_gibbs_to_be_filled /= temperature;
    } // toBoltzmannGibbs( RealType, MultinomialDistribution<ValueType, RealType>& ) const

  template <typename ValueType,
            typename ProbabilityType>
    template <typename RealType>
  GALOSH_INLINE_MULTINOMIALDISTRIBUTION_FROM_BOLTZMANNGIBBS
    void
  MultinomialDistribution<ValueType, ProbabilityType>::
    fromBoltzmannGibbs (
      RealType const temperature,
      MultinomialDistribution<ValueType, RealType> from_dist
    )
    {
      from_dist *= temperature;
      for( uint32_t i = 0; i < m_elementCount; i++ ) {
        m_probs[ i ] =
          exp( toDouble( from_dist.m_probs[ i ] ) );
        if( isinf( m_probs[ i ] ) ) {
          m_probs[ i ] = 1.0;
        }
      }
      normalize();
    } // fromBoltzmannGibbs( RealType, MultinomialDistribution<ValueType,RealType> )

} // End namespace galosh

#endif // __GALOSH_MULTINOMIALDISTRIBUTION_HPP__
