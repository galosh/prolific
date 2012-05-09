/**
 * \file DynamicProgramming.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org with some additions
 * by Ted Holzman
 * \par Library:
 * galosh::prolific
 * \brief Class definition(s) for the Galosh DynamicProgramming class.  This is a
 *  super (duper) class with many subclasses and static member functions
 *  for doing useful things.
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

#ifndef __GALOSH_DYNAMICPROGRAMMING_HPP__
#define __GALOSH_DYNAMICPROGRAMMING_HPP__

 // \todo REMOVE
//#include <cassert>

#include "Prolific.hpp"

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "Profile.hpp"

#include "ProfileTree.hpp"
using galosh::ProfileTree;

#include "Sequence.hpp"
using galosh::Sequence;

#include "Fasta.hpp"
using galosh::Fasta;

//#include "BooleanState.hpp"
//using galosh::BooleanState;

#include "MultinomialDistribution.hpp"
using galosh::MultinomialDistribution;

#include "Random.hpp"
using galosh::Random;

/*
#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include <set>
using std::set;
#include <map>
using std::map;
#include <exception>
using std::exception;
*/
#include <vector>
using std::vector;
#include <list>
using std::list;
#include <iterator>
using std::ostream_iterator;
#include <algorithm>
using std::copy;

#include <string>
#include <iostream>
#include <sstream>
//#include <limits> // for numeric_limits<double>
                  /// TAH 3/12 for converting between FILE* and istream
#include <map>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/stream.hpp>
#include <cstdio>

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

/**
 * For scanning and parsing alignment profiles.  \see AlignmentProfileAccessor
 */
#include "ProlificIOHelpers.hpp"
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/xpressive/regex_actions.hpp>

// For conversion between muscle's MSA class and our MultipleAlignment class
#ifdef __HAVE_MUSCLE
#include "muscle/muscle.h"
#include "muscle/msa.h"
#endif //__HAVE_MUSCLE

// TODO: REMOVE.  TESTING.
//#include "Algebra.hpp"

namespace galosh {

  struct dynamicprogramming_subcell_tag
  {
    // Identifying tag for subcell types
  };
  struct dynamicprogramming_Match_subcell_tag : public dynamicprogramming_subcell_tag
  {
    // Identifying tag for the match subcell
  };
  struct dynamicprogramming_Insertion_subcell_tag : public dynamicprogramming_subcell_tag
  {
    // Identifying tag for the insertion subcell
  };
  struct dynamicprogramming_Deletion_subcell_tag : public dynamicprogramming_subcell_tag
  {
    // Identifying tag for the deletion subcell
  };
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  struct dynamicprogramming_DeletionIn_subcell_tag : public dynamicprogramming_subcell_tag
  {
    // Identifying tag for the deletion-in subcell
  };
  struct dynamicprogramming_DeletionOut_subcell_tag : public dynamicprogramming_subcell_tag
  {
    // Identifying tag for the deletion-in subcell
  };
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

  class Subcell : public dynamicprogramming_subcell_tag
  {
  protected:
    enum Type {
      Match = 0,
      Insertion,
      Deletion
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      ,
      DeletionIn,
      DeletionOut
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
    } m_type;
  public:
    Subcell () :
      m_type( Match )
    {
      // Do nothing else
    } // <init>()

    Subcell ( Subcell const& other ) :
      m_type( other.m_type )
    {
      // Do nothing else
    } // <init>( Subcell const& )

    Subcell ( dynamicprogramming_Match_subcell_tag const unused ) :
      m_type( Match )
    {
      // Do nothing else
    } // <init>( match )
    Subcell ( dynamicprogramming_Insertion_subcell_tag const unused ) :
      m_type( Insertion )
    {
      // Do nothing else
    } // <init>( insertion )
    Subcell ( dynamicprogramming_Deletion_subcell_tag const unused ) :
      m_type( Deletion )
    {
      // Do nothing else
    } // <init>( deletion )
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    Subcell ( dynamicprogramming_DeletionIn_subcell_tag const unused ) :
      m_type( DeletionIn )
    {
      // Do nothing else
    } // <init>( DeletionIn )
    Subcell ( dynamicprogramming_DeletionOut_subcell_tag const unused ) :
      m_type( DeletionOut )
    {
      // Do nothing else
    } // <init>( DeletionOut )
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    Subcell & operator= ( Subcell const& other )
    {
      m_type = other.m_type;
      return *this;
    } // operator=( Subcell const & )
    Subcell & operator= ( dynamicprogramming_Match_subcell_tag const )
    {
      m_type = Match;
      return *this;
    } // operator=( match )
    Subcell & operator= ( dynamicprogramming_Insertion_subcell_tag const )
    {
      m_type = Insertion;
      return *this;
    } // operator=( insertion )
    Subcell & operator= ( dynamicprogramming_Deletion_subcell_tag const )
    {
      m_type = Deletion;
      return *this;
    } // operator=( deletion )
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    Subcell & operator= ( dynamicprogramming_DeletionIn_subcell_tag const )
    {
      m_type = DeletionIn;
      return *this;
    } // operator=( DeletionIn )
    Subcell & operator= ( dynamicprogramming_DeletionOut_subcell_tag const )
    {
      m_type = DeletionOut;
      return *this;
    } // operator=( DeletionOut )
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    bool operator== ( dynamicprogramming_Match_subcell_tag const ) const
    {
      return( m_type == Match );
    } // operator==( match ) const
    bool operator== ( dynamicprogramming_Insertion_subcell_tag const ) const
    {
      return( m_type == Insertion );
    } // operator==( insertion ) const
    bool operator== ( dynamicprogramming_Deletion_subcell_tag const ) const
    {
      return( m_type == Deletion );
    } // operator==( deletion ) const
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    bool operator== ( dynamicprogramming_DeletionIn_subcell_tag const ) const
    {
      return( m_type == DeletionIn );
    } // operator==( DeletionIn ) const
    bool operator== ( dynamicprogramming_DeletionOut_subcell_tag const ) const
    {
      return( m_type == DeletionOut );
    } // operator==( DeletionOut ) const
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    void
    setToMatch ()
    {
      m_type = Match;
    } // setToMatch()
    void
    setToInsertion ()
    {
      m_type = Insertion;
    } // setToInsertion()
    void
    setToDeletion ()
    {
      m_type = Deletion;
    } // setToDeletion()
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    void
    setToDeletionIn ()
    {
      m_type = DeletionIn;
    } // setToDeletionIn()
    void
    setToDeletionOut ()
    {
      m_type = DeletionOut;
    } // setToDeletionOut()
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    bool
    isMatch () const
    {
      return( m_type == Match );
    } // isMatch() const
    bool
    isInsertion () const
    {
      return( m_type == Insertion );
    } // isInsertion() const
    bool
    isDeletion () const
    {
      return( m_type == Deletion );
    } // isDeletion() const
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    bool
    isDeletionIn () const
    {
      return( m_type == DeletionIn );
    } // isDeletionIn() const
    bool
    isDeletionOut () const
    {
      return( m_type == DeletionOut );
    } // isDeletionOut() const
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      Subcell const &subcell )
    {
      if( subcell.m_type == Match ) {
        os << "Match";
      } else if( subcell.m_type == Insertion ) {
        os << "Insertion";
      } else if( subcell.m_type == Deletion ) {
        os << "Deletion";
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      } else if( subcell.m_type == DeletionIn ) {
        os << "DeletionIn";
      } else if( subcell.m_type == DeletionOut ) {
        os << "DeletionOut";
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      } else {
        /// \todo ?
        os << "ERROR!";
      }

      return os;
    } // friend operator<< ( basic_ostream, Subcell const & )
  }; // End class Subcell

static dynamicprogramming_Match_subcell_tag const Match =
              dynamicprogramming_Match_subcell_tag();
static dynamicprogramming_Insertion_subcell_tag const Insertion =
              dynamicprogramming_Insertion_subcell_tag();
static dynamicprogramming_Deletion_subcell_tag const Deletion =
              dynamicprogramming_Deletion_subcell_tag();
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
static dynamicprogramming_DeletionIn_subcell_tag const DeletionIn =
              dynamicprogramming_DeletionIn_subcell_tag();
static dynamicprogramming_DeletionOut_subcell_tag const DeletionOut =
              dynamicprogramming_DeletionOut_subcell_tag();
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef ALLOW_BOLTZMANN_GIBBS
    /**
     * If we use Baldi-style gradient ascent, we convert the profile to
     * Boltzmann-Gibbs scale.  Here is where we store it.  It should be a type
     * that can represent the entire Real range.
     */
    /// \todo PUT BACK!  TESTING.
    ///typedef ProbabilityType BoltzmannGibbsValueType;
    ///typedef bfloat BoltzmannGibbsValueType;
    ///typedef logspace BoltzmannGibbsValueType;
    typedef realspace BoltzmannGibbsValueType;
    ///typedef double BoltzmannGibbsValueType;

    // For Baldi-style gradient ascent:
    template <typename ResidueType, typename ScoreType>
    class PositionBoltzmannGibbs :
    public ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, BoltzmannGibbsValueType> >
    {
      // Same as the PositionEntente but with a temperature...
    protected:
      BoltzmannGibbsValueType m_temperature;

    public:
      PositionBoltzmannGibbs () :
        ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, BoltzmannGibbsValueType> >(),
        m_temperature( 1.0 )
      {
        // Do nothing else.
      } // <init>()

      //template <typename scalar_type>
      //PositionBoltzmannGibbs (
      //  scalar_type const & temperature
      //) :
      //  ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, BoltzmannGibbsValueType> >(),
      //  m_temperature( temperature )
      //{
      //  // Do nothing else.
      //} // <init>( scalar_type const & )

      template <typename AnyProbabilityType>
      PositionBoltzmannGibbs (
        PositionSpecificParameters<ResidueType, AnyProbabilityType> const & profile_position
      ) :
        ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, BoltzmannGibbsValueType> >(),
        m_temperature( 1.0 )
      {
        fromProfilePosition( profile_position );
      } // <init>( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& )

      template <typename AnyProbabilityType,
                typename scalar_type>
      PositionBoltzmannGibbs (
        PositionSpecificParameters<ResidueType, AnyProbabilityType> const & profile_position,
        scalar_type const & temperature
      ) :
        ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, BoltzmannGibbsValueType> >(),
        m_temperature( temperature )
      {
        fromProfilePosition( profile_position );
      } // <init>( PositionSpecificParameters<ResidueType, AnyProbabilityType> const&, scalar_type const & )

      void
      reinitialize ()
      {
        reinitialize( 1.0 );
      }

      template <typename scalar_type>
      void
      reinitialize (
        scalar_type const & temperature
      )
      {
        ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, BoltzmannGibbsValueType> >::reinitialize();
        m_temperature = temperature;
      } // reinitialize( scalar_type const & )

      /**
       * Make this PositionBoltzmannGibbs the Boltzmann-Gibbs equivalent of the
       * given profile position (using the current m_temperature).  It also
       * resets the m_scalar to 1.
       */
      template <typename AnyProbabilityType>
      void
      fromProfilePosition (
        PositionSpecificParameters<ResidueType, AnyProbabilityType> const & profile_position
      )
      {
        this->m_scalar = 1.0;
        profile_position[ Emission::Match ].toBoltzmannGibbs(
          m_temperature,
          ( *this )[ Emission::Match ]
        );
      } // fromProfilePosition( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& )
      
      /**
       * Fill the given profile position with the non-BoltzmannGibbs equivalent
       * of this PositionBoltzmannGibbs.
       */
      template <typename AnyProbabilityType>
      void
      toProfilePosition (
        PositionSpecificParameters<ResidueType, AnyProbabilityType> & profile_position
      )
      {
        profile_position[ Emission::Match ].fromBoltzmannGibbs(
          m_temperature,
          ( *this )[ Emission::Match ]
        );
      } // toProfilePosition( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& )

      PositionBoltzmannGibbs &
      operator= ( PositionBoltzmannGibbs const& other_pos )
      {
        ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, BoltzmannGibbsValueType> >::operator=( other_pos );
        m_temperature = other_pos.m_temperature;
        return *this;
      } // operator=( PositionBoltzmannGibbs const& )

      // Overridden to not actually rescale.
      void
      rescale (
         BoltzmannGibbsValueType const & max
      )
      {
        // Do nothing.
      } // rescale( BoltzmannGibbsValueType const & )

    }; // End class PositionBoltzmannGibbs

  // Explicitly confirm that PositionBoltzmannGibbs and
  // PositionBoltzmannGibbsChange are scalable.
  template <typename ResidueType, typename ScoreType>
  struct has_m_scalar<PositionBoltzmannGibbs<ResidueType, ScoreType> > : boost::true_type
  {}; // PositionBoltzmannGibbs has_m_scalar
#endif // ALLOW_BOLTZMANN_GIBBS

    /**
     * An AlignmentProfile stores the total probability of all paths that flow
     * through each state (for a particular profile, sequence pair).  This
     * class represents one position of the alignment profile (corresponding to
     * one position of the dp matrix)
     */
    template <typename ResidueType, typename ParameterType>
    class AlignmentProfilePositionParameters :
      public PositionSpecificParameters<ResidueType, ParameterType>,
      public GlobalParameters<ResidueType, ParameterType>
    {
      // Boost serialization
    private:
      typedef PositionSpecificParameters<ResidueType, ParameterType> PositionSpecificParameters_ParameterType;
      typedef GlobalParameters<ResidueType, ParameterType> GlobalParameters_ParameterType;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PositionSpecificParameters_ParameterType );
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( GlobalParameters_ParameterType );
      } // serialize( Archive &, const unsigned int )

    public:
      typedef ResidueType APPPResidueType;

      AlignmentProfilePositionParameters ();

      void
      reinitialize ();

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.
       */
      template <typename other_type>
      void
      normalize ( other_type const & min );

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.
       */
      void
      normalize ( ParameterType const& min );

      /**
       * Set all values to 0.
       */
      void
      zero ();

      AlignmentProfilePositionParameters &
      operator= ( AlignmentProfilePositionParameters const& other_pos );

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        AlignmentProfilePositionParameters const& node
      )
      {
        os << "[ ";
        node.writeAlignmentProfilePositionParameters( os );
        os << " ]";
        return os;
      } // friend operator<< ( basic_ostream &, AlignmentProfilePositionParameters const& )
    
      /**
       * Change the probabilities of the contained Multinomials to values drawn
       * from dirichlet distributions with the given counts.  AnyCountType can be
       * a Real type (anything coercible to a double using toDouble( count )).
       */
      template <typename AnyAlignmentProfilePositionParametersType>
      void
      dirichlet (
        AnyAlignmentProfilePositionParametersType const & counts,
        Random & random
      );

      /**
       * Write a comma-separated list of the parameters to the stream.
       */
      template<class CharT, class Traits>
      void
      writeAlignmentProfilePositionParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const;

      /**
       * Write a comma-separated list of the parameters to the stream.
       */
      template<class CharT, class Traits>
      void
      writeParameterCollection (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        writeAlignmentProfilePositionParameters( os );
      } // writeParameterCollection( basic_ostream & ) const

      // operators
      GALOSH_OPERATORS_PositionSpecificParameters( ParameterType )
      GALOSH_OPERATORS_GlobalParameters( ParameterType )

      /**
       * Divide each contained distribution value by denominator.
       */
      template <typename AnyValueType>
      AlignmentProfilePositionParameters &
      operator/= ( AnyValueType const& denominator );
  
      /**
       * Multiply each contained distribution value by scalar.
       */
      template <typename AnyValueType>
      AlignmentProfilePositionParameters &
      operator*= ( AnyValueType const& scalar );
  
      /**
       * Add to each contained distribution the values in the given other
       * AlignmentProfilePositionParameters.  Note that this may violate the rule that
       * the probabilities sum to 1.
       */
      AlignmentProfilePositionParameters &
      operator+= ( AlignmentProfilePositionParameters const& other_pos );
  
      /**
       * Set all values such that each distrubution is randomly distributed.
       */
      void
      uniform ( Random & random );
  
      /**
       * Set all values such that each distrubution is evenly distributed.
       */
      void
      even ();
  
      /**
       * Return the largest value in the contained distributions.
       */
      ParameterType
      maximumValue () const;
  
      /**
       * Calculate and return the Euclidean distance between this position and
       * another position (treating every probability as an orthogonal
       * dimension).
       */
      double
      euclideanDistance (
        AlignmentProfilePositionParameters const& other_pos
      ) const;

      /**
       * Calculate and return the square of the Euclidean distance between this
       * position and another position (treating every probability as an
       * orthogonal dimension).
       */
      double
      euclideanDistanceSquared (
        AlignmentProfilePositionParameters const& other_node
      ) const;

      /**
       * Calculate and return the probability that the Match state of this
       * alignment profile position is identical to the Match state of the
       * given alignment profile position, using the Match state emission
       * values of the two alignment profile positions as the relevant
       * probabilities.  Note that these are (usually) not conditioned on there
       * being a Match at all (that is, the probs don't sum to 1; instead they
       * sum to 1-P(Deletion)).  Note also that this is symmetric, so calling
       * the other position's calculateMatchEmissionCooccurrenceProbability
       * with this pos as an argument will return the same thing.
       */
      ParameterType
      calculateMatchEmissionCooccurrenceProbability (
        AlignmentProfilePositionParameters const & other_pos,
        const bool include_nonemission_prob = false
      ) const;

      /**
       * Calculate and return the expected number of times that the Match state
       * of this alignment profile position is identical to the Match state of
       * the given alignment profile position, using the Match state emission
       * values of the two alignment profile positions as the relevant
       * probabilities.  Note that these are (usually) not conditioned on there
       * being a Match at all (that is, the probs don't sum to 1; instead they
       * sum to 1-P(Deletion)).  Note also that this is symmetric, so calling
       * the other position's calculateMatchEmissionCooccurrenceExpectedCount
       * with this pos as an argument will return the same thing.
       */
      ParameterType
      calculateMatchEmissionCooccurrenceExpectedCount (
        AlignmentProfilePositionParameters const & other_pos,
        const bool include_nonemission_prob = false
      ) const;

      /**
       * Calculate and return the cross entropy E(-log(other_node)).  Note that
       * the cross entropy is non-symmetric (calling other_node.crossEntropy(
       * *this ) will return a different value).
       *
       * This can be used to calculate the (self-)entropy by calling
       * this->crossEntropy( *this ).  It can also be used to calculate the KL
       * divergence by taking the difference of the cross-entropy and the
       * self-entropy, or the symmetrized KL divergence by summing the KL
       * divergences computed both ways.
       *
       * See also the other crossEntropy method, which accepts a weights argument.
       */
      template <typename AnyParameterType>
      double
      crossEntropy (
        AlignmentProfilePositionParameters<ResidueType, AnyParameterType> const& other_node
      ) const;

      /**
       * Calculate and return the (possibly weighted) cross entropy
       * E(-log(other_pos)).  The weights (if non-null) may be of any type
       * convertible to a double, and the cross entropy will be computed as
       * E(-log(weights*other_pos)).  Note that the cross entropy is
       * non-symmetric (calling other_pos.crossEntropy( *this, weights ) will
       * return a different value).
       *
       * This can be used to calculate the (possibly weighted) (self-)entropy by
       * calling this->crossEntropy( *this, weights ).  It can also be used to
       * calculate the KL divergence by taking the difference of the
       * cross-entropy and the self-entropy, or the symmeterized KL divergence by
       * summing the KL divergences computed both ways.  Note, though, that
       * weights should all be the same within any Multinomial distribution if
       * this is to be used to calculate weighted entropies or KL divergences.
       * Otherwise the usual properties of these metrics will be violated.
       */
      template <typename AnyParameterType,
                typename AnyAlignmentProfilePositionParametersType>
      double
      crossEntropy (
        AlignmentProfilePositionParameters<ResidueType, AnyParameterType> const& other_pos,
        AnyAlignmentProfilePositionParametersType const * const weights
      ) const;

    }; // End class AlignmentProfilePositionParameters

  ////// Class galosh::AlignmentProfilePositionParameters ////
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_INIT
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      AlignmentProfilePositionParameters () :
        PositionSpecificParameters<ResidueType, ParameterType>(),
        GlobalParameters<ResidueType, ParameterType>()
      {
        // Do nothing else.
      } // <init>()

  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_REINITIALIZE
  void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
  reinitialize ()
  {
    PositionSpecificParameters<ResidueType, ParameterType>::reinitialize();
    GlobalParameters<ResidueType, ParameterType>::reinitialize();
  } // reinitialize()

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.
       */
  template <typename ResidueType, typename ParameterType>
      template <typename other_type>
  GALOSH_INLINE_TRIVIAL
      void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      normalize ( other_type const & min )
      {
        normalize( static_cast<ParameterType>( min ) );
      } // normalize( other_type const & )

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.
       */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_TRIVIAL
      void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      normalize ( ParameterType const& min )
      {
        PositionSpecificParameters<ResidueType, ParameterType>::normalize( min );
        GlobalParameters<ResidueType, ParameterType>::normalize( min );
      } // normalize( ParameterType const& )

      /**
       * Set all values to 0.
       */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_TRIVIAL
      void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      zero ()
      {
        PositionSpecificParameters<ResidueType, ParameterType>::zero();
        GlobalParameters<ResidueType, ParameterType>::zero();
      } // zero()

  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_TRIVIAL
      void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      even ()
      {
        PositionSpecificParameters<ResidueType, ParameterType>::even();
        GlobalParameters<ResidueType, ParameterType>::even();
      } // even()

  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_TRIVIAL
  AlignmentProfilePositionParameters<ResidueType, ParameterType> &
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      operator= ( AlignmentProfilePositionParameters const& other_pos )
      {
        PositionSpecificParameters<ResidueType, ParameterType>::operator=( other_pos );
        GlobalParameters<ResidueType, ParameterType>::operator=( other_pos );

        return *this;
      } // operator=( AlignmentProfilePositionParameters const& )

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
  template <typename ResidueType, typename ParameterType>
  template <typename AnyAlignmentProfilePositionParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
    dirichlet (
      AnyAlignmentProfilePositionParametersType const & counts,
      Random & random
    )
    {
      PreAlignParameters<ResidueType, ParameterType>::dirichlet( counts, random );
      PostAlignParameters<ResidueType, ParameterType>::dirichlet( counts, random );
      InsertionEmissionParameters<ResidueType, ParameterType>::dirichlet( counts, random );
      PositionTransitionParameters<ResidueType, ParameterType>::dirichlet( counts, random );
    } // dirichlet( AnyAlignmentProfilePositionParametersType const &, Random & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ParameterType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
    void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
    writeAlignmentProfilePositionParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      this->writePositionSpecificParameters( os );
      os << ", ";
      this->writeGlobalParameters( os );
    } // writeAlignmentProfilePositionParameters( basic_ostream & )

    /**
     * Divide each contained distrubution value by denominator.
     */
  template <typename ResidueType, typename ParameterType>
    template <typename AnyValueType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  AlignmentProfilePositionParameters<ResidueType, ParameterType> &
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
    operator/= ( AnyValueType const& denominator )
    {
      PositionSpecificParameters<ResidueType, ParameterType>::operator/=( denominator );
      GlobalParameters<ResidueType, ParameterType>::operator/=( denominator );

      return *this;
    } // operator/=( AnyValueType const& )

      /**
       * Multiply each contained distrubution value by scalar.
       */
  template <typename ResidueType, typename ParameterType>
    template <typename AnyValueType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  AlignmentProfilePositionParameters<ResidueType, ParameterType> &
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
    operator*= ( AnyValueType const& scalar )
    {
      PositionSpecificParameters<ResidueType, ParameterType>::operator*=( scalar );
      GlobalParameters<ResidueType, ParameterType>::operator*=( scalar );

      return *this;
    } // operator*=( AnyValueType const& )

    /**
     * Add to each contained distribution the values in the given other
     * ProfilePosition.  Note that this may violate the rule that the
     * probabilities sum to 1.
     */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  AlignmentProfilePositionParameters<ResidueType, ParameterType> &
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
    operator+= ( AlignmentProfilePositionParameters const& other_pos )
    {
      PositionSpecificParameters<ResidueType, ParameterType>::operator+=( other_pos );
      GlobalParameters<ResidueType, ParameterType>::operator+=( other_pos );

      return *this;
    } // operator+=( PositionSpecificParameters const& )

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_REINITIALIZE
    void
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
    uniform ( Random & random )
    {
      PositionSpecificParameters<ResidueType, ParameterType>::uniform( random );
      GlobalParameters<ResidueType, ParameterType>::uniform( random );
    } // uniform( Random & )

    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ParameterType
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
    maximumValue () const
    {
      ParameterType largest_value =
        PositionSpecificParameters<ResidueType, ParameterType>::maximumValue();
      ParameterType largest_value_tmp =
        GlobalParameters<ResidueType, ParameterType>::maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
      return largest_value;
    } // maximumValue() const

  /**
   * Calculate and return the Euclidean distance between this set of
   * transition parameters and another set (treating every probability as an
   * orthogonal dimension).
   */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_TRIVIAL
  double
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
  euclideanDistance (
    AlignmentProfilePositionParameters const& other_pos
  ) const
  {
    return sqrt( euclideanDistanceSquared( other_pos ) );
  } // euclideanDistance( AlignmentProfilePositionParameters const & )
  
  /**
   * Calculate and return the square of the Euclidean distance between this
   * position and another position (treating every probability as an
   * orthogonal dimension).
   */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
  double
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
  euclideanDistanceSquared (
    AlignmentProfilePositionParameters const& other_pos
  ) const
  {
      double squared_euclidean_distance = 0.0;

      squared_euclidean_distance +=
        PositionSpecificParameters<ResidueType, ParameterType>::euclideanDistanceSquared( other_pos );
      squared_euclidean_distance +=
        GlobalParameters<ResidueType, ParameterType>::euclideanDistanceSquared( other_pos );

      return squared_euclidean_distance;
  } // euclideanDistanceSquared( AlignmentProfilePositionParameters const& )

      /**
       * Calculate and return the probability that the Match state of this
       * alignment profile position is identical to the Match state of the
       * given alignment profile position, using the Match state emission
       * values of the two alignment profile positions as the relevant
       * probabilities.  Note that these are (usually) not conditioned on there
       * being a Match at all (that is, the probs don't sum to 1; instead they
       * sum to 1-P(Deletion)).  If the include_nonemission_prob argument is
       * true, the probability of a co-occurring deletion will be included in the
       * calculation. Note also that this is symmetric, so calling
       * the other position's calculateMatchEmissionCooccurrenceProbability
       * with this pos as an argument will return the same thing.
       */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_TRIVIAL
      ParameterType
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      calculateMatchEmissionCooccurrenceProbability (
        AlignmentProfilePositionParameters const & other_pos,
        const bool include_nonemission_prob
      ) const
      {
        ParameterType prob =
          galosh::calculateEmissionCooccurrenceProbability<ParameterType>( *this, other_pos, Emission::Match );

        if( include_nonemission_prob ) {
          /// Include the probability that there's no match at all.
          /// \todo PUT BACK
          prob += ( ( 1.0 - ( *this )[ Emission::Match ].total() ) * ( 1.0 - other_pos[ Emission::Match ].total() ) );
          /// \todo REMOVE.  TESTING.
          ///ParameterType tmp1 = ( *this )[ Emission::Match ].total();
          ///if( tmp1 == 1.0 ) {
          ///  return prob;
          ///}
          ///if( tmp1 > 1.0 ) {
          ///  // Is the > operator busted for Algebra?
          ///  if( tmp1 == 1 ) {
          ///    cerr << "IT IS 1" << endl;
          ///  }
          ///  cerr << "Uh oh: this has a prob(Match) = " << tmp1 << " > 1.0!?: it is " << (*this) << endl;
          ///}
          ///assert( tmp1 <= 1.0 );
          ///ParameterType tmp2 = other_pos[ Emission::Match ].total();
          ///if( tmp2 == 1.0 ) {
          ///  return prob;
          ///}
          //if( tmp2 > 1.0 ) {
          //  cerr << "Uh oh: other_pos has a prob(Match) = " << tmp2 << " > 1.0!?: it is " << other_pos << endl;
          //}
          //assert( tmp2 <= 1.0 );
          //prob += ( ( 1.0 - tmp1 ) * ( 1.0 - tmp2 ) );
        } // End if include_nonemission_prob

        return prob;
      } // calculateMatchEmissionCooccurrenceProbability ( AlignmentProfilePositionParameters const & ) const

      /**
       * Calculate and return the expected number of times that the Match state
       * of this alignment profile position is identical to the Match state of
       * the given alignment profile position, using the Match state emission
       * values of the two alignment profile positions as the relevant
       * probabilities.  Note that these are (usually) not conditioned on there
       * being a Match at all (that is, the probs don't sum to 1; instead they
       * sum to 1-P(Deletion)). If the include_nonemission_counts argument is
       * true, the probability of a cooccurring deletion will be included in the
       * calculation. Note also that this is symmetric, so calling
       * the other position's calculateMatchEmissionCooccurrenceExpectedCount
       * with this pos as an argument will return the same thing.
       */
  template <typename ResidueType, typename ParameterType>
  GALOSH_INLINE_TRIVIAL
      ParameterType
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      calculateMatchEmissionCooccurrenceExpectedCount (
        AlignmentProfilePositionParameters<ResidueType, ParameterType> const & other_pos,
        const bool include_nonemission_counts
      ) const
      {
        ParameterType ec =
          // defined in Profile.hpp
          galosh::calculateEmissionCooccurrenceProbability<ParameterType>( *this, other_pos, Emission::Match );

        if( include_nonemission_counts ) {
          // Include the probability that there's no match at all.
          ec += ( ( 1.0 - ( *this )[ Emission::Match ].total() ) * ( 1.0 - other_pos[ Emission::Match ].total() ) );
        } // End if include_nonemission_counts

        return ec;
      } // calculateMatchEmissionCooccurrenceExpectedCount ( AlignmentProfilePositionParameters const & ) const

      /**
       * Calculate and return the cross entropy
       * E(-log(other_node)).  Note that the cross entropy is
       * non-symmetric (calling other_node.crossEntropy( *this ) will
       * return a different value).
       *
       * This can be used to calculate the (self-)entropy by
       * calling this->crossEntropy( *this ).  It can also be used to
       * calculate the KL divergence by taking the difference of the
       * cross-entropy and the self-entropy, or the symmeterized KL divergence by
       * summing the KL divergences computed both ways.
       */
  template <typename ResidueType, typename ParameterType>
  template <typename AnyParameterType>
  GALOSH_INLINE_TRIVIAL
      double
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      crossEntropy (
        AlignmentProfilePositionParameters<ResidueType, AnyParameterType> const& other_pos
      ) const
    {
      return
        crossEntropy(
          other_pos,
          ( AlignmentProfilePositionParameters const * )0
        );
    } // crossEntropy( AlignmentProfilePositionParameters const & )

      /**
       * Calculate and return the (possibly weighted) cross entropy
       * E(-log(other_pos)).  The weights (if non-null) may be of any type
       * convertible to a double, and the cross entropy will be computed as
       * E(-log(weights*other_pos)).  Note that the cross entropy is
       * non-symmetric (calling other_pos.crossEntropy( *this, weights ) will
       * return a different value).
       *
       * This can be used to calculate the (possibly weighted) (self-)entropy by
       * calling this->crossEntropy( *this, weights ).  It can also be used to
       * calculate the KL divergence by taking the difference of the
       * cross-entropy and the self-entropy, or the symmeterized KL divergence by
       * summing the KL divergences computed both ways.  Note, though, that
       * weights should all be the same within any Multinomial distribution if
       * this is to be used to calculate weighted entropies or KL divergences.
       * Otherwise the usual properties of these metrics will be violated.
       */
  template <typename ResidueType, typename ParameterType>
  template <typename AnyParameterType,
            typename AnyAlignmentProfilePositionParametersType>
  GALOSH_INLINE_TRIVIAL
      double
  AlignmentProfilePositionParameters<ResidueType, ParameterType>::
      crossEntropy (
        AlignmentProfilePositionParameters<ResidueType, AnyParameterType> const& other_pos,
        AnyAlignmentProfilePositionParametersType const * const weights
      ) const
    {
      double cross_entropy = 0.0;
      cross_entropy +=
        PositionSpecificParameters<ResidueType, ParameterType>::crossEntropy(
          other_pos,
          weights
        );
      cross_entropy +=
        GlobalParameters<ResidueType, ParameterType>::crossEntropy(
          other_pos,
          weights
        );

      return cross_entropy;
    } // crossEntropy( AlignmentProfilePositionParameters<ResidueType, AnyParameterType> const &, AnyAlignmentProfilePositionParametersType const * ) const


  template <typename ResidueType, typename ParameterType>
  struct parameter_collection_traits<AlignmentProfilePositionParameters<ResidueType, ParameterType> >
  {
    typedef ParameterType ProbabilityType;
  }; // AlignmentProfilePositionParameters parameter_collection_traits

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  class DynamicProgramming
  {
  public:
    class Parameters :
      public galosh::Parameters
    {
      // Boost serialization
    private:
      typedef galosh::Parameters galosh_parameters_t;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( galosh_parameters_t );

        ar & BOOST_SERIALIZATION_NVP( useRabinerScaling );
        ar & BOOST_SERIALIZATION_NVP( rabinerScaling_useMaximumValue );
        ar & BOOST_SERIALIZATION_NVP( matrixRowScaleFactor );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// PARAMETERS
      /**
       * Use (our rotated, state-specific) Rabiner-style scaling to keep
       * forward matrix values from underflowing?  Note that if this is false,
       * you should use a log type for the ScoreType (unless the profile and
       * sequences are very short and there are few sequences, in which case
       * you might not encounter underflow).
       */
      bool useRabinerScaling;
      /// NOTE THAT RIGHT NOW, RABINER SCALING SEEMS TO BE BROKEN.  \todo FIX!
  #define DEFAULT_useRabinerScaling false

      /**
       * It turns out that rabinerScaling can use any scale.  Rabiner
       * suggested scaling each time separately, by the total of the matrix
       * values.  We rotate, so we scale each state separately (not each time)
       * -- but actually we group together the states in each position, so we
       * scale the positions.  We can either scale them by the total of the
       * values there (or by the total of just the Match and Deletion values)
       * or by the maximum value, which lets us use more of the range of the
       * MatrixValueType.  We use the total of all values in the position
       * unless rabinerScaling_useMaximumValue is true.
       */
      bool rabinerScaling_useMaximumValue;
      /// \todo Make this depend on the MatrixValueType?  If using bfloats or logs, we don't need to do this.  We can avoid a few calculations if we have this false.
//#define DEFAULT_rabinerScaling_useMaximumValue true
/// \todo put back true.  Testing.
#define DEFAULT_rabinerScaling_useMaximumValue false

      /**
       * Vanilla rabiner scaling would divide all matrix row values by the
       * total of all values in that row, but when there is a wide discrepency
       * between the largest and smallest values, the smallest ones will
       * underflow (with a non-log MatrixValueType).  So we can multiply
       * everything in the row by this scale factor to use more of the range of
       * the MatrixValueType.
       */
      double matrixRowScaleFactor;
      /// \todo Make this depend on the MatrixValueType?  If using bfloats or logs, we don't need to scale it.
//#define DEFAULT_matrixRowScaleFactor ( pow( numeric_limits<float>::max(), .25f ) - 1.0f )
//#define DEFAULT_matrixRowScaleFactor ( pow( numeric_limits<double>::max(), .0625 ) - 1.0 )
#define DEFAULT_matrixRowScaleFactor 1

      Parameters ();
      virtual ~Parameters () {};
    
      // Copy constructor
      template <class AnyParameters>
      Parameters ( const AnyParameters & copy_from );
    
      // Copy constructor/operator
      template <class AnyParameters>
      Parameters & operator= (
        const AnyParameters & copy_from
      );
    
      template <class AnyParameters>
      void
      copyFromNonVirtual (
        AnyParameters const & copy_from
      );

      template <class AnyParameters>
      void
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      );

      virtual void
      copyFrom ( const Parameters & copy_from );
    
      virtual void
      resetToDefaults ();

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        Parameters const& parameters
      )
      {
        parameters.writeParameters( os );
        return os;
      } // friend operator<< ( basic_ostream &, Parameters const& )

      template<class CharT, class Traits>
      void
      writeParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const;

    }; // End inner class Parameters

    template <class ParametersType>
    class ParametersModifierTemplate :
      public galosh::ParametersModifierTemplate<ParametersType>
    {
      typedef typename galosh::ParametersModifierTemplate<ParametersType> base_parameters_modifier_t;

      // Boost serialization
    private:
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information.  This will serialize the
        // parameters too.
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( base_parameters_modifier_t );

        // Serialize the new isModified_ stuff
        ar & BOOST_SERIALIZATION_NVP( isModified_useRabinerScaling );
        ar & BOOST_SERIALIZATION_NVP( isModified_rabinerScaling_useMaximumValue );
        ar & BOOST_SERIALIZATION_NVP( isModified_matrixRowScaleFactor );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// isModified flags for Parameters
      /**
       * Use (our rotated, state-specific) Rabiner-style scaling to keep
       * forward matrix values from underflowing?  Note that if this is false,
       * you should use a log type for the ScoreType (unless the profile and
       * sequences are very short and there are few sequences, in which case
       * you might not encounter underflow).
       */
      bool isModified_useRabinerScaling;

      /**
       * It turns out that rabinerScaling can use any scale.  Rabiner
       * suggested scaling each time separately, by the total of the matrix
       * values.  We rotate, so we scale each state separately (not each time)
       * -- but actually we group together the states in each position, so we
       * scale the positions.  We can either scale them by the total of the
       * values there (or by the total of just the Match and Deletion values)
       * or by the maximum value, which lets us use more of the range of the
       * MatrixValueType.  We use the total of all values in the position
       * unless rabinerScaling_useMaximumValue is true.
       */
      bool isModified_rabinerScaling_useMaximumValue;

      /**
       * Vanilla rabiner scaling would divide all matrix row values by the
       * total of all values in that row, but when there is a wide discrepency
       * between the largest and smallest values, the smallest ones will
       * underflow (with a non-log MatrixValueType).  So we can multiply
       * everything in the row by this scale factor to use more of the range of
       * the MatrixValueType.
       */
      bool isModified_matrixRowScaleFactor;

      //using galosh::ParametersModifierTemplate<ParametersType>::parameters;
      //using base_parameters_modifier_t::parameters;

      ParametersModifierTemplate ();
    
      // Copy constructor
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from );
    
      // Copy constructor/operator
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate & operator= (
        const AnyParametersModifierTemplate & copy_from
      );
    
      template <class AnyParametersModifierTemplate>
      void
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );

      template <class AnyParametersModifierTemplate>
      void
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );

      void
      reset ();

      void
      isModified_reset ();

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        ParametersModifierTemplate const& parameters_modifier
      )
      {
        parameters_modifier.writeParametersModifier( os );

        return os;
      } // friend operator<< ( basic_ostream &, ParametersModifierTemplate const& )

      template<class CharT, class Traits>
      void
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      ) const;

      template<class AnyParameters>
      void
      applyModifications ( AnyParameters & target_parameters ) const;

    }; // End inner class ParametersModifierTemplate

    typedef ParametersModifierTemplate<typename DynamicProgramming::Parameters> ParametersModifier;

    ///// Static utility fns
    /**
     * Utility fn for getting match emission probabilities from a profile
     * position.
     */
    template <typename ProfileType, typename SequenceResidueType>
    static
    ProbabilityType
    getEmissionProbability (
      ProfilePosition<ResidueType, ProfileType> const& pos,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      return( pos[ Emission::Match ][ sequence[ seq_pos_i ] ] );
    } // static getEmissionProbability ( ProfilePosition const&, Sequence<SequenceResidueType> const&, uint32_t const )

    /**
     * Utility fn for getting insertion emission probabilities from a profile position.
     */
    template <typename ProfileType, typename SequenceResidueType>
    static
    ProbabilityType
    getInsertionEmissionProbability (
      ProfilePosition<ResidueType, ProfileType> const& pos,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      return( pos[ Emission::Insertion ][ sequence[ seq_pos_i ] ] );
    } // static getInsertionEmissionProbability ( ProfilePosition const&, Sequence<SequenceResidueType> const&, uint32_t const )

    /**
     * Utility fn for getting pre-align emission probabilities from a profile.
     */
    template <typename ProfileType, typename SequenceResidueType>
    static
    ProbabilityType
    getPreAlignEmissionProbability (
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      return( profile[ Emission::PreAlignInsertion ][ sequence[ seq_pos_i ] ] );
    } // static getPreAlignEmissionProbability ( Profile const&, Sequence<SequenceResidueType> const&, uint32_t const )

    /**
     * Utility fn for getting post-align emission probabilities from a profile.
     */
    template <typename ProfileType, typename SequenceResidueType>
    static
    ProbabilityType
    getPostAlignEmissionProbability (
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      return( profile[ Emission::PostAlignInsertion ][ sequence[ seq_pos_i ] ] );
    } // static getPostAlignEmissionProbability ( Profile const&, Sequence<SequenceResidueType> const&, uint32_t const )

    /**
     * Utility fn for updating match emission counts for possibly-
     * polymorphic sequences.
     */
    template <typename PositionCountsType, typename SequenceResidueType>
    static
    void
    incrementEmissionCounts (
      PositionCountsType & position_counts,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      position_counts[ Emission::Match ][ sequence[ seq_pos_i ] ] += 1;
    } // static incrementEmissionCounts ( PositionCountsType &, Sequence<SequenceResidueType> const&, uint32_t const )

    /**
     * Utility fn for incrementing insertion emission probabilities..
     */
    template <typename PositionCountsType, typename SequenceResidueType>
    static
    void
    incrementInsertionEmissionCounts (
      PositionCountsType & position_counts,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      position_counts[ Emission::Insertion ][ sequence[ seq_pos_i ] ] += 1;
    } // static incrementInsertionEmissionCounts ( PositionCountsType &, Sequence<SequenceResidueType> const&, uint32_t const )

    /**
     * Utility fn for incrementing pre-align emission counts..
     */
    template <typename GlobalCountsType, typename SequenceResidueType>
    static
    void
    incrementPreAlignEmissionCounts (
      GlobalCountsType & global_counts,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      global_counts[ Emission::PreAlignInsertion ][ sequence[ seq_pos_i ] ] += 1;
    } // static incrementPreAlignEmissionCounts ( GlobalCountsType &, Sequence<SequenceResidueType> const&, uint32_t const )

    /**
     * Utility fn for incrementing post-align emission counts.
     */
    template <typename GlobalCountsType, typename SequenceResidueType>
    static
    void
    incrementPostAlignEmissionCounts (
      GlobalCountsType & global_counts,
      Sequence<SequenceResidueType> const& sequence,
      uint32_t const seq_pos_i
    )
    {
      global_counts[ Emission::PostAlignInsertion ][ sequence[ seq_pos_i ] ] += 1;
    } // static incrementPostAlignEmissionCounts ( GlobalCountsType &, Sequence<SequenceResidueType> const&, uint32_t const )
    ///// End static utility fns

    /**
     * An entente stores the sum of relative path likelihoods for paths
     * involving each parameter of the profile.  Path likelihoods are made
     * relative by dividing by the sequence score.  Ententes store sums over
     * sequences of relative path likelihoods.  For parameters that are used
     * multiple times on a path, the path will be counted multiple times.  So,
     * for instance, if a path has 4 C->C transitions, then the
     * postAlign->postAlign position of the entente will include 4 times the
     * relative likelihood of that path.
     *
     * The GlobalEntente is the Entente for the global parameters.
     */
    typedef ScalableParameterCollection<ScoreType, GlobalParameters<ResidueType, MatrixValueType> > GlobalEntente;

    /**
     * An entente stores the sum of relative path likelihoods for paths
     * involving each parameter of the profile.  Path likelihoods are made
     * relative by dividing by the sequence score.  Ententes store sums over
     * sequences of relative path likelihoods.  For parameters that are used
     * multiple times on a path, the path will be counted multiple times.  So,
     * for instance, if a path has 4 C->C transitions, then the
     * postAlign->postAlign position of the entente will include 4 times the
     * relative likelihood of that path.
     *
     * The PositionEntente is the Entente for the parameters of one position.
     */
    typedef ScalableParameterCollection<ScoreType, PositionSpecificParameters<ResidueType, MatrixValueType> > PositionEntente;

    // \todo  PUT BACK.. TESTING.
    typedef ScalableParameterCollection<ScoreType, AlignmentProfilePositionParameters<ResidueType, MatrixValueType> > AlignmentProfilePosition;
    // \todo  REMOVE.  TESTING.
    //typedef ScalableParameterCollection<ScoreType, AlignmentProfilePositionParameters<ResidueType, doublerealspace> > AlignmentProfilePosition;

    /**
     * An AlignmentProfile stores the total probability of all paths that flow
     * through each state (for a particular profile, sequence pair).
     */
    class AlignmentProfile :
      public vector<AlignmentProfilePosition>
    {
      // Boost serialization
    private:
      typedef PositionSpecificParameters<ResidueType, MatrixValueType> PositionSpecificParameters_MatrixValueType;
      typedef GlobalParameters<ResidueType, MatrixValueType> GlobalParameters_MatrixValueType;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        uint32_t length = this->size();
        ar & BOOST_SERIALIZATION_NVP( length );
        if( length != this->size() ) {
          this->resize( length );
        }
      
        // Serialize positions
        if( length > 0 ) {
         uint32_t last_pos = length - 1;
         uint32_t pos_i;
         string tmp_str;
         for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
           tmp_str = "pos_" + boost::lexical_cast<std::string>(pos_i);
           ar & boost::serialization::make_nvp(tmp_str.c_str(), ( this )[ pos_i ]);
         }
        } // End if length > 0, serialize positions
      } // serialize( Archive &, const unsigned int )

    public:
      AlignmentProfile ();

      AlignmentProfile (
        uint32_t length
      );

      void
      reinitialize (
        uint32_t length
      );

      uint32_t
      length () const;

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.  Also sets each
       * contained AlignmentProfilePosition's m_scalar to 1.
       */
      template <typename other_type>
      void
      normalize ( other_type const & min );

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.  Also sets each
       * contained AlignmentProfilePosition's m_scalar to 1.
       */
      void
      normalize ( MatrixValueType const& min );

      /**
       * Set all values to 0, and the scalars to 1.
       */
      void
      zero ();

      AlignmentProfile &
      operator= ( AlignmentProfile const& other_profile );

      /**
       * Adjust each distribution's values such that the largest value is 1.0,
       * utilizing the m_scalar values.
       */
      void
      rescale ();

      /**
       * Adjust each distribution's values such that the largest value is the
       * given maximum, utilizing the m_scalar values.
       */
      void
      rescale ( MatrixValueType const & max );

      /**
       * Divide each value by m_scalar and set the m_scalars to 1.0.
       */
      void
      unscale ();

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        AlignmentProfile const& alignment_profile
      )
      {
        if( alignment_profile.length() == 0 ) {
          return os;
        }
        uint32_t last_pos = alignment_profile.length() - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          os << alignment_profile[ pos_i ] << endl;
        }

        return os;
      } // friend operator<< ( basic_ostream &, AlignmentProfilePosition const& )

      /**
       * Divide each contained distrubution value by denominator.
       */
      template <typename AnyValueType>
      AlignmentProfile &
      operator/= ( AnyValueType const& denominator );
  
      /**
       * Add to each contained distribution the values in the given other
       * AlignmentProfile.  Note that this may violate the rule that
       * the probabilities sum to 1.
       */
      AlignmentProfile &
      operator+= ( AlignmentProfile const& other_profile );
  
      /**
       * Set all values such that each distrubution is randomly distributed.
       */
      void
      uniform ( Random & random );
  
      /**
       * Return the largest value in the contained distributions.
       */
      MatrixValueType
      maximumValue () const;
  
      /**
       * Calculate and return the Euclidean distance between this alignment
       * profile and another one (treating every probability as an orthogonal
       * dimension).
       */
      double
      euclideanDistance (
        AlignmentProfile const& other_profile
      ) const;

      /**
       * Calculate and return the square of the Euclidean distance between this
       * alignment profile and another one (treating every probability as an
       * orthogonal dimension).
       */
      double
      euclideanDistanceSquared (
        AlignmentProfile const& other_profile
      ) const;

      /**
       * Calculate and return the product of the probabilities that the Match
       * states of this alignment profile are identical to the Match states of
       * the given alignment profile, using the Match state emission values of
       * the internal positions of the two alignment profiles as the relevant
       * probabilities.  Note that these are (usually) not conditioned on there
       * being a Match at all (that is, the probs don't sum to 1; instead they
       * sum to 1-P(Deletion)).   If the include_nonemission_prob argument is
       * true, the probability of a cooccurring deletion will be included in the
       * calculation.  Note also that this is symmetric, so calling
       * the other profile's calculateMatchEmissionCooccurrenceProbability
       * with this profile as an argument will return the same thing.
       * NOTE: We assume the profiles are the same length!
       */
      MatrixValueType
      calculateMatchEmissionCooccurrenceProbability (
        AlignmentProfile const & other_profile,
        const bool include_nonemission_prob = false
      ) const;

      /**
       * Calculate and return the expected number of times that the Match state
       * of this alignment profile position is identical to the Match state of
       * the given alignment profile position, using the Match state emission
       * values of the two alignment profile positions as the relevant
       * probabilities.  Note that these are (usually) not conditioned on there
       * being a Match at all (that is, the probs don't sum to 1; instead they
       * sum to 1-P(Deletion)).  If the include_nonemission_counts argument is
       * true, the probability of a cooccurring deletion will be included in the
       * calculation.  Note also that this is symmetric, so calling
       * the other position's calculateMatchEmissionCooccurrenceExpectedCount
       * with this pos as an argument will return the same thing.
       */
      MatrixValueType
      calculateMatchEmissionCooccurrenceExpectedCount (
        AlignmentProfile const & other_profile,
        const bool include_nonemission_counts = false
      ) const;

      /**
       * Calculate and return the cross entropy
       * E(-log(other_node)).  Note that the cross entropy is
       * non-symmetric (calling other_node.crossEntropy( *this ) will
       * return a different value).
       *
       * This can be used to calculate the (self-)entropy by
       * calling this->crossEntropy( *this ).  It can also be used to
       * calculate the KL divergence by taking the difference of the
       * cross-entropy and the self-entropy, or the symmeterized KL divergence by
       * summing the KL divergences computed both ways.
       */
      double
      crossEntropy (
        AlignmentProfile const& other_node
      ) const;

      /**
       * Calculate and return the (possibly weighted) cross entropy
       * E(-log(other_node)).  The weights (if non-null) may be of any type
       * convertible to a double, and the cross entropy will be computed as
       * E(-log(weights*other_node)).  Note that the cross entropy is
       * non-symmetric (calling other_node.crossEntropy( *this, weights ) will
       * return a different value).
       *
       * This can be used to calculate the (possibly weighted) (self-)entropy by
       * calling this->crossEntropy( *this, weights ).  It can also be used to
       * calculate the KL divergence by taking the difference of the
       * cross-entropy and the self-entropy, or the symmeterized KL divergence by
       * summing the KL divergences computed both ways.  Note, though, that
       * weights should all be the same within any Multinomial distribution if
       * this is to be used to calculate weighted entropies or KL divergences.
       * Otherwise the usual properties of these metrics will be violated.
       */
      double
      crossEntropy (
        AlignmentProfile const& other_node,
        AlignmentProfile const * const weights
      ) const;

    }; // End inner class DynamicProgramming::AlignmentProfile

    class DirichletStuff
    {
    public:

      /**
       * Calculates ln( P( observed_counts | dirichlet ) ), the log probability
       * of a count vector given a Dirichlet distribution with the given
       * Dirichlet parameters ("prior_pseudocounts").  Pass the m_scalar from
       * the entente as the count_scalar argument if the entente values have
       * been scaled.  Adapted from Logp_cvec() in hmmer's mathsupport.c, which
       * in turn was "adapted from an implementation by Graeme Mitchison".
       *
       * Both arrays must be of size at least vectors_size.
       */
      template <typename CountType, typename DirichletAlphaType>
      static double
      calculateLogProbability (
        CountType const * observed_counts,
        ScoreType const count_scalar,
        DirichletAlphaType const * prior_pseudocounts,
        uint32_t const vectors_size
      )
      {
        double lnp = 0.0;//log( 1.0 );// log likelihood of P(observed_counts | Dirichlet)
        ProbabilityType sum1, sum2, sum3, tmp;
        sum1 = sum2 = sum3 = 0.0;
        ProbabilityType unscaled_observed_count;
        for( uint32_t i = 0; i < vectors_size; i++ ) {
          unscaled_observed_count =
            observed_counts[ i ];
          unscaled_observed_count /=
            count_scalar;
          sum1 += unscaled_observed_count;
          sum1 += prior_pseudocounts[ i ];
          sum2 += prior_pseudocounts[ i ];
          sum3 += unscaled_observed_count;
          tmp = prior_pseudocounts[ i ];
          tmp += unscaled_observed_count;
          lnp  += gammln( tmp );
          tmp = 1;
          tmp += unscaled_observed_count;
          lnp  -= gammln( tmp );
          lnp  -= gammln( prior_pseudocounts[ i ] );
        }
        lnp -= gammln( sum1 );
        lnp += gammln( sum2 );
        tmp = sum3;
        tmp += 1;
        lnp += gammln( tmp );
        return lnp;
      } // static calculateLogProbability( CountType[] const, CountType, DirichletAlphaType[] const, uint32_t )

      /**
       * Return the natural log of the gamma function of x, where x is > 0.0.
       *
       * Adapted from Gammln() in hmmer's squid's sre_math.c, which in turn was
       * "Adapted from a public domain implementation in the NCBI core math
       * library. Thanks to John Spouge and the NCBI."
       */
      static double
      gammln ( ProbabilityType const & x )
      {
        // Sean Eddy's code had this:
        //
        //   Protect against x=0. We see this in Dirichlet code,
        //   for terms alpha = 0. This is a severe hack but it is effective
        //   and (we think?) safe. (due to GJM)
        // if (x <= 0.0) return 999999.; 
        assert( x > 0 ); // Paul's replacement for now.

        static double cof[ 11 ] = {
          4.694580336184385e+04,
          -1.560605207784446e+05,
          2.065049568014106e+05,
          -1.388934775095388e+05,
          5.031796415085709e+04,
          -9.601592329182778e+03,
          8.785855930895250e+02,
          -3.155153906098611e+01,
          2.908143421162229e-01,
          -2.319827630494973e-04,
          1.251639670050933e-10
        };
      
        double xx = toDouble( x ) - 1.0;
        double tx, tmp;
        tx = tmp = xx + 11.0;
        double value = 1.0;
        for( int i = 10; i >= 0; i-- ) { // sum least significant terms first
          value += cof[ i ] / tmp;
          tmp   -= 1.0;
        }
        value  = log( value );
        tx    += 0.5;
        value += 0.918938533 + ( xx + 0.5 ) * log( tx ) - tx;
        return value;
      } // static gammln( double )

      /**
       * Normalize a vector of log-probabilities (as doubles), putting the
       * result in the given ProbabilityType vector. Be careful of overflowing
       * exp().  Implementation adapted from hmmer's LogNorm() in
       * mathsupport.c, which was in turn adapted from Graeme Mitchison.
       */
      static void
      exponentiateAndNormalize (
        vector<double> const & lp_vec,
        vector<ProbabilityType> & prob_vec
      )
      {
        uint32_t n = prob_vec.size();
        assert( n == lp_vec.size() );

        double max    = log( 0.0 ); // - numeric_limits<double>::infinity();
        static double const fifty = 50.0;

        double denom = 0.0;
        uint32_t i;
        for( i = 0; i < n; i++) {
          if( lp_vec[ i ] > max ) {
            max = lp_vec[ i ];
          }
        }
        for( i = 0; i < n; i++) {
          if( lp_vec[ i ] > ( max - fifty ) ) {
            denom += lp_vec[ i ] - max;
          }
        }
        for( i = 0; i < n; i++ ) {
          if( lp_vec[ i ] > ( max - fifty ) ) {
            prob_vec[ i ] = exp( ( lp_vec[ i ] - max ) - denom );
          } else {
            prob_vec[ i ] = 0.0; // \todo  Use parameter for profile minimumValue?
          }
        }
      } // static exponentiateAndNormalize( vector<double> const&, vector<ProbabilityType> & )

    }; // End inner class DynamicProgramming::DirichletStuff

    template <typename DirichletParameterType>
    class DirichletMixtureMatchEmissionPrior :
      public DirichletStuff,
      public vector<MatchEmissionParameters<ResidueType, DirichletParameterType> >
    {
    public:
      // \todo  Protect access, ensure it sums to 1 (and is non-negative), etc.
      vector<ProbabilityType> m_mixingProbs;
    
      DirichletMixtureMatchEmissionPrior () :
        m_mixingProbs()
      {
        // Do nothing else.
      }; // <init>()
    
      DirichletMixtureMatchEmissionPrior ( DirichletMixtureMatchEmissionPrior const & copy_from )
      {
        *this = copy_from;
      }; // <init>( DirichletMixtureMatchEmissionPrior const & )
    
      DirichletMixtureMatchEmissionPrior ( uint32_t length ) :
        vector<MatchEmissionParameters<ResidueType, DirichletParameterType> >( length ),
        m_mixingProbs( length )
      {
        // Do nothing else.
        reinitialize( length );
      } // <init>( uint32_t )
    
      void
      reinitialize ()
      {
        reinitialize( 0 );
      }

      /**
       * Reset all to their defaults, and mixing probs to 1/length.
       */
      void
      reinitialize (
        uint32_t length
      )
      {
        if( this->size() != length ) {
          this->resize( length );
        }
        if( m_mixingProbs.size() != length ) {
          m_mixingProbs.resize( length );
        }
    
        if( length > 0 ) {
          uint32_t last_pos = length - 1;
          uint32_t pos_i;
          for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
            this->vector<MatchEmissionParameters<ResidueType, DirichletParameterType> >::operator[]( pos_i ).reinitialize();
            m_mixingProbs[ pos_i ] = ( 1.0 / length );
          }
        }
      } // reinitialize( uint32_t )
    
      /**
       * Reinitialize with size 1 and set all counts to 1.0.
       */
      void
      reinitializeToLaplace ()
      {
        reinitializeToEven( 1.0 );
      } // reinitializeToLaplace()
    
      /**
       * Reinitialize with size 1 and set all counts to count.
       */
      void
      reinitializeToEven (
        DirichletParameterType const count
      )
      {
        reinitialize( 1 );
        ( *this )[ 0 ][ Emission::Match ] = count;
      } // reinitializeToEven( DirichletParameterType )
    
      DirichletMixtureMatchEmissionPrior &
      operator= ( DirichletMixtureMatchEmissionPrior const& other_prior )
      {
        if( vector<MatchEmissionParameters<ResidueType, DirichletParameterType> >::size() != other_prior.size() ) {
          vector<MatchEmissionParameters<ResidueType, DirichletParameterType> >::resize( other_prior.size() );
          m_mixingProbs.resize( other_prior.size() );
        }

        uint32_t last_component = this->size() - 1;
        uint32_t component_i;
        for( component_i = 0; component_i <= last_component ; component_i++ ) {
          this->operator[]( component_i ) = other_prior[ component_i ];
          m_mixingProbs[ component_i ] = other_prior.m_mixingProbs[ component_i ];
        }
        return *this;
      } // operator=( DirichletMixtureMatchEmissionPrior const& )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        DirichletMixtureMatchEmissionPrior<DirichletParameterType> const& prior
      )
      {
        uint32_t last_component = prior.size() - 1;
        uint32_t component_i;
        for( component_i = 0; component_i <= last_component ; component_i++ ) {
          if( last_component != 0 ) {
            os << "Component " << ( component_i - 1 ) << ", probability " << prior.m_mixingProbs[ component_i ] << ": " << endl;
          }
          os << prior[ component_i ] << endl;
        }
        return os;
      } // friend operator<< ( basic_ostream, DirichletMixtureMatchEmissionPrior const & )

      /**
       * Calculate the posterior probabilities of the mixing components, given
       * the observed expected counts for one sequence (at one position).
       */
      template <typename MatchEmissionParametersType>
      void
      calculateMixingPosteriors (
        MatchEmissionParametersType const & entente,
        vector<ProbabilityType> & mixing_posteriors
      ) const
      {
        uint32_t num_mixing_components = this->size();
        if( mixing_posteriors.size() != num_mixing_components ) {
          mixing_posteriors.resize( num_mixing_components );
        }
    
        uint32_t mixing_component_i;
        if( num_mixing_components == 1 ) {
          mixing_posteriors[ 0 ] = 1.0;
        } else { // if num_mixing_components == 1 .. else ..
          // Do it in log space, so for now mixing_posteriors is a vector of
          // logs.  We will convert it back later.
          vector<double> log_mixing_posteriors( num_mixing_components );
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            if( m_mixingProbs[ mixing_component_i ] == 0 ) {
              log_mixing_posteriors[ mixing_component_i ] = log( 0.0 ); // - numeric_limits<double>::infinity();
              continue; // Don't try to get posterior prob if prior prob is 0.
            }
            log_mixing_posteriors[ mixing_component_i ] =
              log( m_mixingProbs[ mixing_component_i ] );
            log_mixing_posteriors[ mixing_component_i ] +=
              calculateLogPosterior( entente, mixing_component_i );
            // \todo  REMOVE
            //cout << "calculateMixingPosteriors: log mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          } // End foreach mixing_component_i
          // Now exponentiate and normalize
          DynamicProgramming::DirichletStuff::exponentiateAndNormalize(
            log_mixing_posteriors,
            mixing_posteriors
          );
          // \todo  REMOVE
          //for( mixing_component_i = 0;
          //     mixing_component_i < num_mixing_components;
          //     mixing_component_i++ ) {
          //  cout << "calculateMixingPosteriors: mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          //} // End foreach mixing_component_i
        } // End if num_mixing_components == 1 .. else ..
    
      } // calculateMixingPosteriors( MatchEmissionParametersType const &, vector<ProbabilityType> & )
    
      /**
       * Add prior pseudocounts to an observed emission count vector (an
       * Entente).  Modified from hmmer's P7PriorifyEmissionVector in prior.c.
       */
      template <typename MatchEmissionParametersType>
      void
      incorporatePrior (
        MatchEmissionParametersType & entente
      ) const
      {
        // Note that as a side effect, the entente will be normalized (if/when
        // unscaled), unless there is only 1 mixing component.
  
        uint32_t num_mixing_components = this->size();
        if( num_mixing_components == 1 ) {
          // Then we don't have to mix, we just have to add in the prior
          // pseudocounts.  This is faster than using the mixing code (since we
          // can perform the operations in type Probabilty, which is not a
          // log type).
          MatchEmissionParametersType prior_counts( entente );
          prior_counts.operator=( ( *this )[ 0 ] );
          // \todo  REMOVE
          //cout << "In incorporatePrior(..): prior_counts before scaling is " << prior_counts << endl;
          if( has_m_scalar<MatchEmissionParametersType>() ) {
            getScalar<ScoreType, MatchEmissionParametersType>( prior_counts ) /=
              getScalar<ScoreType, MatchEmissionParametersType>( entente );
          }
          //cout << "In incorporatePrior(..): prior_counts after scaling is " << prior_counts << endl;
          //cout << "In incorporatePrior(..): unscaled entente before is " << entente << endl;
          entente += prior_counts;
          //cout << "In incorporatePrior(..): unscaled entente after is " << entente << endl;
    
          return;
        } // End if num_mixing_components == 1
    
        // \todo  We could further speed this up by doing everything in
        // ProbabilityType type, and instead of dividing the entente by the scalar,
        // multiply the prior by the scalar.  If it ends up being 0, then just
        // return.  Also, why *is* the entente scalar in scoretype anyway?
    
        // We first need posterior probs P( q | entente ) for each mixture
        // component q.
        vector<ProbabilityType> mixing_posteriors( this->size() );
        calculateMixingPosteriors( entente, mixing_posteriors );
    
        // Incorporate the priors into the entente counts.
        // Sean Eddy's code says "following Sjolander (1996)"
        MatrixValueType xi, tmp_mvt;
        ScoreType unscaled_entente_value;
        ScoreType totc, tot, tmp_st;
        uint32_t mixing_component_i;
        uint32_t alphabet_size;
        uint32_t i;

        ScoreType entente_scalar;
        if( has_m_scalar<MatchEmissionParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, MatchEmissionParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // \todo  Code is copied for each type.  Could separate it out into a
        // function of the type tag.
        alphabet_size =
          entente[ Emission::Match ].size();
        // \todo  REMOVE
        //cout << "entente[ Match ].total() returns " << entente[ Emission::Match ].total() << endl;
        totc =
          entente[ Emission::Match ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < alphabet_size; i++ ) {
          unscaled_entente_value =
            entente[ Emission::Match ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Emission::Match ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Emission::Match ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Emission::Match ].total();
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Emission::Match ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
                
            xi += tmp_mvt;

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Emission::Match ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Emission::Match ][ i ] = tmp_st;

          // \todo  REMOVE
          //cout << "[" << i << "] new entente value is " << entente[ Emission::Match ][ i ] << endl;
        } // End foreach residue i
      } // incorporatePrior( MatchEmissionParametersType & ) const
    
      /**
       * Calculate and return ln( P( entente | dirichlet ) ) for component
       * mixing_component_i.  Sums up over the different alphabet types.
       */
      // \todo  Is adding the logposteriors the right thing to do?
      template <typename MatchEmissionParametersType>
      double
      calculateLogPosterior (
        MatchEmissionParametersType const & entente,
        uint32_t mixing_component_i
      ) const
      {
        ScoreType entente_scalar;
        if( has_m_scalar<MatchEmissionParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, MatchEmissionParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // Note that we directly access the m_probs members of the
        // MultinomialDistributions, which are public.
        return calculateLogProbability( entente[ Emission::Match ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Emission::Match ].m_probs, ( *this )[ mixing_component_i ][ Emission::Match ].size() );
      } // calculateLogPosterior( MatchEmissionParametersType const &, uint32_t ) const
    
    }; // End inner class DynamicProgramming::DirichletMixtureMatchEmissionPrior

    // Note that I could have made this subclass all of the (below) sub-priors,
    // but that would have derived from many vectors, and I'm not feeling so
    // brave at the moment.  Unfortunately this results in some code
    // duplication.
    template <typename DirichletParameterType>
    class DirichletMixtureGlobalPrior :
      public DirichletStuff,
      public vector<GlobalParameters<ResidueType, DirichletParameterType> >
    {
    public:
      // \todo  Protect access, ensure it sums to 1 (and is non-negative), etc.
      vector<ProbabilityType> m_mixingProbs;
    
      DirichletMixtureGlobalPrior () :
        m_mixingProbs()
      {
        // Do nothing else.
      }; // <init>()
    
      DirichletMixtureGlobalPrior ( uint32_t length ) :
        vector<GlobalParameters<ResidueType, DirichletParameterType> >( length ),
        m_mixingProbs( length )
      {
        reinitialize( length );
      } // <init>( uint32_t )
    
      DirichletMixtureGlobalPrior ( DirichletMixtureGlobalPrior const & copy_from )
      {
        *this = copy_from;
      }; // <init>( DirichletMixtureGlobalPrior const & )
    
      void
      reinitialize ()
      {
        reinitialize( 0 );
      } // reinitialize()

      /**
       * Reset all to their defaults, and mixing probs to 1/length.
       */
      void
      reinitialize (
        uint32_t length
      )
      {
        if( this->size() != length ) {
          this->resize( length );
        }
        if( m_mixingProbs.size() != length ) {
          m_mixingProbs.resize( length );
        }
    
        if( length > 0 ) {
          uint32_t last_pos = length - 1;
          uint32_t pos_i;
          for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
            this->vector<GlobalParameters<ResidueType, DirichletParameterType> >::operator[]( pos_i ).reinitialize();
            m_mixingProbs[ pos_i ] = 1.0;
            m_mixingProbs[ pos_i ] /= length;
          }
        }
      } // reinitialize( uint32_t )
    
      /**
       * Reinitialize with size 1 and set all counts to 1.0.
       */
      void
      reinitializeToLaplace ()
      {
        reinitializeToEven( 1.0 );
      } // reinitializeToLaplace()
    
      /**
       * Reinitialize with size 1 and set all counts to count.
       */
      void
      reinitializeToEven (
        DirichletParameterType const count
      )
      {
        reinitialize( 1 );
        ( *this )[ 0 ][ Emission::Insertion ] = count;
        ( *this )[ 0 ][ Transition::fromMatch ] = count;
        ( *this )[ 0 ][ Transition::fromInsertion ] = count;
        ( *this )[ 0 ][ Transition::fromDeletion ] = count;
        ( *this )[ 0 ][ Transition::fromPreAlign ] = count;
        ( *this )[ 0 ][ Transition::fromBegin ] = count;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ( *this )[ 0 ][ Transition::fromDeletionIn ] = count;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        ( *this )[ 0 ][ Emission::PreAlignInsertion ] = count;
#ifdef USE_END_DISTRIBUTION
        // For now we we never allow loops, so the End distribution is constant.
        // We don't use End.  Make sure its counts are all 0
        //( *this )[ 0 ][ Transition::fromEnd ] = count;
        ( *this )[ 0 ][ Transition::fromEnd ] = 0;
#endif // USE_END_DISTRIBUTION
        ( *this )[ 0 ][ Transition::fromPostAlign ] = count;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ( *this )[ 0 ][ Transition::fromDeletionOut ] = count;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        ( *this )[ 0 ][ Emission::PostAlignInsertion ] = count;
      } // reinitializeToEven( DirichletParameterType )
    
      DirichletMixtureGlobalPrior &
      operator= ( DirichletMixtureGlobalPrior const& other_prior )
      {
        if( vector<GlobalParameters<ResidueType, DirichletParameterType> >::size() != other_prior.size() ) {
          vector<GlobalParameters<ResidueType, DirichletParameterType> >::resize( other_prior.size() );
          m_mixingProbs.resize( other_prior.size() );
        }

        uint32_t last_component = this->size() - 1;
        uint32_t component_i;
        for( component_i = 0; component_i <= last_component ; component_i++ ) {
          this->operator[]( component_i ) = other_prior[ component_i ];
          m_mixingProbs[ component_i ] = other_prior.m_mixingProbs[ component_i ];
        }
        return *this;
      } // operator=( DirichletMixtureGlobalPrior const& )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        DirichletMixtureGlobalPrior<DirichletParameterType> const& prior
      )
      {
        uint32_t last_component = prior.size() - 1;
        uint32_t component_i;
        for( component_i = 0; component_i <= last_component ; component_i++ ) {
          if( last_component != 0 ) {
            os << "Component " << ( component_i - 1 ) << ", probability " << prior.m_mixingProbs[ component_i ] << ": " << endl;
          }
          os << prior[ component_i ] << endl;
        }
        return os;
      } // friend operator<< ( basic_ostream, DirichletMixtureGlobalPrior const & )

      /**
       * Calculate the posterior probabilities of the mixing components, given
       * the observed expected counts for one sequence (at one position).
       */
      template <typename GlobalParametersType>
      void
      calculateMixingPosteriors (
        GlobalParametersType const & entente,
        vector<ProbabilityType> & mixing_posteriors
      ) const
      {
        uint32_t num_mixing_components = this->size();
        if( mixing_posteriors.size() != num_mixing_components ) {
          mixing_posteriors.resize( num_mixing_components );
        }
    
        uint32_t mixing_component_i;
        if( num_mixing_components == 1 ) {
          mixing_posteriors[ 0 ] = 1.0;
        } else { // if num_mixing_components == 1 .. else ..
          // Do it in log space, so for now mixing_posteriors is a vector of
          // logs.  We will convert it back later.
          vector<double> log_mixing_posteriors( num_mixing_components );
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            if( m_mixingProbs[ mixing_component_i ] == 0 ) {
              log_mixing_posteriors[ mixing_component_i ] = log( 0.0 ); // - numeric_limits<double>::infinity();
              continue; // Don't try to get posterior prob if prior prob is 0.
            }
            log_mixing_posteriors[ mixing_component_i ] =
              log( m_mixingProbs[ mixing_component_i ] );
            log_mixing_posteriors[ mixing_component_i ] +=
              calculateLogPosterior( entente, mixing_component_i );
            // \todo  REMOVE
            //cout << "calculateMixingPosteriors: log mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          } // End foreach mixing_component_i
          // Now exponentiate and normalize
          DynamicProgramming::DirichletStuff::exponentiateAndNormalize(
            log_mixing_posteriors,
            mixing_posteriors
          );
          // \todo  REMOVE
          //for( mixing_component_i = 0;
          //     mixing_component_i < num_mixing_components;
          //     mixing_component_i++ ) {
          //  cout << "calculateMixingPosteriors: mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          //} // End foreach mixing_component_i
        } // End if num_mixing_components == 1 .. else ..
    
      } // calculateMixingPosteriors( GlobalParametersType const &, vector<ProbabilityType> & )
    
      /**
       * Add prior pseudocounts to an observed emission count vector (an
       * Entente).  Modified from hmmer's P7PriorifyEmissionVector in prior.c.
       */
      template <typename GlobalParametersType>
      void
      incorporatePrior (
        GlobalParametersType & entente
      ) const
      {
        // Note that as a side effect, the entente will be normalized (if/when
        // unscaled), unless there is only 1 mixing component.
  
        uint32_t num_mixing_components = this->size();
        if( num_mixing_components == 1 ) {
          // Then we don't have to mix, we just have to add in the prior
          // pseudocounts.  This is faster than using the mixing code (since we
          // can perform the operations in type Probabilty, which is not a
          // log type).
          GlobalParametersType prior_counts( entente );
          //ScalableParameterCollection<ScoreType, GlobalParameters<ResidueType, DirichletParameterType> prior_counts =
          //  ScalableParameterCollection<ScoreType, GlobalParameters<ResidueType, DirichletParameterType> >();
          prior_counts.operator=( ( *this )[ 0 ] );
          // \todo  REMOVE
          //cout << "In incorporatePrior(..): prior_counts before scaling is " << prior_counts << endl;
          if( has_m_scalar<GlobalParametersType>() ) {
             getScalar<ScoreType, GlobalParametersType>( prior_counts ) /=
              getScalar<ScoreType, GlobalParametersType>( entente );
          }
          //cout << "In incorporatePrior(..): prior_counts after scaling is " << prior_counts << endl;
          //cout << "In incorporatePrior(..): unscaled entente before is " << entente << endl;
          entente += prior_counts;
          //cout << "In incorporatePrior(..): unscaled entente after is " << entente << endl;
    
          return;
        } // End if num_mixing_components == 1
    
        // \todo  We could further speed this up by doing everything in
        // ProbabilityType type, and instead of dividing the entente by the scalar,
        // multiply the prior by the scalar.  If it ends up being 0, then just
        // return.  Also, why *is* the entente scalar in scoretype anyway?
    
        // We first need posterior probs P( q | entente ) for each mixture
        // component q.
        vector<ProbabilityType> mixing_posteriors( this->size() );
        calculateMixingPosteriors( entente, mixing_posteriors );
    
        // Incorporate the priors into the entente counts.
        // Sean Eddy's code says "following Sjolander (1996)"
        MatrixValueType xi, tmp_mvt;
        ScoreType unscaled_entente_value;
        ScoreType totc, tot, tmp_st;
        uint32_t mixing_component_i;
        uint32_t alphabet_size;
        uint32_t i;

        ScoreType entente_scalar;
        if( has_m_scalar<GlobalParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, GlobalParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // \todo  Code is copied for each type.  Could separate it out into a
        // function of the type tag.
        alphabet_size =
          entente[ Emission::Insertion ].size();
        // \todo  REMOVE
        //cout << "entente[ Match ].total() returns " << entente[ Emission::Insertion ].total() << endl;
        totc =
          entente[ Emission::Insertion ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < alphabet_size; i++ ) {
          unscaled_entente_value =
            entente[ Emission::Insertion ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Insertion ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Emission::Insertion ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Emission::Insertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Emission::Insertion ].total();
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Emission::Insertion ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
                
            xi += tmp_mvt;

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Emission::Insertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Emission::Insertion ][ i ] = tmp_st;

          // \todo  REMOVE
          //cout << "[" << i << "] new entente value is " << entente[ Emission::Insertion ][ i ] << endl;
        } // End foreach residue i
    
        ////////
        // Match
        totc =
          entente[ Transition::fromMatch ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromMatch>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromMatch ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromMatch ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromMatch ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromMatch ].total();

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromMatch ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromMatch ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromMatch ][ i ] = tmp_st;
          // \todo  REMOVE
          //cout << "[" << i << "] new Match entente value is " << entente[ Transition::fromMatch ][ i ] << endl;
        } // End foreach Match transition target i

        ////////
        // Insertion transitions
        totc =
          entente[ Transition::fromInsertion ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromInsertion>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromInsertion ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromInsertion ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromInsertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromInsertion ].total();

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromInsertion ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromInsertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromInsertion ][ i ] = tmp_st;
          // \todo  REMOVE
          //cout << "[" << i << "] new Insertion entente value is " << entente[ Transition::fromInsertion ][ i ] << endl;
        } // End foreach Insertion transition target i

        ////////
        // Deletion
        totc =
          entente[ Transition::fromDeletion ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromDeletion>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromDeletion ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromDeletion ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromDeletion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromDeletion ].total();

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromDeletion ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromDeletion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromDeletion ][ i ] = tmp_st;
          // \todo  REMOVE
          //cout << "[" << i << "] new Deletion entente value is " << entente[ Transition::fromDeletion ][ i ] << endl;
        } // End foreach Deletion transition target i

        ////////
        // PreAlign
        totc =
          entente[ Transition::fromPreAlign ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromPreAlign>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromPreAlign ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].total();

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromPreAlign ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromPreAlign ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromPreAlign ][ i ] = tmp_st;
          // \todo  REMOVE
          //cout << "[" << i << "] new PreAlign entente value is " << entente[ Transition::fromPreAlign ][ i ] << endl;
        } // End foreach PreAlign transition target i

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        alphabet_size =
          entente[ Emission::PreAlignInsertion ].size();
        // \todo  REMOVE
        //cout << "entente[ Match ].total() returns " << entente[ Emission::PreAlignInsertion ].total() << endl;
        totc =
          entente[ Emission::PreAlignInsertion ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < alphabet_size; i++ ) {
          unscaled_entente_value =
            entente[ Emission::PreAlignInsertion ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::PreAlignInsertion ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].total();
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
                
            xi += tmp_mvt;
        
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Emission::PreAlignInsertion ][ i ] = tmp_st;
        
          // \todo  REMOVE
          //cout << "[" << i << "] new entente value is " << entente[ Emission::PreAlignInsertion ][ i ] << endl;
        } // End foreach residue i
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS


        ////////
        // Begin
        totc =
          entente[ Transition::fromBegin ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromBegin>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromBegin ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // \todo  REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // \todo  REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromBegin ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromBegin ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromBegin ].total();

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromBegin ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // \todo  REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromBegin ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromBegin ][ i ] = tmp_st;
          // \todo  REMOVE
          //cout << "[" << i << "] new Begin entente value is " << entente[ Transition::fromBegin ][ i ] << endl;
        } // End foreach Begin transition target i

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ////////
        // DeletionIn
        totc =
          entente[ Transition::fromDeletionIn ].total();
        // \todo  REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // \todo  REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromDeletionIn>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromDeletionIn ][ i ];
          // \todo  REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].total();

            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromDeletionIn ][ i ] = tmp_st;
          // TODO: REMOVE
          //cout << "[" << i << "] new DeletionIn entente value is " << entente[ Transition::fromDeletionIn ][ i ] << endl;
        } // End foreach DeletionIn transition target i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT


#ifdef USE_END_DISTRIBUTION
        // For now we we never allow loops, so the End distribution is constant.
        //////////
        //// End
        //totc =
        //  entente[ Transition::fromEnd ].total();
        //// TODO: REMOVE
        ////cout << "totc before scaling is " << totc << endl;
        //totc /= entente_scalar;
        //// TODO: REMOVE
        ////cout << "totc after scaling is " << totc << endl;
        ////cout << "as a double, totc is " << toDouble( totc ) << endl;
        //for( i = 0; i < seqan::ValueSize<TransitionFromEnd>::VALUE; i++ ) {
        //  unscaled_entente_value =
        //    entente[ Transition::fromEnd ][ i ];
        //  // TODO: REMOVE
        //  //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
        //  // TODO: REMOVE
        //  //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
        //  unscaled_entente_value /= entente_scalar;
        //  // TODO: REMOVE
        //  //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
        //  // TODO: REMOVE
        //  //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
        //  xi = 0.0;
        //  for( mixing_component_i = 0;
        //       mixing_component_i < num_mixing_components;
        //       mixing_component_i++ ) {
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromEnd ].total() << endl;
        //    //tot = ( *this )[ mixing_component_i ][ Transition::fromEnd ].total();
        //    //cout << "\t as a ScoreType, this is " << tot << endl;
        //    tot = totc;
        //    tot +=
        //      ( *this )[ mixing_component_i ][ Transition::fromEnd ].total();
        //
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
        //    tmp_mvt = unscaled_entente_value;
        //    tmp_mvt +=
        //      ( *this )[ mixing_component_i ][ Transition::fromEnd ][ i ];
        //    tmp_mvt /= tot;
        //    tmp_mvt *=
        //      mixing_posteriors[ mixing_component_i ];
        //    
        //    xi += tmp_mvt;
        //
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromEnd ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
        //    //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
        //  }
        //  tmp_st = xi;
        //  tmp_st *= entente_scalar;
        //  entente[ Transition::fromEnd ][ i ] = tmp_st;
        //  // TODO: REMOVE
        //  //cout << "[" << i << "] new End entente value is " << entente[ Transition::fromEnd ][ i ] << endl;
        //} // End foreach End transition target i
#endif // USE_END_DISTRIBUTION

        ////////
        // PostAlign
        totc =
          entente[ Transition::fromPostAlign ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromPostAlign>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromPostAlign ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].total();

            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromPostAlign ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromPostAlign ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromPostAlign ][ i ] = tmp_st;
          // TODO: REMOVE
          //cout << "[" << i << "] new PostAlign entente value is " << entente[ Transition::fromPostAlign ][ i ] << endl;
        } // End foreach PostAlign transition target i

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ////////
        // DeletionOut
        totc =
          entente[ Transition::fromDeletionOut ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromDeletionOut>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromDeletionOut ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].total();

            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
            
            xi += tmp_mvt;

            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Transition::fromDeletionOut ][ i ] = tmp_st;
          // TODO: REMOVE
          //cout << "[" << i << "] new DeletionOut entente value is " << entente[ Transition::fromDeletionOut ][ i ] << endl;
        } // End foreach DeletionOut transition target i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        alphabet_size =
          entente[ Emission::PostAlignInsertion ].size();
        // TODO: REMOVE
        //cout << "entente[ Match ].total() returns " << entente[ Emission::PostAlignInsertion ].total() << endl;
        totc =
          entente[ Emission::PostAlignInsertion ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < alphabet_size; i++ ) {
          unscaled_entente_value =
            entente[ Emission::PostAlignInsertion ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::PostAlignInsertion ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].total() << endl;
            //tot = ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot = totc;
            tot +=
              ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            tmp_mvt = unscaled_entente_value;
            tmp_mvt +=
              ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ][ i ];
            tmp_mvt /= tot;
            tmp_mvt *=
              mixing_posteriors[ mixing_component_i ];
                
            xi += tmp_mvt;
        
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          tmp_st = xi;
          tmp_st *= entente_scalar;
          entente[ Emission::PostAlignInsertion ][ i ] = tmp_st;
        
          // TODO: REMOVE
          //cout << "[" << i << "] new entente value is " << entente[ Emission::PostAlignInsertion ][ i ] << endl;
        } // End foreach residue i
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      } // incorporatePrior( GlobalParametersType & ) const
    
      /**
       * Calculate and return ln( P( entente | dirichlet ) ) for component
       * mixing_component_i.  Sums up over the different alphabet types.
       */
      // TODO: Is adding the logposteriors the right thing to do?
      template <typename GlobalParametersType>
      double
      calculateLogPosterior (
        GlobalParametersType const & entente,
        uint32_t mixing_component_i
      ) const
      {
        ScoreType entente_scalar;
        if( has_m_scalar<GlobalParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, GlobalParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // Note that we directly access the m_probs members of the
        // MultinomialDistributions, which are public.
        double prod = 0.0; // log( 1.0 );
        prod += calculateLogProbability( entente[ Emission::Insertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Emission::Insertion ].m_probs, ( *this )[ mixing_component_i ][ Emission::Insertion ].size() );
        prod += calculateLogProbability( entente[ Transition::fromMatch ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromMatch ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromMatch ].size() );
        prod += calculateLogProbability( entente[ Transition::fromInsertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromInsertion ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromInsertion ].size() );
        prod += calculateLogProbability( entente[ Transition::fromDeletion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromDeletion ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromDeletion ].size() );
        prod += calculateLogProbability( entente[ Transition::fromPreAlign ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].size() );
        prod += calculateLogProbability( entente[ Transition::fromBegin ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromBegin ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromBegin ].size() );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        prod += calculateLogProbability( entente[ Transition::fromDeletionIn ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].size() );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        prod += calculateLogProbability( entente[ Emission::PreAlignInsertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].m_probs, ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].size() );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
#ifdef USE_END_DISTRIBUTION
         // For now we we never allow loops, so the End distribution is constant.
         //prod += calculateLogProbability( entente[ Transition::fromEnd ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromEnd ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromEnd ].size() );
#endif // USE_END_DISTRIBUTION
        prod += calculateLogProbability( entente[ Transition::fromPostAlign ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].size() );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        prod += calculateLogProbability( entente[ Transition::fromDeletionOut ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].size() );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        prod += calculateLogProbability( entente[ Emission::PostAlignInsertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].m_probs, ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].size() );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

        return prod;
      } // calculateLogPosterior( GlobalParametersType const &, uint32_t ) const
    
    }; // End inner class DynamicProgramming::DirichletMixtureGlobalPrior

    template <typename DirichletParameterType>
    class DirichletMixtureInsertionEmissionPrior :
      public DirichletStuff,
      public vector<InsertionEmissionParameters<ResidueType, DirichletParameterType> >
    {
    public:
      // TODO: Protect access, ensure it sums to 1 (and is non-negative), etc.
      vector<ProbabilityType> m_mixingProbs;
    
      DirichletMixtureInsertionEmissionPrior () :
        m_mixingProbs()
      {
        // Do nothing else.
      }; // <init>()
    
      DirichletMixtureInsertionEmissionPrior (
        uint32_t length
      ) :
        vector<InsertionEmissionParameters<ResidueType, DirichletParameterType> >( length ),
        m_mixingProbs()
      {
        reinitialize( length );
      } // <init>( uint32_t )
    
      /**
       * Reset all to their defaults, and mixing probs to 1/length.
       */
      void
      reinitialize (
        uint32_t length
      )
      {
        if( this->size() != length ) {
          this->resize( length );
        }
        if( m_mixingProbs.size() != length ) {
          m_mixingProbs.resize( length );
        }
    
        if( length > 0 ) {
          uint32_t last_pos = length - 1;
          uint32_t pos_i;
          for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
            this->vector<InsertionEmissionParameters<ResidueType, DirichletParameterType> >::operator[]( pos_i ).reinitialize();
            m_mixingProbs[ pos_i ] = 1.0 / length;
          }
        }
      } // reinitialize( uint32_t )
    
      /**
       * Reinitialize with size 1 and set all counts to 1.0.
       */
      void
      reinitializeToLaplace ()
      {
        reinitializeToEven( 1.0 );
      } // reinitializeToLaplace()
    
      /**
       * Reinitialize with size 1 and set all counts to count.
       */
      void
      reinitializeToEven (
        DirichletParameterType const count
      )
      {
        reinitialize( 1 );
        ( *this )[ 0 ][ Emission::Insertion ] = count;
      } // reinitializeToEven( DirichletParameterType )
    
      /**
       * Calculate the posterior probabilities of the mixing components, given
       * the observed expected counts for one sequence (at one position).
       */
      template <typename InsertionEmissionParametersType>
      void
      calculateMixingPosteriors (
        InsertionEmissionParametersType const & entente,
        vector<ProbabilityType> & mixing_posteriors
      ) const
      {
        uint32_t num_mixing_components = this->size();
        if( mixing_posteriors.size() != num_mixing_components ) {
          mixing_posteriors.resize( num_mixing_components );
        }
    
        uint32_t mixing_component_i;
        if( num_mixing_components == 1 ) {
          mixing_posteriors[ 0 ] = 1.0;
        } else { // if num_mixing_components == 1 .. else ..
          // Do it in log space, so for now mixing_posteriors is a vector of
          // logs.  We will convert it back later.
          vector<double> log_mixing_posteriors( num_mixing_components );
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            if( m_mixingProbs[ mixing_component_i ] == 0 ) {
              log_mixing_posteriors[ mixing_component_i ] = log( 0.0 ); // - numeric_limits<double>::infinity();
              continue; // Don't try to get posterior prob if prior prob is 0.
            }
            log_mixing_posteriors[ mixing_component_i ] =
              (
                log( m_mixingProbs[ mixing_component_i ] ) +
                calculateLogPosterior( entente, mixing_component_i )
              );
            // TODO: REMOVE
            //cout << "calculateMixingPosteriors: log mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          } // End foreach mixing_component_i
          // Now exponentiate and normalize
          DynamicProgramming::DirichletStuff::exponentiateAndNormalize(
            log_mixing_posteriors,
            mixing_posteriors
          );
          // TODO: REMOVE
          //for( mixing_component_i = 0;
          //     mixing_component_i < num_mixing_components;
          //     mixing_component_i++ ) {
          //  cout << "calculateMixingPosteriors: mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          //} // End foreach mixing_component_i
        } // End if num_mixing_components == 1 .. else ..
    
      } // calculateMixingPosteriors( InsertionEmissionParametersType const &, vector<ProbabilityType> & )
    
      /**
       * Add prior pseudocounts to an observed emission count vector (an
       * Entente).  Modified from hmmer's P7PriorifyEmissionVector in prior.c.
       */
      template <typename InsertionEmissionParametersType>
      void
      incorporatePrior (
        InsertionEmissionParametersType & entente
      ) const
      {
        // Note that as a side effect, the entente will be normalized (if/when
        // unscaled), unless there is only 1 mixing component.
  
        uint32_t num_mixing_components = this->size();
        if( num_mixing_components == 1 ) {
          // Then we don't have to mix, we just have to add in the prior
          // pseudocounts.  This is faster than using the mixing code (since we
          // can perform the operations in type Probabilty, which is not a
          // log type).
          InsertionEmissionParametersType prior_counts( entente );
          //ScalableParameterCollection<ScoreType, InsertionEmissionParameters<ResidueType, DirichletParameterType> prior_counts =
          //  ScalableParameterCollection<ScoreType, InsertionEmissionParameters<ResidueType, DirichletParameterType> >();
          prior_counts.operator=( ( *this )[ 0 ] );
          // TODO: REMOVE
          //cout << "In incorporatePrior(..): prior_counts before scaling is " << prior_counts << endl;
          if( has_m_scalar<InsertionEmissionParametersType>() ) {
             getScalar<ScoreType, InsertionEmissionParametersType>( prior_counts ) /=
               getScalar<ScoreType, InsertionEmissionParametersType>( entente );
          }
          //cout << "In incorporatePrior(..): prior_counts after scaling is " << prior_counts << endl;
          //cout << "In incorporatePrior(..): unscaled entente before is " << entente << endl;
          entente += prior_counts;
          //cout << "In incorporatePrior(..): unscaled entente after is " << entente << endl;
    
          return;
        } // End if num_mixing_components == 1
    
        // TODO: We could further speed this up by doing everything in
        // ProbabilityType type, and instead of dividing the entente by the scalar,
        // multiply the prior by the scalar.  If it ends up being 0, then just
        // return.  Also, why *is* the entente scalar in scoretype anyway?
    
        // We first need posterior probs P( q | entente ) for each mixture
        // component q.
        vector<ProbabilityType> mixing_posteriors( this->size() );
        calculateMixingPosteriors( entente, mixing_posteriors );
    
        // Incorporate the priors into the entente counts.
        // Sean Eddy's code says "following Sjolander (1996)"
        MatrixValueType xi;
        ScoreType unscaled_entente_value;
        ScoreType totc, tot;
        uint32_t mixing_component_i;
        uint32_t alphabet_size;
        uint32_t i;

        ScoreType entente_scalar;
        if( has_m_scalar<InsertionEmissionParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, InsertionEmissionParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        alphabet_size =
          entente[ Emission::Insertion ].size();
        // TODO: REMOVE
        //cout << "entente[ Match ].total() returns " << entente[ Emission::Insertion ].total() << endl;
        totc =
          entente[ Emission::Insertion ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < alphabet_size; i++ ) {
          unscaled_entente_value =
            entente[ Emission::Insertion ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Insertion ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Emission::Insertion ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Emission::Insertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Emission::Insertion ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Emission::Insertion ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Emission::Insertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          entente[ Emission::Insertion ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new entente value is " << entente[ Emission::Insertion ][ i ] << endl;
        } // End foreach residue i
    
      } // incorporatePrior( InsertionEmissionParametersType & ) const
    
      /**
       * Calculate and return ln( P( entente | dirichlet ) ) for component
       * mixing_component_i.  Sums up over the different alphabet types.
       */
      // TODO: Is adding the logposteriors the right thing to do?
      template <typename InsertionEmissionParametersType>
      double
      calculateLogPosterior (
        InsertionEmissionParametersType const & entente,
        uint32_t mixing_component_i
      ) const
      {
        ScoreType entente_scalar;
        if( has_m_scalar<InsertionEmissionParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, InsertionEmissionParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // Note that we directly access the m_probs members of the
        // MultinomialDistributions, which are public.
        return calculateLogProbability( entente[ Emission::Insertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Emission::Insertion ].m_probs, ( *this )[ mixing_component_i ][ Emission::Insertion ].size() );
      } // calculateLogPosterior( InsertionEmissionParametersType const &, uint32_t ) const
    
    }; // End inner class DynamicProgramming::DirichletMixtureInsertionEmissionPrior

    template <typename DirichletParameterType>
    class DirichletMixturePositionTransitionPrior :
      public DirichletStuff,
      public vector<PositionTransitionParameters<ResidueType, DirichletParameterType> >
    {
    public:
      /// \todo  Protect access, ensure it sums to 1 (and is non-negative), etc.
      vector<ProbabilityType> m_mixingProbs;
    
      DirichletMixturePositionTransitionPrior () :
        m_mixingProbs()
      {
        // Do nothing else.
      }; // <init>()
    
      DirichletMixturePositionTransitionPrior ( uint32_t length ) :
        vector<PositionTransitionParameters<ResidueType, DirichletParameterType> >( length ),
        m_mixingProbs( length )
      {
        reinitialize( length );
      } // <init>( uint32_t )
    
      /**
       * Reset all to their defaults, and mixing probs to 1/length.
       */
      void
      reinitialize (
        uint32_t length
      )
      {
        if( this->size() != length ) {
          this->resize( length );
        }
        if( m_mixingProbs.size() != length ) {
          m_mixingProbs.resize( length );
        }
    
        if( length > 0 ) {
          uint32_t last_pos = length - 1;
          uint32_t pos_i;
          for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
            this->vector<PositionTransitionParameters<ResidueType, DirichletParameterType> >::operator[]( pos_i ).reinitialize();
            m_mixingProbs[ pos_i ] = 1.0 / length;
          }
        }
      } // reinitialize( uint32_t )
    
      /**
       * Reinitialize with size 1 and set all counts to 1.0.
       */
      void
      reinitializeToLaplace ()
      {
        reinitializeToEven( 1.0 );
      } // reinitializeToLaplace()
    
      /**
       * Reinitialize with size 1 and set all counts to count.
       */
      void
      reinitializeToEven (
        DirichletParameterType const count
      )
      {
        reinitialize( 1 );
        uint32_t i;
        for( i = 0; i < seqan::ValueSize<TransitionFromMatch>::VALUE; i++ ) {
          ( *this )[ 0 ][ Transition::fromMatch ][ i ] = count;
        } // End foreach transition i
        for( i = 0; i < seqan::ValueSize<TransitionFromInsertion>::VALUE; i++ ) {
          ( *this )[ 0 ][ Transition::fromInsertion ][ i ] = count;
        } // End foreach transition i
        for( i = 0; i < seqan::ValueSize<TransitionFromDeletion>::VALUE; i++ ) {
          ( *this )[ 0 ][ Transition::fromDeletion ][ i ] = count;
        } // End foreach transition i
      } // reinitializeToEven( DirichletParameterType )
    
      /**
       * Calculate the posterior probabilities of the mixing components, given
       * the observed expected counts for one sequence.
       */
      template <typename PositionTransitionParametersType>
      void
      calculateMixingPosteriors (
        PositionTransitionParametersType const & entente,
        vector<ProbabilityType> & mixing_posteriors
      ) const
      {
        uint32_t num_mixing_components = this->size();
        if( mixing_posteriors.size() != num_mixing_components ) {
          mixing_posteriors.resize( num_mixing_components );
        }
    
        uint32_t mixing_component_i;
        if( num_mixing_components == 1 ) {
          mixing_posteriors[ 0 ] = 1.0;
        } else { // if num_mixing_components == 1 .. else ..
          // Do it in log space, so for now mixing_posteriors is a vector of
          // logs.  We will convert it back later.
          vector<double> log_mixing_posteriors( num_mixing_components );
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            if( m_mixingProbs[ mixing_component_i ] == 0 ) {
              log_mixing_posteriors[ mixing_component_i ] = log( 0.0 ); // - numeric_limits<double>::infinity();
              continue; // Don't try to get posterior prob if prior prob is 0.
            }
            log_mixing_posteriors[ mixing_component_i ] =
              ( 
                log( m_mixingProbs[ mixing_component_i ] ) +
                calculateLogPosterior( entente, mixing_component_i )
              );
            // TODO: REMOVE
            //cout << "calculateMixingPosteriors: log mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          } // End foreach mixing_component_i
          // Now exponentiate and normalize
          DynamicProgramming::DirichletStuff::exponentiateAndNormalize(
            log_mixing_posteriors,
            mixing_posteriors
          );
          // TODO: REMOVE
          //for( mixing_component_i = 0;
          //     mixing_component_i < num_mixing_components;
          //     mixing_component_i++ ) {
          //  cout << "calculateMixingPosteriors: mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          //} // End foreach mixing_component_i
        } // End if num_mixing_components == 1 .. else ..
    
      } // calculateMixingPosteriors( PositionTransitionParametersType const &, vector<ProbabilityType> & )
    
      /**
       * Add prior pseudocounts to an observed Position transition count vector
       * (an Entente).  Modified from hmmer's P7PriorifyTransitionVector in
       * prior.c.
       */
      template <typename PositionTransitionParametersType>
      void
      incorporatePrior (
        PositionTransitionParametersType & entente
      ) const
      {
        uint32_t num_mixing_components = this->size();
        if( num_mixing_components == 1 ) {
          // Then we don't have to mix, we just have to add in the prior
          // pseudocounts.  This is faster than using the mixing code (since we
          // can perform the operations in type Probabilty, which is not a
          // log type).
          PositionTransitionParametersType prior_counts( entente );
          //ScalableParameterCollection<ScoreType, PositionTransitionParameters<ResidueType, DirichletParameterType> prior_counts =
          //  ScalableParameterCollection<ScoreType, PositionTransitionParameters<ResidueType, DirichletParameterType> >();
          prior_counts.operator=( ( *this )[ 0 ] );
          // TODO: REMOVE
          //cout << "In incorporatePrior(..): prior_counts before scaling is " << prior_counts << endl;
          if( has_m_scalar<PositionTransitionParametersType>() ) {
             getScalar<ScoreType, PositionTransitionParametersType>( prior_counts ) /=
              getScalar<ScoreType, PositionTransitionParametersType>( entente );
          }
          //cout << "In incorporatePrior(..): prior_counts after scaling is " << prior_counts << endl;
          //cout << "In incorporatePrior(..): unscaled entente before is " << entente << endl;
          entente += prior_counts;
          //cout << "In incorporatePrior(..): unscaled entente after is " << entente << endl;
    
          return;
        } // End if num_mixing_components == 1
    
        // TODO: We could further speed this up by doing everything in
        // ProbabilityType type, and instead of dividing the entente by the scalar,
        // multiply the prior by the scalar.  If it ends up being 0, then just
        // return.  Also, why *is* the entente scalar in scoretype anyway?
    
        // We first need posterior probs P( q | entente ) for each mixture
        // component q.
        vector<ProbabilityType> mixing_posteriors( this->size() );
        calculateMixingPosteriors( entente, mixing_posteriors );
    
        // Incorporate the priors into the entente counts.
        // Sean Eddy's code says "following Sjolander (1996)"
        MatrixValueType xi;
        ScoreType unscaled_entente_value;
        ScoreType totc, tot;
        uint32_t mixing_component_i;
        uint32_t i;

        ScoreType entente_scalar;
        if( has_m_scalar<PositionTransitionParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, PositionTransitionParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        ////////
        // Match
        totc =
          entente[ Transition::fromMatch ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromMatch>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromMatch ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromMatch ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromMatch ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromMatch ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromMatch ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromMatch ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromMatch ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new Match entente value is " << entente[ Transition::fromMatch ][ i ] << endl;

        } // End foreach Match transition target i

        ////////
        // Insertion
        totc =
          entente[ Transition::fromInsertion ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromInsertion>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromInsertion ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromInsertion ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromInsertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromInsertion ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromInsertion ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromInsertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromInsertion ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new Insertion entente value is " << entente[ Transition::fromInsertion ][ i ] << endl;

        } // End foreach Insertion transition target i

        ////////
        // Deletion
        totc =
          entente[ Transition::fromDeletion ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromDeletion>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromDeletion ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromDeletion ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromDeletion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromDeletion ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromDeletion ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromDeletion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromDeletion ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new Deletion entente value is " << entente[ Transition::fromDeletion ][ i ] << endl;

        } // End foreach Deletion transition target i
    
      } // incorporatePrior( PositionTransitionParametersType & ) const
    
      /**
       * Calculate and return ln( P( entente | dirichlet ) ) for component
       * mixing_component_i
       */
      template <typename PositionTransitionParametersType>
      double
      calculateLogPosterior (
        PositionTransitionParametersType const & entente,
        uint32_t mixing_component_i
      ) const
      {
        ScoreType entente_scalar;
        if( has_m_scalar<PositionTransitionParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, PositionTransitionParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // Note that we directly access the m_probs members of the
        // MultinomialDistributions, which are public.
        double prod = 0.0; // log( 1.0 );
        prod += calculateLogProbability( entente[ Transition::fromMatch ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromMatch ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromMatch ].size() );
        prod += calculateLogProbability( entente[ Transition::fromInsertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromInsertion ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromInsertion ].size() );
        prod += calculateLogProbability( entente[ Transition::fromDeletion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromDeletion ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromDeletion ].size() );
        return prod;
      } // calculateLogPosterior( PositionTransitionParametersType const &, uint32_t ) const
    
    }; // End inner class DynamicProgramming::DirichletMixturePositionTransitionPrior

    template <typename DirichletParameterType>
    class DirichletMixturePreAlignPrior :
      public DirichletStuff,
      public vector<PreAlignParameters<ResidueType, DirichletParameterType> >
    {
    public:
      // TODO: Protect access, ensure it sums to 1 (and is non-negative), etc.
      vector<ProbabilityType> m_mixingProbs;
    
      DirichletMixturePreAlignPrior () :
        m_mixingProbs()
      {
        // Do nothing else.
      }; // <init>()
    
      DirichletMixturePreAlignPrior (
        uint32_t length
      ) :
        vector<PreAlignParameters<ResidueType, DirichletParameterType> >( length ),
        m_mixingProbs()
      {
        reinitialize( length );
      } // <init>( uint32_t )
    
      /**
       * Reset all to their defaults, and mixing probs to 1/length.
       */
      void
      reinitialize (
        uint32_t length
      )
      {
        if( this->size() != length ) {
          this->resize( length );
        }
        if( m_mixingProbs.size() != length ) {
          m_mixingProbs.resize( length );
        }
    
        if( length > 0 ) {
          uint32_t last_pos = length - 1;
          uint32_t pos_i;
          for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
            this->vector<PreAlignParameters<ResidueType, DirichletParameterType> >::operator[]( pos_i ).reinitialize();
            m_mixingProbs[ pos_i ] = 1.0 / length;
          }
        }
      } // reinitialize( uint32_t )
    
      /**
       * Reinitialize with size 1 and set all counts to 1.0.
       */
      void
      reinitializeToLaplace ()
      {
        reinitializeToEven( 1.0 );
      } // reinitializeToLaplace()
    
      /**
       * Reinitialize with size 1 and set all counts to count.
       */
      void
      reinitializeToEven (
        DirichletParameterType const count
      )
      {
        reinitialize( 1 );
        ( *this )[ 0 ][ Transition::fromPreAlign ] = count;
        ( *this )[ 0 ][ Transition::fromBegin ] = count;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ( *this )[ 0 ][ Transition::fromDeletionIn ] = count;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        ( *this )[ 0 ][ Emission::PreAlignInsertion ] = count;
      } // reinitializeToEven( DirichletParameterType )
    
      /**
       * Calculate the posterior probabilities of the mixing components, given
       * the observed expected counts for one sequence (at one position).
       */
      template <typename PreAlignParametersType>
      void
      calculateMixingPosteriors (
        PreAlignParametersType const & entente,
        vector<ProbabilityType> & mixing_posteriors
      ) const
      {
        uint32_t num_mixing_components = this->size();
        if( mixing_posteriors.size() != num_mixing_components ) {
          mixing_posteriors.resize( num_mixing_components );
        }
    
        uint32_t mixing_component_i;
        if( num_mixing_components == 1 ) {
          mixing_posteriors[ 0 ] = 1.0;
        } else { // if num_mixing_components == 1 .. else ..
          // Do it in log space, so for now mixing_posteriors is a vector of
          // logs.  We will convert it back later.
          vector<double> log_mixing_posteriors( num_mixing_components );
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            if( m_mixingProbs[ mixing_component_i ] == 0 ) {
              log_mixing_posteriors[ mixing_component_i ] = log( 0.0 ); // - numeric_limits<double>::infinity();
              continue; // Don't try to get posterior prob if prior prob is 0.
            }
            log_mixing_posteriors[ mixing_component_i ] =
              (
                log( m_mixingProbs[ mixing_component_i ] ) +
                calculateLogPosterior( entente, mixing_component_i )
              );
            // TODO: REMOVE
            //cout << "calculateMixingPosteriors: log mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          } // End foreach mixing_component_i
          // Now exponentiate and normalize
          DynamicProgramming::DirichletStuff::exponentiateAndNormalize(
            log_mixing_posteriors,
            mixing_posteriors
          );
          // TODO: REMOVE
          //for( mixing_component_i = 0;
          //     mixing_component_i < num_mixing_components;
          //     mixing_component_i++ ) {
          //  cout << "calculateMixingPosteriors: mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          //} // End foreach mixing_component_i
        } // End if num_mixing_components == 1 .. else ..
    
      } // calculateMixingPosteriors( PreAlignParametersType const &, vector<ProbabilityType> & )
    
      /**
       * Add prior pseudocounts to an observed emission count vector (an
       * Entente).  Modified from hmmer's P7PriorifyEmissionVector in prior.c.
       */
      template <typename PreAlignParametersType>
      void
      incorporatePrior ( PreAlignParametersType & entente ) const
      {
        // Note that as a side effect, the entente will be normalized (if/when
        // unscaled), unless there is only 1 mixing component.
  
        uint32_t num_mixing_components = this->size();
        if( num_mixing_components == 1 ) {
          // Then we don't have to mix, we just have to add in the prior
          // pseudocounts.  This is faster than using the mixing code (since we
          // can perform the operations in type Probabilty, which is not a
          // log type).
          PreAlignParametersType prior_counts( entente );
          //ScalableParameterCollection<ScoreType, PreAlignParameters<ResidueType, DirichletParameterType> prior_counts =
          //  ScalableParameterCollection<ScoreType, PreAlignParameters<ResidueType, DirichletParameterType> >();
          prior_counts.operator=( ( *this )[ 0 ] );
          // TODO: REMOVE
          //cout << "In incorporatePrior(..): prior_counts before scaling is " << prior_counts << endl;
          if( has_m_scalar<PreAlignParametersType>() ) {
             getScalar<ScoreType, PreAlignParametersType>( prior_counts ) /=
              getScalar<ScoreType, PreAlignParametersType>( entente );
          }
          //cout << "In incorporatePrior(..): prior_counts after scaling is " << prior_counts << endl;
          //cout << "In incorporatePrior(..): unscaled entente before is " << entente << endl;
          entente += prior_counts;
          //cout << "In incorporatePrior(..): unscaled entente after is " << entente << endl;
    
          return;
        } // End if num_mixing_components == 1
    
        // We first need posterior probs P( q | entente ) for each mixture
        // component q.
        vector<ProbabilityType> mixing_posteriors( this->size() );
        calculateMixingPosteriors( entente, mixing_posteriors );
    
        // Incorporate the priors into the entente counts.
        // Sean Eddy's code says "following Sjolander (1996)"
        MatrixValueType xi;
        ScoreType unscaled_entente_value;
        ScoreType totc, tot;
        uint32_t mixing_component_i;
        uint32_t alphabet_size;
        uint32_t i;

        ScoreType entente_scalar;
        if( has_m_scalar<PreAlignParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, PreAlignParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        ////////
        // PreAlign
        totc =
          entente[ Transition::fromPreAlign ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromPreAlign>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromPreAlign ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromPreAlign ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromPreAlign ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromPreAlign ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new PreAlign entente value is " << entente[ Transition::fromPreAlign ][ i ] << endl;

        } // End foreach PreAlign transition target i

        ////////
        // Begin
        totc =
          entente[ Transition::fromBegin ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromBegin>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromBegin ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromBegin ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromBegin ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromBegin ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromBegin ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromBegin ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromBegin ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new Begin entente value is " << entente[ Transition::fromBegin ][ i ] << endl;

        } // End foreach Begin transition target i

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ////////
        // DeletionIn
        totc =
          entente[ Transition::fromDeletionIn ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromDeletionIn>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromDeletionIn ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromDeletionIn ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new DeletionIn entente value is " << entente[ Transition::fromDeletionIn ][ i ] << endl;

        } // End foreach DeletionIn transition target i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        alphabet_size =
          entente[ Emission::PreAlignInsertion ].size();
        // TODO: REMOVE
        //cout << "entente[ Match ].total() returns " << entente[ Emission::PreAlignInsertion ].total() << endl;
        totc =
          entente[ Emission::PreAlignInsertion ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < alphabet_size; i++ ) {
          unscaled_entente_value =
            entente[ Emission::PreAlignInsertion ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::PreAlignInsertion ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          entente[ Emission::PreAlignInsertion ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new entente value is " << entente[ Emission::PreAlignInsertion ][ i ] << endl;
        } // End foreach residue i
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      } // incorporatePrior( PreAlignParametersType & ) const
    
      /**
       * Calculate and return ln( P( entente | dirichlet ) ) for component
       * mixing_component_i.  Sums up over the different alphabet types.
       */
      // TODO: Is adding the logposteriors the right thing to do?
      template <typename PreAlignParametersType>
      double
      calculateLogPosterior (
        PreAlignParametersType const & entente,
        uint32_t mixing_component_i
      ) const
      {
        ScoreType entente_scalar;
        if( has_m_scalar<PreAlignParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, PreAlignParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // Note that we directly access the m_probs members of the
        // MultinomialDistributions, which are public.
        double prod = 0.0; // log( 1.0 );
        prod += calculateLogProbability( entente[ Transition::fromPreAlign ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromPreAlign ].size() );
        prod += calculateLogProbability( entente[ Transition::fromBegin ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromBegin ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromBegin ].size() );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        prod += calculateLogProbability( entente[ Transition::fromDeletionIn ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromDeletionIn ].size() );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        prod += calculateLogProbability( entente[ Emission::PreAlignInsertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].m_probs, ( *this )[ mixing_component_i ][ Emission::PreAlignInsertion ].size() );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
        return prod;
      } // calculateLogPosterior( PreAlignParametersType const &, uint32_t ) const
    
    }; // End inner class DynamicProgramming::DirichletMixturePreAlignPrior

    template <typename DirichletParameterType>
    class DirichletMixturePostAlignPrior :
      public DirichletStuff,
      public vector<PostAlignParameters<ResidueType, DirichletParameterType> >
    {
    public:
      // TODO: Protect access, ensure it sums to 1 (and is non-negative), etc.
      vector<ProbabilityType> m_mixingProbs;
    
      DirichletMixturePostAlignPrior () :
        m_mixingProbs()
      {
        // Do nothing else.
      }; // <init>()
    
      DirichletMixturePostAlignPrior (
        uint32_t length
      ) :
        vector<PostAlignParameters<ResidueType, DirichletParameterType> >( length ),
        m_mixingProbs()
      {
        reinitialize( length );
      } // <init>( uint32_t )
    
      /**
       * Reset all to their defaults, and mixing probs to 1/length.
       */
      void
      reinitialize (
        uint32_t length
      )
      {
        if( this->size() != length ) {
          this->resize( length );
        }
        if( m_mixingProbs.size() != length ) {
          m_mixingProbs.resize( length );
        }
    
        if( length > 0 ) {
          uint32_t last_pos = length - 1;
          uint32_t pos_i;
          for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
            this->vector<PostAlignParameters<ResidueType, DirichletParameterType> >::operator[]( pos_i ).reinitialize();
            m_mixingProbs[ pos_i ] = 1.0 / length;
          }
        }
      } // reinitialize( uint32_t )
    
      /**
       * Reinitialize with size 1 and set all counts to 1.0.
       */
      void
      reinitializeToLaplace ()
      {
        reinitializeToEven( 1.0 );
      } // reinitializeToLaplace()
    
      /**
       * Reinitialize with size 1 and set all counts to count.
       */
      void
      reinitializeToEven (
        DirichletParameterType const count
      )
      {
        reinitialize( 1 );
#ifdef USE_END_DISTRIBUTION
        // For now we we never allow loops, so the End distribution is constant.
        //( *this )[ 0 ][ Transition::fromEnd ] = count;
#endif // USE_END_DISTRIBUTION
        ( *this )[ 0 ][ Transition::fromPostAlign ] = count;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ( *this )[ 0 ][ Transition::fromDeletionOut ] = count;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        ( *this )[ 0 ][ Emission::PostAlignInsertion ] = count;
      } // reinitializeToEven( DirichletParameterType )
    
      /**
       * Calculate the posterior probabilities of the mixing components, given
       * the observed expected counts for one sequence (at one position).
       */
      template <typename PostAlignParametersType>
      void
      calculateMixingPosteriors (
        PostAlignParametersType const & entente,
        vector<ProbabilityType> & mixing_posteriors
      ) const
      {
        uint32_t num_mixing_components = this->size();
        if( mixing_posteriors.size() != num_mixing_components ) {
          mixing_posteriors.resize( num_mixing_components );
        }
    
        uint32_t mixing_component_i;
        if( num_mixing_components == 1 ) {
          mixing_posteriors[ 0 ] = 1.0;
        } else { // if num_mixing_components == 1 .. else ..
          // Do it in log space, so for now mixing_posteriors is a vector of
          // logs.  We will convert it back later.
          vector<double> log_mixing_posteriors( num_mixing_components );
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            if( m_mixingProbs[ mixing_component_i ] == 0 ) {
              log_mixing_posteriors[ mixing_component_i ] = log( 0.0 ); // - numeric_limits<double>::infinity();
              continue; // Don't try to get posterior prob if prior prob is 0.
            }
            log_mixing_posteriors[ mixing_component_i ] =
              (
                log( m_mixingProbs[ mixing_component_i ] ) *
                calculateLogPosterior( entente, mixing_component_i )
              );
            // TODO: REMOVE
            //cout << "calculateMixingPosteriors: log mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          } // End foreach mixing_component_i
          // Now exponentiate and normalize
          DynamicProgramming::DirichletStuff::exponentiateAndNormalize(
            log_mixing_posteriors,
            mixing_posteriors
          );
          // TODO: REMOVE
          //for( mixing_component_i = 0;
          //     mixing_component_i < num_mixing_components;
          //     mixing_component_i++ ) {
          //  cout << "calculateMixingPosteriors: mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
          //} // End foreach mixing_component_i
        } // End if num_mixing_components == 1 .. else ..
    
      } // calculateMixingPosteriors( PostAlignParametersType const &, vector<ProbabilityType> & )
    
      /**
       * Add prior pseudocounts to an observed emission count vector (an
       * Entente).  Modified from hmmer's P7PriorifyEmissionVector in prior.c.
       */
      template <typename PostAlignParametersType>
      void
      incorporatePrior (
        PostAlignParametersType & entente
      ) const
      {
        // Note that as a side effect, the entente will be normalized (if/when
        // unscaled), unless there is only 1 mixing component.
  
        uint32_t num_mixing_components = this->size();
        if( num_mixing_components == 1 ) {
          // Then we don't have to mix, we just have to add in the prior
          // pseudocounts.  This is faster than using the mixing code (since we
          // can perform the operations in type Probabilty, which is not a
          // log type).
          PostAlignParametersType prior_counts( entente );
          //ScalableParameterCollection<ScoreType, PostAlignParameters<ResidueType, DirichletParameterType> > prior_counts =
          //  ScalableParameterCollection<ScoreType, PostAlignParameters<ResidueType, DirichletParameterType> >();
          prior_counts.operator=( ( *this )[ 0 ] );
          // TODO: REMOVE
          //cout << "In incorporatePrior(..): prior_counts before scaling is " << prior_counts << endl;
          if( has_m_scalar<PostAlignParametersType>() ) {
             getScalar<ScoreType, PostAlignParametersType>( prior_counts ) /=
              getScalar<ScoreType, PostAlignParametersType>( entente );
          }
          //cout << "In incorporatePrior(..): prior_counts after scaling is " << prior_counts << endl;
          //cout << "In incorporatePrior(..): unscaled entente before is " << entente << endl;
          entente += prior_counts;
          //cout << "In incorporatePrior(..): unscaled entente after is " << entente << endl;
    
          return;
        } // End if num_mixing_components == 1
    
        // We first need posterior probs P( q | entente ) for each mixture
        // component q.
        vector<ProbabilityType> mixing_posteriors( this->size() );
        calculateMixingPosteriors( entente, mixing_posteriors );

        // Incorporate the priors into the entente counts.
        // Sean Eddy's code says "following Sjolander (1996)"
        MatrixValueType xi;
        ScoreType unscaled_entente_value;
        ScoreType totc, tot;
        uint32_t mixing_component_i;
        uint32_t alphabet_size;
        uint32_t i;

        ScoreType entente_scalar;
        if( has_m_scalar<PostAlignParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, PostAlignParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

#ifdef USE_END_DISTRIBUTION
        // For now we we never allow loops, so the End distribution is constant.
        //////////
        //// End
        //totc =
        //  entente[ Transition::fromEnd ].total();
        //// TODO: REMOVE
        ////cout << "totc before scaling is " << totc << endl;
        //totc /= entente_scalar;
        //// TODO: REMOVE
        ////cout << "totc after scaling is " << totc << endl;
        ////cout << "as a double, totc is " << toDouble( totc ) << endl;
        //for( i = 0; i < seqan::ValueSize<TransitionFromEnd>::VALUE; i++ ) {
        //  unscaled_entente_value =
        //    entente[ Transition::fromEnd ][ i ];
        //  // TODO: REMOVE
        //  //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
        //  // TODO: REMOVE
        //  //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
        //  unscaled_entente_value /= entente_scalar;
        //  // TODO: REMOVE
        //  //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
        //  // TODO: REMOVE
        //  //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
        //  xi = 0.0;
        //  for( mixing_component_i = 0;
        //       mixing_component_i < num_mixing_components;
        //       mixing_component_i++ ) {
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromEnd ].total() << endl;
        //    tot = ( *this )[ mixing_component_i ][ Transition::fromEnd ].total();
        //    //cout << "\t as a ScoreType, this is " << tot << endl;
        //    tot =
        //      totc +
        //      ( *this )[ mixing_component_i ][ Transition::fromEnd ].total();
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
        //    xi +=
        //      mixing_posteriors[ mixing_component_i ] *
        //      ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromEnd ][ i ] ) / tot;
        //    // TODO: REMOVE
        //    //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromEnd ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
        //    //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
        //  } // End foreach mixing component ..
        //  entente[ Transition::fromEnd ][ i ] =
        //    ( entente_scalar * xi );
        //  // TODO: REMOVE
        //  //cout << "[" << i << "] new End entente value is " << entente[ Transition::fromEnd ][ i ] << endl;
        //
        //} // End foreach End transition target i
#endif // USE_END_DISTRIBUTION

        ////////
        // PostAlign
        totc =
          entente[ Transition::fromPostAlign ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromPostAlign>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromPostAlign ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromPostAlign ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromPostAlign ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromPostAlign ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new PostAlign entente value is " << entente[ Transition::fromPostAlign ][ i ] << endl;

        } // End foreach PostAlign transition target i

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        ////////
        // DeletionOut
        totc =
          entente[ Transition::fromDeletionOut ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < seqan::ValueSize<TransitionFromDeletionOut>::VALUE; i++ ) {
          unscaled_entente_value =
            entente[ Transition::fromDeletionOut ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::Match ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          } // End foreach mixing component ..
          entente[ Transition::fromDeletionOut ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new DeletionOut entente value is " << entente[ Transition::fromDeletionOut ][ i ] << endl;

        } // End foreach DeletionOut transition target i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT


#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        alphabet_size =
          entente[ Emission::PostAlignInsertion ].size();
        // TODO: REMOVE
        //cout << "entente[ Match ].total() returns " << entente[ Emission::PostAlignInsertion ].total() << endl;
        totc =
          entente[ Emission::PostAlignInsertion ].total();
        // TODO: REMOVE
        //cout << "totc before scaling is " << totc << endl;
        totc /= entente_scalar;
        // TODO: REMOVE
        //cout << "totc after scaling is " << totc << endl;
        //cout << "as a double, totc is " << toDouble( totc ) << endl;
        for( i = 0; i < alphabet_size; i++ ) {
          unscaled_entente_value =
            entente[ Emission::PostAlignInsertion ][ i ];
          // TODO: REMOVE
          //cout << "[" << i << "] scaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t (before converting to ScoreType, this was " << entente[ Emission::PostAlignInsertion ][ i ] << ")" << endl;
          unscaled_entente_value /= entente_scalar;
          // TODO: REMOVE
          //cout << "[" << i << "] unscaled_entente_value is " << unscaled_entente_value << endl;
          // TODO: REMOVE
          //cout << "\t as double, this is " << toDouble( unscaled_entente_value ) << endl;
          xi = 0.0;
          for( mixing_component_i = 0;
               mixing_component_i < num_mixing_components;
               mixing_component_i++ ) {
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] component prior total is " << ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].total() << endl;
            tot = ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].total();
            //cout << "\t as a ScoreType, this is " << tot << endl;
            tot =
              totc +
              ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].total();
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] tot is " << tot << endl;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] mixing_posteriors[ " << mixing_component_i << " ] is " << mixing_posteriors[ mixing_component_i ] << endl;
            xi +=
              mixing_posteriors[ mixing_component_i ] *
              ( unscaled_entente_value + ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ][ i ] ) / tot;
            // TODO: REMOVE
            //cout << "[" << i << ", " << mixing_component_i << "] Just added to xi " << mixing_posteriors[ mixing_component_i ] << " * " << "( ( " << toDouble( unscaled_entente_value ) << " + " << toDouble( ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ][ i ] ) << " ) / " << toDouble( tot ) << " )" << endl;
            //cout << "[" << i << ", " << mixing_component_i << "] xi is now " << xi << endl;
          }
          entente[ Emission::PostAlignInsertion ][ i ] =
            ( entente_scalar * xi );
          // TODO: REMOVE
          //cout << "[" << i << "] new entente value is " << entente[ Emission::PostAlignInsertion ][ i ] << endl;
        } // End foreach residue i
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      } // incorporatePrior( PostAlignParametersType & ) const
    
      /**
       * Calculate and return ln( P( entente | dirichlet ) ) for component
       * mixing_component_i.  Sums up over the different alphabet types.
       */
      // TODO: Is adding the logposteriors the right thing to do?
      template <typename PostAlignParametersType>
      double
      calculateLogPosterior (
        PostAlignParametersType const & entente,
        uint32_t mixing_component_i
      ) const
      {
        ScoreType entente_scalar;
        if( has_m_scalar<PostAlignParametersType>() ) {
          entente_scalar =
            getScalar<ScoreType, PostAlignParametersType>( entente );
        } else {
          entente_scalar = 1.0;
        }

        // Note that we directly access the m_probs members of the
        // MultinomialDistributions, which are public.
        double prod = 0.0; // log( 1.0 );
#ifdef USE_END_DISTRIBUTION
        // For now we we never allow loops, so the End distribution is constant.
        //prod += calculateLogProbability( entente[ Transition::fromEnd ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromEnd ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromEnd ].size() );
#endif // USE_END_DISTRIBUTION
        prod += calculateLogProbability( entente[ Transition::fromPostAlign ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromPostAlign ].size() );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        prod += calculateLogProbability( entente[ Transition::fromDeletionOut ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].m_probs, ( *this )[ mixing_component_i ][ Transition::fromDeletionOut ].size() );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
        prod += calculateLogProbability( entente[ Emission::PostAlignInsertion ].m_probs, entente_scalar, ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].m_probs, ( *this )[ mixing_component_i ][ Emission::PostAlignInsertion ].size() );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
        return prod;
      } // calculateLogPosterior( PostAlignParametersType const &, uint32_t ) const
    
    }; // End inner class DynamicProgramming::DirichletMixturePostAlignPrior

    /**
     * A vector of distances from a particular sequence to all of the sequences
     * for a particular position of the Entente.
     */
    class SequencePositionEntenteDistances : public vector<double>
    {
    public:
      SequencePositionEntenteDistances () :
        vector<double>()
      {
        // Do nothing else.
      } // <init>()

      SequencePositionEntenteDistances ( uint32_t const & sequence_count ) :
        vector<double>( sequence_count )
      {
        reset();
      } // <init>( uint32_t const & )

      void reinitialize (
        uint32_t const & sequence_count
      )
      {
        if( sequence_count != size() ) {
          resize( sequence_count );
        }
        reset();
      } // reinitialize( uint32_t const& )

      void reset ()
      {
        uint32_t length = size();
        // Initialize everything to 0.
        uint32_t seq_i;
        for( seq_i = 0; seq_i < length; seq_i++ ) {
          this->operator[]( seq_i ) = 0.0;
        }
      } // reset()

      SequencePositionEntenteDistances &
      operator= ( SequencePositionEntenteDistances const& other_seq_ped )
      {
        if( size() != other_seq_ped.size() ) {
          resize( other_seq_ped.size() );
        }
        uint32_t last_seq = other_seq_ped.size() - 1;
        uint32_t seq_i;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          this->operator[]( seq_i ) = other_seq_ped[ seq_i ];
        }
        return *this;
      } // operator=( SequencePositionEntenteDistances const& )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        SequencePositionEntenteDistances const& seq_ped )
      {
        uint32_t last_seq = seq_ped.size() - 1;
        uint32_t seq_i;
        for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
          os << "\t\tTo Sequence " << seq_i << ":\t" << seq_ped[ seq_i ] << endl;
        }
        return os;
      } // friend operator<< ( basic_ostream, SequencePositionEntenteDistances const & )

    }; // End inner class DynamicProgramming::SequencePositionEntenteDistances

    /**
     * A matrix of distances, represented as a vector of
     * SequencePositionEntenteDistances.  Its index is a sequence number.
     */
    class PositionEntenteDistances : public vector<SequencePositionEntenteDistances>
    {
    public:
      PositionEntenteDistances () :
        vector<SequencePositionEntenteDistances>()
      {
        // Do nothing else.
      } // <init>()

      PositionEntenteDistances (
        uint32_t const & sequence_count
      ) :
        vector<SequencePositionEntenteDistances>( sequence_count )
      {
        // TODO: For symmetric distances, we need only allocate the upper triangle of the matrix (or the lower triangle).

        uint32_t last_seq = sequence_count - 1;
        uint32_t seq_i;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          this->operator[]( seq_i ).reinitialize( sequence_count );
        }
      } // <init>( uint32_t const& )

      void reinitialize (
        uint32_t const & sequence_count
      )
      {
        // TODO: For symmetric distances, we need only allocate the upper triangle of the matrix (or the lower triangle).

        if( sequence_count != vector<SequencePositionEntenteDistances>::size() ) {
          vector<SequencePositionEntenteDistances>::resize( sequence_count );
        }
        uint32_t last_seq = sequence_count - 1;
        uint32_t seq_i;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          this->operator[]( seq_i ).reinitialize( sequence_count );
        }
      } // reinitialize( uint32_t const& )

      void reset ()
      {
        uint32_t length = vector<SequencePositionEntenteDistances>::size();
        // Initialize everything to 0.
        uint32_t seq_i;
        for( seq_i = 0; seq_i < length; seq_i++ ) {
          ( this->operator[]( seq_i ) ).reset();
        }
      } // reset()


      PositionEntenteDistances &
      operator= ( PositionEntenteDistances const& other_ped )
      {
        if( vector<SequencePositionEntenteDistances>::size() != other_ped.size() ) {
          vector<SequencePositionEntenteDistances>::resize( other_ped.size() );
        }
        uint32_t last_seq = other_ped.size() - 1;
        uint32_t seq_i;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          this->operator[]( seq_i ) = other_ped[ seq_i ];
        }
        return *this;
      } // operator=( PositionEntenteDistances const& )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        PositionEntenteDistances const& ped )
      {
        uint32_t last_seq = ped.size() - 1;
        uint32_t seq_i;
        for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
          os << ped[ seq_i ] << endl;
        }
        return os;
      } // friend operator<< ( basic_ostream, PositionEntenteDistances const & )

    }; // End inner class DynamicProgramming::PositionEntenteDistances

    /**
     * A vector of PositionEntenteDistances, one for each position.
     */
    class EntenteDistances : public vector<PositionEntenteDistances>
    {
    public:
      EntenteDistances () :
        vector<PositionEntenteDistances>()
      {
        // Do nothing else.
      } // <init>()

      EntenteDistances (
        uint32_t const & entente_length,
        uint32_t const & sequence_count
      ) :
        vector<PositionEntenteDistances>( entente_length )
      {
        // Initialize each contained PositionEntenteDistances
        uint32_t last_row = entente_length - 1;
        uint32_t row_i;
        for( row_i = 0; row_i <= last_row; row_i++ ) {
          this->operator[]( row_i ).reinitialize( sequence_count );
        } // End foreach row
      } // <init>( uint32_t const &, uint32_t const & )

      template <typename ProfileType>
      EntenteDistances (
        ProfileType const& profile,
        uint32_t const & sequence_count
      ) :
        vector<PositionEntenteDistances>( profile.length() )
      {
        // Initialize each contained PositionEntenteDistances
        uint32_t last_row = profile.length() - 1;
        uint32_t row_i;
        for( row_i = 0; row_i <= last_row; row_i++ ) {
          this->operator[]( row_i ).reinitialize( sequence_count );
        } // End foreach row
      } // <init>( Profile const&, uint32_t const& )

      template <typename ProfileType>
      void reinitialize (
        ProfileType const& profile,
        uint32_t const & sequence_count
      )
      {
        reinitialize( profile.length(), sequence_count );
      } // reinitialize( Profile const&, uint32_t const& )

      void reinitialize (
        uint32_t const & entente_length,
        uint32_t const & sequence_count
      )
      {
        if( entente_length != vector<PositionEntenteDistances>::size() ) {
          vector<PositionEntenteDistances>::resize( entente_length );
        }
        // Initialize each contained PositionEntenteDistances
        uint32_t last_row = entente_length - 1;
        uint32_t row_i;
        for( row_i = 0; row_i <= last_row; row_i++ ) {
          this->operator[]( row_i ).reinitialize( sequence_count );
        } // End foreach row
      } // reinitialize( uint32_t const&, uint32_t const& )

      void reset ()
      {
        uint32_t length = vector<PositionEntenteDistances>::size();
        // Initialize everything to 0.
        uint32_t seq_i;
        for( seq_i = 0; seq_i < length; seq_i++ ) {
          ( this->operator[]( seq_i ) ).reset();
        }
      } // reset()

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        EntenteDistances const& distances )
      {
        uint32_t last_pos = distances.size() - 1;
        uint32_t pos_i;
        uint32_t last_seq = ((std::vector<PositionEntenteDistances>)(distances[ 0 ])).size() - 1;
        uint32_t seq_i;
        for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
          os << "From Sequence " << seq_i << ":" << endl;
          for( pos_i = 0; pos_i <= last_pos ; pos_i++ ) {
            os << "\tAt Position " << pos_i << ": " << endl;
            os << distances[ pos_i ][ seq_i ] << endl;
          }
          os << endl;
        }
        return os;
      } // friend operator<< ( basic_ostream, EntenteDistances const & )

    }; // End inner class DynamicProgramming::EntenteDistances

    /**
     * The idea of this is that the forward score of a sequence can be
     * calculated by multiplying each profile position parameter by a
     * coefficient and adding.  Or, looked at another way, the sequence score
     * is a product of a row vector (this) times a column vector (the profile
     * position's emission parameters, with an extra 1 for the constant that
     * accounts for deletions of the position).
     */
    template <class CoefficientsType>
    class PositionSpecificSequenceScoreCoefficientsTemplate :
      public PositionSpecificParameters<ResidueType, CoefficientsType>
    {
    public:
      /**
       * For efficiency we can scale all of the coefficients.  When
       * useRabinerScaling is true, m_inverseScalar will typically be the score
       * of the sequences calculated using the values of the corresponding
       * profile position at the time that the coefficients were last updated.
       *
       * You can return to the unscaled values by multiplying each scaled value
       * by m_inverseScalar.
       */
      ScoreType m_inverseScalar;

      CoefficientsType m_constant;
      PositionSpecificSequenceScoreCoefficientsTemplate () :
        PositionSpecificParameters<ResidueType, CoefficientsType>(), // set the missingValueElement.
        m_inverseScalar( 1.0 ),
        m_constant( 0.0 )
      {
        // Do nothing else.
      } // <init>()

      void
      reinitialize ()
      {
        PositionSpecificParameters<ResidueType, CoefficientsType>::reinitialize();
        m_inverseScalar = 1.0;
        m_constant = 0.0;
      } // reinitialize()

      /**
       * Return an unscaled version of this.  Note that its CoefficientsType
       * will be ScoreType.
       */
      PositionSpecificSequenceScoreCoefficientsTemplate<ScoreType>
      createUnscaledCopy ()
      {
        PositionSpecificSequenceScoreCoefficientsTemplate<ScoreType> unscaled_copy;
        unscaled_copy = *this;
        unscaled_copy *= unscaled_copy.m_inverseScalar;
        unscaled_copy.m_inverseScalar = 1.0;
        return unscaled_copy;
      } // createUnscaledCopy()

      /**
       * Calculate the (forward/backward) score of the sequence with this
       * PositionSpecificSequenceScoreCoefficientsTemplate using that and the
       * corresponding ProfilePosition.
       */
      template <typename ProfileType>
      CoefficientsType calculateScore (
        ProfilePosition<ResidueType, ProfileType> const & profile_position
      ) const {
        return calculateScore( profile_position, NULL );
      } // calculateScore( ProfilePosition const & )

      /**
       * Calculate the (forward/backward) score of the sequence with this
       * PositionSpecificSequenceScoreCoefficientsTemplate using that and the
       * corresponding ProfilePosition.  If the given pointer to a
       * PositionEntente object is non-null, fill that object with the
       * expected number of emissions of each type, scaled by
       * 1/(its maximum value*m_inverseScalar).
       */
      template <typename ProfileType>
      CoefficientsType calculateScore (
        ProfilePosition<ResidueType, ProfileType> const & profile_position,
        PositionEntente * expected_emission_counts
      ) const {

        // TODO: REMOVE
        //cout << "in calculateScore(): profile_position is " << profile_position << endl;

        if( expected_emission_counts != NULL ) {
          expected_emission_counts->zero();
        }

        CoefficientsType coefficients_type_score = m_constant;

        MultinomialDistribution<ResidueType, CoefficientsType> products =
          this->operator[]( Emission::Match );
        products *= profile_position[ Emission::Match ];
        // TODO: REMOVE
        //cout << "in calculateScore(): products is now " << products << endl;
        coefficients_type_score += products.total();

        if( expected_emission_counts != NULL ) {
          expected_emission_counts->m_scalar =
            1.0 / products.maximumValue();
          products *= expected_emission_counts->m_scalar;
            
          ( *expected_emission_counts )[ Emission::Match ] =
            products;
          // TODO: REMOVE
          //cout << "in calculateScore(): expected_emission_counts->m_scalar is now " << expected_emission_counts->m_scalar << ", and ( *expected_emission_counts )[ Emission::Match ] is now " << ( *expected_emission_counts )[ Emission::Match ] << endl;
        } // End if expected_emission_counts != NULL

        ScoreType score = coefficients_type_score;
        score *= m_inverseScalar;

        if( expected_emission_counts != NULL ) {
          // Also update the scalar for the expected emission counts.
          expected_emission_counts->m_scalar /= m_inverseScalar;
        }

        return score;
      } // calculateScore( ProfilePosition const&, PositionEntente * )

      /**
       * Write a comma-separated list of the parameters to the stream.
       */
      template<class CharT, class Traits>
      void
      writePositionSpecificParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        os << this->operator[]( Emission::Match );
      } // writePositionSpecificParameters( basic_ostream & )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        PositionSpecificSequenceScoreCoefficientsTemplate const& pos
      )
      {
        os << "[ ";
        pos.writePositionSpecificParameters( os );
        os << ", constant:" << pos.m_constant;
        os << ", inverseScalar:" << pos.m_inverseScalar;
        os << " ]";
        return os;
      } // friend operator<< ( basic_ostream &, PositionSpecificSequenceScoreCoefficientsTemplate const& )

      /**
       * Copy the values of the other coefficients into this one.
       */
      template <typename AnyCoefficientsType>
      PositionSpecificSequenceScoreCoefficientsTemplate &
      operator= ( PositionSpecificSequenceScoreCoefficientsTemplate<AnyCoefficientsType> const& other_pos )
      {
        PositionSpecificParameters<ResidueType, CoefficientsType>::operator=( other_pos );
        m_constant = other_pos.m_constant;
        m_inverseScalar = other_pos.m_inverseScalar;
  
        return *this;
      } // operator=( PositionSpecificSequenceScoreCoefficientsTemplate<AnyCoefficientsType> const& )
  
      /**
       * Multiply each contained "distrubution" value by scalar.
       */
      PositionSpecificSequenceScoreCoefficientsTemplate &
      operator*= ( CoefficientsType const& scalar )
      {
        PositionSpecificParameters<ResidueType, CoefficientsType>::operator*=( scalar );
        m_constant *= scalar;
  
        return *this;
      } // operator*=( CoefficientsType& )
  
      /**
       * Divide each contained "distrubution" value by denominator.
       */
      PositionSpecificSequenceScoreCoefficientsTemplate &
      operator/= ( CoefficientsType const& denominator )
      {
        PositionSpecificParameters<ResidueType, CoefficientsType>::operator/=( denominator );
        m_constant /= denominator;
  
        return *this;
      } // operator/=( CoefficientsType& )
  
      /**
       * Calculate and return the Euclidean distance between these coefficients
       * and another vector of coefficients (treating every coefficient as an
       * orthogonal dimension).
       */
      double
      euclideanDistance (
        PositionSpecificSequenceScoreCoefficientsTemplate const& other_pos
      ) const {
        return sqrt( euclideanDistanceSquared( other_pos ) );
      } // euclideanDistance( PositionSpecificSequenceScoreCoefficientsTemplate const& other_pos )
  
      /**
       * Calculate and return the square of the Euclidean distance between
       * these coefficients and another vector of coefficients (treating every
       * coefficient as an orthogonal dimension).  Note that if the two
       * coefficients use different m_inverseScalars, this will be quite a bit
       * slower than if they use the same m_inverseScalars.
       */
      double
      euclideanDistanceSquared (
        PositionSpecificSequenceScoreCoefficientsTemplate const& other_pos
      )  {
        // Since the coefficients class uses scaling, we need to be careful
        // here.  Fortunately we don't need this function to be extremely fast.
        if( m_inverseScalar != other_pos->inverseScalar ) {
          return createUnscaledCopy().euclideanDistanceSquared( other_pos->createUnscaledCopy()) ;
        }
        double squared_euclidean_distance =
            PositionSpecificParameters<ResidueType, CoefficientsType>::euclideanDistanceSquared( other_pos );
        squared_euclidean_distance +=
          ( ( m_constant - other_pos.m_constant ) *
            ( m_constant - other_pos.m_constant ) );
            
        return squared_euclidean_distance;
      } // euclideanDistanceSquared( PositionSpecificSequenceScoreCoefficientsTemplate const& )
  
      /**
       * Set all values to 0, and the inverseScalar to 1.
       */
      void zero ()
      {
        PositionSpecificParameters<ResidueType, CoefficientsType>::zero();
        m_inverseScalar = 1.0;
        m_constant = 0;
      } // zero()

    }; // End inner class DynamicProgramming::PositionSpecificSequenceScoreCoefficientsTemplate

    typedef PositionSpecificSequenceScoreCoefficientsTemplate<ScoreType> PositionSpecificSequenceScoreCoefficients;

    /**
     * Generally we will want to store one coefficients object for each
     * sequence.  This is a handydandy way to store a vector of them.
     */
    class PositionSpecificSequenceScoreCoefficientsVector :
      public vector<PositionSpecificSequenceScoreCoefficients>
    {
    public:
      PositionSpecificSequenceScoreCoefficientsVector () :
        vector<PositionSpecificSequenceScoreCoefficients>()
      {
        // Do nothing else.
      } // <init>()

      PositionSpecificSequenceScoreCoefficientsVector ( uint32_t size ) :
        vector<PositionSpecificSequenceScoreCoefficients>( size )
      {
        // Do nothing else.
      } // <init>( uint32_t )

      void
      reinitialize ()
      {
        reinitialize( 0 );
      } // reinitialize()

      void
      reinitialize ( 
        const uint32_t size
      )
      {
        if( vector<PositionSpecificSequenceScoreCoefficients>::size() != size ) {
          vector<PositionSpecificSequenceScoreCoefficients>::resize( size );
        }
        for( uint32_t seq_i = 0; seq_i < size; seq_i++ ) {
          this->operator[]( seq_i ).reinitialize();
        }
      } // reinitialize( uint32_t )

    }; // End inner class DynamicProgramming::PositionSpecificSequenceScoreCoefficientsVector

    /**
     * For (space) efficiency, this doesn't store del-in or del-out stuff in
     * the case in which we are using del-in and del-out but are disallowing
     * pre-aligns and post-aligns.  They can be stored on a per-row basis
     * instead.
     */
    template <typename ValueType>
    class MatrixCell
    {
    protected:
      ValueType m_match;
      ValueType m_insertion;
      ValueType m_deletion;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
      ValueType m_deletionIn;
      ValueType m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
    
    public:
      MatrixCell () :
        m_match( 0 ),
        m_insertion( 0 ),
        m_deletion( 0 )
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
       ,
       m_deletionIn( 0 ),
       m_deletionOut( 0 )
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
      {
        // Do nothing else.
      } // <init>()
    
      MatrixCell ( MatrixCell const & other_cell ) :
        m_match( other_cell.m_match ),
        m_insertion( other_cell.m_insertion ),
        m_deletion( other_cell.m_deletion )
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        ,
        m_deletionIn( other_cell.m_deletionIn ),
        m_deletionOut( other_cell.m_deletionOut )
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
      {
        // Do nothing else.
      } // <init>( MatrixCell const & )

      // Zero everything.
      void
      reinitialize ()
      {
        m_match = 0;
        m_insertion = 0;
        m_deletion = 0;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        m_deletionIn = 0;
        m_deletionOut = 0;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
      } // reinitialize()

      ValueType &
      operator[] ( dynamicprogramming_Match_subcell_tag const )
      {
        // TODO: REMOVE
        //cout << "Match subcell!" << endl;
        return m_match;
      } // operator[] ( dynamicprogramming_Match_subcell_tag const )
    
      ValueType const&
      operator[] ( dynamicprogramming_Match_subcell_tag const ) const
      {
        // TODO: REMOVE
        //cout << "Match subcell!" << endl;
        return m_match;
      } // operator[] ( dynamicprogramming_Match_subcell_tag const ) const
    
      ValueType &
      operator[] ( dynamicprogramming_Insertion_subcell_tag const )
      {
        // TODO: REMOVE
        //cout << "Insertion subcell!" << endl;
        return m_insertion;
      } // operator[] ( dynamicprogramming_Deletion_subcell_tag const )
    
      ValueType const&
      operator[] ( dynamicprogramming_Insertion_subcell_tag const ) const
      {
        // TODO: REMOVE
        //cout << "Insertion subcell!" << endl;
        return m_insertion;
      } // operator[] ( dynamicprogramming_Deletion_subcell_tag const ) const
    
      ValueType &
      operator[] ( dynamicprogramming_Deletion_subcell_tag const )
      {
        // TODO: REMOVE
        //cout << "Deletion subcell!" << endl;
        return m_deletion;
      } // operator[] ( dynamicprogramming_Deletion_subcell_tag const )
    
      ValueType const&
      operator[] ( dynamicprogramming_Deletion_subcell_tag const ) const
      {
        // TODO: REMOVE
        //cout << "Deletion subcell!" << endl;
        return m_deletion;
      } // operator[] ( dynamicprogramming_Deletion_subcell_tag const ) const
    
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
     ValueType &
     operator[] ( dynamicprogramming_DeletionIn_subcell_tag const )
     {
       // TODO: REMOVE
       //cout << "DeletionIn subcell!" << endl;
       return m_deletionIn;
     } // operator[] ( dynamicprogramming_DeletionIn_subcell_tag const )
   
     ValueType const&
     operator[] ( dynamicprogramming_DeletionIn_subcell_tag const ) const
     {
       // TODO: REMOVE
       //cout << "DeletionIn subcell!" << endl;
       return m_deletionIn;
     } // operator[] ( dynamicprogramming_DeletionIn_subcell_tag const ) const
   
     ValueType &
     operator[] ( dynamicprogramming_DeletionOut_subcell_tag const )
     {
       // TODO: REMOVE
       //cout << "DeletionOut subcell!" << endl;
       return m_deletionOut;
     } // operator[] ( dynamicprogramming_DeletionOut_subcell_tag const )
   
     ValueType const&
     operator[] ( dynamicprogramming_DeletionOut_subcell_tag const ) const
     {
       // TODO: REMOVE
       //cout << "DeletionOut subcell!" << endl;
       return m_deletionOut;
     } // operator[] ( dynamicprogramming_DeletionOut_subcell_tag const ) const
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )

      ValueType &
      operator[] ( Subcell const & subcell )
      {
        if( subcell.isMatch() ) {
          return m_match;
        } else if( subcell.isInsertion() ) {
          return m_insertion;
        } else if( subcell.isDeletion() ) {
          return m_deletion;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        } else if( subcell.isDeletionIn() ) {
          return m_deletionIn;
        } else if( subcell.isDeletionOut() ) {
          return m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        } else {
          // ?!
          cerr << "ERROR: Unrecognized subcell value!" << endl;
          assert( false );
          return m_deletion; // I have to return *something*..
        }
      } // operator[] ( Subcell const & )
    
      ValueType const &
      operator[] ( Subcell const & subcell ) const
      {
        if( subcell.isMatch() ) {
          return m_match;
        } else if( subcell.isInsertion() ) {
          return m_insertion;
        } else if( subcell.isDeletion() ) {
          return m_deletion;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        } else if( subcell.isDeletionIn() ) {
          return m_deletionIn;
        } else if( subcell.isDeletionOut() ) {
          return m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        } else {
          // ?!
          cerr << "ERROR: Unrecognized subcell value!" << endl;
          assert( false );  //TAH 2/12 copied from Paul's earlier code
          return m_deletion; // I have to return *something*..

        }
      } // operator[] ( Subcell const & ) const
    
      MatrixCell &
      operator+= ( MatrixCell const& other_cell )
      {
        m_match += other_cell.m_match;
        m_insertion += other_cell.m_insertion;
        m_deletion += other_cell.m_deletion;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        m_deletionIn += other_cell.m_deletionIn;
        m_deletionOut += other_cell.m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
    
        return *this;
      } // operator+=( MatrixCell const& )
    
      MatrixCell &
      operator*= ( ValueType const& scalar )
      {
        m_match *= scalar;
        m_insertion *= scalar;
        m_deletion *= scalar;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        m_deletionIn *= scalar;
        m_deletionOut *= scalar;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
    
        return *this;
      } // operator*=( MatrixCell const& )
    
      MatrixCell &
      operator/= ( ValueType const& denominator )
      {
        m_match /= denominator;
        m_insertion /= denominator;
        m_deletion /= denominator;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        m_deletionIn /= denominator;
        m_deletionOut /= denominator;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
    
        return *this;
      } // operator/=( MatrixCell const& )
    
      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        MatrixCell const& cell
      )
      {
        os << "[ " <<
          cell.m_match << ", " <<
          cell.m_insertion << ", " <<
          cell.m_deletion;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        os <<
          ", " << cell.m_deletionIn <<
          ", " << cell.m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        os << " ]";
        return os;
      } // friend operator<< ( basic_ostream, MatrixCell const & )
    
    }; // End inner class MatrixCell

    /**
     * For Profile-Profile alignment only (can't handle del-in and del-out)
     */
    template <typename ValueType>
    class BacktraceablePPMatrixCell : public MatrixCell<ValueType>
    {
    public:
      bool m_matchFromMatch;
      bool m_matchFromInsertion;
      bool m_matchFromDeletion;
      bool m_insertionFromMatch;
      bool m_insertionFromInsertion;
      bool m_insertionFromDeletion;
      bool m_deletionFromMatch;
      bool m_deletionFromInsertion;
      bool m_deletionFromDeletion;
//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//        ,
//        m_deletionInFromDeletionIn( false ),
//        m_deletionInFromMatch( false ),
//        m_matchFromDeletionIn( false ),
//        m_deletionOutFromDeletionOut( false ),
//        m_deletionOutFromMatch( false )
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
    
    public:
      BacktraceablePPMatrixCell () :
        MatrixCell<ValueType>(),
        m_matchFromMatch( false ),
        m_matchFromInsertion( false ),
        m_matchFromDeletion( false ),
        m_insertionFromMatch( false ),
        m_insertionFromInsertion( false ),
        m_insertionFromDeletion( false ),
        m_deletionFromMatch( false ),
        m_deletionFromInsertion( false ),
        m_deletionFromDeletion( false )
      {
        // Do nothing else.
      } // <init>()
    
      BacktraceablePPMatrixCell ( BacktraceablePPMatrixCell const & other_cell ) :
        MatrixCell<ValueType>( other_cell ),
        m_matchFromMatch( other_cell.m_matchFromMatch ),
        m_matchFromInsertion( other_cell.m_matchFromInsertion ),
        m_matchFromDeletion( other_cell.m_matchFromDeletion ),
        m_insertionFromMatch( other_cell.m_insertionFromMatch ),
        m_insertionFromInsertion( other_cell.m_insertionFromInsertion ),
        m_insertionFromDeletion( other_cell.m_insertionFromDeletion ),
        m_deletionFromMatch( other_cell.m_deletionFromMatch ),
        m_deletionFromInsertion( other_cell.m_deletionFromInsertion ),
        m_deletionFromDeletion( other_cell.m_deletionFromDeletion )
//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//        ,
//        m_deletionInFromDeletionIn( other_cell.m_deletionInFromDeletionIn ),
//        m_deletionInFromMatch( other_cell.m_deletionInFromMatch ),
//        m_matchFromDeletionIn( other_cell.m_matchFromDeletionIn ),
//        m_deletionOutFromDeletionOut( other_cell.m_deletionOutFromDeletionOut ),
//        m_deletionOutFromMatch( other_cell.m_deletionOutFromMatch )
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      {
        // Do nothing else.
      } // <init>( BacktraceablePPMatrixCell const & )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        BacktraceablePPMatrixCell const& cell
      )
      {
        os << "[ " <<
          cell.m_match <<
          ( cell.m_matchFromMatch ? " m" : "" ) <<
          ( cell.m_matchFromInsertion ? " i" : "" ) <<
          ( cell.m_matchFromDeletion ? " d" : "" ) <<
//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//          ( cell.m_matchFromDeletionIn ? " D" : "" ) <<
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          ", " <<
          cell.m_insertion <<
          ( cell.m_insertionFromMatch ? " m" : "" ) <<
          ( cell.m_insertionFromInsertion ? " i" : "" ) <<
          ( cell.m_insertionFromDeletion ? " d" : "" ) <<
          ", " <<
          cell.m_deletion <<
          ( cell.m_deletionFromMatch ? " m" : "" ) <<
          ( cell.m_deletionFromInsertion ? " i" : "" ) <<
          ( cell.m_deletionFromDeletion ? " d" : "" ) <<
//#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
//          ", " <<
//          cell.m_deletionIn <<
//          ( cell.m_deletionInFromMatch ? " m" : "" ) <<
//          ( cell.m_deletionInFromDeletionIn ? " D" : "" ) <<
//          ", " <<
//          cell.m_deletionOut <<
//          ( cell.m_deletionOutFromMatch ? " m" : "" ) <<
//          ( cell.m_deletionOutFromDeletionOut ? " D" : "" ) <<
//#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          " ]";
        return os;
      } // friend operator<< ( basic_ostream, BacktraceablePPMatrixCell const & )

    }; // End inner class BacktraceablePPMatrixCell    

    /**
     * The RabinerScalableMatrix class defines a SequentialAccessContainer to
     * access a set of matrices intended to be indexed first by row, then by
     * sequence, then by column.  For the row-by-row algorithms, this promotes
     * wise paging behavior (keeping in memory just the local rows).
     */
    template <typename ValueType>
    class RabinerScalableMatrix
    {
    public:
      class Row : public vector<MatrixCell<ValueType> >
      {
        public:
        /**
         * If parameters.useRabinerScaling is true, we scale the dynamic
         * programming matrix rows by matrixRowScaleFactor/rabinerInverseScalar.
         * Each row has its own scalar.  Actually to get back to what it would
         * be without scaling, you would need to multiply by
         * rabinerInverseScalar/matrixRowScaleFactor.
         *
         * This is that scalar.
         */
        ScoreType m_rabinerInverseScalar;

        /**
         * If parameters.useRabinerScaling is true, we scale the dynamic
         * programming matrix rows by matrixRowScaleFactor/rabinerInverseScalar.
         * Each row has its own scalar.  Actually to get back to what it would
         * be without scaling, you would need to multiply by
         * rabinerInverseScalar/matrixRowScaleFactor.
         *
         * This is the product of those scalars, up to and including this row
         * (for forward rows, they accumulate from the first row to this one,
         * and for the backward rows, they accumulate from the last row to this
         * one).
         */
        ScoreType m_rabinerCumulativeInverseScalar;

#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
        ValueType m_deletionIn;
        ValueType m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )

        Row () :
          m_rabinerInverseScalar( 1.0 ),
          m_rabinerCumulativeInverseScalar( 1.0 )
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          ,
          m_deletionIn( 0 ),
          m_deletionOut( 0 )
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        {
          // Do nothing.
        } // <init>()

        Row ( uint32_t row_length ) :
          vector<MatrixCell<ValueType> >( row_length ),
          m_rabinerInverseScalar( 1.0 ),
          m_rabinerCumulativeInverseScalar( 1.0 )
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          ,
          m_deletionIn( 0 ),
          m_deletionOut( 0 )
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        {
          // Do nothing else.
        } // <init>( uint32_t )

        Row ( Row const & other_row ) :
          m_rabinerInverseScalar( other_row.m_rabinerInverseScalar ),
          m_rabinerCumulativeInverseScalar( other_row.m_rabinerCumulativeInverseScalar )
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          ,
          m_deletionIn( other_row.m_deletionIn ),
          m_deletionOut( other_row.m_deletionOut )
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        {
          assign( other_row.begin(), other_row.end() );
        } // <init>( Row const & )

        // Note: this does not reset the Cells to 0...
        void
        reinitialize ( uint32_t row_length )
        {
          if( vector<MatrixCell<ValueType> >::size() != row_length ) {
            vector<MatrixCell<ValueType> >::resize( row_length );
          }
          // For efficiency, don't bother zeroing the values.
          //uint32_t last_cell = length() - 1;
          //uint32_t cell_i;
          //for( cell_i = 0; cell_i <= last_cell ; cell_i++ ) {
          //  ( *this )[ cell_i ].reinitialize();
          //}
          m_rabinerInverseScalar = 1.0;
          m_rabinerCumulativeInverseScalar = 1.0;
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          m_deletionIn = 0;
          m_deletionOut = 0;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        } // reinitialize ( uint32_t )

        uint32_t
        length () const
        {
          return vector<MatrixCell<ValueType> >::size();
        } // length() const

        /**
         * This is used (when useRabinerScaling is true) to copy the forward
         * matrix rabinerInverseScalar values to the backward matrices.
         *
         * Assumption: the argument has the same dimensions as this.
         */
        //void
        //copyRabinerInverseScalarsFrom (
        //  Row const & other_row
        //)
        //{
        //  m_rabinerInverseScalar = other_row.m_rabinerInverseScalar;
        //} // copyRabinerInverseScalarsFrom( Row const & )

        template<class CharT, class Traits>
        friend std::basic_ostream<CharT,Traits>&
        operator<< (
          std::basic_ostream<CharT,Traits>& os,
          typename RabinerScalableMatrix::Row const& row
        )
        {
          row.writeCells( os );
          os << ", inverseScalar: " << row.m_rabinerInverseScalar;
          os << ", cumulativeInverseScalar: " << row.m_rabinerCumulativeInverseScalar;
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          os << ", deletionIn: " << row.m_deletionIn;
          os << ", deletionOut: " << row.m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
          return os;
        } // friend operator<< ( basic_ostream, Row const & )

        template<class CharT, class Traits>
        void
        writeCells (
          std::basic_ostream<CharT,Traits>& os
        ) const
        {
          // TODO: Why won't this work?
          //copy( row.begin(), row.end(),
          //ostream_iterator< typename
          //RabinerScalableMatrix::Row >( os, ":" ) );
          // TODO: Why won't this work?
//           for( MatrixCell<ValueType>  & c = row.begin(); c < row.end(); c++ ) {
//             os << c << ", ";
//           }
//           os << row.end();
          uint32_t last_cell = length() - 1;
          uint32_t cell_i;
          for( cell_i = 0; cell_i <= last_cell ; cell_i++ ) {
            os << ( *this )[ cell_i ];
            if( cell_i != last_cell ) {
              os << ", ";
            }
          }
        } // writeCells( basic_ostream ) const

        /**
         * Save a subset of the columns in the given Row to this one.  Any
         * column with an index that divides evenly by store_every_Nth_column
         * will be stored, including row 0, and additionally the final column
         * will be stored (it may even be stored twice if it is also divisible
         * by store_every_Nth_column).  The del-in and del-out values will be
         * stored, too.
         *
         * @see restoreColumns(..)
         */
        void
        storeColumns (
          Row const & copy_columns_from,
          uint32_t const & store_every_Nth_column
        )
        {
          const uint32_t last_col = ( copy_columns_from.size() - 1 );
          register uint32_t col_i;
          uint32_t stored_col_i = 0;

          for( col_i = 0; col_i < last_col; col_i++ ) {
            if( ( col_i % store_every_Nth_column ) > 0 ) {
              continue;
            }
            ( *this )[ stored_col_i ] = copy_columns_from[ col_i ];
            stored_col_i += 1;
          } // End foreach col_i
          // And always store the final column in the last place.
          ( *this )[ stored_col_i ] = copy_columns_from[ last_col ];

#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          m_deletionIn = copy_columns_from.m_deletionIn;
          m_deletionOut = copy_columns_from.m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        } // storeColumns( Row const &, uint32_t const & )

        /**
         * Restore a subset of the columns into the given Row from this one.
         * Any column with an index that divides evenly by
         * store_every_Nth_column will be restored, including row 0, and
         * additionally the final column will be restored.  del-in and del-out
         * values will be restored, too.
         *
         * @see storeColumns(..)
         */
        void
        restoreColumns (
          uint32_t const & store_every_Nth_column,
          Row & copy_columns_to
        ) const
        {
          const uint32_t last_col = ( copy_columns_to.size() - 1 );
          register uint32_t col_i;
          uint32_t stored_col_i = 0;

          for( col_i = 0; col_i < last_col; col_i++ ) {
            if( ( col_i % store_every_Nth_column ) > 0 ) {
              continue;
            }
            copy_columns_to[ col_i ] = ( *this )[ stored_col_i ];
            stored_col_i += 1;
          } // End foreach col_i
          // And we always store the final column in the last place.
          copy_columns_to[ last_col ] = ( *this )[ stored_col_i ];
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          copy_columns_to.m_deletionIn = m_deletionIn;
          copy_columns_to.m_deletionOut = m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
        } // restoreColumns( uint32_t const &, Row & ) const

        Row &
        operator= ( Row const& other_row )
        {
          if( vector<MatrixCell<ValueType> >::size() != other_row.size() ) {
            vector<MatrixCell<ValueType> >::resize( other_row.size() );
          }
          uint32_t last_cell = length() - 1;
          uint32_t cell_i;
          for( cell_i = 0; cell_i <= last_cell ; cell_i++ ) {
            this->operator[]( cell_i ) = other_row[ cell_i ];
          }
          m_rabinerInverseScalar = other_row.m_rabinerInverseScalar;
          m_rabinerCumulativeInverseScalar = other_row.m_rabinerCumulativeInverseScalar;
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          m_deletionIn = other_row.m_deletionIn;
          m_deletionOut = other_row.m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
          return *this;
        } // operator=( Row const& )

        // Note that this does not account for the rabiner scalars, nor does it
        // alter them.
        Row &
        operator+= ( Row const& other_row )
        {
          //if( vector<MatrixCell<ValueType> >::size() != other_row.size() ) {
          //  // TODO: ?
          //  vector<MatrixCell<ValueType> >::resize( other_row.size() );
          //}
          assert( vector<MatrixCell<ValueType> >::size() == other_row.size() );
          uint32_t last_cell = length() - 1;
          uint32_t cell_i;
          for( cell_i = 0; cell_i <= last_cell ; cell_i++ ) {
            this->operator[]( cell_i ) += other_row[ cell_i ];
          }
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          m_deletionIn += other_row.m_deletionIn;
          m_deletionOut += other_row.m_deletionOut;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
          return *this;
        } // operator+=( Row const& )

        // Note that this does not account for the rabiner scalars, nor does it
        // alter them.
        Row &
        operator*= ( ValueType const& scalar )
        {
          uint32_t last_cell = length() - 1;
          uint32_t cell_i;
          for( cell_i = 0; cell_i <= last_cell ; cell_i++ ) {
            this->operator[]( cell_i ) *= scalar;
          }
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          m_deletionIn *= scalar;
          m_deletionOut *= scalar;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
          return *this;
        } // operator*=( Row const& )

        // Note that this does not account for the rabiner scalars, nor does it
        // alter them.
        Row &
        operator/= ( ValueType const& denominator )
        {
          uint32_t last_cell = length() - 1;
          uint32_t cell_i;
          for( cell_i = 0; cell_i <= last_cell ; cell_i++ ) {
            this->operator[]( cell_i ) /= denominator;
          }
#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
          m_deletionIn /= denominator;
          m_deletionOut /= denominator;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
          return *this;
        } // operator/=( Row const& )

      }; // End inner class DynamicProgramming::RabinerScalableMatrix::Row

      /**
       * A vector of Rows.  Its index is usually a sequence number.
       */
      class RowVector : public vector<Row>
      {
      public:
        RowVector () :
          vector<Row>()
        {
          // Do nothing else.
        } // <init>()

        RowVector (
          RowVector const& other_row_vector
        )
        {
          assign( other_row_vector.begin(), other_row_vector.end() );
        } // <init>( RowVector const& )

        //RowVector (
        //  vector<Sequence<SequenceResidueType> > const& sequences
        //) :
        //  vector<Row>( sequences.size() )
        //{
        //  // TODO: REMOVE
        //  //cout << "Initializing RowVector with " << sequences.size() << " elements." << endl;
        //  uint32_t last_seq = sequences.size() - 1;
        //  uint32_t seq_i;
        //  for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
        //    // TODO: REMOVE
        //    //cout << "Initializing Sequence " << seq_i << " of RowVector.  This one has length " << sequences[ seq_i ].length() << "." << endl;
        //    this->operator[]( seq_i ).reinitialize(
        //      sequences[ seq_i ].length() + 1
        //    );
        //  }
        //  
        //} // <init>( vector<Sequence<SequenceResidueType> > const& )

        template <typename SequenceResidueType>
        RowVector (
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count
        ) :
          vector<Row>( sequence_count )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
          // TODO: REMOVE
          //cout << "<init> Initializing RowVector with " << sequence_count << " elements." << endl;
          uint32_t last_seq = sequence_count - 1;
          uint32_t seq_i;
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
            //if( seq_i == 0 ) {
              // TODO: REMOVE
              //cout << "Initializing Sequence " << seq_i << " of RowVector.  This one has length " << sequences[ seq_i ].length() << "." << endl;
              //}
            this->operator[]( seq_i ).reinitialize(
              sequences[ seq_i ].length() + 1
            );
          }
          
        } // <init>( vector<Sequence<SequenceResidueType> > const&, uint32_t )

        template <typename SequenceResidueType>
        RowVector (
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count,
          uint32_t const & store_every_Nth_column
        ) :
          vector<Row>( sequence_count )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
          // TODO: REMOVE
          //cout << "<init> Initializing RowVector with " << sequence_count << " elements." << endl;
          uint32_t last_seq = sequence_count - 1;
          uint32_t seq_i;
          uint32_t columns_needed;
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
            //if( seq_i == 0 ) {
              // TODO: REMOVE
              //cout << "Initializing Sequence " << seq_i << " of RowVector.  This one has length " << sequences[ seq_i ].length() << "." << endl;
              //}
            columns_needed = sequences[ seq_i ].length() + 1;
            if(
              ( store_every_Nth_column != numeric_limits<uint32_t>::max() ) &&
              ( store_every_Nth_column != 0 )
            ) {
              columns_needed /= store_every_Nth_column;
              columns_needed += 1;
              // Plus an extra one for the last column, which we always store.
              columns_needed += 1;
            }
            this->operator[]( seq_i ).reinitialize(
              columns_needed
            );
          } // End foreach sequence...
          
        } // <init>( vector<Sequence<SequenceResidueType> > const&, uint32_t, uint32_t const & )

        //template <typename SequenceResidueType>
        //void reinitialize (
        //  vector<Sequence<SequenceResidueType> > const& sequences
        //)
        //{
        //  reinitialize( sequences, sequences.size() );
        //} // reinitialize( vector<Sequence<SequenceResidueType> > const & )

        template <typename SequenceResidueType>
        void
        reinitialize (
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t const & sequence_count
        )
        {
          reinitialize(
            sequences,
            sequence_count,
            0
          );
        } // reinitialize( vector<Sequence<SequenceResidueType> > const &, uint32_t const & )

        template <typename SequenceResidueType>
        void
        reinitialize (
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count,
          uint32_t const & store_every_Nth_column
        )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
          if( vector<Row>::size() != sequence_count ) {
            vector<Row>::resize( sequence_count );
          }
          
          // TODO: REMOVE
          //cout << "reInitializing RowVector with " << sequence_count << " elements." << endl;
          uint32_t last_seq = sequence_count - 1;
          uint32_t seq_i;
          uint32_t columns_needed;
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
            //if( seq_i == 0 ) {
              // TODO: REMOVE
              //cout << "Initializing Sequence " << seq_i << " of RowVector.  This one has length " << sequences[ seq_i ].length() << "." << endl;
              //}
            columns_needed = sequences[ seq_i ].length() + 1;
            if(
              ( store_every_Nth_column != numeric_limits<uint32_t>::max() ) &&
              ( store_every_Nth_column != 0 )
            ) {
              columns_needed /= store_every_Nth_column;
              columns_needed += 1;
              // Plus an extra one for the last column, which we always store.
              columns_needed += 1;
            }
            this->operator[]( seq_i ).reinitialize(
              columns_needed
            );
          } // End foreach sequence...
        } // reinitialize( vector<Sequence<SequenceResidueType> > const &, uint32_t, uint32_t const & )

        void
        reinitialize ()
        {
          reinitialize( 0, 0 );
        } // reinitialize()

        void
        reinitialize (
          uint32_t longest_sequence_length,
          uint32_t sequence_count
        )
        {
          if( vector<Row>::size() != sequence_count ) {
            vector<Row>::resize( sequence_count );
          }
          if( sequence_count == 0 ) {
            return;
          }
          
          // TODO: REMOVE
          //cout << "reInitializing RowVector with " << sequence_count << " elements, all of length " << longest_sequence_length << endl;
          uint32_t last_seq = sequence_count - 1;
          uint32_t seq_i;
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
            this->operator[]( seq_i ).reinitialize(
              longest_sequence_length + 1
            );
          }
        } // reinitialize( uint32_t, uint32_t )

        RowVector &
        operator= ( RowVector const& other_row_vector )
        {
          // TODO: REMOVE
          //cout << "operator= RowVector with " << vector<Row>::size() << " elements." << endl;
          //cout << "argument is a RowVector with " << other_row_vector.size() << " elements." << endl;
          if( vector<Row>::size() != other_row_vector.size() ) {
            vector<Row>::resize( other_row_vector.size() );
          }
          uint32_t last_seq = other_row_vector.size() - 1;
          uint32_t seq_i;
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
            // TODO: REMOVE
            //if( seq_i == 0 ) {
            //cout << "Copying Sequence " << seq_i << " of other RowVector.  This one has length " << other_row_vector[ seq_i ].length() << "." << endl;
            //}
            this->operator[]( seq_i ) = other_row_vector[ seq_i ];
          }
          return *this;
        } // operator=( RowVector const& )

        RowVector &
        operator+= ( RowVector const& other_row_vector )
        {
          if( vector<Row>::size() != other_row_vector.size() ) {
            // TODO: ?
            vector<Row>::resize( other_row_vector.size() );
          }
          uint32_t last_seq = other_row_vector.size() - 1;
          uint32_t seq_i;
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
            this->operator[]( seq_i ) += other_row_vector[ seq_i ];
          }
          return *this;
        } // operator+=( RowVector const& )

        RowVector &
        operator*= ( ValueType const& scalar )
        {
          uint32_t last_seq = vector<Row>::size() - 1;
          uint32_t seq_i;
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
            this->operator[]( seq_i ) *= scalar;
          }
          return *this;
        } // operator*=( ValueType const& )

        template<class CharT, class Traits>
        friend std::basic_ostream<CharT,Traits>&
        operator<< (
          std::basic_ostream<CharT,Traits>& os,
          typename RabinerScalableMatrix::RowVector const& row_vec )
        {
          uint32_t last_seq = row_vec.size() - 1;
          uint32_t seq_i;
          for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
            os << row_vec[ seq_i ] << endl;
          }
          return os;
        } // friend operator<< ( basic_ostream, RowVector const &)

      }; // End inner class DynamicProgramming::RabinerScalableMatrix::RowVector

      /**
       * A vector of RowVectors.
       */
      /*
       * TODO: If m_compressed is true, then there will be
       * some compression going on under the hood.  This compression absolutely
       * requires that the access be sequential, and in reverse order (from
       * profile.length() downto 0).
       *    NOTE: This compression is presently accomplished through anchor rows and cols in ProfileTrainer.hpp, which are SequentialAccessContainers with options for thinning (storing a subset of the rows or of the cols -- but it's not automagic.
       */
      class SequentialAccessContainer : public list<RowVector>
      //public vector<RowVector>
      {
        // TODO: implement compression

      public:
        typedef typename list<RowVector>::iterator iterator;
        typedef typename list<RowVector>::const_iterator const_iterator;
        typedef typename list<RowVector>::reverse_iterator reverse_iterator;
        typedef typename list<RowVector>::const_reverse_iterator const_reverse_iterator;

        uint32_t m_storeEveryNthColumn;
        uint32_t m_storeEveryNthRow;

        SequentialAccessContainer () :
          list<RowVector>(),
          m_storeEveryNthColumn( 0 ),
          m_storeEveryNthRow( 0 )
        {
          // Do nothing else.
        } // <init>()

        template <typename SequenceResidueType>
        SequentialAccessContainer (
          uint32_t const & profile_length,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count
        ) :
          list<RowVector>( profile_length + 1 ),
          m_storeEveryNthColumn( numeric_limits<double>::max() ),
          m_storeEveryNthRow( numeric_limits<double>::max() )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
          // Initialize each contained RowVector
          iterator iter;
          for( iter = this->begin(); iter != this->end(); iter++ ) {
            ( *iter ).reinitialize( sequences, sequence_count );
          } // End foreach row
        } // <init>( uint32_t const &, vector<Sequence<SequenceResidueType> > const&, uint32_t const & )

        template <typename SequenceResidueType>
        SequentialAccessContainer (
          uint32_t const & profile_length,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count,
          uint32_t const & store_every_Nth_column
        ) :
          list<RowVector>( profile_length + 1 ),
          m_storeEveryNthColumn( store_every_Nth_column ),
          m_storeEveryNthRow( numeric_limits<double>::max() )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
          // Initialize each contained RowVector
          iterator iter;
          for( iter = this->begin(); iter != this->end(); iter++ ) {
            ( *iter ).reinitialize(
              sequences,
              min( ( size_t )sequence_count, sequences.size() ),
              store_every_Nth_column
            );
          } // End foreach row
        } // <init>( uint32_t const &, vector<Sequence<SequenceResidueType> > const&, uint32_t const &, uint32_t const & )

        template <typename ProfileType, typename SequenceResidueType>
        SequentialAccessContainer (
          ProfileType const& profile,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count
        ) :
          list<RowVector>( profile.length() + 1 ),
          m_storeEveryNthColumn( numeric_limits<double>::max() ),
          m_storeEveryNthRow( numeric_limits<double>::max() )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
          // Initialize each contained RowVector
          iterator iter;
          for( iter = this->begin(); iter != this->end(); iter++ ) {
            ( *iter ).reinitialize( sequences, sequence_count );
          } // End foreach row
        } // <init>( ProfileType&, vector<Sequence<SequenceResidueType> > const&, uint32_t const& )

        template <typename ProfileType, typename SequenceResidueType>
        SequentialAccessContainer (
          ProfileType const& profile,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count,
          uint32_t const & store_every_Nth_column
        ) :
          list<RowVector>( profile.length() + 1 ),
          m_storeEveryNthColumn( store_every_Nth_column ),
          m_storeEveryNthRow( numeric_limits<double>::max() )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
          // Initialize each contained RowVector
          iterator iter;
          for( iter = this->begin(); iter != this->end(); iter++ ) {
            ( *iter ).reinitialize(
              sequences,
              sequence_count,
              store_every_Nth_column
            );
          } // End foreach row
        } // <init>( ProfileType&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, uint32_t const& )

        void
        reinitialize ()
        {
          this->resize( 0 );
          m_storeEveryNthColumn = 0;
          m_storeEveryNthRow = 0;
        } // reinitialize()

        template <typename SequenceResidueType>
        void
        reinitialize (
          uint32_t const & profile_length,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t const & sequence_count
        )
        {
          reinitialize(
            profile_length,
            sequences,
            sequence_count,
            numeric_limits<uint32_t>::max(),
            numeric_limits<uint32_t>::max()
          );
        } // reinitialize( uint32_t const &, vector<Sequence<SequenceResidueType> > const&, uint32_t const & )

        template <typename SequenceResidueType>
        void
        reinitialize (
          uint32_t const & profile_length,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count,
          uint32_t const & store_every_Nth_column,
          uint32_t const & store_every_Nth_row
        )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );

          uint32_t length = profile_length + 1;
          m_storeEveryNthRow = store_every_Nth_row;
          m_storeEveryNthColumn = store_every_Nth_column;

          if(
            ( store_every_Nth_row != 0 ) &&
            ( store_every_Nth_row != numeric_limits<uint32_t>::max() )
          ) {
            // For example, if length is 10, and store_every_Nth_row is 2, we
            // want 5 anchor rows: 0, 2, 4, 6, and 8.  If length is 11, we want
            // a 6th at 10.
            length -= 1; // now it is the index of the last row (eg. 9 or 10 in the above examples)
            length /= store_every_Nth_row; // (eg. 4 or 5)
            length += 1; // (eg. 5 or 6)
          }

          // Special: if they are both 0, don't store anything at all.
          if(
            ( store_every_Nth_row == 0 ) &&
            ( store_every_Nth_column == 0 )
          ) {
            length = 0;
          }
          
          if( list<RowVector>::size() != length ) {
            list<RowVector>::resize( length );
          }

          // Initialize each contained RowVector to have rows of varying lengths
          if( length > 0 ) {
            size_t size = sequence_count;
            // Initialize each contained RowVector
            iterator iter;
            for( iter = this->begin(); iter != this->end(); iter++ ) {
              ( *iter ).reinitialize(
                sequences,
                size,
                store_every_Nth_column
              );
            } // End foreach row
          } // End if length > 0
        } // reinitialize( uint32_t const &, vector<Sequence<SequenceResidueType> > const&, uint32_t const &, uint32_t const & )

        template <typename SequenceResidueType>
        void
        reinitialize (
          uint32_t const & profile_length,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t sequence_count,
          uint32_t const & longest_sequence_length
        )
        {
          sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );

          if( list<RowVector>::size() != ( profile_length + 1 ) ) {
            list<RowVector>::resize( profile_length + 1 );
          }

          // Initialize each contained RowVector to have rows of length longest_sequence_length
          size_t size = sequence_count;
          iterator iter;
          for( iter = this->begin(); iter != this->end(); iter++ ) {
            ( *iter ).reinitialize( longest_sequence_length, size );
          } // End foreach row
        } // reinitialize( uint32_t const &, vector<Sequence<SequenceResidueType> > const&, uint32_t const &, uint32_t const & )

        template <typename ProfileType, typename SequenceResidueType>
        void
        reinitialize (
          ProfileType const& profile,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t const & sequence_count
        )
        {
          reinitialize(
            profile.length(),
            sequences,
            sequence_count
          );
        } // reinitialize( ProfileType const &, template<Sequence<SequenceResidueType>> const&, uint32_t const & )

        template <typename ProfileType, typename SequenceResidueType>
        void
        reinitialize (
          ProfileType const& profile,
          vector<Sequence<SequenceResidueType> > const& sequences,
          uint32_t const & sequence_count,
          uint32_t const & row_lengths
        )
        {
          reinitialize(
            profile.length(),
            sequences,
            sequence_count,
            row_lengths
          );
        } // reinitialize( ProfileType const &, template<Sequence<SequenceResidueType>> const&, uint32_t const &, uint32_t const & )

        /**
         * This is used (when useRabinerScaling is true) to copy the forward
         * matrix rabinerInverseScalar values to the backward matrices.  In
         * addition to copying the rabinerInverseScalar values, it also
         * accumulates their products in the corresponding
         * rabinerCumulativeInverseScalar variables.  It does so from the last
         * row, backwards, of course.
         *
         * Assumption: the argument has the same dimensions as this.
         */
        void
        copyRabinerInverseScalarsFrom (
          SequentialAccessContainer const & other_matrices
        )
        {
          assert( other_matrices.size() == this->size() );

          reverse_iterator iter = this->rbegin();
          reverse_iterator iter_prev = this->rbegin();
          const_reverse_iterator other_iter = other_matrices.rbegin();

          uint32_t last_seq = other_matrices.begin()->size() - 1;
          assert( last_seq == ( this->begin()->size() - 1 ) );

          uint32_t seq_i;
          for(
            ;//iter = this->rbegin(), other_iter = other_matrices.rbegin();
            ( iter != this->rend() ); // redundant: && ( other_iter < other_matrices.rend() ) );
            iter++, other_iter++
          ) {
            for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
              ( *iter )[ seq_i ].m_rabinerInverseScalar =
                ( *other_iter )[ seq_i ].m_rabinerInverseScalar;
              // Also update the rabinerCumulativeInverseScalars...
              if( iter == this->rbegin() ) {
                ( *iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                  ( *iter )[ seq_i ].m_rabinerInverseScalar;
              } else {
                ( *iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                  ( // Remember here that -1 refers to the next row, since we're iterating in reverse.
                    ( *( iter_prev ) )[ seq_i ].m_rabinerCumulativeInverseScalar *
                    ( *( iter ) )[ seq_i ].m_rabinerInverseScalar
                  );
              }
            } // End foreach seq_i
            if( iter != this->rbegin() ) {
              iter_prev++;
            }
          } // End foreach row, in reverse order
        } // copyRabinerInverseScalarsFrom( SequentialAccessContainer const & )

        template<class CharT, class Traits>
        friend std::basic_ostream<CharT,Traits>&
        operator<< (
          std::basic_ostream<CharT,Traits>& os,
          typename RabinerScalableMatrix::SequentialAccessContainer const& matrices )
        {
          const_iterator iter;
          uint32_t last_seq = matrices.begin()->size() - 1;
          uint32_t seq_i;
          for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
            os << "Sequence " << seq_i << ":" << endl;
            for( iter = matrices.begin(); iter != matrices.end(); iter++ ) {
              os << ( *iter )[ seq_i ] << endl;
            }
            os << endl;
          }
          return os;
        } // friend operator<< ( basic_ostream, SequentialAccessContainer const & )

     }; // End inner class DynamicProgramming::RabinerScalableMatrix::SequentialAccessContainer
    }; // End inner class DynamicProgramming::RabinerScalableMatrix

    typedef RabinerScalableMatrix<MatrixValueType> Matrix;

    /**
     * A MultipleAlignment stores one path per sequence for a vector of
     * Sequences and a profile model.  Typically it is the most
     * likely path, calculated via the viterbi algorithm.  Internally we also
     * store sufficient information to display these paths as a pileup multiple
     * alignment.
     */
    template <typename ProfileType, typename SequenceResidueType>
    class MultipleAlignment
    {
    public:
      ProfileType const * m_profile;
      vector<Sequence<SequenceResidueType> > const * m_sequences;
      /**
       * The number of sequences to use (the index one greater than that of the
       * last one to be used; must be <= m_sequences.size().
       */
      uint32_t m_sequence_count;

    public:


      // A vector of length m_sequence_count that contains each number in the
      // range 0 .. ( m_sequence_count - 1 ) exactly once, specifying a
      // reordering of the sequences.
      vector<uint32_t > m_sequenceOrder;

      /**
       * m_sequenceAdvances is a vector of int-vectors, one int-vector per
       * sequence.  Each int-vector has length four greater than that of the
       * profile, and the integer at index i indicates the number of sequence
       * positions advanced at profile position (i-2) by the most probable path
       * in generating the sequence.  So if a profile position's integer is 0,
       * that profile position was a deletion in the path generating the
       * sequence.  If it is 1, then it was a match.  Any number n greater than
       * 1 indicates that there were n-1 insertions following a match at that
       * position.  If the profile has length M, index 0 indicates the number
       * of pre-align insertions, and the index M-1 indicates the number of
       * post-align insertions.  Index 0 indicates the number of deletion-ins,
       * and index M-2 indicates the number of deletion-outs.  Note that the
       * value at index M-3 must be either 0 or 1 (since post-align insertions
       * are stored in the final position).
       */
      vector<vector<uint32_t> > m_sequenceAdvances;

      //MultipleAlignment (
      //  ProfileType const& profile,
      //  vector<Sequence<SequenceResidueType> > const& sequences
      //) :
      //  m_profile( profile ),
      //  m_sequences( sequences ),
      //  m_sequence_count( sequences.size() ),
      //  m_sequenceAdvances( sequences.size() )
      //{
      //  // Also initialize the m_sequenceAdvances rows.
      //  uint32_t last_seq = sequences.size() - 1;
      //  uint32_t seq_i;
      //  uint32_t profile_length = profile.length();
      //  for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
      //    m_sequenceAdvances[ seq_i ].resize( profile_length + 4 );
      //  }
      //} // <init>( ProfileType const&, vector<Sequence<SequenceResidueType> > )

      MultipleAlignment () :
        m_profile( NULL ),
        m_sequences( NULL ),
        m_sequence_count( 0 ),
        m_sequenceOrder(),
        m_sequenceAdvances()
      {
      } // <init>()

      /**
       * Only use the first sequence_count sequences of the sequences
       * vector.
       */
      MultipleAlignment (
        ProfileType const * const profile,
        vector<Sequence<SequenceResidueType> > const * const sequences,
        uint32_t const& arg_sequence_count_arg
      ) :
        m_profile( profile ),
        m_sequences( sequences ),
        m_sequence_count( ( arg_sequence_count_arg == 0 ) ? sequences->size() : min( static_cast<size_t>( arg_sequence_count_arg ), sequences->size() ) ),
        m_sequenceOrder( m_sequence_count ),
        m_sequenceAdvances( m_sequence_count )
      {
        // Also initialize the m_sequenceAdvances rows.
        uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        uint32_t profile_length = m_profile->length();
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          m_sequenceAdvances[ seq_i ].resize( profile_length + 4 );
          m_sequenceOrder[ seq_i ] = seq_i;
        }
      } // <init>( ProfileType const * const, vector<Sequence<SequenceResidueType> > const * const, uint32_t const & )

      /**
       * Only use the first sequence_count sequences of the sequences
       * vector.
       */
      void
      reinitialize (
        ProfileType const * const profile,
        vector<Sequence<SequenceResidueType> > const * const sequences,
        uint32_t const& sequence_count
      )
      {
        m_profile = profile;
        m_sequences = sequences;
        m_sequence_count = ( ( sequence_count == 0 ) ? sequences->size() : min( static_cast<size_t>( sequence_count ), sequences->size() ) );
        m_sequenceOrder.resize( sequence_count );
        m_sequenceAdvances.resize( m_sequence_count );

        if( m_sequence_count == 0 ) {
          return;
        }
        // Also initialize the m_sequenceAdvances rows.
        uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        uint32_t profile_length = m_profile->length();
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          m_sequenceOrder[ seq_i ] = seq_i;
          m_sequenceAdvances[ seq_i ].resize( profile_length + 4 );
        }
      } // reinitialize( ProfileType const&, vector<Sequence<SequenceResidueType> >, uint32_t const & )

      /**
       * Calculate and return the score of the given path.
       */
      static
      ScoreType
      scorePath (
        ProfileType const & profile,
        Sequence<SequenceResidueType> const & sequence,
        vector<uint32_t> const & sequence_advances
      ) {
        static const bool do_extra_debugging = false;//true;

        assert( sequence_advances.size() == ( profile.length() + 4 ) );
        uint32_t last_seq_adv_i = ( profile.length() + 3 );
        ScoreType score = 1.0;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> > scaled_match_distributions( profile.length() - 1 );
        for( uint32_t pos_i = 0; pos_i < ( profile.length() - 1 ); pos_i++ ) {
          profile.createScaledMatchDistributionForPosition(
            pos_i,
            scaled_match_distributions[ pos_i ]
          );
        } // End foreach pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        register uint32_t i, seq_pos_i, seq_adv_i;
        seq_pos_i = 0;
        seq_adv_i = 0;
        for( i = 0; i < sequence_advances[ 0 ]; i++ ) {
          // Pre-aligns.
          score *=
            profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
          score *=
            getPreAlignEmissionProbability(
              profile,
              sequence,
              seq_pos_i
            );
          ++seq_pos_i;
        } // End foreach pre-align insertion i
        score *= profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ];
        if( do_extra_debugging ) {
          cout << "ma.calculateScore(..): score after pre-aligns: " << score << endl;
        }

        seq_adv_i = 2;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        if( do_extra_debugging ) {
          cout << "ma.calculateScore(..): deletion-in count: " << sequence_advances[ 1 ] << endl;
        }
       if( sequence_advances[ 1 ] > 0 ) {
         score *=
           profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletionIn ];
         // There will also be 0s at these positions..
         assert( sequence_advances[ seq_adv_i ] == 0 );
         // This is also the deletion of the corresponding profile position..
         seq_adv_i++;
         for( i = 1; i < sequence_advances[ 1 ]; i++ ) {
           // There will also be 0s at these positions..
           assert( sequence_advances[ seq_adv_i ] == 0 );
           // Del-ins
           score *=
             profile[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toDeletionIn ];
           // This is also the deletion of the corresponding profile position..
           seq_adv_i++;
         } // End foreach del-in
         if( do_extra_debugging ) {
           cout << "ma.calculateScore(..): score after deletion-ins: " << score << endl;
         }
         // The position after the end of the del-ins must be a Match.
         assert( sequence_advances[ seq_adv_i ] > 0 );
         score *=
           profile[ Transition::fromDeletionIn ][ TransitionFromDeletionIn::toMatch ];
         score *=
           getEmissionProbability(
             profile[ seq_adv_i - 2 ],
             sequence,
             seq_pos_i
           );
         ++seq_pos_i;
         for( i = 1; i < sequence_advances[ seq_adv_i ]; i++ ) {
           // Insertions
           if( i == 1 ) {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
             score *=
               scaled_match_distributions[ seq_adv_i - 2 ][ TransitionFromMatch::toInsertion ];
#else
             score *=
               profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
             // Also multiply in the closing prob now.
             score *=
               profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];

           } else {
             score *=
               profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
           }
           score *=
             getInsertionEmissionProbability(
               profile[ seq_adv_i - 2 ],
               sequence,
               seq_pos_i
             );
           ++seq_pos_i;
         } // End foreach insertion i
         if( do_extra_debugging ) {
           cout << "ma.calculateScore(..): score (after deletion-ins) after seq_adv_i==" << seq_adv_i << ": " << score << endl;
         }
         seq_adv_i++;
       } // End if there's any deletion ins

#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        const uint32_t deletion_outs = sequence_advances[ last_seq_adv_i - 1 ];
        if( do_extra_debugging ) {
          cout << "ma.calculateScore(..): deletion-out count: " << deletion_outs << endl;
        }
        for( ; seq_adv_i <= ( last_seq_adv_i - 2 - deletion_outs ); seq_adv_i++ ) {
          if( sequence_advances[ seq_adv_i ] == 0 ) {
            if( seq_adv_i == 2 ) { // Corresponding to first profile pos.
              score *= profile[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
            } else if( sequence_advances[ seq_adv_i - 1 ] == 0 ) {
              score *= profile[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
            } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              score *=
                scaled_match_distributions[ seq_adv_i - 3 ][ TransitionFromMatch::toDeletion ];
#else
              score *=
                profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            }
          } else { // if sequence_advances[ seq_adv_i ] == 0 .. else ..
            if( seq_adv_i == 2 ) { // Corresponding to first profile pos.
              score *= profile[ Transition::fromBegin ][ TransitionFromBegin::toMatch ];
            } else if( sequence_advances[ seq_adv_i - 1 ] == 0 ) {
              score *= profile[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ];
            } else if( sequence_advances[ seq_adv_i - 1 ] == 1 ) {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              score *=
                scaled_match_distributions[ seq_adv_i - 3 ][ TransitionFromMatch::toMatch ];
#else
              score *=
                profile[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            } else {
              // Already incorporated.
              //score *= profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
            }
            score *=
              getEmissionProbability(
                profile[ seq_adv_i - 2 ],
                sequence,
                seq_pos_i
              );
            ++seq_pos_i;

            for( i = 1; i < sequence_advances[ seq_adv_i ]; i++ ) {
              if( do_extra_debugging ) {
                cout << "ma.calculateScore(..): score before insertions at seq_adv_i==" << seq_adv_i << ": " << score << endl;
              }
              // Insertions
              if( i == 1 ) {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                score *=
                  scaled_match_distributions[ seq_adv_i - 2 ][ TransitionFromMatch::toInsertion ];
#else
                score *=
                  profile[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                // Also multiply in the closing prob now.
                score *=
                  profile[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];

              } else {
                score *=
                  profile[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
              }
              score *=
                getInsertionEmissionProbability(
                  profile[ seq_adv_i - 2 ],
                  sequence,
                  seq_pos_i
                );
              ++seq_pos_i;
              if( do_extra_debugging ) {
                cout << "ma.calculateScore(..): score after insertion " << i << " at seq_adv_i==" << seq_adv_i << ": " << score << endl;
              }
            } // End foreach insertion i
          } // End if sequence_advances[ seq_adv_i ] == 0 .. else ..
          if( do_extra_debugging ) {
            cout << "ma.calculateScore(..): score after seq_adv_i==" << seq_adv_i << ": " << score << endl;
          }
        } // End foreach internal seq_adv_i

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
       if( deletion_outs > 0 ) {
         assert( seq_adv_i == ( last_seq_adv_i - 1 - deletion_outs ) );
         // Note we don't use the scaled version here.  We want the del-out
         // OPEN prob.
         score *=
           profile[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ];
         // There will also be 0s at these positions..
         assert( sequence_advances[ seq_adv_i ] == 0 );
         for( i = 1; i < sequence_advances[ last_seq_adv_i - 1 ]; i++ ) {
           // Del-outs
           score *=
             profile[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ];
           // This is also the deletion of the corresponding profile position..
           seq_adv_i++;
           // There will also be 0s at these positions..
           assert( sequence_advances[ seq_adv_i ] == 0 );
         } // End foreach del-out
         score *=
           profile[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ];
         // The position after the end of the del-ins ought to be beyond all
         // of the Profile positions.
         seq_adv_i++;

        if( do_extra_debugging ) {
          cout << "ma.calculateScore(..): score after deletion-outs: " << score << endl;
        }
       } // End if there's any deletion outs
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        assert( seq_adv_i == ( last_seq_adv_i - 1 ) );
        // Now skip beyond the del-outs pos
        seq_adv_i++;

#ifdef USE_END_DISTRIBUTION
        score *=
          profile[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ];
#endif // USE_END_DISTRIBUTION
        for( i = 0; i < sequence_advances[ last_seq_adv_i ]; i++ ) {
          // Post-aligns.
          score *=
            profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
          score *=
            getPostAlignEmissionProbability(
              profile,
              sequence,
              seq_pos_i
            );
          ++seq_pos_i;
        } // End foreach post-align insertion i
        score *= profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ];

        if( do_extra_debugging ) {
          cout << "ma.calculateScore(..): score after post-aligns: " << score << endl;
        }
        assert( seq_pos_i == sequence.length() );

        return score;
      } // static scorePath( ProfileType const &, Sequence<SequenceResidueType> const &, vector<uint32_t> const & )

      ScoreType
      calculateScore () const
      {
        return calculateScore( *m_profile );
      } // calculateScore () const

      ScoreType
      calculateScore (
        ProfileType const & profile
      ) const
      {
        ScoreType score = 1.0;
        
        for( uint32_t seq_i = 0; seq_i < m_sequence_count; seq_i++ ) {
          // TODO: REMOVE
          //cout << "calculateScore: seq_i is " << seq_i << ".  m_sequenceAdvances[ seq_i ] is\n( " << m_sequenceAdvances[ seq_i ][ 0 ];
          //for( int i = 0; i < m_sequenceAdvances[ seq_i ].size(); i++ ) {
          //  cout << ", " << m_sequenceAdvances[ seq_i ][ i ];
          //}
          //cout << " )\n";

          score *=
            scorePath(
              profile,
              ( *m_sequences )[ seq_i ],
              m_sequenceAdvances[ seq_i ]
            );
        }
        return score;
      } // calculateScore ( ProfileType const & ) const

      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toStream ( std::basic_ostream<CharT,Traits>& os ) const
      {
        uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        uint32_t last_align_pos = m_profile->length() + 3;
        uint32_t align_pos_i;
        for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
          os << "Sequence " << m_sequenceOrder[ seq_i ] << ":" << endl;
          os << "( ";
          for( align_pos_i = 0; align_pos_i <= last_align_pos ; align_pos_i++ ) {
            if( align_pos_i != 0 ) {
              os << ", ";
            }
            os << m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ align_pos_i ];
          }
          os << " )" << endl;
        }
        return os;
      } // toStream( ostream& )

      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toPairwiseStream (
        std::basic_ostream<CharT,Traits>& os,
        vector<string> const* seq_labels_ptr = NULL,
        uint32_t const wrap_column = 80
      ) const
      {
        // Show flanking deletions or no?
        static const bool suppress_flanking_deletions = false;

        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_ref_pos = m_profile->length();

        uint32_t block_width;

        vector<string> reference_strings;
        vector<string> alignment_strings;
        vector<string> instance_strings;

        uint32_t alignment_i;
        uint32_t alignment_length = 0;

        const uint32_t last_seq_adv_i = m_profile->length() + 3;
        uint32_t ref_i; // in coords of 1..last_ref_pos
        uint32_t ref_i_width; // space req'd to print ref_i
        uint32_t inst_i; // in coords of 1..last_inst_pos
        uint32_t inst_i_width; // space req'd to print inst_i
        uint32_t last_inst_pos; // in coords of 1..last_inst_pos
        uint32_t pos_width; // max space req'd to print inst_i or ref_i
        uint32_t padding_i;
        uint32_t substring_i;
        size_t find_pos;

        string substring;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          if( seq_labels_ptr == NULL ) {
            os << "> Sequence " << m_sequenceOrder[ seq_i ] << endl;
          } else {
            os << "> " << ( *seq_labels_ptr )[ m_sequenceOrder[ seq_i ] ] << endl;
          }
          os << endl;

          last_inst_pos = ( *m_sequences )[ m_sequenceOrder[ seq_i ] ].length();
          // TODO: REMOVE
          //cout << "last_ref_pos is " << last_ref_pos << "; last_inst_pos is " << last_inst_pos << endl;

          pos_width =
            max(
              ( log10( last_ref_pos ) + 1 ),
              ( log10( last_inst_pos ) + 1 )
            );
          block_width = // "Sequence" is 8 long
            max( 0, ( int )wrap_column  -  8 - 1  -  ( int )pos_width - 1  -  1 - ( int )pos_width );
          if( block_width == 0 ) {
            // TODO: Prevent this from happening!
            os << "ERROR: You must use a larger wrap_column argument; " << wrap_column << " is not large enough to accommodate the sequence numbers on either side of the output." << endl;
            return os;
          }
          // TODO: REMOVE
          //cout << "pos_width is " << pos_width << "; block_width is " << block_width << "; wrap_column is " << wrap_column << endl;

          if( reference_strings.size() == 0 ) {
            getPairwiseStrings(
              reference_strings,
              alignment_strings,
              instance_strings
            );
          }
          alignment_length =
            reference_strings[ m_sequenceOrder[ seq_i ] ].length(); // all vecs same len
          // TODO: REMOVE
          //cout << "alignment_length is " << alignment_length << endl;

          // Blocks

          alignment_i = 0;
          // Start ref count after deletion-ins and other starting deletions,
          // unless there's pre-align insertions.
          ref_i = 1;
          if(
            suppress_flanking_deletions &&
            ( m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ 0 ] == 0 )
          ) {
            ref_i = 1 + m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ 1 ];
            // Also count any other flanking deletions
            for( ; ref_i <= last_ref_pos ; ref_i++ ) {
              if( m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ ref_i + 1 ] != 0 ) {
                break;
              }
            }
            alignment_i = ref_i - 1;
          } // End if suppress_flanking_deletions and there's no pre-align insertions ..
          // Truncate the alignment if there's any deletion-outs or other
          // post-align deletions (unless there's post-align insertions)
          // TODO: REMOVE
          //cout << "There are " << m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ last_seq_adv_i ] << " post-align insertions." << endl;
          if(
            suppress_flanking_deletions &&
            ( m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ last_seq_adv_i ] == 0 )
          ) {
            // TODO: REMOVE
            //cout << "Truncating aligment length by the " << m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ last_seq_adv_i - 1 ] << " del-outs." << endl;
            alignment_length -= 
              m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ last_seq_adv_i - 1 ];
            // Also count any other flanking deletions
            for( uint32_t tmp_i = last_seq_adv_i - 2 - m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ last_seq_adv_i - 1 ]; tmp_i >= 2 ; tmp_i-- ) {
              if( m_sequenceAdvances[ m_sequenceOrder[ seq_i ] ][ tmp_i ] != 0 ) {
                break;
              }
              alignment_length -= 1;
            }
          } // End if suppress_flanking_deletions and there's no pos-align insertions ..
          inst_i = 1;
          for( ; alignment_i < alignment_length; alignment_i += block_width ) {
            // TODO: REMOVE
            //cout << "alignment_i is " << alignment_i << endl;

            // On top is the Profile line.
            substring =
              reference_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                min( block_width, ( alignment_length - alignment_i ) )
              );
            os << "Model    " << ref_i;
            ref_i_width = log10( ref_i ) + 1;
            for( padding_i = 0; padding_i < ( pos_width - ref_i_width ); padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            // Count ref advances
            ref_i += substring.length();
            // (Note ref_i--)
            for( substring_i = 0; substring_i < substring.length(); ref_i-- ) {
              find_pos = substring.find( '-', substring_i );
              if( find_pos == string::npos ) {
                break;
              }
              substring_i = ( find_pos + 1 );
              //cout << "substring_i is " << substring_i;
              //cout << "ref_i is " << ref_i;
            } // End foreach substring_i
            os << " " << ( ref_i - 1 );
            os << endl;

            // In the middle is the Alignment line.
            substring =
              alignment_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                min( block_width, ( alignment_length - alignment_i ) )
              );
            os << "         "; // 9 long
            for( padding_i = 0; padding_i < pos_width; padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            os << endl;
    
            // On bottom is the Sequence line.
            substring =
              instance_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                min( block_width, ( alignment_length - alignment_i ) )
              );
            os << "Sequence " << inst_i;
            inst_i_width = log10( inst_i ) + 1;
            for( padding_i = 0; padding_i < ( pos_width - inst_i_width ); padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            // Count inst advances
            inst_i += substring.length();
            // (Note inst_i--)
            for( substring_i = 0; substring_i < substring.length(); inst_i-- ) {
              find_pos = substring.find_first_of( "-. ", substring_i );
              if( find_pos == string::npos ) {
                break;
              }
              substring_i = ( find_pos + 1 );
              //cout << "substring_i is " << substring_i;
              //cout << "inst_i is " << inst_i;
            } // End foreach substring_i
            os << " " << ( inst_i - 1 );

            os << endl;
            // An extra line between blocks:
            os << endl;
          } // End foreach block
  
          os << endl;
          // An extra line between sequences:
          os << endl;
        } // End foreach sequence

        return os;
      } // toPairwiseStream( ostream&, ... )

      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toAlignedFastaStream (
        std::basic_ostream<CharT,Traits>& os,
        vector<string> const* seq_labels_ptr = NULL,
        uint32_t const wrap_column = 60,
        bool const & include_model_seq = false,
        string const & model_seq_label = "Model"
      ) const
      {
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;

        string reference_string;
        vector<string> instance_strings;

        uint32_t alignment_i;
        uint32_t alignment_length = 0;

        if( instance_strings.size() == 0 ) {
          reference_string =
            getPileupStrings(
              instance_strings
            );
        }
        alignment_length =
          reference_string.length();

        if( include_model_seq ) {
          // The model sequence
          os << ">" << model_seq_label << endl;
          for( alignment_i = 0; alignment_i < alignment_length; alignment_i += wrap_column ) {
            os <<
              reference_string.substr(
                alignment_i,
                min( wrap_column, ( alignment_length - alignment_i ) )
              );
            os << endl;
          } // End foreach line ...
          os << endl; // An extra line between seqs
        } // End if include_model_seq

        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          // The model sequence
          if( seq_labels_ptr == NULL ) {
            os << ">" << "Seq" << m_sequenceOrder[ seq_i ] << endl;
          } else {
            os << ">" << ( *seq_labels_ptr )[ m_sequenceOrder[ seq_i ] ] << endl;
          }
          for( alignment_i = 0; alignment_i < alignment_length; alignment_i += wrap_column ) {
            os <<
              instance_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                min( wrap_column, ( alignment_length - alignment_i ) )
              );
            os << endl;
          } // End foreach line ...
          if( seq_i != last_seq ) {
            os << endl; // An extra line between seqs
          }
        } // End foreach seq_i

        return os;
      } // toAlignedFastaStream( ostream&, ... )

      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toPileupStream (
        std::basic_ostream<CharT,Traits>& os,
        vector<string> const* seq_labels_ptr = NULL,
        uint32_t const wrap_column = 80,
        bool const aligned_fasta_format = false
      ) const
      {
        //static const double log10 = log( 10 );

        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_ref_pos = m_profile->length();

        uint32_t block_width;

        string reference_string;
        vector<string> instance_strings;

        uint32_t alignment_i;
        uint32_t alignment_length = 0;

        uint32_t ref_i; // in coords of 1..last_ref_pos
        uint32_t ref_i_width; // space req'd to print ref_i
        vector<uint32_t> inst_i( m_sequence_count, 1 ); // in coords of 1..last_inst_pos
        uint32_t inst_i_width; // space req'd to print inst_i
        uint32_t last_inst_pos; // in coords of 1..last_inst_pos
        uint32_t pos_width; // max space req'd to print inst_i or ref_i
        //uint32_t seq_i_width; // space req'd to print seq_i
        uint32_t padding_i;
        uint32_t substring_i;
        size_t find_pos;

        string substring;

        string ref_label = "=RF";
        const uint32_t ref_label_width = 3;
        vector<string> local_seq_labels( m_sequence_count );
        uint32_t seq_label_width;
        uint32_t max_seq_label_width = 0;
        // TODO: Make this a parameter:
        static const uint32_t hard_max_seq_label_width = 15;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          //if( seq_labels_ptr == NULL ) {
          local_seq_labels[ seq_i ] = "Seq ";
          local_seq_labels[ seq_i ].append(
            lexical_cast<string>( seq_i )
          ); // use original ordering
          seq_label_width =
            local_seq_labels[ seq_i ].length();
          if( seq_labels_ptr != NULL ) {
            if( seq_label_width < ( hard_max_seq_label_width - 2 ) ) {
              local_seq_labels[ seq_i ].append( ": " );
              local_seq_labels[ seq_i ].append(
                ( *seq_labels_ptr )[ seq_i ],
                0,
                ( hard_max_seq_label_width - seq_label_width )
              );
            }
            seq_label_width =
              local_seq_labels[ seq_i ].length();
          } // End if seq_labels_ptr != NULL
          if( seq_label_width > max_seq_label_width ) {
            max_seq_label_width = seq_label_width;
          }
        } // End foreach seq_i
        uint32_t max_label_width =
          max( max_seq_label_width, ref_label_width ); 

        uint32_t largest_last_inst_pos = 0;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          last_inst_pos = ( *m_sequences )[ seq_i ].length();
          if( last_inst_pos > largest_last_inst_pos ) {
            largest_last_inst_pos = last_inst_pos;
          }
        } // End foreach seq_i..
        uint32_t largest_inst_i_width =
          log10( largest_last_inst_pos ) + 1;
        //( uint32_t )( log( largest_last_inst_pos ) / log10 ) + 1;
        uint32_t max_ref_i_width =
          log10( last_ref_pos ) + 1;
        //( uint32_t )( log( last_ref_pos ) / log10 ) + 1;
        uint32_t max_seq_i_width =
          ( ( last_seq == 0 ) ? 1 : ( log10( last_seq ) + 1 ) );
        //( ( last_seq == 0 ) ? 1 : ( ( uint32_t )( log( last_seq ) / log10 ) + 1 ) );

        // TODO: REMOVE
        //cout << "last_ref_pos is " << last_ref_pos << "; largest_last_inst_pos is " << largest_last_inst_pos << endl;

        pos_width =
          max( max_ref_i_width, largest_inst_i_width );
        block_width =
          max( 0, ( int )wrap_column - (int)max_label_width - 1  - ( int )max_seq_i_width - 1 - ( int )pos_width - 1  -  1 - ( int )pos_width );
        if( block_width == 0 ) {
          // TODO: Prevent this from happening!
          os << "ERROR: You must use a larger wrap_column argument; " << wrap_column << " is not large enough to accommodate the sequence numbers on either side of the output." << endl;
          return os;
        }
        // TODO: REMOVE
        //cout << "pos_width is " << pos_width << "; block_width is " << block_width << "; wrap_column is " << wrap_column << endl;

        if( instance_strings.size() == 0 ) {
          reference_string =
            getPileupStrings(
              instance_strings
            );
        }
        alignment_length =
          reference_string.length();

        // Blocks
        ref_i = 1;
        for( alignment_i = 0; alignment_i < alignment_length; alignment_i += block_width ) {
          // TODO: REMOVE
          //cout << "alignment_i is " << alignment_i << endl;

          // On top is the Profile line.
          substring =
            reference_string.substr(
              alignment_i,
              block_width
            );
          os << ref_label;
          seq_label_width = ref_label_width;
          for( padding_i = 0; padding_i < ( max_label_width - seq_label_width ); padding_i++ ) {
            os << " ";
          }
          os << " ";
          os << ref_i;
          ref_i_width = log10( ref_i ) + 1; //( uint32_t )( log( ref_i ) / log10 ) + 1;
          //ref_i_width = 0;
          for( padding_i = 0; padding_i < ( pos_width - ref_i_width ); padding_i++ ) {
            os << " ";
          }
          os << " ";
          os << substring;
          // Count ref advances
          ref_i += substring.length();
          // Don't count the insertions..
          // (Note ref_i--)
          for( substring_i = 0; substring_i < substring.length(); ref_i-- ) {
            find_pos = substring.find( '.', substring_i );
            if( find_pos == string::npos ) {
              break;
            }
            substring_i = ( find_pos + 1 );
            //cout << "substring_i is " << substring_i;
            //cout << "ref_i is " << ref_i;
          }
          os << " " << ( ref_i - 1 );
          os << endl;

          // Now come the sequence lines
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {

            substring =
              instance_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                block_width
              );
            os << 
              local_seq_labels[ m_sequenceOrder[ seq_i ] ];
            seq_label_width =
              local_seq_labels[ m_sequenceOrder[ seq_i ] ].length();

            for( padding_i = 0; padding_i < ( max_label_width - seq_label_width ); padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << inst_i[ m_sequenceOrder[ seq_i ] ];
            inst_i_width =
              log10( inst_i[ m_sequenceOrder[ seq_i ] ] ) + 1;
            //( uint32_t )( log( inst_i[ m_sequenceOrder[ seq_i ] ] ) / log10 ) + 1;
            for( padding_i = 0; padding_i < ( pos_width - inst_i_width ); padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            // Count inst advances
            inst_i[ m_sequenceOrder[ seq_i ] ] += substring.length();
            // Don't count the insertions..
            // (Note inst_i[ m_sequenceOrder[ seq_i ] ]--)
            for( substring_i = 0; substring_i < substring.length(); inst_i[ m_sequenceOrder[ seq_i ] ]-- ) {
              find_pos = substring.find_first_of( "-. ", substring_i );
              if( find_pos == string::npos ) {
                break;
              }
              substring_i = ( find_pos + 1 );
              //cout << "substring_i is " << substring_i;
              //cout << "inst_i[ m_sequenceOrder[ seq_i ] ] is " << inst_i[ m_sequenceOrder[ seq_i ] ];
            }
            os << " " << ( inst_i[ m_sequenceOrder[ seq_i ] ] - 1 );
            os << endl;

          } // End foreach sequence

            // An extra line between blocks:
          os << endl;
        } // End foreach block
  
        os << endl;

        return os;
      } // toPileupStream( ostream& )

      /**
       * Represent this multiple alignment in a BLAST-like way.  Get the
       * m_sequence_count (unbroken) strings for the model ("reference") lines,
       * the alignment lines, and the sequence ("instance") lines.  Modifies
       * the given vectors to hold these lines.
       */
      void
      getPairwiseStrings (
        vector<string> & reference_strings,
        vector<string> & alignment_strings,
        vector<string> & instance_strings
      ) const
      {
        // TODO: REMOVE
        //cout << "In getPairwiseStrings(..)" << endl;

        reference_strings.clear();
        reference_strings.resize( m_sequence_count );
        alignment_strings.clear();
        alignment_strings.resize( m_sequence_count );
        instance_strings.clear();
        instance_strings.resize( m_sequence_count )
;
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_seq_adv_i = m_profile->length() + 3;
        uint32_t seq_adv_i;
        uint32_t last_seq_pos;
        uint32_t seq_pos_i;
        uint32_t tmp_seq_pos_i;
        uint32_t pairwise_pos_i;

        stringstream os;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          last_seq_pos = ( *m_sequences )[ seq_i ].length() - 1;
          seq_adv_i = 0;

          while( seq_adv_i <= last_seq_adv_i ) {
          
            // On top is the Profile line.
            os.str( "" ); // Empty the stringstream
            seq_adv_i = 0;
            seq_pos_i = 0;
            pairwise_pos_i = 0;
            if( seq_pos_i == 0 ) {
              // Take care of any pre-align insertions.
              pairwise_pos_i = m_sequenceAdvances[ seq_i ][ 0 ];
              while( seq_pos_i < pairwise_pos_i ) {
                os << "-"; // An insertion
                seq_pos_i++;
              }
              seq_adv_i = 2; // Skip the del-in count, since it's redundant (the deleted positions will have 0s too)
            } // End if it's the start of the very first block
            for( ;
                 ( seq_adv_i <= last_seq_adv_i );
                 seq_adv_i++ ) {
              if( seq_adv_i == ( last_seq_adv_i - 1 ) ) {
                // Skip the del-out count, since it's redundant (the deleted
                // positions will have 0s too).
                continue;
              }
              if(
                ( seq_adv_i != last_seq_adv_i ) &&
                ( m_sequenceAdvances[ seq_i ][ seq_adv_i ] == 0 )
              ) {
                // A deletion
                //os << ( *m_profile )[ seq_adv_i - 2 ][ Emission::Match ].maximumValueType().lowercase();
                os << ( *m_profile )[ seq_adv_i - 2 ][ Emission::Match ].maximumValueType();
                pairwise_pos_i++;
              } else {
                if( seq_adv_i != last_seq_adv_i ) {
                  // A match, with maybe some insertions following.
                    os << ( *m_profile )[ seq_adv_i - 2 ][ Emission::Match ].maximumValueType();
                  tmp_seq_pos_i = seq_pos_i + 1;
                } else {
                  // A post-align insertion
                  tmp_seq_pos_i = seq_pos_i;
                }
                seq_pos_i += m_sequenceAdvances[ seq_i ][ seq_adv_i ];
                while( tmp_seq_pos_i++ < seq_pos_i ) {
                  os << "-"; // An insertion
                }
                pairwise_pos_i += m_sequenceAdvances[ seq_i ][ seq_adv_i ];
              } // End if it's a deletion .. else ..
            } // End foreach seq_adv_i
            reference_strings[ seq_i ] = os.str();
          
            // In the middle is the Alignment line.
            os.str( "" ); // Empty the stringstream
            seq_adv_i = 0;
            seq_pos_i = 0;
            // Take care of any pre-align insertions.
            pairwise_pos_i = m_sequenceAdvances[ seq_i ][ 0 ];
            while( seq_pos_i < pairwise_pos_i ) {
              os << "^"; // An insertion
              seq_pos_i++;
            }
            for( seq_adv_i = 2; seq_adv_i <= last_seq_adv_i; seq_adv_i++ ) {
              if( seq_adv_i == ( last_seq_adv_i - 1 ) ) {
                continue;
              }
              if(
                ( seq_adv_i != last_seq_adv_i ) &&
                ( m_sequenceAdvances[ seq_i ][ seq_adv_i ] == 0 )
              ) {
                // A deletion
                os << "v";
                pairwise_pos_i++;
              } else {
                if( seq_adv_i != last_seq_adv_i ) {
                  // A match, with maybe some insertions following.
                  if(  
                         ( *m_sequences )[ seq_i ][ seq_pos_i ] ==
                         ( *m_profile )[ seq_adv_i - 2 ][ Emission::Match ].maximumValueType()
                       
                  ) {
                    os << "|";
                  } else {
                    os << " ";
                  }
                  tmp_seq_pos_i = seq_pos_i + 1;
                } else {
                  tmp_seq_pos_i = seq_pos_i;
                }
                seq_pos_i += m_sequenceAdvances[ seq_i ][ seq_adv_i ];
                while( tmp_seq_pos_i++ < seq_pos_i ) {
                  os << "^"; // An insertion
                }
                pairwise_pos_i += m_sequenceAdvances[ seq_i ][ seq_adv_i ];
              } // End if it's a deletion .. else ..
            } // End foreach seq_adv_i
            alignment_strings[ seq_i ] = os.str();
          
            // On bottom is the Sequence line.
            os.str( "" ); // Empty the stringstream
            seq_adv_i = 0;
            seq_pos_i = 0;
            // Take care of any pre-align insertions.
            pairwise_pos_i = m_sequenceAdvances[ seq_i ][ 0 ];
            while( seq_pos_i < pairwise_pos_i ) {
              //os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ].lowercase();
              os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ];
            }
            for( seq_adv_i = 2; seq_adv_i <= last_seq_adv_i; seq_adv_i++ ) {
              if( seq_adv_i == ( last_seq_adv_i - 1 ) ) {
                continue;
              }
              if(
                ( seq_adv_i != last_seq_adv_i ) &&
                ( m_sequenceAdvances[ seq_i ][ seq_adv_i ] == 0 )
              ) {
                // Mark del-ins and del-outs differently
                if( seq_adv_i < ( 2 + m_sequenceAdvances[ seq_i ][ 1 ] ) ) {
                  // A DeletionIn
                  os << ".";
                } else if( seq_adv_i > ( last_seq_adv_i - 2 - m_sequenceAdvances[ seq_i ][ last_seq_adv_i - 1 ] ) ) {
                  // A DeletionOut
                  os << ".";
                } else {
                  // A deletion
                  os << "-";
                }
                pairwise_pos_i++;
              } else {
                if( seq_adv_i != last_seq_adv_i ) {
                  // A match, with maybe some insertions following.
                    os << ( *m_sequences )[ seq_i ][ seq_pos_i ];
                  tmp_seq_pos_i = seq_pos_i + 1;
                } else {
                  // Post-align insertion
                  tmp_seq_pos_i = seq_pos_i;
                }
                seq_pos_i += m_sequenceAdvances[ seq_i ][ seq_adv_i ];
                while( tmp_seq_pos_i < seq_pos_i ) {
                  //os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ].lowercase(); // An insertion
                  os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ]; // An insertion
                  tmp_seq_pos_i++;
                }
                pairwise_pos_i += m_sequenceAdvances[ seq_i ][ seq_adv_i ];
              } // End if it's a deletion .. else ..
            } // End foreach seq_adv_i
            instance_strings[ seq_i ] = os.str();
          
          } // End while( seq_adv_i < last_seq_adv_i )

          // TODO: REMOVE
          //cout << reference_strings[ seq_i ] << endl;
          //cout << alignment_strings[ seq_i ] << endl;
          //cout << instance_strings[ seq_i ] << endl;
        } // End foreach sequence

      } // getPairwiseStrings( vector<string> &, vector<string> &, vector<string> &, uint8_t ) const

      // mark
      /**
       * Represent this multiple alignment in a SELEX-like way.  Get the
       * (unbroken) string for the model ("reference") line and the
       * m_sequence_count (unbroken) strings for the sequence ("instance")
       * lines.  Modifies the given vector to hold the sequence lines, and
       * returns the reference line.
       */
      string
      getPileupStrings (
        vector<string> & instance_strings
      ) const
      {
        // TODO: REMOVE
        //cout << "In getPileupStrings(..)" << endl;

        instance_strings.clear();
        instance_strings.resize( m_sequence_count );
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_seq_adv_i = m_profile->length() + 3;
        uint32_t seq_adv_i;
        uint32_t last_seq_pos;
        uint32_t seq_pos_i;
        uint32_t tmp_seq_pos_i;
        uint32_t pileup_pos_i;
        uint32_t extra_insertions_counter_i;
        uint32_t extra_insertions_counter_max;

        // Calculate a sequenceAdvances vector which has the max of all
        // seq's sequenceAdvances vectors at each position, and while
        // we're at it, calculate the length of the longest sequence.
        vector<uint32_t> max_sequenceAdvances( last_seq_adv_i + 1, 0 );
        uint32_t largest_last_seq_pos = 0;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          last_seq_pos = ( *m_sequences )[ seq_i ].length();
          if( last_seq_pos > largest_last_seq_pos ) {
            largest_last_seq_pos = last_seq_pos;
          }
          for( seq_adv_i = 0; seq_adv_i <= last_seq_adv_i; seq_adv_i++ ) {
            if( ( seq_adv_i == 1 ) || ( seq_adv_i == ( last_seq_adv_i - 1 ) ) ) {
              continue;
            }
            if(
              m_sequenceAdvances[ seq_i ][ seq_adv_i ] >
              max_sequenceAdvances[ seq_adv_i ]
            ) {
              max_sequenceAdvances[ seq_adv_i ] =
                m_sequenceAdvances[ seq_i ][ seq_adv_i ];
            }
          } // Foreach position of the sequenceAdvances vectors...
        } // End foreach seq_i..

        stringstream os;

        // the Sequence lines.
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          os.str( "" ); // Empty the stringstream
          seq_adv_i = 0;
          seq_pos_i = 0;
          // Take care of any pre-align insertions.
          pileup_pos_i = m_sequenceAdvances[ seq_i ][ 0 ];
          while( seq_pos_i < pileup_pos_i ) {
            //os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ].lowercase();
            os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ];
          }
          extra_insertions_counter_max =
            (
              max_sequenceAdvances[ 0 ] -
              m_sequenceAdvances[ seq_i ][ 0 ]
            );
          for(
            extra_insertions_counter_i = 0;
            extra_insertions_counter_i < extra_insertions_counter_max;
            extra_insertions_counter_i++
          ) {
            os << "."; //" "; //".";
          }
          pileup_pos_i = max_sequenceAdvances[ 0 ];

          for( seq_adv_i = 2; seq_adv_i <= last_seq_adv_i; seq_adv_i++ ) {
            if( seq_adv_i == ( last_seq_adv_i - 1 ) ) {
              continue;
            }
            if(
              ( seq_adv_i != last_seq_adv_i ) &&
              ( m_sequenceAdvances[ seq_i ][ seq_adv_i ] == 0 )
            ) {
              // A deletion
              os << "-";
              pileup_pos_i += 1;
            } else {
              if( seq_adv_i != last_seq_adv_i ) {
                // A match, with maybe some insertions following.
                os << ( *m_sequences )[ seq_i ][ seq_pos_i ];
                tmp_seq_pos_i = seq_pos_i + 1;
                pileup_pos_i += 1;
              } else {
                // Post-align insertion
                tmp_seq_pos_i = seq_pos_i;
              }
              seq_pos_i += m_sequenceAdvances[ seq_i ][ seq_adv_i ];
              while( tmp_seq_pos_i < seq_pos_i ) {
                //os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ].lowercase(); // An insertion
                os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ]; // An insertion
                tmp_seq_pos_i++;
                pileup_pos_i++;
              } // End while adding insertions..
            } // End if it's a deletion .. else ..

            // Now add any extra insertions that aren't in this seuqence but
            // are in the alignment.
            // OLD:
            //extra_insertions_counter_max =
            //  (
            //    max_sequenceAdvances[ seq_adv_i ] -
            //    max(
            //      ( uint32_t )( ( seq_adv_i == last_seq_adv_i ) ? 0 : 1 ),
            //      m_sequenceAdvances[ seq_i ][ seq_adv_i ]
            //    )
            //  );
            // NEW:
            extra_insertions_counter_max =
              max_sequenceAdvances[ seq_adv_i ];
            if(
              ( seq_adv_i != last_seq_adv_i ) &&
              ( m_sequenceAdvances[ seq_i ][ seq_adv_i ] == 0 )
            ) {
              if( extra_insertions_counter_max > 0 ) {
                extra_insertions_counter_max -= 1;
              } // else leave it at 0
            } else {
              extra_insertions_counter_max -=
                max(
                  ( uint32_t )( ( seq_adv_i == last_seq_adv_i ) ? 0 : 1 ),
                  m_sequenceAdvances[ seq_i ][ seq_adv_i ]
                );
            } // End subtracting from the extra insertions counter...

            for(
              extra_insertions_counter_i = 0;
              extra_insertions_counter_i < extra_insertions_counter_max;
              extra_insertions_counter_i++
            ) {
              os << ".";//" "; //".";
            }
              
            pileup_pos_i += max_sequenceAdvances[ seq_adv_i ];

          } // End foreach seq_adv_i
          instance_strings[ seq_i ] = os.str();
          
        } // End foreach sequence

        // the reference line.
        os.str( "" ); // Empty the stringstream
        seq_adv_i = 0;
        seq_pos_i = 0;
        pileup_pos_i = 0;

        // Take care of any pre-align insertions.
        pileup_pos_i = max_sequenceAdvances[ 0 ];
        while( seq_pos_i < pileup_pos_i ) {
          os << "."; // An insertion
          seq_pos_i++;
        }
        for( seq_adv_i = 2;
             ( seq_adv_i <= last_seq_adv_i );
             seq_adv_i++
        ) {
          if( seq_adv_i == ( last_seq_adv_i - 1 ) ) {
            continue;
          }
          if( seq_adv_i != last_seq_adv_i ) {
            // A match, with maybe some insertions following.
              os << ( *m_profile )[ seq_adv_i - 2 ][ Emission::Match ].maximumValueType();
            pileup_pos_i++;
            tmp_seq_pos_i = seq_pos_i + 1;
          } else {
            // A post-align insertion
            tmp_seq_pos_i = seq_pos_i;
          }
          seq_pos_i += max_sequenceAdvances[ seq_adv_i ];
          while( tmp_seq_pos_i++ < seq_pos_i ) {
            os << "."; // An insertion
          }
          pileup_pos_i += max_sequenceAdvances[ seq_adv_i ];
        } // End foreach seq_adv_i
        // return the reference string.
        return os.str();
      } // getPileupStrings( vector<string> &, uint8_t ) const
    // endmark

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        MultipleAlignment<ProfileType, SequenceResidueType> const &ma )
      {
        return ma.toStream( os );
      } // friend operator<< ( basic_ostream, MultipleAlignment const & )

      // An interface for handling named fasta-style alignment strings (gapped
      // sequences).  Used by toGappedSequenceContainer(..).
      class GappedSequenceContainer
      {
      public:
        // This will always be called first.
        virtual void setNumSequences ( uint32_t num_seqs ) = 0;
        // This will always be called second.
        // numColumns is the length of each gapped sequence string
        virtual void setNumColumns ( uint32_t num_cols ) = 0;
        virtual void addGappedSequence ( string const & seq_name, string const & gapped_sequence ) = 0;
      }; // End inner class DynamicProgramming::MultipleAlignment::GappedSequenceContainer

#ifdef __HAVE_MUSCLE
      class MuscleGappedSequenceContainer : public GappedSequenceContainer
      {
      public:
        MSA * m_msa;

        MuscleGappedSequenceContainer ( MSA * msa ) :
          m_msa( msa )
        {
          // Do nothing else.
        } // <init>( MSA & )

        // This will always be called first.
        void setNumSequences ( uint32_t num_seqs )
        {
          // Do nothing.
        } // setNumSequences( uint32_t )

        // This will always be called second.
        // numColumns is the length of each gapped sequence string
        void setNumColumns ( uint32_t num_cols )
        {
          // Do nothing.
        } // setNumColumns ( uint32_t )

        void addGappedSequence (
          string const & seq_label_string,
          string const & gapped_sequence_string
        )
        {
          char * gapped_sequence = new char[ gapped_sequence_string.length() + 1 ];
          memcpy( gapped_sequence, gapped_sequence_string.c_str(), gapped_sequence_string.length() + 1 );
          char * seq_label = new char[ seq_label_string.length() + 1 ];
          memcpy( seq_label, seq_label_string.c_str(), seq_label_string.length() + 1 );
          // TODO: REMOVE
          //cout << "APPENDING \"" << seq_label << "\"" << ", length " << gapped_sequence_string.length() << endl;
          m_msa->AppendSeq(
            gapped_sequence,
            gapped_sequence_string.length(),
            seq_label
          );
        } // addGappedSequence( string const &, string const & )

      }; // End inner class DynamicProgramming::MultipleAlignment::MuscleGappedSequenceContainer

      /**
       * For conversion between MultipleAlignment objects and muscle's MSA
       * objects.  Appends this MultipleAlignment's data to the end of the
       * given MSA.
       */
      void
      appendToMSA (
        MSA * msa,
        vector<string> const* seq_labels_ptr = NULL
      ) const
      {
        MuscleGappedSequenceContainer gsc =
          MuscleGappedSequenceContainer( msa );

        toGappedSequenceContainer( gsc, seq_labels_ptr );
      } // appendToMSA( MSA &, vector<string> * ) const
#endif //__HAVE_MUSCLE
    
      /**
       * For conversion between MultipleAlignment objects and anything
       * implementing the GappedSequenceContainer interface.
       */
      void
      toGappedSequenceContainer (
        GappedSequenceContainer & gsc,
        vector<string> const* seq_labels_ptr = NULL
      ) const
      {
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;

        string reference_string;
        vector<string> instance_strings;

        uint32_t alignment_i;
        uint32_t alignment_length = 0;

        if( instance_strings.size() == 0 ) {
          reference_string =
            getPileupStrings(
              instance_strings
            );
        }
        alignment_length =
          reference_string.length();

        gsc.setNumSequences( last_seq + 1 );
        gsc.setNumColumns( alignment_length );

        string seq_label_string;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          // Now append it to the GSC.
          if( seq_labels_ptr == NULL ) {
            seq_label_string = "Sequence " + seq_i;
          } else {
            seq_label_string = ( *seq_labels_ptr )[ seq_i ];
          }
          gsc.addGappedSequence( seq_label_string, instance_strings[ seq_i ] );
        } // End foreach seq_i

      } // toGappedSequenceContainer( GappedSequenceContainer &, vector<string> * ) const
    
    }; // End inner class DynamicProgramming::MultipleAlignment

    /**
     * A TreeMultipleAlignment stores one path per sequence for a vector of
     * Sequences and a ProfileTree.  Internally we also store
     * sufficient information to display these paths as a pileup multiple
     * alignment.
     */
    template <typename ProfileTreeType, typename SequenceResidueType>
    class TreeMultipleAlignment
    {
      //protected:
    public:
      // The tree
      ProfileTreeType const * m_profileTree;
      vector<Sequence<SequenceResidueType> > const * m_sequences;
      /**
       * The number of sequences to use (the index one greater than that of the
       * last one to be used; must be <= m_sequences.size().
       */
      uint32_t m_sequence_count;

      // A pointer to a vector of the same length as the number of nodes in the
      // tree, indexed by vertex, with a vector of indices at each position --
      // indicating assignments of sequences to nodes of the tree.  Indices are
      // into the m_sequences vector.  NOTE: For now we assume in
      // alignLeaves(..) that all internal nodes have no sequences assigned to
      // them, and that all leaf nodes have at least one sequence.
      vector<vector<uint32_t> > const * m_sequenceIndices;

      // A pointer to a vector of the same length as the number of nodes in the
      // tree, indexed by parent vertex, with a vector of indices at each
      // position -- profile-profile alignment advances for the alignment
      // between the parent's two children.  Zero-length for leaves.
      vector<vector<uint32_t> > const * m_profileProfileAlignments;

      // A vector of length m_sequence_count that contains each number in the
      // range 0 .. ( m_sequence_count - 1 ) exactly once, specifying a
      // reordering of the sequences.
      vector<uint32_t > m_sequenceOrder;

      /**
       * m_insertionsAfterPosition is a vector of int-vectors, one int-vector
       * per sequence.  Each int-vector has length one greater than that of the
       * root profile, and the integer at index i indicates the number of
       * sequence positions inserted by root position (i-1).  Index 0 indicates
       * the number of pre-align insertions.
       */
      vector<vector<uint32_t> > m_insertionsAfterPosition;

      /**
       * m_matchIndicators is a vector of bool-vectors, one bool-vector per
       * sequence (actually for some reason I can't get gcc to compile if I use
       * bools, so I use uint8_t instead -- meaning you must explicitly test
       * (== true)).  Each vector has the same length as the root profile, and
       * the value at index i indicates whether there was a Match (as opposed
       * to a Deletion) at root position i.
       */
      vector<vector<uint8_t> > m_matchIndicators;

      TreeMultipleAlignment () :
        m_profileTree( NULL ),
        m_sequences( NULL ),
        m_sequence_count( 0 ),
        m_sequenceIndices( NULL ),
        m_profileProfileAlignments( NULL ),
        m_sequenceOrder(),
        m_insertionsAfterPosition(),
        m_matchIndicators()
      {
        // Do nothing else.
      } // <init>()

      /**
       * Only use the first sequence_count sequences of the sequences
       * vector.
       */
      TreeMultipleAlignment (
        ProfileTreeType const * const profile_tree,
        vector<Sequence<SequenceResidueType> > const * const sequences,
        uint32_t const& arg_sequence_count_arg,
        vector<vector<uint32_t> > const * sequence_indices,
        vector<vector<uint32_t> > const * profile_profile_alignments
      ) :
        m_profileTree( profile_tree ),
        m_sequences( sequences ),
        m_sequence_count( ( arg_sequence_count_arg == 0 ) ? sequences->size() : min( arg_sequence_count_arg, sequences->size() ) ),
        m_sequenceOrder( m_sequence_count ),
        m_sequenceIndices( sequence_indices ),
        m_profileProfileAlignments( profile_profile_alignments ),
        m_insertionsAfterPosition( m_sequence_count ),
        m_matchIndicators( m_sequence_count )
      {
        // Also initialize the vectors.
        uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        uint32_t root_length = profile_tree->getProfileTreeRoot()->length();
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          m_insertionsAfterPosition[ seq_i ].resize( root_length + 1, 0 );
          m_matchIndicators[ seq_i ].resize( root_length, false );
          m_sequenceOrder[ seq_i ] = seq_i;
        }
      } // <init>( ProfileTreeType const * const, vector<Sequence<SequenceResidueType> > const * const, uint32_t const &, vector<vector<uint32_t> > const * )

      /**
       * Only use the first sequence_count sequences of the sequences
       * vector.
       */
      void
      reinitialize (
        ProfileTreeType const * const profile_tree,
        vector<Sequence<SequenceResidueType> > const * const sequences,
        uint32_t const& sequence_count,
        vector<vector<uint32_t> > const * sequence_indices,
        vector<vector<uint32_t> > const * profile_profile_alignments
      )
      {
        m_profileTree = profile_tree;
        m_sequences = sequences;
        m_sequence_count = ( ( sequence_count == 0 ) ? sequences->size() : min( (size_t &)sequence_count, sequences->size() ) );
        m_sequenceIndices = sequence_indices;
        m_profileProfileAlignments = profile_profile_alignments;

        m_sequenceOrder.resize( m_sequence_count );
        m_insertionsAfterPosition.resize( m_sequence_count );
        m_matchIndicators.resize( m_sequence_count );

        if( ( m_sequence_count == 0 ) || ( profile_tree == 0 ) ) {
          return;
        }
        // Also initialize the vectors.
        uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        uint32_t root_length = profile_tree->getProfileTreeRoot()->length();
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          m_insertionsAfterPosition[ seq_i ].resize( root_length + 1, 0 );
          m_matchIndicators[ seq_i ].resize( root_length, false );
          m_sequenceOrder[ seq_i ] = seq_i;
        }

      } // reinitialize( ProfileTreeType const&, vector<Sequence<SequenceResidueType> >, uint32_t const &, vector<vector<uint32_t> > const * )

      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toStream ( std::basic_ostream<CharT,Traits>& os ) const
      {
        uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        uint32_t last_align_pos = m_profileTree->getProfileTreeRoot()->length();
        uint32_t align_pos_i;
        for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
          os << "Sequence " << m_sequenceOrder[ seq_i ] << ":" << endl;
          os << "( ";
          for( align_pos_i = 0; align_pos_i <= last_align_pos ; align_pos_i++ ) {
            if( align_pos_i != 0 ) {
              os << ", ";
              os << ( ( m_matchIndicators[ m_sequenceOrder[ seq_i ] ][ align_pos_i - 1 ] == true ) ? "T " : "F " );
            }
            os << m_insertionsAfterPosition[ m_sequenceOrder[ seq_i ] ][ align_pos_i ];
          }
          os << " )" << endl;
        }
        return os;
      } // toStream( ostream& )

      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toPairwiseStream (
        std::basic_ostream<CharT,Traits>& os,
        vector<string> const* seq_labels_ptr = NULL,
        uint32_t wrap_column = 80
      ) const
      {
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_ref_pos =
          m_profileTree->getProfileTreeRoot()->length();

        uint32_t block_width;

        vector<string> reference_strings;
        vector<string> alignment_strings;
        vector<string> instance_strings;

        uint32_t alignment_i;
        uint32_t alignment_length;

        uint32_t ref_i; // in coords of 1..last_ref_pos
        uint32_t ref_i_width; // space req'd to print ref_i
        uint32_t inst_i; // in coords of 1..last_inst_pos
        uint32_t inst_i_width; // space req'd to print inst_i
        uint32_t last_inst_pos; // in coords of 1..last_inst_pos
        uint32_t pos_width; // max space req'd to print inst_i or ref_i
        uint32_t padding_i;
        uint32_t substring_i;
        size_t find_pos;

        string substring;

        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          if( seq_labels_ptr == NULL ) {
            os << "> Sequence " << m_sequenceOrder[ seq_i ] << endl;
          } else {
            os << "> " << ( *seq_labels_ptr )[ m_sequenceOrder[ seq_i ] ] << endl;
          }
          os << endl;

          last_inst_pos = ( *m_sequences )[ m_sequenceOrder[ seq_i ] ].length();
          // TODO: REMOVE
          //cout << "last_ref_pos is " << last_ref_pos << "; last_inst_pos is " << last_inst_pos << endl;

          pos_width =
            max(
              ( log10( last_ref_pos ) + 1 ),
              ( log10( last_inst_pos ) + 1 )
            );
          block_width = // "Sequence" is 8 long
            max( 0, ( int )wrap_column  -  8 - 1  -  ( int )pos_width - 1  -  1 - ( int )pos_width );
          if( block_width == 0 ) {
            // TODO: Prevent this from happening!
            os << "ERROR: You must use a larger wrap_column argument; " << wrap_column << " is not large enough to accommodate the sequence numbers on either side of the output." << endl;
            return os;
          }
          // TODO: REMOVE
          //cout << "pos_width is " << pos_width << "; block_width is " << block_width << "; wrap_column is " << wrap_column << endl;

          if( reference_strings.size() == 0 ) {
            getPairwiseStrings(
              reference_strings,
              alignment_strings,
              instance_strings
            );
          }
          alignment_length =
            reference_strings[ 0 ].length(); // all vecs same len

          // TODO: REMOVE
          //cout << "alignment_length is " << alignment_length << endl;

          // Blocks
          ref_i = 1;
          inst_i = 1;
          for( alignment_i = 0; alignment_i < alignment_length; alignment_i += block_width ) {
            // TODO: REMOVE
            //cout << "alignment_i is " << alignment_i << endl;

            // On top is the Profile line.
            substring =
              reference_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                block_width
              );
            os << "Model    " << ref_i;
            ref_i_width = log10( ref_i ) + 1;
            for( padding_i = 0; padding_i < ( pos_width - ref_i_width ); padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            // Count ref advances
            ref_i += substring.length();
            // (Note ref_i--)
            for( substring_i = 0; substring_i < substring.length(); ref_i-- ) {
              find_pos = substring.find( '-', substring_i );
              if( find_pos == string::npos ) {
                break;
              }
              substring_i = ( find_pos + 1 );
              //cout << "substring_i is " << substring_i;
              //cout << "ref_i is " << ref_i;
            }
            os << " " << ( ref_i - 1 );
            os << endl;
  
            // In the middle is the Alignment line.
            substring =
              alignment_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                block_width
              );
            os << "         "; // 9 long
            for( padding_i = 0; padding_i < pos_width; padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            os << endl;
    
            // On bottom is the Sequence line.
            substring =
              instance_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                block_width
              );
            os << "Sequence " << inst_i;
            inst_i_width = log10( inst_i ) + 1;
            for( padding_i = 0; padding_i < ( pos_width - inst_i_width ); padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            // Count inst advances
            inst_i += substring.length();
            // (Note inst_i--)
            for( substring_i = 0; substring_i < substring.length(); inst_i-- ) {
              find_pos = substring.find( '-', substring_i );
              if( find_pos == string::npos ) {
                break;
              }
              substring_i = ( find_pos + 1 );
              //cout << "substring_i is " << substring_i;
              //cout << "inst_i is " << inst_i;
            }
            os << " " << ( inst_i - 1 );
            os << endl;

            // An extra line between blocks:
            os << endl;
          } // End foreach block
  
          os << endl;

          // An extra line between sequences:
          os << endl;

        } // End foreach sequence

        return os;
      } // toPairwiseStream( ostream& )

      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toAlignedFastaStream (
        std::basic_ostream<CharT,Traits>& os,
        vector<string> const* seq_labels_ptr = NULL,
        uint32_t const wrap_column = 60,
        bool const & include_model_seq = false,
        string const & model_seq_label = "Model"
      ) const
      {
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;

        string reference_string;
        vector<string> instance_strings;

        uint32_t alignment_i;
        uint32_t alignment_length = 0;

        if( instance_strings.size() == 0 ) {
          reference_string =
            getPileupStrings(
              instance_strings
            );
        }
        alignment_length =
          reference_string.length();

        if( include_model_seq ) {
          // The model sequence
          os << ">" << model_seq_label << endl;
          for( alignment_i = 0; alignment_i < alignment_length; alignment_i += wrap_column ) {
            os <<
              reference_string.substr(
                alignment_i,
                min( wrap_column, ( alignment_length - alignment_i ) )
              );
            os << endl;
          } // End foreach line ...
          os << endl; // An extra line between seqs
        } // End if include_model_seq

        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          // The model sequence
          if( seq_labels_ptr == NULL ) {
            os << ">" << "Seq" << m_sequenceOrder[ seq_i ] << endl;
          } else {
            os << ">" << ( *seq_labels_ptr )[ m_sequenceOrder[ seq_i ] ] << endl;
          }
          for( alignment_i = 0; alignment_i < alignment_length; alignment_i += wrap_column ) {
            os <<
              instance_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                min( wrap_column, ( alignment_length - alignment_i ) )
              );
            os << endl;
          } // End foreach line ...
          if( seq_i != last_seq ) {
            os << endl; // An extra line between seqs
          }
        } // End foreach seq_i

        return os;
      } // toAlignedFastaStream( ostream&, ... )

      // mark
      template<class CharT, class Traits>
      std::basic_ostream<CharT,Traits>&
      toPileupStream (
        std::basic_ostream<CharT,Traits>& os,
        vector<string> const* seq_labels_ptr = NULL,
        uint32_t wrap_column = 80
      ) const
      {
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_ref_pos =
          m_profileTree->getProfileTreeRoot()->length();

        uint32_t block_width;

        string reference_string;
        vector<string> instance_strings;

        uint32_t alignment_i;
        uint32_t alignment_length = 0;

        uint32_t ref_i; // in coords of 1..last_ref_pos
        uint32_t ref_i_width; // space req'd to print ref_i
        vector<uint32_t> inst_i( m_sequence_count, 1 ); // in coords of 1..last_inst_pos
        uint32_t inst_i_width; // space req'd to print inst_i
        uint32_t last_inst_pos; // in coords of 1..last_inst_pos
        uint32_t pos_width; // max space req'd to print inst_i or ref_i
        //uint32_t seq_i_width; // space req'd to print seq_i
        uint32_t padding_i;
        uint32_t substring_i;
        size_t find_pos;

        string substring;

        string ref_label = "=RF";
        const uint32_t ref_label_width = 3;
        vector<string> local_seq_labels( m_sequence_count );
        uint32_t seq_label_width;
        uint32_t max_seq_label_width = 0;
        // TODO: Make this a parameter:
        static const uint32_t hard_max_seq_label_width = 15;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          //if( seq_labels_ptr == NULL ) {
          local_seq_labels[ seq_i ] = "Seq ";
          local_seq_labels[ seq_i ].append(
            lexical_cast<string>( seq_i )
          ); // use original ordering
          seq_label_width =
            local_seq_labels[ seq_i ].length();
          if( seq_labels_ptr != NULL ) {
            if( seq_label_width < ( hard_max_seq_label_width - 2 ) ) {
              local_seq_labels[ seq_i ].append( ": " );
              local_seq_labels[ seq_i ].append(
                ( *seq_labels_ptr )[ seq_i ],
                0,
                ( hard_max_seq_label_width - seq_label_width )
              );
            }
            seq_label_width =
              local_seq_labels[ seq_i ].length();
          } // End if seq_labels_ptr != NULL
          if( seq_label_width > max_seq_label_width ) {
            max_seq_label_width = seq_label_width;
          }
        } // End foreach seq_i
        uint32_t max_label_width =
          max( max_seq_label_width, ref_label_width ); 

        uint32_t largest_last_inst_pos = 0;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          last_inst_pos = ( *m_sequences )[ seq_i ].length();
          if( last_inst_pos > largest_last_inst_pos ) {
            largest_last_inst_pos = last_inst_pos;
          }
        } // End foreach seq_i..
        uint32_t largest_inst_i_width =
          log10( largest_last_inst_pos ) + 1;
        uint32_t max_ref_i_width =
          log10( last_ref_pos ) + 1;
        uint32_t max_seq_i_width =
          ( ( last_seq == 0 ) ? 1 : ( log10( last_seq ) + 1 ) );

        // TODO: REMOVE
        //cout << "last_ref_pos is " << last_ref_pos << "; largest_last_inst_pos is " << largest_last_inst_pos << endl;

        pos_width =
          max( max_ref_i_width, largest_inst_i_width );
        block_width =
          max( 0, ( int )wrap_column - (int)max_label_width - 1  - ( int )max_seq_i_width - 1 - ( int )pos_width - 1  -  1 - ( int )pos_width );
        if( block_width == 0 ) {
          // TODO: Prevent this from happening!
          os << "ERROR: You must use a larger wrap_column argument; " << wrap_column << " is not large enough to accommodate the sequence numbers on either side of the output." << endl;
          return os;
        }
        // TODO: REMOVE
        //cout << "pos_width is " << pos_width << "; block_width is " << block_width << "; wrap_column is " << wrap_column << endl;

        if( instance_strings.size() == 0 ) {
          reference_string =
            getPileupStrings(
              instance_strings
            );
        }
        alignment_length =
          reference_string.length();

        // Blocks
        ref_i = 1;
        for( alignment_i = 0; alignment_i < alignment_length; alignment_i += block_width ) {
          // TODO: REMOVE
          //cout << "alignment_i is " << alignment_i << endl;

          // On top is the Profile line.
          substring =
            reference_string.substr(
              alignment_i,
              block_width
            );
          os << ref_label;
          seq_label_width = ref_label_width;
          for( padding_i = 0; padding_i < ( max_label_width - seq_label_width ); padding_i++ ) {
            os << " ";
          }
          os << " ";
          ref_i_width = 0;
          for( padding_i = 0; padding_i < ( pos_width - ref_i_width ); padding_i++ ) {
            os << " ";
          }
          os << " ";
          os << substring;
          // Count ref advances
          ref_i += substring.length();
          // Don't count the insertions..
          // (Note ref_i--)
          for( substring_i = 0; substring_i < substring.length(); ref_i-- ) {
            find_pos = substring.find( '.', substring_i );
            if( find_pos == string::npos ) {
              break;
            }
            substring_i = ( find_pos + 1 );
            //cout << "substring_i is " << substring_i;
            //cout << "ref_i is " << ref_i;
          }
          os << " " << ( ref_i - 1 );
          os << endl;

          // Now come the sequence lines
          for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {

            substring =
              instance_strings[ m_sequenceOrder[ seq_i ] ].substr(
                alignment_i,
                block_width
              );
            os << 
              local_seq_labels[ m_sequenceOrder[ seq_i ] ];
            seq_label_width =
              local_seq_labels[ m_sequenceOrder[ seq_i ] ].length();

            for( padding_i = 0; padding_i < ( max_label_width - seq_label_width ); padding_i++ ) {
              os << " ";
            }
            os << " " << inst_i[ m_sequenceOrder[ seq_i ] ];
            inst_i_width = log10( inst_i[ m_sequenceOrder[ seq_i ] ] ) + 1;
            for( padding_i = 0; padding_i < ( pos_width - inst_i_width ); padding_i++ ) {
              os << " ";
            }
            os << " ";
            os << substring;
            // Count inst advances
            inst_i[ m_sequenceOrder[ seq_i ] ] += substring.length();
            // Don't count the insertions..
            // (Note inst_i[ m_sequenceOrder[ seq_i ] ]--)
            for( substring_i = 0; substring_i < substring.length(); inst_i[ m_sequenceOrder[ seq_i ] ]-- ) {
              find_pos = substring.find_first_of( "-. ", substring_i );
              if( find_pos == string::npos ) {
                break;
              }
              substring_i = ( find_pos + 1 );
              //cout << "substring_i is " << substring_i;
              //cout << "inst_i[ m_sequenceOrder[ seq_i ] ] is " << inst_i[ m_sequenceOrder[ seq_i ] ];
            }
            os << " " << ( inst_i[ m_sequenceOrder[ seq_i ] ] - 1 );
            os << endl;

          } // End foreach sequence

            // An extra line between blocks:
          os << endl;
        } // End foreach block
  
        os << endl;

        return os;
      } // toPileupStream( ostream& )
      // endmark

      /**
       * Represent this multiple alignment in a BLAST-like way.  Get the
       * m_sequence_count (unbroken) strings for the model ("reference") lines,
       * the alignment lines, and the sequence ("instance") lines.  Modifies
       * the given vectors to hold these lines.
       */
      void
      getPairwiseStrings (
        vector<string> & reference_strings,
        vector<string> & alignment_strings,
        vector<string> & instance_strings
      ) const
      {
        // TODO: REMOVE
        //cout << "In getPairwiseStrings(..)" << endl;

        reference_strings.clear();
        reference_strings.resize( m_sequence_count );
        alignment_strings.clear();
        alignment_strings.resize( m_sequence_count );
        instance_strings.clear();
        instance_strings.resize( m_sequence_count );
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_ins_pos =
          m_profileTree->getProfileTreeRoot()->length();
        uint32_t ins_pos_i;
        uint32_t last_seq_pos;
        uint32_t seq_pos_i;
        uint32_t tmp_seq_pos_i;
        uint32_t pairwise_pos_i;

        stringstream os;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {

          last_seq_pos = ( *m_sequences )[ seq_i ].length() - 1;
          ins_pos_i = 0;

          while( ins_pos_i <= last_ins_pos ) {
          
            // On top is the Profile line.
            os.str( "" ); // Empty the stringstream
            ins_pos_i = 0;
            seq_pos_i = 0;
            pairwise_pos_i = 0;
            if( seq_pos_i == 0 ) {
              // Take care of any pre-align insertions.
              pairwise_pos_i = m_insertionsAfterPosition[ seq_i ][ 0 ];
              while( seq_pos_i < pairwise_pos_i ) {
                os << "-"; // An insertion
                seq_pos_i++;
              }
              ins_pos_i = 1;
            } // End if it's the start of the very first block
            for( ;
                 ( ins_pos_i <= last_ins_pos );
                 ins_pos_i++ ) {
              if(
                ( m_matchIndicators[ seq_i ][ ins_pos_i - 1 ] == false )
              ) {
                // A deletion
                //os << ( *m_profileTree->getProfileTreeRoot() )[ ins_pos_i - 1 ][ Emission::Match ].maximumValueType().lowercase();
                os << ( *m_profileTree->getProfileTreeRoot() )[ ins_pos_i - 1 ][ Emission::Match ].maximumValueType();
              } else {
                // A match
                  os << ( *m_profileTree->getProfileTreeRoot() )[ ins_pos_i - 1 ][ Emission::Match ].maximumValueType();
                seq_pos_i++;
              } // End if it's a deletion .. else ..
              pairwise_pos_i++;
              tmp_seq_pos_i = seq_pos_i;
              seq_pos_i += m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];
              while( tmp_seq_pos_i++ < seq_pos_i ) {
                os << "-"; // An insertion
              }
              pairwise_pos_i += m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];
            } // End foreach ins_pos_i
            reference_strings[ seq_i ] = os.str();
          
            // In the middle is the Alignment line.
            os.str( "" ); // Empty the stringstream
            ins_pos_i = 0;
            seq_pos_i = 0;
            // Take care of any pre-align insertions.
            pairwise_pos_i = m_insertionsAfterPosition[ seq_i ][ 0 ];
            while( seq_pos_i < pairwise_pos_i ) {
              os << "^"; // An insertion
              seq_pos_i++;
            }
            for( ins_pos_i = 1; ins_pos_i <= last_ins_pos; ins_pos_i++ ) {
              if(
                ( m_matchIndicators[ seq_i ][ ins_pos_i - 1 ] == false )
              ) {
                // A deletion
                os << "v";
              } else {
                // A match
                if( 
                  ( *m_sequences )[ seq_i ][ seq_pos_i ] ==
                  ( *m_profileTree->getProfileTreeRoot() )[ ins_pos_i - 1 ][ Emission::Match ].maximumValueType()
                ) {
                  os << "|";
                } else {
                  os << " ";
                }
                seq_pos_i++;
              } // End if it's a deletion .. else ..
              pairwise_pos_i++;
              tmp_seq_pos_i = seq_pos_i;
              seq_pos_i += m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];
              while( tmp_seq_pos_i++ < seq_pos_i ) {
                os << "^"; // An insertion
              }
              pairwise_pos_i += m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];
            } // End foreach ins_pos_i
            alignment_strings[ seq_i ] = os.str();
          
            // On bottom is the Sequence line.
            os.str( "" ); // Empty the stringstream
            ins_pos_i = 0;
            seq_pos_i = 0;
            // Take care of any pre-align insertions.
            pairwise_pos_i = m_insertionsAfterPosition[ seq_i ][ 0 ];
            while( seq_pos_i < pairwise_pos_i ) {
              //os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ].lowercase();
              os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ];
            }
            for( ins_pos_i = 1; ins_pos_i <= last_ins_pos; ins_pos_i++ ) {
              if(
                ( m_matchIndicators[ seq_i ][ ins_pos_i - 1 ] == false )
              ) {
                // A deletion
                os << "-";
              } else {
                // A match
                  os << ( *m_sequences )[ seq_i ][ seq_pos_i ];
                seq_pos_i++;
              } // End if it's a deletion .. else ..
              pairwise_pos_i++;
              tmp_seq_pos_i = seq_pos_i;

              seq_pos_i += m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];
              while( tmp_seq_pos_i < seq_pos_i ) {
                //os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ].lowercase(); // An insertion
                os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ]; // An insertion
                tmp_seq_pos_i++;
              }
              pairwise_pos_i += m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];

            } // End foreach ins_pos_i
            instance_strings[ seq_i ] = os.str();
          
          } // End while( ins_pos_i < last_ins_pos )

          // TODO: REMOVE
          //cout << reference_strings[ seq_i ] << endl;
          //cout << alignment_strings[ seq_i ] << endl;
          //cout << instance_strings[ seq_i ] << endl;
        } // End foreach sequence

        return;
      } // getPairwiseStrings( vector<string> &, vector<string> &, vector<string> &, uint8_t ) const

      // mark
      /**
       * Represent this multiple alignment in a SELEX-like way.  Get the
       * (unbroken) string for the model ("reference") line and the
       * m_sequence_count (unbroken) strings for the sequence ("instance")
       * lines.  Modifies the given vector to hold the sequence lines, and
       * returns the reference line.
       */
      string
      getPileupStrings (
        vector<string> & instance_strings
      ) const
      {
        // TODO: REMOVE
        //cout << "In getPileupStrings(..)" << endl;

        instance_strings.clear();
        instance_strings.resize( m_sequence_count );
        const uint32_t last_seq = m_sequence_count - 1;
        uint32_t seq_i;
        const uint32_t last_ins_pos =
          m_profileTree->getProfileTreeRoot()->length();
        uint32_t ins_pos_i;
        uint32_t last_seq_pos;
        uint32_t seq_pos_i;
        uint32_t tmp_seq_pos_i;
        uint32_t pileup_pos_i;
        uint32_t extra_insertions_counter_i;
        uint32_t extra_insertions_counter_max;

        // Calculate an insertionsAfterPosition vector which has the max of all
        // seq's insertionsAfterPosition vectors at each position, and while
        // we're at it, calculate the length of the longest sequence.
        vector<uint32_t> max_insertionsAfterPosition( last_ins_pos + 1, 0 );
        uint32_t largest_last_seq_pos = 0;
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          last_seq_pos = ( *m_sequences )[ seq_i ].length();
          if( last_seq_pos > largest_last_seq_pos ) {
            largest_last_seq_pos = last_seq_pos;
          }
          for( ins_pos_i = 0; ins_pos_i <= last_ins_pos; ins_pos_i++ ) {
            if(
              m_insertionsAfterPosition[ seq_i ][ ins_pos_i ] >
              max_insertionsAfterPosition[ ins_pos_i ]
            ) {
              max_insertionsAfterPosition[ ins_pos_i ] =
                m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];
            }
          } // Foreach position of the insertionsAfterPosition vectors...
        } // End foreach seq_i..

        stringstream os;

        // the Sequence lines.
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          os.str( "" ); // Empty the stringstream
          ins_pos_i = 0;
          seq_pos_i = 0;
          // Take care of any pre-align insertions.
          pileup_pos_i = m_insertionsAfterPosition[ seq_i ][ 0 ];
          while( seq_pos_i < pileup_pos_i ) {
            //os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ].lowercase();
            os << ( *m_sequences )[ seq_i ][ seq_pos_i++ ];
          }
          extra_insertions_counter_max =
            (
              max_insertionsAfterPosition[ 0 ] -
              m_insertionsAfterPosition[ seq_i ][ 0 ]
            );
          for(
            extra_insertions_counter_i = 0;
            extra_insertions_counter_i < extra_insertions_counter_max;
            extra_insertions_counter_i++
          ) {
            os << " "; //".";
          }
          pileup_pos_i = max_insertionsAfterPosition[ 0 ];

          for( ins_pos_i = 1; ins_pos_i <= last_ins_pos; ins_pos_i++ ) {
            if(
              ( m_matchIndicators[ seq_i ][ ins_pos_i - 1 ] == false )
            ) {
              // A deletion
              os << "-";
            } else {
              // A match
              os << ( *m_sequences )[ seq_i ][ seq_pos_i ];
              seq_pos_i++;
            } // End if it's a deletion .. else ..
            pileup_pos_i++;
            tmp_seq_pos_i = seq_pos_i;

            seq_pos_i += m_insertionsAfterPosition[ seq_i ][ ins_pos_i ];
            while( tmp_seq_pos_i < seq_pos_i ) {
              //os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ].lowercase(); // An insertion
              os << ( *m_sequences )[ seq_i ][ tmp_seq_pos_i ]; // An insertion
              tmp_seq_pos_i++;
            }
            // Now add any extra insertions that aren't in this seuqence but
            // are in the alignment.
            extra_insertions_counter_max =
              (
                max_insertionsAfterPosition[ ins_pos_i ] -
                m_insertionsAfterPosition[ seq_i ][ ins_pos_i ]
              );
            for(
              extra_insertions_counter_i = 0;
              extra_insertions_counter_i < extra_insertions_counter_max;
              extra_insertions_counter_i++
            ) {
              os << " "; //".";
            }

            pileup_pos_i += max_insertionsAfterPosition[ ins_pos_i ];

          } // End foreach ins_pos_i
          instance_strings[ seq_i ] = os.str();
          
        } // End foreach sequence

        // the reference line.
        os.str( "" ); // Empty the stringstream
        ins_pos_i = 0;
        seq_pos_i = 0;
        pileup_pos_i = 0;

        if( seq_pos_i == 0 ) {
          // Take care of any pre-align insertions.
          pileup_pos_i = max_insertionsAfterPosition[ 0 ];
          while( seq_pos_i < pileup_pos_i ) {
            os << "."; // An insertion
            seq_pos_i++;
          }
          ins_pos_i = 1;
        } // End if it's the start of the very first block
        for( ;
             ( ins_pos_i <= last_ins_pos );
             ins_pos_i++
        ) {
          os << ( *m_profileTree->getProfileTreeRoot() )[ ins_pos_i - 1 ][ Emission::Match ].maximumValueType();
          seq_pos_i++;
          pileup_pos_i++;
          tmp_seq_pos_i = seq_pos_i;
          seq_pos_i += max_insertionsAfterPosition[ ins_pos_i ];
          while( tmp_seq_pos_i++ < seq_pos_i ) {
            os << "."; // An insertion
          }
          pileup_pos_i += max_insertionsAfterPosition[ ins_pos_i ];
        } // End foreach ins_pos_i
        // return the reference string.
        return os.str();
      } // getPileupStrings( vector<string> &, uint8_t ) const
    // endmark

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        TreeMultipleAlignment<ProfileTreeType, SequenceResidueType> const &ma )
      {
        return ma.toStream( os );
      } // friend operator<< ( basic_ostream, TreeMultipleAlignment const & )

      // An interface for handling named fasta-style alignment strings (gapped
      // sequences).  Used by toGappedSequenceContainer(..).
      class GappedSequenceContainer
      {
      public:
        // This will always be called first.
        virtual void setNumSequences ( uint32_t num_seqs ) = 0;
        // This will always be called second.
        // numColumns is the length of each gapped sequence string
        virtual void setNumColumns ( uint32_t num_cols ) = 0;
        virtual void addGappedSequence ( string const & seq_name, string const & gapped_sequence ) = 0;
      }; // End inner class DynamicProgramming::TreeMultipleAlignment::GappedSequenceContainer

#ifdef __HAVE_MUSCLE
      class MuscleGappedSequenceContainer : public GappedSequenceContainer
      {
      public:
        MSA * m_msa;

        MuscleGappedSequenceContainer ( MSA * msa ) :
          m_msa( msa )
        {
          // Do nothing else.
        } // <init>( MSA & )

        // This will always be called first.
        void setNumSequences ( uint32_t num_seqs )
        {
          // Do nothing.
        } // setNumSequences( uint32_t )

        // This will always be called second.
        // numColumns is the length of each gapped sequence string
        void setNumColumns ( uint32_t num_cols )
        {
          // Do nothing.
        } // setNumColumns ( uint32_t )

        void addGappedSequence (
          string const & seq_label_string,
          string const & gapped_sequence_string
        )
        {
          char * gapped_sequence = new char[ gapped_sequence_string.length() + 1 ];
          memcpy( gapped_sequence, gapped_sequence_string.c_str(), gapped_sequence_string.length() + 1 );
          char * seq_label = new char[ seq_label_string.length() + 1 ];
          memcpy( seq_label, seq_label_string.c_str(), seq_label_string.length() + 1 );
          // TODO: REMOVE
          //cout << "APPENDING \"" << seq_label << "\"" << endl;
          m_msa->AppendSeq(
            gapped_sequence,
            gapped_sequence_string.length(),
            seq_label
          );
        } // addGappedSequence( string const &, string const & )

      }; // End inner class DynamicProgramming::TreeMultipleAlignment::MuscleGappedSequenceContainer

      /**
       * For conversion between TreeMultipleAlignment objects and muscle's MSA
       * objects.  Appends this TreeMultipleAlignment's data to the end of the
       * given MSA.
       */
      void
      appendToMSA (
        MSA * msa,
        vector<string> const* seq_labels_ptr = NULL
      ) const
      {
        MuscleGappedSequenceContainer gsc =
          MuscleGappedSequenceContainer( msa );

        toGappedSequenceContainer( gsc, seq_labels_ptr );
      } // appendToMSA( MSA &, vector<string> * ) const
#endif //__HAVE_MUSCLE
    
    }; // End inner class DynamicProgramming::TreeMultipleAlignment

    // Represents the current row, column, and subcell of a random walker
    // through the forward-backward matrix.
    class SamplerState
    {
    public:
      uint32_t column;
      Subcell subcell;

      SamplerState () :
        column( 0 ),
        subcell( Match ) 
      {
        // Do nothing else
      } // <init>()

      SamplerState ( SamplerState const & other_state ) :
        column( other_state.column ),
        subcell( other_state.subcell ) 
      {
        // Do nothing else
      } // <init>( SamplerState const & )

      SamplerState &
      operator= ( SamplerState const & other_state )
      {
        column = other_state.column;
        subcell = other_state.subcell;
        return *this;
      } // operator=( SamplerState const & )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        SamplerState const & sampler_state )
      {
        os << "{column=" << sampler_state.column << ", ";
        os << "subcell=" << sampler_state.subcell << "}";

        return os;
      } // friend operator<< ( basic_ostream, SamplerState const & )

    }; // End inner class DynamicProgramming::SamplerState

    /**
     * This class represents a matrix of distances (or proximities) between
     * e.g. sequences.  It might be symmetric (public bool isSymmetric will be
     * true if so); it might be a proximity matrix rather than a distance
     * matrix (isProximityMatrix will be true if so).
     *
     * See also MultinomialDistribution::DistanceMatrix for distances among
     * Residues.
     */
    template <typename DistanceType>
    class DistanceMatrix
    {
    protected:
      uint32_t m_size;
      vector<DistanceType> m_backingVector;
      bool m_isSymmetric;
      bool m_isProximityMatrix; // versus a distance matrix

    public:

      DistanceMatrix (
        uint32_t const size = 0,
        bool const is_symmetric = false,
        bool const is_proximity_matrix = false
      ) :
        m_size( size ),
        m_backingVector(),
        m_isSymmetric( is_symmetric ),
        m_isProximityMatrix( is_proximity_matrix )
      {
        if( size == 0 ) {
          return;
        }
        reinitialize( size, is_symmetric, is_proximity_matrix );
      } // <init>( .. )

      /**
       * Reinitialize with the given parameters and set all distances to 0.
       */
      void
      reinitialize (
        uint32_t const size,
        bool const is_symmetric = false,
        bool const is_proximity_matrix = false
      )
      {
        m_size = size;
        m_isSymmetric = is_symmetric;
        m_isProximityMatrix = is_proximity_matrix;

        uint32_t backing_vector_size =
          ( m_isSymmetric ? ( ( ( size + 1 ) * size ) / 2 ) : ( size * size ) );
        if( m_backingVector.size() != backing_vector_size ) {
          m_backingVector.resize( backing_vector_size );
        }
        // Set all values to 0
        *this = static_cast<DistanceType>( 0 );
      } // reinitialize( uint32_t, [ bool, bool ] )

      uint32_t
      size () const
      {
        return m_size;
      } // size ()

      DistanceType
      maximumValue () const
      {
        DistanceType maximum_value = m_backingVector[ 0 ];
        for( uint32_t i = 1; i < m_backingVector.size(); i++ ) {
          if( m_backingVector[ i ] > maximum_value ) {
            maximum_value = m_backingVector[ i ];
          }
        }
        return maximum_value;
      } // maximumValue () const

      bool
      isSymmetric () const
      {
        return m_isSymmetric;
      } // isSymmetric () const

      bool
      isProximityMatrix () const
      {
        return m_isProximityMatrix;
      } // isProximityMatrix () const

      DistanceType const &
      operator () ( uint32_t const & from, uint32_t const & to ) const
      {
        return m_backingVector[ getBackingVectorIndex( from, to ) ];
      } // operator () ( uint32_t const &, uint32_t const & ) const

      DistanceType &
      operator () ( uint32_t const & from, uint32_t const & to )
      {
        return m_backingVector[ getBackingVectorIndex( from, to ) ];
      } // operator () ( uint32_t const &, uint32_t const & )

      /**
       * Set all distances in this matrix to the given value.
       */
      DistanceMatrix<DistanceType> &
      operator= ( DistanceType const & dist )
      {
        uint32_t backing_vector_size = m_backingVector.size();
        for( uint32_t i = 0; i < backing_vector_size; i++ ) {
          m_backingVector[ i ] = dist;
        }

        return *this;
      } // operator= ( DistanceType const & )

      /**
       * Set the distances in this matrix to be equivalent to those in the
       * given matrix.  Also reinitializes the size of this matrix.  If the
       * other matrix is not symmetric and this matrix is symmetric, this
       * matrix will lose its symmetry.
       */
      template <typename AnyDistanceType>
      DistanceMatrix<DistanceType> &
      operator= ( DistanceMatrix<AnyDistanceType> const & other_matrix )
      {
        uint32_t backing_vector_size = m_backingVector.size();
        if( backing_vector_size != other_matrix.m_backingVector.size() ) {
          m_backingVector.resize( other_matrix.m_backingVector.size() );
          backing_vector_size = m_backingVector.size();
        }
        if( !other_matrix.isSymmetric() ) {
          // Lose symmetry.
          m_isSymmetric = false;
        }
        for( uint32_t i = 0; i < backing_vector_size; i++ ) {
          m_backingVector[ i ] =
            static_cast<DistanceType>( other_matrix.m_backingVector[ i ] );
        }

        return *this;
      } // operator= ( DistanceMatrix<AnyDistanceType> const & )

      template<class CharT, class Traits>
      friend std::basic_ostream<CharT,Traits>&
      operator<< (
        std::basic_ostream<CharT,Traits>& os,
        DistanceMatrix const & matrix )
      {
        uint32_t size = matrix.size();
        for( uint32_t from = 0; from < size; from++ ) {
          os << "[ ";
          for( uint32_t to = 0; to < size; to++ ) {
            if( to > 0 ) {
              os << ", ";
            }
            os << matrix( from, to );
          } // End foreach to
          os << " ]" << endl;
        } // End foreach from

        return os;
      } // friend operator<< ( basic_ostream, DistanceMatrix const & )

      /**
       * @brief Reset the proximities to the
       * MatchEmissionCooccurrenceProbabilities among the given alignment
       * profiles.
       */
      void
      setToMatchEmissionCooccurrenceProbabilities (
        vector<AlignmentProfile> const & alignment_profiles,
        const bool include_nonemission_counts = false
      )
      {
        // Reinitialize to the size of the given vector, YES symmetric, YES a
        // proximity (not a distance)
        this->reinitialize( alignment_profiles.size(), true, true );

        uint32_t size = this->size();
        for( uint32_t to = 0; to < size; to++ ) {
          for( uint32_t from = 0; from <= to; from++ ) {
            this->operator()( from, to ) =
              static_cast<DistanceType>( alignment_profiles[ from ].calculateMatchEmissionCooccurrenceProbability( alignment_profiles[ to ], include_nonemission_counts ) );
          } // End foreach from
        } // End foreach to
      } // setToMatchEmissionCooccurrenceProbabilities( vector<AlignmentProfile> const & )

      /**
       * @brief Reset the proximities to the
       * MatchEmissionCooccurrenceExpectedCounts among the given alignment
       * profiles.
       */
      void
      setToMatchEmissionCooccurrenceExpectedCounts (
        vector<AlignmentProfile> const & alignment_profiles,
        const bool include_nonemission_counts = false
      )
      {
        // Reinitialize to the size of the given vector, YES symmetric, YES a
        // proximity (not a distance)
        this->reinitialize( alignment_profiles.size(), true, true );

        uint32_t size = this->size();
        for( uint32_t to = 0; to < size; to++ ) {
          for( uint32_t from = 0; from <= to; from++ ) {
            this->operator()( from, to ) =
              static_cast<DistanceType>( alignment_profiles[ from ].calculateMatchEmissionCooccurrenceExpectedCount( alignment_profiles[ to ], include_nonemission_counts ) );
          } // End foreach from
        } // End foreach to
      } // setToMatchEmissionCooccurrenceExpectedCounts( vector<AlignmentProfile> const &, bool )

      /**
       * @brief Reset the distances to the crossEntropies among the given
       * alignment profiles.
       */
      void
      setToCrossEntropies (
        vector<AlignmentProfile> const & alignment_profiles,
        AlignmentProfile const * const weights = NULL
      )
      {
        // Reinitialize to the size of the given vector, NO, not symmetric,
        // NOT, not a proximity (yes, a distance)
        this->reinitialize( alignment_profiles.size(), false, false );

        uint32_t size = this->size();
        for( uint32_t to = 0; to < size; to++ ) {
          for( uint32_t from = 0; from < size; from++ ) {
            this->operator()( from, to ) =
              static_cast<DistanceType>( toDouble( alignment_profiles[ from ].crossEntropy( alignment_profiles[ to ], weights ) ) );
          } // End foreach from
        } // End foreach to
      } // setToCrossEntropies( vector<AlignmentProfile> const &, AlignmentProfile const * )

      /**
       * @brief Reset the distances to the crossEntropies among the given
       * alignment profiles, just among the Match Emission distributions at position pos_i.
       */
      void
      setToCrossEntropies (
        vector<AlignmentProfile> const & alignment_profiles,
        uint32_t const & pos_i,
        AlignmentProfile const * const weights = NULL
      )
      {
        // Reinitialize to the size of the given vector, NO, not symmetric,
        // NOT, not a proximity (yes, a distance)
        this->reinitialize( alignment_profiles.size(), false, false );

        uint32_t size = this->size();
        for( uint32_t to = 0; to < size; to++ ) {
          for( uint32_t from = 0; from < size; from++ ) {
            this->operator()( from, to ) =
              static_cast<DistanceType>( toDouble( static_cast<AlignmentProfilePositionParameters<ResidueType,MatrixValueType> >(alignment_profiles[ from ][ pos_i ]).crossEntropy( alignment_profiles[ to ][ pos_i ],
                                                                                                                        ( weights ? &( *weights )[ pos_i ] : ( AlignmentProfilePosition const * )0 )

) ) );
          } // End foreach from
        } // End foreach to
      } // setToCrossEntropies( vector<AlignmentProfile> const &, uint32_t const &, AlignmentProfile const * )

      /**
       * @brief Reset the distances to the Kullback-Leibler divergences.
       * @see setToCrossEntropies
       * @see setToKullbackLeiblerDivergences()
       */
      void
      setToKullbackLeiblerDivergences (
        vector<AlignmentProfile> const & alignment_profiles,
        AlignmentProfile const * const weights = NULL
      )
      {
        setToCrossEntropies( alignment_profiles, weights );

        setToKullbackLeiblerDivergences();
      } // setToKullbackLeiblerDivergences( vector<AlignmentProfile> const &, AlignmentProfile const * )

      /**
       * @brief Reset the distances to the Kullback-Leibler divergences among
       * Match Emission distributions at position pos_i.
       *
       * @see setToCrossEntropies
       * @see setToKullbackLeiblerDivergences()
       */
      void
      setToKullbackLeiblerDivergences (
        vector<AlignmentProfile> const & alignment_profiles,
        uint32_t const & pos_i,
        AlignmentProfile const * const weights = NULL
      )
      {
        setToCrossEntropies( alignment_profiles, pos_i, weights );

        setToKullbackLeiblerDivergences();
      } // setToKullbackLeiblerDivergences( vector<AlignmentProfile> const &, uint32_t const &, AlignmentProfile const * )

      /**
       * @brief Reset the distances to the Kullback-Leibler divergences.
       * Assumes that setToCrossEntropies(..) has already been called (so that
       * the current values are those cross-entropies.
       */
      void
      setToKullbackLeiblerDivergences ()
      {
        uint32_t size = this->size();
        for( uint32_t from = 0; from < size; from++ ) {
          for( uint32_t to = 0; to < size; to++ ) {
            if( to == from ) {
              continue;
            }
            this->operator()( from, to ) -= this->operator()( from, from );
          } // End foreach to
          this->operator()( from, from ) = 0.0;
        } // End foreach from
      } // setToKullbackLeiblerDivergences()

      /**
       * @brief Reset the distances to the symmeterized Kullback-Leibler
       * divergences.
       * @see setToCrossEntropies
       * @see setToKullbackLeiblerDivergences()
       * @see setToSymmeterizedKullbackLeiblerDivergences()
       *
       * This will make the distance matrix symmetric as a side-effect.
       */
      void
      setToSymmeterizedKullbackLeiblerDivergences (
        vector<AlignmentProfile> const & alignment_profiles,
        AlignmentProfile const * const weights = NULL
      )
      {
        setToKullbackLeiblerDivergences( alignment_profiles, weights );
        setToSymmeterizedKullbackLeiblerDivergences();
      } // setToSymmeterizedKullbackLeiblerDivergences( vector<AlignmentProfile> const &, AlignmentProfile const * )

      /**
       * @brief Reset the distances to the symmeterized Kullback-Leibler
       * divergences among the Match Emission distributions at position pos_i.
       *
       * @see setToCrossEntropies
       * @see setToKullbackLeiblerDivergences()
       * @see setToSymmeterizedKullbackLeiblerDivergences()
       *
       * This will make the distance matrix symmetric as a side-effect.
       */
      void
      setToSymmeterizedKullbackLeiblerDivergences (
        vector<AlignmentProfile> const & alignment_profiles,
        uint32_t const & pos_i,
        AlignmentProfile const * const weights = NULL
      )
      {
        setToKullbackLeiblerDivergences( alignment_profiles, pos_i, weights );

        setToSymmeterizedKullbackLeiblerDivergences();
      } // setToSymmeterizedKullbackLeiblerDivergences( vector<AlignmentProfile> const &, uint32_t const &, AlignmentProfile const * )

      /**
       * @brief Reset the distances to the symmeterized Kullback-Leibler divergences.
       * Assumes that setToKullbackLeiblerDivergences(..) has already been
       * called (so that the current values are those (asymmetric) divergences.
       *
       * This will make the distance matrix symmetric as a side-effect.
       */
      void
      setToSymmeterizedKullbackLeiblerDivergences ()
      {
        uint32_t size = this->size();
        for( uint32_t to = 0; to < size; to++ ) {
          for( uint32_t from = 0; from < to; from++ ) {
            this->operator()( to, from ) +=
              this->operator()( from, to );
            this->operator()( from, to ) =
              this->operator()( to, from );
          } // End foreach to
        } // End foreach from
        // It is symmetric, but things are still stored in the full matrix. The
        // storage schemes are different, so I can't make it symmetrized
        // inline -- though one could create a copy of this with is_symmetric
        // true.
      } // setToSymmeterizedKullbackLeiblerDivergences()

      /**
       * @brief Reset the distances to the euclidean distances among the given
       * alignment profiles.
       *
       * This will make the distance matrix symmetric as a side-effect.
       */
      void
      setToEuclideanDistances (
        vector<AlignmentProfile> const & alignment_profiles
      )
      {
        // Reinitialize to the size of the given vector, YES, symmetric,
        // NO, not a proximity (yes, a distance)
        this->reinitialize( alignment_profiles.size(), true, false );

        uint32_t size = this->size();
        for( uint32_t from = 0; from < size; from++ ) {
          for( uint32_t to = 0; to < from; to++ ) {
            this->operator()( from, to ) =
              static_cast<DistanceType>( alignment_profiles[ from ].euclideanDistance( alignment_profiles[ to ] ) );
          } // End foreach to
          this->operator()( from, from ) = 0;
        } // End foreach from
      } // setToEuclideanDistances( vector<AlignmentProfile> const & )

      /**
       * @brief Reset the distances to the euclidean distances among the given
       * alignment profiles, just among the Match Emission distributions at position pos_i.
       *
       * This will make the distance matrix symmetric as a side-effect.
       */
      void
      setToEuclideanDistances (
        vector<AlignmentProfile> const & alignment_profiles,
        uint32_t const & pos_i
      )
      {
        // Reinitialize to the size of the given vector, YES, symmetric,
        // NO, not a proximity (yes, a distance)
        this->reinitialize( alignment_profiles.size(), true, false );

        uint32_t size = this->size();
        for( uint32_t from = 0; from < size; from++ ) {
          for( uint32_t to = 0; to < from; to++ ) {
            this->operator()( from, to ) =
              // TODO: Include something other than just the Match Emission dist?
              static_cast<DistanceType>(static_cast<MultinomialDistribution<ResidueType,ProbabilityType> >( alignment_profiles[ from ][ pos_i ][ Emission::Match ]).euclideanDistance( alignment_profiles[ to ][ pos_i ][ Emission::Match ] ) );
          } // End foreach to
          this->operator()( from, from ) = 0;
        } // End foreach from
      } // setToEuclideanDistances( vector<AlignmentProfile> const &, uint32_t const & )

      /**
       * Calculate and return the sum of the off-diagonal values of the matrix.
       *
       * NOTE that this will count BOTH i->k and k->i, even if the matrix is
       * symmetric.
       */
      DistanceType
      totalOffDiagonal ()
      {
        return totalOffDiagonal( 0, this->size() );
      } // totalOffDiagonal()

      /**
       * Calculate and return the sum of the off-diagonal values of the matrix.
       *
       * NOTE that this will count BOTH i->k and k->i, even if the matrix is
       * symmetric.
       */
      DistanceType
      totalOffDiagonal (
        uint32_t const & first_element,
        uint32_t const & num_elements
      )
      {
        DistanceType total = 0;
        if( num_elements == 1 ) {
          //  ?
          return 0;
        }
        assert( ( first_element + num_elements - 1 ) < this->size() );
        DistanceType max = this->operator()( first_element, ( first_element + 1 ) );
        DistanceType val;
        for( uint32_t from = first_element; from < ( first_element + num_elements ); from++ ) {
          for( uint32_t to = first_element; to < ( m_isSymmetric ? from : ( first_element + num_elements ) ); to++ ) {
            if( to == from ) {
              continue;
            }
            total += this->operator()( from, to );
          } // End foreach to
        } // End foreach from
        if( m_isSymmetric ) {
          total *= 2;
        }
        return total;
      } // totalOffDiagonal( uint32_t const &, uint32_t const & )
      
      /**
       * Calculate and return the average of the off-diagonal values of the matrix.
       */
      DistanceType
      averageOffDiagonal ()
      {
        return averageOffDiagonal( 0, this->size() );
      } // averageOffDiagonal()

      /**
       * Calculate and return the average of the off-diagonal values of the matrix.
       */
      DistanceType
      averageOffDiagonal (
        uint32_t const & first_element,
        uint32_t const & num_elements
      )
      {
        DistanceType average = this->totalOffDiagonal( first_element, num_elements );
        average /= ( num_elements * ( num_elements - 1 ) );

        return average;
      } // averageOffDiagonal( uint32_t const &, uint32_t const & )

      /**
       * Calculate and return the variance of the off-diagonal values of the
       * matrix.  Note that this estimates the average ( using
       * averageOffDiagonal() ) so we use the -1 adjustment for an unbiased
       * estimate of the variance.
       */
      DistanceType
      varianceOffDiagonal ()
      {
        return varianceOffDiagonal( 0, this->size() );
      } // varianceOffDiagonal()

      /**
       * Calculate and return the variance of the off-diagonal values of the
       * matrix.  Note that this estimates the average ( using
       * averageOffDiagonal() ) so we use the -1 adjustment for an unbiased
       * estimate of the variance.
       */
      DistanceType
      varianceOffDiagonal (
        uint32_t const & first_element,
        uint32_t const & num_elements
      )
      {
        DistanceType const average = this->averageOffDiagonal( first_element, num_elements );
        return varianceOffDiagonal( first_element, num_elements, average );
      } // varianceOffDiagonal( uint32_t const &, uint32_t const & )

      /**
       * Calculate and return the variance of the off-diagonal values of the
       * matrix.  Note that this estimates the average ( using
       * averageOffDiagonal() ) so we use the -1 adjustment for an unbiased
       * estimate of the variance.
       */
      DistanceType
      varianceOffDiagonal (
         DistanceType const & average
      )
      {
        return varianceOffDiagonal( 0, this->size(), average );
      } // varianceOffDiagonal( DistanceType const & )

      /**
       * Calculate and return the variance of the off-diagonal values of the
       * matrix.  Note that this estimates the average ( using
       * averageOffDiagonal() ) so we use the -1 adjustment for an unbiased
       * estimate of the variance.
       */
      DistanceType
      varianceOffDiagonal (
        uint32_t const & first_element,
        uint32_t const & num_elements,
        DistanceType const & average
      )
      {
        DistanceType var = 0;
        if( num_elements == 1 ) {
          //  ?
          return 0;
        }
        assert( ( first_element + num_elements - 1 ) < this->size() );
        DistanceType max = this->operator()( first_element, ( first_element + 1 ) );
        DistanceType val;
        for( uint32_t from = first_element; from < ( first_element + num_elements ); from++ ) {
          for( uint32_t to = first_element; to < ( m_isSymmetric ? from : ( first_element + num_elements ) ); to++ ) {
            if( to == from ) {
              continue;
            }
            val = this->operator()( from, to );
            val -= average;
            var += ( val * val );
          } // End foreach to
        } // End foreach from
        if( m_isSymmetric ) {
          var *= 2;
        }
        var /= ( ( num_elements * ( num_elements - 1 ) ) - 1 ); // -1 adjustment ...
        //var /= ( num_elements - 1 ); // -1 adjustment ...

        return var;
      } // varianceOffDiagonal( uint32_t const &, uint32_t const &, DistanceType const & )

      /**
       * Return the largest value in the off-diagonal elements of this matrix.
       * Note that this does not depend on whether this is a proximity or a
       * distance matrix; it just returns the largest value.
       */
      DistanceType
      maximumOffDiagonalValue ()
      {
        return maximumOffDiagonalValue( 0, this->size() );
      } // maximumOffDiagonalValue()

      /**
       * Return the largest value in the off-diagonal elements of this matrix,
       * including only those elements with indices between first_index and ( first_index +
       * num_elements - 1 ), inclusive.
       *
       * Note that this does not depend on whether this is a proximity or a
       * distance matrix; it just returns the largest value.
       */
      DistanceType
      maximumOffDiagonalValue (
        uint32_t const & first_element,
        uint32_t const & num_elements
      )
      {
        if( num_elements == 1 ) {
          //  ?
          return 0;
        }
        assert( ( first_element + num_elements - 1 ) < this->size() );
        DistanceType max = this->operator()( first_element, ( first_element + 1 ) );
        DistanceType val;
        for( uint32_t from = first_element; from < ( first_element + num_elements ); from++ ) {
          for( uint32_t to = first_element; to < ( m_isSymmetric ? from : ( first_element + num_elements ) ); to++ ) {
            if( to == from ) {
              continue;
            }
            val = this->operator()( from, to );
            if( val > max ) {
              max = val;
            }
          } // End foreach to
        } // End foreach from
        return max;
      } // maximumOffDiagonalValue( uint32_t const &, uint32_t const & )

      /**
       * Return the smallest value in the off-diagonal elements of this matrix.
       * Note that this does not depend on whether this is a proximity or a
       * distance matrix; it just returns the smallest value.
       */
      DistanceType
      minimumOffDiagonalValue ()
      {
        return minimumOffDiagonalValue( 0, this->size() );
      } // minimumOffDiagonalValue()

      /**
       * Return the smallest value in the off-diagonal elements of this matrix,
       * including only those elements with indices between first_index and ( first_index +
       * num_elements - 1 ), inclusive.
       *
       * Note that this does not depend on whether this is a proximity or a
       * distance matrix; it just returns the smallest value.
       */
      DistanceType
      minimumOffDiagonalValue (
        uint32_t const & first_element,
        uint32_t const & num_elements
      )
      {
        if( num_elements == 1 ) {
          //  ?
          return 0;
        }
        assert( ( first_element + num_elements - 1 ) < this->size() );
        DistanceType min = this->operator()( first_element, ( first_element + 1 ) );
        DistanceType val;
        for( uint32_t from = first_element; from < ( first_element + num_elements ); from++ ) {
          for( uint32_t to = first_element; to < ( m_isSymmetric ? from : ( first_element + num_elements ) ); to++ ) {
            if( to == from ) {
              continue;
            }
            val = this->operator()( from, to );
            if( val < min ) {
              min = val;
            }
          } // End foreach to
        } // End foreach from
        return min;
      } // minimumOffDiagonalValue( uint32_t const &, uint32_t const & )

    protected:
      uint32_t
      getBackingVectorIndex ( uint32_t const & from, uint32_t const & to ) const
      {
        assert( from < m_size );
        assert( to < m_size );
        if( m_isSymmetric ) {
          if( from < to ) {
            // swap.  Lower-triangular so we need to <= from.
            return ( ( ( ( to + 1 ) * to ) / 2 ) + from );
          } else { // from >= to
            return ( ( ( ( from + 1 ) * from ) / 2 ) + to );
          }
        } else { // if m_isSymmetric .. else ..
          return ( ( from * m_size ) + to );
        } // End if m_isSymmetric .. else ..
      } // getBackingVectorIndex( uint32_t const &, uint32_t const & ) const

    }; // End inner class DynamicProgramming::DistanceMatrix

  public:
    DynamicProgramming ()
    {
      // Do nothing else.
    } // <init>()

    /**
     * Given the forward viterbi matrices (Calculated as by the
     * forward_score_viterbi(..) method), calculate and return the
     * MultipleAlignment of the sequences with the profile model.
     */
    template <typename ProfileType, typename SequenceResidueType>
    void
    forward_viterbiAlign (
      Parameters const& parameters,
      typename Matrix::SequentialAccessContainer const & viterbi_matrices,
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<ProfileType, SequenceResidueType> & ma
    ) const;

    /**
     * Calculate the minimum-cost alignment between the two profiles, using the
     * symmeterized KL divergence as the distance metric, and the given indel
     * costs.  Fill the given uint32_t vector (which will be resized to one
     * greater than the length of profile a) with the number of advances of
     * profile b for each position of profile a (index i corresponds to
     * position i-1 of profile a, with index 0 indicating the number of
     * positions of profile b that align before all positions of profile a).
     * After calling this method, that vector's values will sum to the length
     * of profile b.  The return value is the cost of the computed alignment.
     *
     * Note: we do not allow I->D or D->I transitions (that is, a match must be
     * a match, not an insertion by both profiles simultaneously).
     */
    template <typename ProfileTypeA, typename ProfileTypeB>
    double
    profileProfile_align_SKL (
      Parameters const& parameters,
      ProfileTypeA const& profile_a,
      ProfileTypeB const& profile_b,
      double const indel_open_cost,
      double const indel_extension_cost,
      vector<uint32_t> & profile_b_advances,
      double const maximum_match_cost = 0 // no maximum
    ) const;

    template <typename ProfileTypeA, typename ProfileTypeB>
    ScoreType
    calculatePathCooccurrenceProbability (
      Parameters const& parameters,
      ProfileTypeA const& profile_a,
      ProfileTypeB const& profile_b
    ) const;

    /**
     * Given the forward row for row_i and the current SamplerState (with the
     * current column and substate of the sampler for the given sequence,
     * representing the beginning of the partial path from model positions
     * row_i through (profile.length() - 1), which has already been drawn),
     * draw more of the partial path, corresponding to position (row_i-1) of
     * the given profile.  The given target SamplerState reference will be
     * modified to hold the sampler state after drawing the partial path (which
     * is sufficient information, along with the 'current' SamplerState, to
     * reconstruct that partial path).
     *
     * If the GlobalCountsType and/or PositionCountsType argument
     * pointers are non-null, then they should be a Profile type and a
     * ProfilePosition type, storing unsigned integer counts.  These
     * counts will be updated with the number of the observations on the drawn
     * partial path corresponding to each parameter.
     *
     * If the final argument (the sequence_advances) pointer is non-NULL, then
     * it should be a vector of length ( profile.length() + 4 ).  It will
     * be modified to reflect the drawn partial path.  Note that Matches into
     * row row_i won't be counted until this method is called for the previous
     * row (row_i - 1).
     */
    template <typename ProfileType,
              typename GlobalCountsType,
              typename PositionCountsType,
              typename SequenceResidueType>
    void
    drawPartialPathForPosition (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& forward_row,
      SamplerState const& current_sampler_state,
      SamplerState & target_sampler_state,
      Random & random,
      GlobalCountsType * global_counts,
      PositionCountsType * position_counts,
      vector<uint32_t> * sequence_advances
    ) const;

    /**
     * Given the forward matrices corresponding to the given profile and
     * sequences, draw paths from the distribution of the paths, given the
     * sequences and profiles; the results are stored in either the given
     * counts arguments or the given multiple alignment argument (one of which
     * should be non-null, or this method is a waste of time).
     *
     * If the CountsType argument pointer is non-null, then it should be a
     * Profile type with the same length as the profile argument, used for
     * storing unsigned integer counts.  These counts will be updated with the
     * number of the observations on the drawn partial path corresponding to
     * each parameter.
     *
     * If the given MultipleAlignment pointer is non-null, then it will be
     * reinitialized and modified to hold the drawn paths.
     */
    template <typename ProfileType,
              typename CountsType,
              typename SequenceResidueType>
    void
    drawPaths (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const & sequences,
      const uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer const& forward_matrices,
      Random & random,
      CountsType * counts,
      MultipleAlignment<ProfileType, SequenceResidueType> * multiple_alignment
    ) const;

    
    /**
     * Given the backward row for row row_i, and the forward row for row_i - 1,
     * calculate the PositionSpecificSequenceScoreCoefficients for position
     * (row_i-1) of the given profile for the given sequence.  Results go in
     * the given PositionSpecificSequenceScoreCoefficients reference.  If
     * parameters.useRabinerScaling is true, the result will be scaled (by the
     * inverse of the product of the rabinerCumulativeInverseScalar values in
     * the given Row references, which is the sequence score).  The inverse of
     * that scalar will be stored in the m_inverseScalar field of the
     * coefficients.  An unscaled version can be retrieved by calling
     * createUnscaledCopy() on the coefficients.
     *
     * NOTE: row_i must be > 0.
     */
    // Note we could one day change the pos.spec. coeffs to include the deletions (easily) and even insertions (but we'd need a vector of the length of the sequence, (times the emissions!) to store every possible combination of insertion extension and match emission.  This is not considering the insertion emissions, but that's usually fine to consider as fixed.
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
    void
    calculatePositionSpecificSequenceScoreCoefficients (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& backward_row,
      PositionSpecificSequenceScoreCoefficients & coefficients
    ) const;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    // Delegating method that creates the prev_row_match_distribution...
    template <typename ProfileType,
              typename SequenceResidueType>
    void
    calculatePositionSpecificSequenceScoreCoefficients (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& backward_row,
      PositionSpecificSequenceScoreCoefficients & coefficients
    ) const
    {
      // We need to use scaled versions of the TransitionFromMatch distributions.
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> prev_row_match_distribution;
      if( ( row_i >= 2 ) && ( ( row_i - 2 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 2 ),
          prev_row_match_distribution
        );
      } // End if row_i >= 2
      calculatePositionSpecificSequenceScoreCoefficients (
        parameters,
        profile,
        prev_row_match_distribution,
        sequence,
        row_i,
        prev_forward_row,
        backward_row,
        coefficients
      );

      return;
    }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    /**
     * Use the given PositionSpecificSequenceScoreCoefficients reference to
     * update the entente position corresponding to the given profile position.
     * The entente position will be updated by adding the expected number of
     * emissions of each type (divided by the given inverse scalar, if it is
     * non-NULL, or by the calculated sequence score otherwise).  The return
     * value is the (not scaled) sequence score.  Don't forget to zero() the
     * position entente (before calling this for all sequences).  If the
     * inverse_scalar pointer is non-NULL but points to the value 0, the score
     * will be calculated and returned but the position entente will not be
     * updated (it's like having a 0 scalar, not like having a 0 inverse
     * scalar, which would result in a divide-by-0 error).
     */
    template <typename ProfileType>
    ScoreType
    updatePositionEntenteForSequence (
      Parameters const& parameters,
      ProfilePosition<ResidueType, ProfileType> const& profile_position,
      PositionSpecificSequenceScoreCoefficients const & coefficients,
      ScoreType const * inverse_scalar,
      PositionEntente & position_entente
    ) const;

    /**
     * Use the given PositionSpecificSequenceScoreCoefficients reference to
     * calculate the update for the entente position corresponding to the given
     * profile position.  The calculated update should be added to the entente
     * position.  This will add the expected number of emissions of each type
     * (divided by the given inverse scalar, if it is non-NULL, or by the
     * calculated sequence score otherwise).  If the inverse_scalar pointer is
     * non-NULL but points to the value 0, the score will be calculated and
     * returned but the position_entente_update will be zero()ed (it's like
     * having a 0 scalar, not like having a 0 inverse scalar, which would
     * result in a divide-by-0 error).  The return value is the (not scaled)
     * sequence score.
     */
    template <typename ProfileType>
    ScoreType
    calculatePositionEntenteUpdate (
      Parameters const& parameters,
      ProfilePosition<ResidueType, ProfileType> const& profile_position,
      PositionSpecificSequenceScoreCoefficients const & coefficients,
      ScoreType const * inverse_scalar,
      PositionEntente & position_entente_update
    ) const;

    /**
     * Given the forward and backward rows for row row_i, and the forward row
     * for row_i - 1, update the given GlobalEntente for the given sequence,
     * profile pair.  Results go in the given GlobalEntente reference.  Don't
     * forget to zero() the entente (before calling this for all rows and all
     * sequences).
     * NOTE: If row_i is 0, prev_forward_row is ignored.
     */
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
    ScoreType
    updateGlobalEntenteForSequence (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& forward_row,
      typename Matrix::Row const& backward_row,
      ScoreType const * inverse_scalar,
      GlobalEntente & global_entente
    ) const;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    // creates prev_row_match_distribution and row_match_distribution and
    // delegates.
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    updateGlobalEntenteForSequence (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& forward_row,
      typename Matrix::Row const& backward_row,
      ScoreType const * inverse_scalar,
      GlobalEntente & global_entente
    ) const
    {
      // We need to use scaled versions of the TransitionFromMatch distributions.
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> prev_row_match_distribution;
      if( ( row_i >= 2 ) && ( ( row_i - 2 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 2 ),
          prev_row_match_distribution
        );
      } // End if row_i >= 2
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> row_match_distribution;
      if( ( row_i >= 1 ) && ( ( row_i - 1 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 1 ),
          row_match_distribution
        );
      } // End if row_i >= 1

      return
        updateGlobalEntenteForSequence (
          parameters,
          profile,
          prev_row_match_distribution,
          row_match_distribution,
          sequence,
          row_i,
          prev_forward_row,
          forward_row,
          backward_row,
          inverse_scalar,
          global_entente
        );
    } // updateGlobalEntenteForSequence(..)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    /**
     * Given the forward matrices for a set of Sequences, calculate their
     * AlignmentProfiles, and place the result in the given vector ref.
     */
    template <typename ProfileType, typename SequenceType>
    void
    calculateAlignmentProfiles (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<SequenceType> const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer const& forward_matrices,
      vector<AlignmentProfile> & alignment_profiles
    ) const;

    /**
     * Given the forward matrices for a set of Sequences, calculate their
     * AlignmentProfiles, and place the result in the given vector ref.  The
     * given RowVectors are for temporary storage -- they will be clobbered.
     * If the scores are not provided, and if useRabinerScaling is true and
     * rabinerScaling_useMaximumValue is also true, then they will be
     * calculated (each is the cumulativeInverseScalar of the last row of the
     * forward matrix) -- otherwise they are not needed.
     */
    template <typename ProfileType, typename SequenceType>
    void
    calculateAlignmentProfiles (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<SequenceType> const& sequences,
      uint32_t const& sequence_count,
      vector<ScoreType> const * sequence_scores,
      typename Matrix::SequentialAccessContainer const& forward_matrices,
      typename Matrix::RowVector & backward_rows_1,
      typename Matrix::RowVector & backward_rows_2,
      vector<AlignmentProfile> & alignment_profiles
    ) const;

    /**
     * Given the forward and backward rows for row row_i, and the forward row
     * for ( row_i - 1 ) and the backward row for ( row_i + 1 ), calculate the
     * expected uses of each profile parameter at profile position ( row_i - 1
     * ) for the given sequence, profile pair.  Results will be divided by the
     * sequence score or by the given inverse_scalar (if it is non-null).
     * Results go in the given AlignmentProfilePosition reference.
     * If the coefficients pointer is non-NULL, the coefficients will be
     * calculated too (as if you'd called
     * calculatePositionSpecificSequenceScoreCoefficients(..).
     *
     * NOTE: If useRabinerScaling is true but rabinerScaling_useMaximumValue is
     * also true, then you *MUST* supply a value for inverse_scalar.  We can't
     * calculate it here.  Likewise when the backward_row is not up-to-date, as
     * happens in ProfileTrainer.hpp.
     * 
     * NOTE: If row_i is 0, prev_forward_row is ignored.  If row_i == last_row,
     * next_backward_row is ignored.
     */
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
    ScoreType
    calculateAlignmentProfilePosition (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& forward_row,
      typename Matrix::Row const& backward_row,
      typename Matrix::Row const& next_backward_row,
      ScoreType const * inverse_scalar,
      AlignmentProfilePosition & pos,
      PositionSpecificSequenceScoreCoefficients * coefficients
    ) const;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    // creates prev_row_match_distribution and row_match_distribution and delegates.
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    calculateAlignmentProfilePosition (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& forward_row,
      typename Matrix::Row const& backward_row,
      typename Matrix::Row const& next_backward_row,
      ScoreType const * inverse_scalar,
      AlignmentProfilePosition & pos,
      PositionSpecificSequenceScoreCoefficients * coefficients
    ) const
    {
      // We need to use scaled versions of the TransitionFromMatch distributions.
      // TODO: REMOVE
      //cout << "Hi from     calculateAlignmentProfilePosition (..) (USE_DEL_IN_DEL_OUT)" << endl;
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> prev_row_match_distribution;
      if( ( row_i >= 2 ) && ( ( row_i - 2 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 2 ),
          prev_row_match_distribution
        );
      } // End if row_i >= 2
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> row_match_distribution;
      if( ( row_i >= 1 ) && ( ( row_i - 1 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 1 ),
          row_match_distribution
        );
      } // End if row_i >= 1

      return
        calculateAlignmentProfilePosition(
          parameters,
          profile,
          prev_row_match_distribution,
          row_match_distribution,
          sequence,
          row_i,
          prev_forward_row,
          forward_row,
          backward_row,
          next_backward_row,
          inverse_scalar,
          pos,
          coefficients
        );
    } // calculateAlignmentProfilePosition(..)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef ALLOW_BOLTZMANN_GIBBS
    /**
     * Use the given ProfilePosition reference
     * and Coefficients reference to update the
     * PositionBoltzmannGibbs corresponding to the given profile position
     * with Baldi's gradient ascent step (equation 3 from Pierre Baldi, Yves
     * Chauvin).  If the given inverse scalar is non-NULL, it will be used
     * instead of the calculated sequence score for scaling the change.  The
     * return value is the (not scaled) sequence score.  Note that the learning
     * rate will *not* be incorporated; the caller is responsible for
     * multiplying the resulting PositionBoltzmannGibbs by the learning
     * rate afterwards (or equivalently dividing its m_scalar by the learning rate). 
     *
     * Don't forget to zero() the PositionBoltzmannGibbs (before calling
     * this for all sequences).  If the inverse_scalar pointer is non-NULL but
     * points to the value 0, the score will be calculated and returned but the
     * position_boltzmann_gibbs_change will not be updated (it's like having a 0 scalar, not
     * like having a 0 inverse scalar, which would result in a divide-by-0
     * error).
     */
    template <typename ProfileType>
    ScoreType
    updatePositionBoltzmannGibbsChangeForSequence (
      Parameters const& parameters,
      ProfilePosition<ResidueType, ProfileType> const& profile_position,
      PositionSpecificSequenceScoreCoefficients const & coefficients,
      ScoreType const * inverse_scalar,
      PositionBoltzmannGibbs<ResidueType, ScoreType> & position_boltzmann_gibbs_change
    ) const;
#endif // ALLOW_BOLTZMANN_GIBBS

    /**
     * Calculate the distances among the given sequence ententes, and replace
     * the values of the given PositionEntenteDistances object with these
     * distances.
     */
    void
    updatePositionEntenteDistances (
      vector<PositionEntente> const& position_ententes,
      PositionEntenteDistances & position_entente_distances
    ) const;

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile tree, using the forward algorithm.  The sequences
     * vector must be of length profile_tree.nodeCount().  Use the given
     * matrices, which must be large enough to accommodate the largest set of
     * sequences.
     */
    template <typename ProfileTreeType, typename SequenceType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileTreeType const& profile_tree,
      vector<SequenceType> const& sequences,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      return
        forward_score(
          parameters,
          false, // Don't use viterbi
          profile_tree,
          sequences,
          matrices
        );
    } // forward_score( Parameters const&, ProfileTree const&, vector<SequenceType> const&, SequentialAccessContainer & ) const

    /**
     * Calculate and return the viterbi score of the given sequences, given
     * the given profile tree, using the viterbi algorithm.  The sequences
     * vector must be of length profile_tree.nodeCount().  Use the given
     * matrices, which must be large enough to accommodate the largest set of
     * sequences.
     */
    template <typename ProfileTreeType, typename SequenceType>
    ScoreType
    forward_score_viterbi (
      Parameters const& parameters,
      ProfileTreeType const& profile_tree,
      vector<SequenceType> const& sequences,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      return
        forward_score(
          parameters,
          true, // Use viterbi
          profile_tree,
          sequences,
          matrices
        );
    } // forward_score_viterbi( Parameters const&, ProfileTree const&, vector<SequenceType> const&, SequentialAccessContainer & ) const

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile tree, using the forward or viterbi algorithm.  The
     * sequences vector must be of length profile_tree.nodeCount().  Use the
     * given matrices, which must be large enough to accommodate the largest
     * set of sequences.
     */
    template <typename ProfileTreeType, typename SequenceType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileTreeType const& profile_tree,
      vector<SequenceType> const& sequences,
      typename Matrix::SequentialAccessContainer & matrices
    ) const;

    /**
     * Calculate and return the likelihood of the given sequence, given the
     * given profile model, using the forward algorithm.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence
    ) const
    {
      typename Matrix::Row prev_forward_row =
        typename Matrix::Row( sequence.length() + 1 );
      typename Matrix::Row forward_row =
        typename Matrix::Row( sequence.length() + 1 );
      ScoreType score;

      forward_score(
        parameters,
        profile,
        sequence,
        prev_forward_row,
        forward_row,
        score
      );
      return score;
    } // forward_score( Parameters const&, Profile const&, Sequence<SequenceResidueType> const& ) const

    /**
     * Calculate the likelihood of the sequences whose forward and backward
     * rows are given.  They are assumed to have been pre-calculated (see
     * forward_calculateRow(..) and backward_calculateRow(..)), for the same
     * profile, sequences, and row number.  As a side effect, the given
     * ScoreType vector (if non-NULL) becomes filled with the scores of the
     * sequences, and the largest_score pointer (if non-null) will be filled
     * with the largest of the scores of the sequences.  We need to know if it
     * is the last row because we store things differently there..
     */
    template <typename ProfileType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      uint32_t const& sequence_count,
      uint32_t const row_i,
      typename Matrix::RowVector const & forward_rows,
      typename Matrix::RowVector const & backward_rows,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score
    ) const;

    /**
     * Calculate the likelihood of the sequence whose forward and backward
     * rows are given.  They are assumed to have been pre-calculated (see
     * forward_calculateRow(..) and backward_calculateRow(..)), for the same
     * profile, sequence, and row number.  We need to know if it is the last
     * row because we store things differently there..
     */
    template <typename ProfileType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      uint32_t const row_i,
      typename Matrix::Row const & forward_row,
      typename Matrix::Row const & backward_row
    ) const;

    /**
     * Calculate the likelihood of the given sequence, given the given profile
     * model, using the given Rows as temporary storage (these will be modified
     * such that afterwards, if the return value is false then forward_row will
     * hold the forward matrix row values for the last row (at index
     * profile.length()), and prev_forward_row will hold the values for the row
     * just before that (at index profile.length()-1)), if the return value is
     * true then these rows will be swapped.  The side-effect, besides changing
     * the forward rows, is to change the value in the score reference
     * argument.  The return value indicates whether the rows are swapped.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    bool
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      typename Matrix::Row & prev_forward_row,
      typename Matrix::Row & forward_row,
      ScoreType & score
    ) const;

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.  As a side effect,
     * the given matrices references becomes filled with the forward matrix
     * values.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      return
        forward_score(
          parameters,
          false, // Don't use viterbi
          profile,
          sequences,
          sequence_count,
          matrices,
          NULL,
          NULL
        );
    } // forward_score( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer & ) const

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.  As a side effect,
     * the given matrices references becomes filled with the forward matrix
     * values, and the the given ScoreType vector (if non-NULL) becomes filled
     * with the scores of the sequences.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices,
      vector<ScoreType> * sequence_scores
    ) const
    {
      return
        forward_score(
          parameters,
          false, // Don't use viterbi
          profile,
          sequences,
          sequence_count,
          matrices,
          sequence_scores,
          NULL
        );
    } // forward_score( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer &, vector<ScoreType> * ) const

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.  As a side effect,
     * the given matrices references becomes filled with the forward matrix
     * values, and the the given ScoreType vector (if non-NULL) becomes filled
     * with the scores of the sequences, and the largest_score pointer (if
     * non-null) will be filled with the largest of the scores of the
     * sequences.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score
    ) const
    {
      return
        forward_score(
          parameters,
          false, // Don't use viterbi
          profile,
          sequences,
          sequence_count,
          matrices,
          sequence_scores,
          largest_score
        );
    } // forward_score( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer &, vector<ScoreType> *, ScoreType * ) const

    /**
     * Calculate and return the product of the likelihoods of the most likely
     * paths for the given sequences, given the given profile model, using the
     * viterbi forward algorithm.  As a side effect, the given matrices
     * references becomes filled with the viterbi forward matrix values.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score_viterbi (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      return
        forward_score(
          parameters,
          true, // Use viterbi
          profile,
          sequences,
          sequence_count,
          matrices
        );
    } // forward_score_viterbi ( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer & ) const

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.  As a side effect,
     * the given matrices references becomes filled with the forward matrix
     * values.  If the use_viterbi argument is true, calculate most-likely-path
     * likelihoods rather than full (all path) likelihoods.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      return forward_score(
        parameters,
        use_viterbi,
        profile,
        sequences,
        sequence_count,
        matrices,
        NULL,
        NULL
      );
    } // forward_score( Parameters const&, bool, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer & ) const

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.  As a side effect,
     * the given matrices references becomes filled with the forward matrix
     * values, and the the given ScoreType vector (if non-NULL) becomes filled
     * with the scores of the sequences.  If the use_viterbi argument is true,
     * calculate most-likely-path likelihoods rather than full (all path)
     * likelihoods.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices,
      vector<ScoreType> * sequence_scores
    ) const
    {
      return forward_score(
        parameters,
        use_viterbi,
        profile,
        sequences,
        sequence_count,
        matrices,
        sequence_scores,
        NULL
      );
    } // forward_score( Parameters const&, bool, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer &, vector<ScoreType> * ) const

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.  As a side effect,
     * the given matrices references becomes filled with the forward matrix
     * values, and the the given ScoreType vector (if non-NULL) becomes filled
     * with the scores of the sequences, and the largest_score pointer (if
     * non-null) will be filled with the largest of the scores of the
     * sequences.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer & matrices,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score
    ) const;

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count
    ) const;

    /**
     * Calculate the total likelihood of the given sequences, given the given
     * profile model, using the given Rows as temporary storage (these will be
     * modified such that afterwards, if the return value is false then
     * forward_rows will hold the forward matrix row values for the last rows
     * (at index profile.length()), and prev_forward_rows will hold the values
     * for the rows just before that (at index profile.length()-1)).  If the
     * return value is true then these rows will be swapped.  The side-effect,
     * besides changing the forward rows, is to change the value in the score
     * reference argument.  The return value indicates whether the rows are
     * swapped. If parameters.useRabinerScaling is
     * true, the rabinerInverseScaling and rabinerCumulativeInverseScalar
     * values of the forward rows will be set, too. If the use_viterbi argument
     * is true, calculate most-likely-path likelihoods rather than full (all
     * path) likelihoods.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    bool
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::RowVector & prev_forward_rows,
      typename Matrix::RowVector & forward_rows,
      ScoreType & score
    ) const
    {
      return forward_score(
        parameters,
        use_viterbi,
        profile,
        sequences,
        sequence_count,
        ( typename Matrix::SequentialAccessContainer * )( NULL ),
        ( typename Matrix::SequentialAccessContainer * )( NULL ),
        prev_forward_rows,
        forward_rows,
        ( vector<ScoreType> * )( NULL ),
        ( ScoreType * )( NULL ),
        score
      );
    } // forward_score( Parameters const&, bool, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, RowVector &, RowVector &, ScoreType & ) const

    /**
     * Calculate the total likelihood of the given sequences, given the given
     * profile model, using the given Rows as temporary storage (these will be
     * modified such that afterwards, if the return value is false then
     * forward_rows will hold the forward matrix row values for the last rows
     * (at index profile.length()), and prev_forward_rows will hold the values
     * for the rows just before that (at index profile.length()-1)).  If the
     * return value is true then these rows will be swapped.  One side-effect,
     * besides changing the forward rows, is to change the value in the score
     * reference argument.  The return value indicates whether the rows are
     * swapped. If the given ScoreType vector is non-NULL, it will become
     * filled with the scores of the sequences, and the largest_score pointer
     * (if non-null) will be filled with the largest of the scores of the
     * sequences.  If parameters.useRabinerScaling is
     * true, the rabinerInverseScaling and rabinerCumulativeInverseScalar
     * values of the forward rows will be set, too.  If the anchor_columns pointer is non-null, then anchor columns (to be
     * used with forward_reverseCalculateRow(..) will be stored also (using the
     * member var anchor_columns->m_storeEveryNthColumn).  If the anchor_rows pointer is non-null, then anchor rows will be
     * stored also (using the member var anchor_rows->m_storeEveryNthRow).  If
     * the use_viterbi argument is true, calculate most-likely-path likelihoods
     * rather than full (all path) likelihoods.
     */
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    template <typename ProfileType,
              typename SequenceResidueType,
              typename ScaledMatchDistributionProbabilityType>
#else
    template <typename ProfileType,
              typename SequenceResidueType>
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
    bool
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> > const & scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer * anchor_columns,
      typename Matrix::SequentialAccessContainer * anchor_rows,
      typename Matrix::RowVector & prev_forward_rows,
      typename Matrix::RowVector & forward_rows,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score,
      ScoreType & score
    ) const;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    // creates the scaled_match_distributions and delegates.
    template <typename ProfileType,
              typename SequenceResidueType>
    bool
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer * anchor_columns,
      typename Matrix::SequentialAccessContainer * anchor_rows,
      typename Matrix::RowVector & prev_forward_rows,
      typename Matrix::RowVector & forward_rows,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score,
      ScoreType & score
    ) const
    {
      vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> > scaled_match_distributions( profile.length() - 1 );
      for( uint32_t pos_i = 0; pos_i < ( profile.length() - 1 ); pos_i++ ) {
        profile.createScaledMatchDistributionForPosition(
          pos_i,
          scaled_match_distributions[ pos_i ]
        );
      } // End foreach pos_i

      return
        forward_score (
          parameters,
          use_viterbi,
          profile,
          scaled_match_distributions,
          sequences,
          sequence_count,
          anchor_columns,
          anchor_rows,
          prev_forward_rows,
          forward_rows,
          sequence_scores,
          largest_score,
          score
        );
    } // forward_score(..)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    /**
     * Using the forward row (of the Forward algorithm matrix) for (row_i - 1),
     * calculate the forward row for row_i.  Return the given forward_row
     * reference.  If parameters.useRabinerScaling is true, then the calculated
     * forward row will be scaled by the sum of all of its values.  That sum
     * will be stored in forward_row.m_rabinerInverseScalar.
     * NOTE: when row_i is last_row, forward_row.m_rabinerInverseScalar is the
     * last Match value times the postalign-to-terminal probability.
     * NOTE: when row_i is 0, prev_forward_row is ignored.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    typename Matrix::Row &
    forward_calculateRow (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row & forward_row
    ) const
    {
      return
        forward_calculateRow(
          parameters,
          false, // Don't use viterbi
          profile,
          sequence,
          row_i,
          prev_forward_row,
          forward_row
        );
    } // forward_calculateRow( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, Row & ) const

    /**
     * Using the forward row (of the Forward algorithm matrix) for (row_i - 1),
     * calculate the forward row for row_i.  Return the given forward_row
     * reference.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     * If parameters.useRabinerScaling is true, then
     * the calculated forward row will be scaled by the sum of all of its
     * values.
     * NOTE: when row_i is last_row, forward_row.m_rabinerInverseScalar is the
     * last Match value times the postalign-to-terminal probability.
     * NOTE: when row_i is 0, prev_forward_row is ignored.
     */
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
    typename Matrix::Row &
    forward_calculateRow (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row & forward_row
    ) const;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    // creates prev_row_match_distribution and row_match_distribution and delegates.
    template <typename ProfileType,
              typename SequenceResidueType>
    typename Matrix::Row &
    forward_calculateRow (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row & forward_row
    ) const
    {
      // We need to use scaled versions of the TransitionFromMatch distributions.
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> prev_row_match_distribution;
      if( ( row_i >= 2 ) && ( ( row_i - 2 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 2 ),
          prev_row_match_distribution
        );
      } // End if row_i >= 2
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> row_match_distribution;
      if( ( row_i >= 1 ) && ( ( row_i - 1 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 1 ),
          row_match_distribution
        );
      } // End if row_i >= 1

      return forward_calculateRow (
        parameters,
        use_viterbi,
        profile,
        prev_row_match_distribution,
        row_match_distribution,
        sequence,
        row_i,
        prev_forward_row,
        forward_row
      );
    } // forward_calculateRow (..)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

     /**
     * Using the forward row (of the Forward algorithm matrix) for (row_i + 1),
     * calculate the forward row for row_i.  Return the given forward_row
     * reference.  row_i must be not be 0, last_row, or last_row - 1.  The
     * values for the last column must be set ahead of time.  Note that this is
     * *not* backward_calculateRow, which computes "backward" values of the
     * "forward-backward" algorithm.  This method computes "forward" values,
     * but computes earlier forward values from later ones (ie. in "reverse").
     */
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
    typename Matrix::Row &
    forward_reverseCalculateRow (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & next_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      uint32_t const & store_every_Nth_column,
      typename Matrix::Row const& next_forward_row,
      typename Matrix::Row & forward_row
    ) const;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    template <typename ProfileType,
              typename SequenceResidueType>
    typename Matrix::Row &
    forward_reverseCalculateRow (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      uint32_t const & store_every_Nth_column,
      typename Matrix::Row const& next_forward_row,
      typename Matrix::Row & forward_row
    ) const
    {
      // We need to use scaled versions of the TransitionFromMatch distributions.
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> row_match_distribution;
      if( ( row_i >= 1 ) && ( ( row_i - 1 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 1 ),
          row_match_distribution
        );
      } // End if row_i >= 1
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> next_row_match_distribution;
      if( ( ( row_i + 1 ) >= 1 ) && ( ( ( row_i + 1 ) - 1 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          row_i,
          next_row_match_distribution
        );
      } // End if ( row_i + 1 ) >= 1

      return
        forward_reverseCalculateRow (
          parameters,
          profile,
          next_row_match_distribution,
          row_match_distribution,
          sequence,
          row_i,
          store_every_Nth_column,
          next_forward_row,
          forward_row
        );
    } // forward_reverseCalculateRow(..)

#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    /**
     * Calculate and return the likelihood of the given sequence, given the
     * given profile model, using the backward algorithm.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    backward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence
    ) const
    {
      return
        backward_score(
          parameters,
          profile,
          sequence,
          NULL
        );
    } // backward_score( Parameters const&, Profile const&, Sequence<SequenceResidueType> const& ) const

    /**
     * Calculate and return the likelihood of the given sequence, given the
     * given profile model, using the backward algorithm.  If the given
     * rabiner_inverse_scalars vector pointer is non-null, it should have
     * length one greater than the length of the profile, and should contain
     * the m_rabinerInverseScalars of the forward rows, calculated using
     * forward_calculateRow(..).
     * NOTE: If parameters.useRabinerScaling is true, you should provide a
     * non-null rabiner_inverse_scalars vector pointer, and vice-versa.  We
     * just assume that if it is non-null, then you mean us to use it.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    backward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const vector<MatrixValueType> * const rabiner_inverse_scalars
    ) const;

    /**
     * Calculate the likelihood of the given sequence, given the given profile
     * model, using the given Rows as temporary storage (these will be modified
     * such that afterwards, if the return value is false then backward_row
     * will hold the backward matrix row values for the first row (at index 0),
     * and prev_backward_row will hold the values for the row just after that
     * (at index 1), if the return value is true then these rows will be
     * swapped.  The side-effect, besides changing the backward rows, is to
     * change the value in the score reference argument.  The return value
     * indicates whether the rows are swapped.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    bool
    backward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      typename Matrix::Row & prev_backward_row,
      typename Matrix::Row & backward_row,
      ScoreType & score
    ) const
    {
      return
        backward_score(
          parameters,
          profile,
          sequence,
          NULL,
          prev_backward_row,
          backward_row,
          score
        );
    } // backward_score( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, Row &, Row &, ScoreType & ) const

    /**
     * Calculate the likelihood of the given sequence, given the given profile
     * model, using the given Rows as temporary storage (these will be modified
     * such that afterwards, if the return value is false then backward_row
     * will hold the backward matrix row values for the first row (at index 0),
     * and prev_backward_row will hold the values for the row just after that
     * (at index 1), if the return value is true then these rows will be
     * swapped.  The side-effect, besides changing the backward rows, is to
     * change the value in the score reference argument.  The return value
     * indicates whether the rows are swapped.  If the given
     * rabiner_inverse_scalars vector pointer is non-null, it should have
     * length one greater than the length of the profile, and should contain
     * the m_rabinerInverseScalars of the forward rows, calculated using
     * forward_calculateRow(..).
     * NOTE: If parameters.useRabinerScaling is true, you should provide a
     * non-null rabiner_inverse_scalars vector pointer, and vice-versa.  We
     * just assume that if it is non-null, then you mean us to use it.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    bool
    backward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const vector<MatrixValueType> * const rabiner_inverse_scalars,
      typename Matrix::Row & prev_backward_row,
      typename Matrix::Row & backward_row,
      ScoreType & score
    ) const;

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the backward algorithm.  As a side effect,
     * the given matrices references becomes filled with the backward matrix
     * values.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    backward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      return
        backward_score(
          parameters,
          false, // Don't use viterbi
          profile,
          sequences,
          sequence_count,
          matrices
        );

    } // backward_score( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer& ) const

    /**
     * Calculate and return the product of the likelihoods of the most likely
     * path of each of the given sequences, given the given profile model,
     * using the backward algorithm.  As a side effect, the given matrices
     * references becomes filled with the viterbi backward matrix values.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    backward_score_viterbi (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t const& sequence_count,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      return
        backward_score (
          parameters,
          true, // Use viterbi
          profile,
          sequences,
          sequence_count,
          matrices
        );

    } // backward_score_viterbi( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer& ) const

    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the backward algorithm.  As a side effect,
     * the given matrices references becomes filled with the backward matrix
     * values.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     * If useRabinerScaling is true, the m_rabinerInverseScalars of the matrices
     * will be used, so you must first compute the forward matrices (using,
     * eg. forward_score(..)), then call
     * matrices.copyRabinerInverseScalarsFrom( forward_matrices ).
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    ScoreType
    backward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer & matrices
    ) const;

    /**
     * Using the backward row (of the Backward algorithm matrix) for (row_i +
     * 1), calculate the backward row for row_i.  Return the given backward_row
     * reference.  If parameters.useRabinerScaling is true, then the calculated
     * forward row will be scaled by backward_row.m_rabinerInverseScalar.
     * backward_row.m_rabinerInverseScalar is the sum of the *forward* row
     * values for row_i, calculated using forward_calculateRow(..).
     * NOTE: when row_i is last_row, backward_row.m_rabinerInverseScalar is the
     * last Match value times the postalign-to-terminal probability.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    typename Matrix::Row &
    backward_calculateRow (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& next_backward_row,
      typename Matrix::Row & backward_row
    ) const
    {
      return
        backward_calculateRow(
          parameters,
          false, // Don't use viterbi
          profile,
          sequence,
          row_i,
          next_backward_row,
          backward_row
        );
    } // backward_calculateRow( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, Row & ) const
        
    /**
     * Using the backward row (of the Backward algorithm matrix) for (row_i +
     * 1), calculate the backward row for row_i.  Return the given backward_row
     * reference.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     * If parameters.useRabinerScaling is true, then
     * the calculated forward row will be scaled by
     * backward_row.m_rabinerInverseScalar.
     * backward_row.m_rabinerInverseScalar is the sum of the *forward* row
     * values for row_i, calculated using forward_calculateRow(..).
     * NOTE: when row_i is last_row, backward_row.m_rabinerInverseScalar is the
     * last Match value times the postalign-to-terminal probability.
     */
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
    typename Matrix::Row &
    backward_calculateRow (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& next_backward_row,
      typename Matrix::Row & backward_row
    ) const;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    // creates row_match_distribution and delegates.
    template <typename ProfileType,
              typename SequenceResidueType>
    typename Matrix::Row &
    backward_calculateRow (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& next_backward_row,
      typename Matrix::Row & backward_row
    ) const
    {
      // We need to use scaled versions of the TransitionFromMatch distributions.
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> row_match_distribution;
      if( ( row_i >= 1 ) && ( ( row_i - 1 ) < ( profile.length() - 1 ) ) ) {
        profile.createScaledMatchDistributionForPosition(
          ( row_i - 1 ),
          row_match_distribution
        );
      } // End if ( row_i >= 1 ) && ( row_i != last_row )

      return
        backward_calculateRow (
          parameters,
          use_viterbi,
          profile,
          row_match_distribution,
          sequence,
          row_i,
          next_backward_row,
          backward_row
        );
    } // backward_calculateRow(..)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

    /**
     * Generate (randomly) sequences from the given Profile model.  Actually,
     * paths are drawn, and in addition to modifying the given Fasta reference
     * to hold the sequences, the given MultipleAlignment reference is modified
     * to hold the paths.  The Fasta type includes sequence description
     * strings.  These descriptions will be set to be ( description_prefix +
     * seq_i ) for each sequence seq_i (so the first sequence would have name
     * "Randomly generated 0" if the description_prefix == "Randomly generated
     * ").
     *
     * NOTE: For now we will only include sequences of length > 0.  We will
     * keep drawing until we get a sequence of non-zero length, if we have to.
     * The return value is the number of 0-length sequences drawn.
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    uint32_t
    drawSequences (
      Parameters const& parameters,
      ProfileType const & profile,
      uint32_t const sequence_count,
      string const & description_prefix,
      Random & random,
      Fasta<SequenceResidueType> & fasta,
      MultipleAlignment<ProfileType, SequenceResidueType> & multiple_alignment
    ) const;

    /**
     * Generate (randomly) a sequence from the given Profile model.  Both the
     * sequence and the particular path drawn are returned to the caller via
     * modifying the provided references.  The sequence_advances int-vector
     * will have length four greater than that of the profile, and the integer
     * at index i will indicate the number of sequence positions advanced at
     * profile position (i-2) by the path.  So if a profile position's integer
     * is 0, that profile position was a deletion in the path generating the
     * sequence.  If it is 1, then it was a match.  Any number n greater than 1
     * indicates that there were n-1 insertions following a match at that
     * position.  Index 0 indicates the number of pre-align insertions.  The
     * final index indicates the number of post-align insertions.  Index 1
     * indicates deletion-ins, and the second-to-last index indicates
     * deletion-outs.  Together, the sequence and the sequence_advances convey
     * all of the information in the path.  The return value is the length of
     * the generated sequence (sequence.length()).
     */
    template <typename ProfileType,
              typename SequenceResidueType>
    uint32_t
    drawSequence (
      Parameters const & parameters,
      ProfileType const & profile,
      Random & random,
      vector<uint32_t> & sequence_advances, // the path to be created
      Sequence<SequenceResidueType> & sequence // the sequence to be created
    ) const;

    template <typename ProfileTreeType,
              typename SequenceResidueType>
    void
    alignLeaves (
      Parameters const& parameters,
      TreeMultipleAlignment<ProfileTreeType, SequenceResidueType> & multiple_alignment
    ) const;

  }; // End class DynamicProgramming

  //======//// traits classes ////========//

// TODO: REMOVE.  Seeing if moving this makes it work...
//  template <typename DPProbabilityType,
//            typename ScoreType,
//            typename MatrixValueType>
//  template <typename ResidueType, typename ParameterType>
//  struct parameter_collection_traits<typename DynamicProgramming<DPProbabilityType, ScoreType, MatrixValueType>::template AlignmentProfilePositionParameters<ResidueType, ParameterType> >
//  {
//    typedef ParameterType ProbabilityType;
//  }; // AlignmentProfilePositionParameters parameter_collection_traits

// TODO: Put back?  GCC 4.3.3 on linux didn't like this.  I've just taken it out -- will that suffice?
//     template <typename ResidueType,
//               typename DPProbabilityType,
//               typename ScoreType,
//               typename MatrixValueType>
//   struct profile_traits<typename DynamicProgramming<ResidueType, DPProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile > 
//   {
//     typedef DPProbabilityType ProbabilityType;
//   }; // AlignmentProfile profile_traits

  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::DynamicProgramming::Parameters ////
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_INIT
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      Parameters ()
      {
        if( DEFAULT_debug >= DEBUG_All ) {
          cout << "[debug] DynamicProgramming::Parameters::<init>()" << endl;
        } // End if DEBUG_All
        resetToDefaults();
      } // <init>()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_INIT
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      // Copy constructor
      Parameters ( const AnyParameters & copy_from )
      {
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] DynamicProgramming::Parameters::<init>( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParameters const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
      // Copy constructor/operator
  typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
  operator= (
        const AnyParameters & copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] DynamicProgramming::Parameters::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParameters const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        galosh::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] DynamicProgramming::Parameters::copyFromNonVirtual( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtualDontDelegate( copy_from );
      } // copyFromNonVirtual( AnyParameters const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
  copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      )
      {
        useRabinerScaling = copy_from.useRabinerScaling;
        rabinerScaling_useMaximumValue = copy_from.rabinerScaling_useMaximumValue;
        matrixRowScaleFactor = copy_from.matrixRowScaleFactor;
      } // copyFromNonVirtualDontDelegate( AnyParameters const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_COPY
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFrom ( const Parameters & copy_from )
      {
        copyFromNonVirtual( copy_from );
      } // copyFrom( Parameters const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_REINITIALIZE
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      resetToDefaults ()
      {
        galosh::Parameters::resetToDefaults();
        // TODO: Why isn't the compiler finding "debug" in galosh::Parameters?
        //if( debug >= DEBUG_All ) {
        //  cout << "[debug] DynamicProgramming::Parameters::resetToDefaults()" << endl;
        //} // End if DEBUG_All
        useRabinerScaling = DEFAULT_useRabinerScaling;
        rabinerScaling_useMaximumValue = DEFAULT_rabinerScaling_useMaximumValue;
        matrixRowScaleFactor = DEFAULT_matrixRowScaleFactor;
      } // resetToDefaults()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      writeParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        //os << static_cast<galosh::Parameters>( parameters ) << endl;
        galosh::Parameters::writeParameters( os );
        os << endl;

        os << "[DynamicProgramming]" << endl;
        os << "useRabinerScaling = " << useRabinerScaling << endl;
        os << "rabinerScaling_useMaximumValue = " << rabinerScaling_useMaximumValue << endl;
        os << "matrixRowScaleFactor = " << matrixRowScaleFactor << endl;
      } // writeParameters( basic_ostream & ) const

  ////// Class galosh::DynamicProgramming::ParametersModifierTemplate ////
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_INIT
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        if( DEFAULT_debug >= DEBUG_All ) {
          cout << "[debug] DynamicProgramming::ParametersModifierTemplate::<init>()" << endl;
        } // End if DEBUG_All
        isModified_reset();
      } // <init>()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_INIT
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] DynamicProgramming::ParametersModifierTemplate::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParametersModifierTemplate const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
      // Copy constructor/operator
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType> &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
  operator= (
        const AnyParametersModifierTemplate & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] DynamicProgramming::ParametersModifierTemplate::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParametersModifierTemplate const & )
    
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] DynamicProgramming::ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
        } // End if DEBUG_All
        isModified_copyFromNonVirtual( copy_from );

        base_parameters_modifier_t::parameters.copyFromNonVirtual( copy_from.parameters );
      } // copyFromNonVirtual( AnyParametersModifierTemplate const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        galosh::ParametersModifierTemplate<ParametersType>::isModified_copyFromNonVirtual( copy_from );
        isModified_useRabinerScaling = copy_from.isModified_useRabinerScaling;
        isModified_rabinerScaling_useMaximumValue = copy_from.isModified_rabinerScaling_useMaximumValue;
        isModified_matrixRowScaleFactor = copy_from.isModified_matrixRowScaleFactor;
      } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
  reset ()
      {
        isModified_reset();
        base_parameters_modifier_t::parameters.resetToDefaults();
      } // reset()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        galosh::ParametersModifierTemplate<ParametersType>::isModified_reset();
        isModified_useRabinerScaling = false;
        isModified_rabinerScaling_useMaximumValue = false;
        isModified_matrixRowScaleFactor = false;
      } // isModified_reset()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        //galosh::ParametersModifierTemplate<ParametersType>::operator<<( os, parameters_modifier );
        galosh::ParametersModifierTemplate<ParametersType>::writeParametersModifier( os );
        os << endl;

        os << "[DynamicProgramming]" << endl;
        if( isModified_useRabinerScaling ) {
          os << "useRabinerScaling = " << base_parameters_modifier_t::parameters.useRabinerScaling << endl;
        }
        if( isModified_rabinerScaling_useMaximumValue ) {
          os << "rabinerScaling_useMaximumValue = " << base_parameters_modifier_t::parameters.rabinerScaling_useMaximumValue << endl;
        }
        if( isModified_matrixRowScaleFactor ) {
          os << "matrixRowScaleFactor = " << base_parameters_modifier_t::parameters.matrixRowScaleFactor << endl;
        }
      } // writeParametersModifier ( basic_ostream & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template <class AnyParameters>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
        applyModifications ( AnyParameters & target_parameters ) const
      {
        galosh::ParametersModifierTemplate<ParametersType>::applyModifications( target_parameters );

        if( isModified_useRabinerScaling ) {
          target_parameters.useRabinerScaling =
            base_parameters_modifier_t::parameters.useRabinerScaling;
        }
        if( isModified_rabinerScaling_useMaximumValue ) {
          target_parameters.rabinerScaling_useMaximumValue =
            base_parameters_modifier_t::parameters.rabinerScaling_useMaximumValue;
        }
        if( isModified_matrixRowScaleFactor ) {
          target_parameters.matrixRowScaleFactor =
            base_parameters_modifier_t::parameters.matrixRowScaleFactor;
        }
      } // applyModifications( AnyParameters & ) const

  ////// Class galosh::DynamicProgramming::AlignmentProfile ////
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_INIT
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      AlignmentProfile () :
        vector<AlignmentProfilePosition>()
      {
        // Do nothing else.
      } // <init>()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_INIT
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      AlignmentProfile (
        uint32_t length
      ) :
        vector<AlignmentProfilePosition>( length )
      {
        this->reinitialize( length );
      } // <init>( uint32_t )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_REINITIALIZE
    void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    reinitialize (
      uint32_t length
    )
    {
      if( vector<AlignmentProfilePosition>::size() != length ) {
        vector<AlignmentProfilePosition>::resize( length );
      }

      if( length != 0 ) {
        uint32_t last_pos = length - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ).reinitialize();
        }
      }

    } // reinitialize( uint32_t )

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_INIT
    uint32_t
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    length () const
    {
      return vector<AlignmentProfilePosition>::size();
    } // length() const

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.  Also sets each
       * contained AlignmentProfilePosition's m_scalar to 1.
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
      template <typename other_type>
  GALOSH_INLINE_TRIVIAL
      void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      normalize ( other_type const & min )
      {
        //normalize( min );
        normalize( static_cast<MatrixValueType>( min ) );
      } // normalize( other_type const & )

      /**
       * Adjust each distribution's values such that they sum to one, ensuring
       * that no value is less than the specified minimum.  Also sets each
       * contained AlignmentProfilePosition's m_scalar to 1.
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_TRIVIAL
      void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      normalize ( MatrixValueType const& min )
      {
        if( length() == 0 ) {
          return;
        }
        uint32_t last_pos = length() - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ).normalize( min );
        }
      } // normalize( MatrixValueType const& )

      /**
       * Set all values to 0, and the scalars to 1.
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_TRIVIAL
      void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      zero ()
      {
        if( length() == 0 ) {
          return;
        }
        uint32_t last_pos = length() - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ).zero();
        }
      } // zero()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_TRIVIAL
  typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      operator= ( AlignmentProfile const& other_profile )
      {
        if( length() == 0 ) {
          return *this;
        }
        uint32_t last_pos = min( length(), other_profile.length() ) - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ) =
            other_profile[ pos_i ];
        }

        return *this;
      } // operator=( AlignmentProfile const& )

      /**
       * Adjust each distribution's values such that the largest value is 1.0,
       * utilizing the m_scalar values.
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_TRIVIAL
      void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      rescale ()
      {
        rescale( 1.0 );
      } // rescale()

      /**
       * Adjust each distribution's values such that the largest value is the
       * given maximum, utilizing the m_scalar values.
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
      void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      rescale ( MatrixValueType const & max )
      {
        if( length() == 0 ) {
          return;
        }
        uint32_t last_pos = length() - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ).rescale( max );
        }
      } // rescale( MatrixValueType const & )

      /**
       * Divide each value by m_scalar and set the m_scalars to 1.0.
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
      void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      unscale ()
      {
        if( length() == 0 ) {
          return;
        }
        uint32_t last_pos = length() - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ).unscale();
        }
      } // unscale()

    /**
     * Divide each contained distribution value by denominator.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
    template <typename AnyValueType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    operator/= ( AnyValueType const& denominator )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<AlignmentProfilePosition>::operator[]( pos_i ) /= denominator;
      }

      return *this;
    } // operator/=( AnyValueType const& )

    /**
     * Add to each contained distribution the values in the given other
     * AlignmentProfile.  Note that this may violate the rule that the
     * probabilities sum to 1.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    operator+= ( AlignmentProfile const& other_profile )
    {
      /*   //TAH  2/12 compiler doesn't like return with no value on non-void funcs.
      if( length() == 0 ) {
        return;
      }
      */
      assert( length() == other_profile.length() );
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<AlignmentProfilePosition>::operator[]( pos_i ) +=
          other_profile[ pos_i ];
      }

      return *this;
    } // operator+=( AlignmentProfile const& )

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_REINITIALIZE
    void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    uniform ( Random & random )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<AlignmentProfilePosition>::operator[]( pos_i ).uniform( random );
      }
    } // uniform( Random & )

    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    MatrixValueType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    maximumValue () const
    {
      if( length() == 0 ) {
        return;
      }
      ProbabilityType largest_value =
        this->vector<AlignmentProfilePosition>::operator[]( 0 ).maximumValue();
      ProbabilityType largest_value_tmp;
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        largest_value_tmp =
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ).maximumValue();
        if( largest_value_tmp > largest_value ) {
          largest_value = largest_value_tmp;
        }
      }
      return largest_value;
    } // maximumValue() const

      /**
       * Calculate and return the Euclidean distance between this alignment
       * profile and another one (treating every probability as an orthogonal
       * dimension).
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_TRIVIAL
  double
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
  euclideanDistance (
    AlignmentProfile const& other_profile
  ) const
  {
    return sqrt( euclideanDistanceSquared( other_profile ) );
  } // euclideanDistance( AlignmentProfile const & )
  
      /**
       * Calculate and return the square of the Euclidean distance between this
       * alignment profile and another one (treating every probability as an
       * orthogonal dimension).
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
  double
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
  euclideanDistanceSquared (
    AlignmentProfile const& other_profile
  ) const
  {
      double squared_euclidean_distance = 0.0;
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        squared_euclidean_distance +=
          this->vector<AlignmentProfilePosition>::operator[]( pos_i ).euclideanDistanceSquared( other_profile[ pos_i ] );
      }

      return squared_euclidean_distance;
  } // euclideanDistanceSquared( AlignmentProfile const& )

      /**
       * Calculate and return the product of the probabilities that the Match
       * states of this alignment profile are identical to the Match states of
       * the given alignment profile, using the Match state emission values of
       * the internal positions of the two alignment profiles as the relevant
       * probabilities.  Note that these are (usually) not conditioned on there
       * being a Match at all (that is, the probs don't sum to 1; instead they
       * sum to 1-P(Deletion)).  Note also that this is symmetric, so calling
       * the other profile's calculateMatchEmissionCooccurrenceProbability
       * with this profile as an argument will return the same thing.
       * NOTE: We assume the profiles are the same length!
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_ALGORITHM
      MatrixValueType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      calculateMatchEmissionCooccurrenceProbability (
        AlignmentProfile const & other_profile,
        const bool include_nonemission_counts
      ) const
      {
        MatrixValueType prob = 1.0;
        uint32_t last_pos = length() - 1;
        uint32_t pos_i;
        // Note we do *not* include pos 0, since it has no Match emissions.
        for( pos_i = 1; pos_i <= last_pos; pos_i++ ) {
          prob *=
            this->vector<AlignmentProfilePosition>::operator[]( pos_i ).calculateMatchEmissionCooccurrenceProbability( other_profile[ pos_i ], include_nonemission_counts );
        }
        return prob;
      } // calculateMatchEmissionCooccurrenceProbability ( AlignmentProfile const & ) const

      /**
       * Calculate and return the sum of the expected number of times that the
       * Match states of this alignment profile are identical to the Match
       * states of the given alignment profile, using the Match state emission
       * values of the internal positions of the two alignment profiles as the
       * relevant probabilities.  Note that these are (usually) not conditioned
       * on there being a Match at all (that is, the probs don't sum to 1;
       * instead they sum to 1-P(Deletion)).  Note also that this is symmetric,
       * so calling the other profile's
       * calculateMatchEmissionCooccurrenceExpectedCount with this profile as an
       * argument will return the same thing.  NOTE: We assume the profiles are
       * the same length!
       */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_ALGORITHM
      MatrixValueType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
      calculateMatchEmissionCooccurrenceExpectedCount (
        AlignmentProfile const & other_profile,
        const bool include_nonemission_counts
      ) const
      {
        MatrixValueType expected_count = 0.0;
        uint32_t last_pos = length() - 1;
        uint32_t pos_i;
        // Note we do *not* include pos 0, since it has no Match emissions.
        for( pos_i = 1; pos_i <= last_pos; pos_i++ ) {
          expected_count +=
            this->vector<AlignmentProfilePosition>::operator[]( pos_i ).calculateMatchEmissionCooccurrenceProbability( other_profile[ pos_i ], include_nonemission_counts );
        }
        return expected_count;
      } // calculateMatchEmissionCooccurrenceExpectedCount ( AlignmentProfile const & ) const

    /**
     * Calculate and return the cross entropy E(-log(other_node)).  Note that
     * the cross entropy is non-symmetric (calling other_node.crossEntropy(
     * *this ) will return a different value).
     *
     * This can be used to calculate the (self-)entropy by calling
     * this->crossEntropy( *this ).  It can also be used to calculate the KL
     * divergence by taking the difference of the cross-entropy and the
     * self-entropy, or the symmeterized KL divergence by summing the KL
     * divergences computed both ways.
     *
     * See also the other crossEntropy method, which accepts a weights argument.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_TRIVIAL
    double
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    crossEntropy (
      AlignmentProfile const& other_node
    ) const
    {
      return crossEntropy( other_node, ( AlignmentProfile const * )0 );
    } // crossEntropy( AnyInternalNodeOrRoot const& )

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_node)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_node)).  Note that the cross entropy is
     * non-symmetric (calling other_node.crossEntropy( *this, weights ) will
     * return a different value).
     *
     * This can be used to calculate the (possibly weighted) (self-)entropy by
     * calling this->crossEntropy( *this, weights ).  It can also be used to
     * calculate the KL divergence by taking the difference of the
     * cross-entropy and the self-entropy, or the symmeterized KL divergence by
     * summing the KL divergences computed both ways.  Note, though, that
     * weights should all be the same within any Multinomial distribution if
     * this is to be used to calculate weighted entropies or KL divergences.
     * Otherwise the usual properties of these metrics will be violated.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_ALGORITHM
    double
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::AlignmentProfile::
    crossEntropy (
      AlignmentProfile const& other_node,
      AlignmentProfile const * const weights
    ) const
    {
      assert( this->length() == other_node.length() );
      assert( ( weights == 0 ) || ( this->length() == weights->length() ) );

      double cross_entropy = 0.0;
      // And each position.
      if( length() == 0 ) {
        return cross_entropy;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      // TODO: REMOVE
      //cout << "crossEntropy(..):" << endl;
      //double ce;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        cross_entropy +=
          vector<AlignmentProfilePosition>::operator[]( pos_i ).crossEntropy(
            other_node[ pos_i ],
            ( weights ? &( *weights )[ pos_i ] : ( AlignmentProfilePosition const * )0 )
        );
        //ce = 
        //  vector<Position>::operator[]( pos_i ).crossEntropy(
        //    other_node[ pos_i ],
        //    ( weights ? &( *weights )[ pos_i ] : ( ProfilePosition<double> const * )0 )
        //);
        //cout << "[ " << pos_i << " ] " << ( ce - vector<Position>::operator[]( pos_i ).crossEntropy(
        //  ( *this )[ pos_i ],
        //    ( weights ? &( *weights )[ pos_i ] : ( ProfilePosition<double> const * )0 )
        //) )
        //<< endl;
        //cross_entropy += ce;
      }
      return cross_entropy;
    } // crossEntropy( AlignmentProfile const &, AlignmentProfile const * ) const


  ////// Class galosh::DynamicProgramming ////
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
    /**
     * Given the forward viterbi matrices (Calculated as by the
     * forward_score_viterbi(..) method), calculate and return the
     * MultipleAlignment of the sequences with the profile model.
     */
  void
    DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    forward_viterbiAlign (
      Parameters const& parameters,
      typename Matrix::SequentialAccessContainer const & viterbi_matrices,
      typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template MultipleAlignment<ProfileType, SequenceResidueType> & ma
    ) const
    {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
      // TODO: Implement for this case.  Since right now the trainer uses forward_reverseCalculateRow(..), which can't handle this case, I'm being lazy and not imlementing it here either.
      assert( false && "TODO: Implement forward_viterbiAlign(..) for the case with USE_DEL_IN_DEL_OUT but not DISALLOW_FLANKING_TRANSITIONS" );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT && !DISALLOW_FLANKING_TRANSITIONS

      static const bool do_extra_debugging = false;//true;

      ProfileType const& profile =
        *( ma.m_profile );
      vector<Sequence<SequenceResidueType> > const& sequences =
        *( ma.m_sequences );
      uint32_t const& sequence_count =
        ma.m_sequence_count;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> > scaled_match_distributions( profile.length() - 1 );
      for( uint32_t pos_i = 0; pos_i < ( profile.length() - 1 ); pos_i++ ) {
        profile.createScaledMatchDistributionForPosition(
          pos_i,
          scaled_match_distributions[ pos_i ]
        );
      } // End foreach pos_i
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      const bool use_rabiner_scaling = parameters.useRabinerScaling;

      uint32_t last_seq = sequence_count - 1;
      uint32_t seq_i;
      uint32_t seq_adv_i;
      uint32_t last_row = profile.length();
      uint32_t row_i;
      typename Matrix::SequentialAccessContainer::const_reverse_iterator viterbi_matrices_rev_iter;
      uint32_t last_col;
      uint32_t col_i;
      Subcell subcell;
      Subcell highest_previous_subcell;
      ScoreType highest_previous_subcell_value = 0;
      ScoreType tmp_previous_subcell_value = 0;
      uint32_t highest_previous_subcell_row_i = 0;
      uint32_t highest_previous_subcell_col_i = 0;

      for( seq_i = 0; seq_i <= last_seq ; seq_i++ ) {
        // Make sure we start with a zeroed m_sequenceAdvances vector.
        // TODO: Is this necessary?
        ma.m_sequenceAdvances[ seq_i ].assign(
          ma.m_sequenceAdvances[ seq_i ].size(),
          0
        );

        last_col = sequences[ seq_i ].length();

        row_i = last_row;
        seq_adv_i = row_i + 1;
        viterbi_matrices_rev_iter = viterbi_matrices.rbegin();
        col_i = last_col;

        subcell = Match;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        // Then maybe we start in the DeletionOut cell.
        if( 
          ( *viterbi_matrices_rev_iter )[ seq_i ].m_deletionOut >
          ( *viterbi_matrices_rev_iter )[ seq_i ][ last_col ][ Match ]
        ) {
          subcell = DeletionOut;
        }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        while( row_i > 0 ) {

          if( do_extra_debugging ) {
            cout << "viterbiAlign [seq " << seq_i << "]: " <<
              "row " << row_i << ", col " << col_i << ", seq_adv_i " << seq_adv_i << ", subcell " << subcell << endl;
          }

          // TODO: REMOVE
          //if( seq_i == 17 ) {
          //  cout << "row_i = " << row_i << ", col_i = " << col_i <<
          //    ", subcell = " << ( ( subcell == Match ) ? "Match" : ( ( subcell == Insertion ) ? "Insertion" : "Deletion" ) ) << endl;
          //}

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          if( subcell == DeletionOut ) {
            assert( col_i == last_col );
            assert( row_i >= 2 );
            // First count the DeletionOut.
            ma.m_sequenceAdvances[ seq_i ][ last_row + 2 ] += 1;

            // Then the two sources are the DeletionOut above (unless this is
            // row 2, in which case there is none above) and the Match above
            // and to the left.
            highest_previous_subcell = Match;
            highest_previous_subcell_row_i = row_i - 1;
            highest_previous_subcell_col_i = col_i;
            // Note that we use the unscaled value here.
            highest_previous_subcell_value =
              profile[
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletionOut
              ];
            if( row_i == last_row ) {
              highest_previous_subcell_value *=
                profile[
                  Transition::fromDeletionOut
                ][
                  TransitionFromDeletionOut::toEnd
                ];
#ifdef USE_END_DISTRIBUTION
              highest_previous_subcell_value *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
#endif // USE_END_DISTRIBUTION
            } // End if row_i == last_row
            // temporarily increment it
            viterbi_matrices_rev_iter++;
            highest_previous_subcell_value *=
              ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i ][ Match ];
            highest_previous_subcell_value /=
              parameters.matrixRowScaleFactor;
            if( use_rabiner_scaling ) {
              highest_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
            } // End if use_rabiner_scaling
            // put it back
            viterbi_matrices_rev_iter--;
            if( row_i > 2 ) {
              tmp_previous_subcell_value =
                profile[
                  Transition::fromDeletionOut
                ][
                  TransitionFromDeletionOut::toDeletionOut
                ];
              if( row_i == last_row ) {
                tmp_previous_subcell_value *=
                  profile[
                    Transition::fromDeletionOut
                  ][
                    TransitionFromDeletionOut::toEnd
                  ];
#ifdef USE_END_DISTRIBUTION
                tmp_previous_subcell_value *=
                  profile[
                    Transition::fromEnd
                  ][
                    TransitionFromEnd::toPostAlign
                  ];
#endif // USE_END_DISTRIBUTION
              } // End if row_i == last_row
              // temporarily increment it
              viterbi_matrices_rev_iter++;
              tmp_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ].m_deletionOut;
              tmp_previous_subcell_value /=
                parameters.matrixRowScaleFactor;
              if( use_rabiner_scaling ) {
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
              } // End if use_rabiner_scaling
              // put it back
              viterbi_matrices_rev_iter--;
              if( tmp_previous_subcell_value >
                  highest_previous_subcell_value ) {
                highest_previous_subcell_value =
                  tmp_previous_subcell_value;
                highest_previous_subcell_row_i = row_i - 1;
                highest_previous_subcell_col_i = col_i;
                highest_previous_subcell = DeletionOut;
              }
            } // End if row_i > 2
            // End if subcell == DeletionOut
          } else if( subcell == DeletionIn ) {
            assert( col_i == 0 );
            assert( row_i != last_row );
            // First count the DeletionIn.
            ma.m_sequenceAdvances[ seq_i ][ 1 ] += 1;

            // There is always only one source.
            if( row_i == 1 ) {
              highest_previous_subcell = Match;
            } else {
              highest_previous_subcell = DeletionIn;
            }
            highest_previous_subcell_row_i = row_i - 1;
            highest_previous_subcell_col_i = col_i;
            // End if subcell == DeletionIn
          } else
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          if( subcell == Match ) {
            // Then the three sources are from the previous row and cell, unless
            // this is the last row (in which case there's an additional source
            // to the left and another two directly above) or the first row (in
            // which case there's only a source to the left).

            // if USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS, then
            // there is an additional source in the previous row when col_i ==
            // 1.
  
            if( col_i != 0 ) {
              highest_previous_subcell = Match;
              highest_previous_subcell_row_i = row_i - 1;
              highest_previous_subcell_col_i = col_i - 1;
              if( row_i == 1 ) {
                highest_previous_subcell_value =
                  profile[
                    Transition::fromPreAlign
                  ][
                    TransitionFromPreAlign::toBegin
                  ];
                highest_previous_subcell_value *=
                  profile[
                    Transition::fromBegin
                  ][
                    TransitionFromBegin::toMatch
                  ];
              } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                highest_previous_subcell_value =
                  scaled_match_distributions[ row_i - 2 ][
                    TransitionFromMatch::toMatch
                  ];
#else
                highest_previous_subcell_value =
                  profile[ row_i - 2 ][
                    Transition::fromMatch
                  ][
                    TransitionFromMatch::toMatch
                  ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
              }
              highest_previous_subcell_value *=
                getEmissionProbability(
                  profile[ row_i - 1 ],
                  sequences[ seq_i ],
                  col_i - 1
                );
              // temporarily increment it
              viterbi_matrices_rev_iter++;
              highest_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i - 1 ][ Match ];
              highest_previous_subcell_value /=
                parameters.matrixRowScaleFactor;
              if( use_rabiner_scaling ) {
                highest_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
              } // End if use_rabiner_scaling
              // put it back
              viterbi_matrices_rev_iter--;
#ifdef USE_END_DISTRIBUTION
              if( row_i == last_row ) {
                highest_previous_subcell_value *=
                  profile[
                    Transition::fromEnd
                  ][
                    TransitionFromEnd::toPostAlign
                  ];
              }
#endif // USE_END_DISTRIBUTION
              if( row_i > 1 ) {  // No Insertion subcell in row 0
                tmp_previous_subcell_value =
                  profile[ row_i - 2 ][
                    Transition::fromInsertion
                  ][
                    TransitionFromInsertion::toMatch
                  ];
                tmp_previous_subcell_value *=
                  getEmissionProbability(
                    profile[ row_i - 1 ],
                    sequences[ seq_i ],
                    col_i - 1
                  );
                // temporarily increment it
                viterbi_matrices_rev_iter++;
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i - 1 ][ Insertion ];
                tmp_previous_subcell_value /=
                  parameters.matrixRowScaleFactor;
                if( use_rabiner_scaling ) {
                  tmp_previous_subcell_value *=
                    ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
                } // End if use_rabiner_scaling
                // put it back
                viterbi_matrices_rev_iter--;
#ifdef USE_END_DISTRIBUTION
                if( row_i == last_row ) {
                  tmp_previous_subcell_value *=
                    profile[
                      Transition::fromEnd
                    ][
                      TransitionFromEnd::toPostAlign
                    ];
                }
#endif // USE_END_DISTRIBUTION
                if( tmp_previous_subcell_value >
                    highest_previous_subcell_value ) {
                  highest_previous_subcell_value =
                    tmp_previous_subcell_value;
                  highest_previous_subcell = Insertion;
                }
              } // End if row_i > 1
              if( row_i > 1 ) {  // No Deletion subcell in row 0
                tmp_previous_subcell_value =
                  profile[ row_i - 2 ][
                    Transition::fromDeletion
                  ][
                    TransitionFromDeletion::toMatch
                  ];
                tmp_previous_subcell_value *=
                  getEmissionProbability(
                    profile[ row_i - 1 ],
                    sequences[ seq_i ],
                    col_i - 1
                  );
                // temporarily increment it
                viterbi_matrices_rev_iter++;
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i - 1 ][ Deletion ];
                tmp_previous_subcell_value /=
                  parameters.matrixRowScaleFactor;
                if( use_rabiner_scaling ) {
                  tmp_previous_subcell_value *=
                    ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
                } // End if use_rabiner_scaling
                // put it back
                viterbi_matrices_rev_iter--;
#ifdef USE_END_DISTRIBUTION
                if( row_i == last_row ) {
                  tmp_previous_subcell_value *=
                    profile[
                      Transition::fromEnd
                    ][
                       TransitionFromEnd::toPostAlign
                    ];
                }
#endif // USE_END_DISTRIBUTION
                if( tmp_previous_subcell_value >
                    highest_previous_subcell_value ) {
                  highest_previous_subcell_value =
                    tmp_previous_subcell_value;
                  highest_previous_subcell = Deletion;
                }
              } // End if row_i > 1
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              // In the first column there's DeletionIn->Match
              if( ( row_i > 1 ) && ( col_i == 1 ) ) {
                tmp_previous_subcell_value =
                  profile[
                    Transition::fromDeletionIn
                  ][
                    TransitionFromDeletionIn::toMatch
                  ];
                tmp_previous_subcell_value *=
                  getEmissionProbability(
                    profile[ row_i - 1 ],
                    sequences[ seq_i ],
                    col_i - 1
                  );
                // temporarily increment it
                viterbi_matrices_rev_iter++;
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ].m_deletionIn;
                tmp_previous_subcell_value /=
                  parameters.matrixRowScaleFactor;
                if( use_rabiner_scaling ) {
                  tmp_previous_subcell_value *=
                    ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
                } // End if use_rabiner_scaling
                // put it back
                viterbi_matrices_rev_iter--;
                if( tmp_previous_subcell_value >
                    highest_previous_subcell_value ) {
                  highest_previous_subcell_value =
                    tmp_previous_subcell_value;
                  highest_previous_subcell = DeletionIn;
                }
              } // End if ( row_i > 1 ) && ( col_i == 1 )
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
              // In the last row there's postAlignInsertion Match->Match.
              if( row_i == last_row ) {
                tmp_previous_subcell_value =
                  profile[
                    Transition::fromPostAlign
                  ][
                    TransitionFromPostAlign::toPostAlign
                  ];
                tmp_previous_subcell_value *=
                  getPostAlignEmissionProbability(
                    profile,
                    sequences[ seq_i ],
                    col_i - 1
                  );
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i - 1 ][ Match ];
                tmp_previous_subcell_value /=
                  parameters.matrixRowScaleFactor;
                if( use_rabiner_scaling ) {
                  tmp_previous_subcell_value *=
                    ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
                } // End if use_rabiner_scaling
                if( tmp_previous_subcell_value >
                    highest_previous_subcell_value ) {
                  highest_previous_subcell_value = tmp_previous_subcell_value;
                  highest_previous_subcell_row_i = row_i;
                  highest_previous_subcell = Match;
                }
              } // End if row_i == last_row
            } // End if col_i != 0
            // In the last row there's Deletions into the Match state.
            if( row_i == last_row ) {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              tmp_previous_subcell_value =
                scaled_match_distributions[ row_i - 2 ][
                  TransitionFromMatch::toDeletion
                ];
#else
              tmp_previous_subcell_value =
                profile[ row_i - 2 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toDeletion
                ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
#ifdef USE_END_DISTRIBUTION
              // And we transition from the end to the postalign.
              tmp_previous_subcell_value *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
#endif // USE_END_DISTRIBUTION
              // temporarily increment it
              viterbi_matrices_rev_iter++;
              tmp_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i ][ Match ];
              tmp_previous_subcell_value /=
                parameters.matrixRowScaleFactor;
              if( use_rabiner_scaling ) {
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
              } // End if use_rabiner_scaling
              // put it back
              viterbi_matrices_rev_iter--;
              if( tmp_previous_subcell_value >
                  highest_previous_subcell_value ) {
                // TODO: REMOVE
                //if( seq_i == 17 ) {
                //  cout << "delOpen into last row" << endl;
                //}
                // deletionOpen into last row.
                highest_previous_subcell_value =
                  tmp_previous_subcell_value;
                highest_previous_subcell_row_i = row_i - 1;
                highest_previous_subcell_col_i = col_i;
                highest_previous_subcell = Match;
              }
              tmp_previous_subcell_value =
                profile[ row_i - 2 ][
                  Transition::fromDeletion
                ][
                  TransitionFromDeletion::toDeletion
                ];
#ifdef USE_END_DISTRIBUTION
              // And we transition from the end to the postalign.
              tmp_previous_subcell_value *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
#endif // USE_END_DISTRIBUTION
              // temporarily increment it
              viterbi_matrices_rev_iter++;
              tmp_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i ][ Deletion ];
              tmp_previous_subcell_value /=
                parameters.matrixRowScaleFactor;
              if( use_rabiner_scaling ) {
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
              } // End if use_rabiner_scaling
              // put it back
              viterbi_matrices_rev_iter--;
              if( tmp_previous_subcell_value >
                  highest_previous_subcell_value ) {
                // deletionExtesion into last row.
                highest_previous_subcell_value =
                  tmp_previous_subcell_value;
                highest_previous_subcell_row_i = row_i - 1;
                highest_previous_subcell_col_i = col_i;
                highest_previous_subcell = Deletion;
              }
            } // End if row_i == last_row (include deletions too)
            // End if subcell == Match
          } else if( subcell == Insertion ) {
            // Insertions come from the left, either from the Match or the
            // Insertion slot.  There's no use of the Insertion subcell in the
            // last or first rows.  Note that we can assert that col_i != 0,
            // since there's no Insertion into that column.
            highest_previous_subcell = Match;
            highest_previous_subcell_row_i = row_i;
            highest_previous_subcell_col_i = col_i - 1;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            highest_previous_subcell_value =
              scaled_match_distributions[ row_i - 1 ][
                TransitionFromMatch::toInsertion
              ];
#else
            highest_previous_subcell_value =
              profile[ row_i - 1 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toInsertion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            highest_previous_subcell_value *=
              getInsertionEmissionProbability(
                profile[ row_i - 1 ],
                sequences[ seq_i ],
                col_i - 1
              );
            highest_previous_subcell_value *=
              ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i - 1 ][ Match ];
            highest_previous_subcell_value /=
              parameters.matrixRowScaleFactor;
            if( use_rabiner_scaling ) {
              highest_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
            } // End if use_rabiner_scaling
            tmp_previous_subcell_value =
              profile[ row_i - 1 ][
                Transition::fromInsertion
              ][
                TransitionFromInsertion::toInsertion
              ];
            tmp_previous_subcell_value *=
              getInsertionEmissionProbability(
                profile[ row_i - 1 ],
                sequences[ seq_i ],
                col_i - 1
              );
            tmp_previous_subcell_value *=
              ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i - 1 ][ Insertion ];
            tmp_previous_subcell_value /=
              parameters.matrixRowScaleFactor;
            if( use_rabiner_scaling ) {
              tmp_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
            } // End if use_rabiner_scaling
            if( tmp_previous_subcell_value >
                highest_previous_subcell_value ) {
              highest_previous_subcell_value =
                tmp_previous_subcell_value;
              highest_previous_subcell_row_i = row_i;
              highest_previous_subcell_col_i = col_i - 1;
              highest_previous_subcell = Insertion;
            }
            // End if subcell == Insertion
          } else { // subcell == Deletion
            // Deletions come from the left, either from the Match or the
            // Deletion slot.  There's no use of the Deletion subcell in the
            // last or first rows.
            // Assume row_i is not 0.  Also we can assert that row_i won't be
            // last_row when subcell is Deletion, since things are stored at the
            // Match position in the last row.
            highest_previous_subcell = Match;
            highest_previous_subcell_row_i = row_i - 1;
            highest_previous_subcell_col_i = col_i;
            if( row_i == 1 ) { 
              highest_previous_subcell_value =
                profile[
                  Transition::fromPreAlign
                ][
                  TransitionFromPreAlign::toBegin
                ];
              highest_previous_subcell_value *=
                profile[
                  Transition::fromBegin
                ][
                  TransitionFromBegin::toDeletion
                ];
            } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              highest_previous_subcell_value =
                scaled_match_distributions[ row_i - 2 ][
                  TransitionFromMatch::toDeletion
                ];
#else
              highest_previous_subcell_value =
                profile[ row_i - 2 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toDeletion
                ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            }
            // temporarily increment it
            viterbi_matrices_rev_iter++;
            highest_previous_subcell_value *=
              ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i ][ Match ];
            highest_previous_subcell_value /=
              parameters.matrixRowScaleFactor;
            if( use_rabiner_scaling ) {
              highest_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
            } // End if use_rabiner_scaling
            // put it back
            viterbi_matrices_rev_iter--;
            if( row_i > 1 ) {
              tmp_previous_subcell_value =
                profile[ row_i - 2 ][
                  Transition::fromDeletion
                ][
                  TransitionFromDeletion::toDeletion
                ];
              // temporarily increment it
              viterbi_matrices_rev_iter++;
              tmp_previous_subcell_value *=
                ( *viterbi_matrices_rev_iter )[ seq_i ][ col_i ][ Deletion ];
              tmp_previous_subcell_value /=
                parameters.matrixRowScaleFactor;
              if( use_rabiner_scaling ) {
                tmp_previous_subcell_value *=
                  ( *viterbi_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
              } // End if use_rabiner_scaling
              // put it back
              viterbi_matrices_rev_iter--;
              if( tmp_previous_subcell_value >
                  highest_previous_subcell_value ) {
                highest_previous_subcell_value =
                  tmp_previous_subcell_value;
                highest_previous_subcell_row_i = row_i - 1;
                highest_previous_subcell_col_i = col_i;
                highest_previous_subcell = Deletion;
              }
            } // End if row_i > 1
            // End if subcell == Deletion
          } // End switch subcell...
  
          if( col_i == highest_previous_subcell_col_i ) {
            // Then this was a deletion.

            // This was causing a bug:  into the last row you can delete & still insert afterwards.  It is unnecessary anyway, since the vector values start at 0.
            //ma.m_sequenceAdvances[ seq_i ][ seq_adv_i ] = 0;

          } else if(
            ( row_i == last_row ) &&
            ( highest_previous_subcell_row_i == row_i )
          ) {
            // A post-align insertion
            ma.m_sequenceAdvances[ seq_i ][ last_row + 3  ]+= 1;
            col_i = highest_previous_subcell_col_i;
          } else {
            // TODO: REMOVE
            //if( seq_i == 17 ) {
            //  cout << "<<an advance>>" << endl;
            //}
            // An insertion or match.
            ma.m_sequenceAdvances[ seq_i ][ seq_adv_i ] += 1;
            col_i = highest_previous_subcell_col_i;
          }
          // TODO: REMOVE
          //if( seq_i == 17 ) {
          //  cout << "ma.m_sequenceAdvances[ " << row_i << " ] is now " << ma.m_sequenceAdvances[ 17 ][ row_i ] << endl;
          //}
          subcell = highest_previous_subcell;
          if( highest_previous_subcell_row_i == ( row_i - 1 ) ) {
            viterbi_matrices_rev_iter++; // Increment matrices iterator, which is going in reverse through the RowVectors
          } else {
            // TODO: REMOVE
            assert( highest_previous_subcell_row_i == row_i );
          }
          row_i = highest_previous_subcell_row_i;
          seq_adv_i = ( row_i + 1 );
        } // End while row_i > 0
        
        if( do_extra_debugging ) {
          cout << "viterbiAlign [seq " << seq_i << "]: " <<
            "row " << row_i << ", col " << col_i << ", seq_adv_i " << seq_adv_i << ", subcell " << subcell << endl;
        }

        // Now if we got to row 0 before column 0, there's pre-align insertions.
        ma.m_sequenceAdvances[ seq_i ][ 0 ] = col_i;

        // TODO: REMOVE
        //if( seq_i == 17 ) {
        //  cout << "ma.m_sequenceAdvances[ 17 ] is " << endl;
        //  for( uint32_t row_i = 0; row_i < ma.m_sequenceAdvances[ 17 ].size(); row_i++ ) {
        //    cout << "ma.m_sequenceAdvances[ " << row_i << " ] is now " << ma.m_sequenceAdvances[ 17 ][ row_i ] << endl;
        //  }
        //}
      } // End for each seq_i

      return;
    } // forward_viterbiAlign ( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileTypeA, typename ProfileTypeB>
  GALOSH_INLINE_ALGORITHM
    /**
     * Calculate the minimum-cost alignment between the two profiles, using the
     * symmeterized KL divergence as the distance metric, and the given indel
     * costs.  Fill the given uint32_t vector (which will be resized to one
     * greater than the length of profile a) with the number of advances of
     * profile b for each position of profile a (index i corresponds to
     * position i-1 of profile a, with index 0 indicating the number of
     * positions of profile b that align before all positions of profile a).
     * After calling this method, that vector's values will sum to the length
     * of profile b.  The return value is the cost of the computed alignment.
     *
     * Note that we do not allow I->D or D->I transitions (that is, a match
     * must be a match, not an insertion by both profiles simultaneously), but
     * if you set the maximum_match_cost argument to something non-zero, all
     * match costs will be truncated to that value -- so you can effectively
     * allow gap matches by setting maximum_match_cost to ( indel_open_cost * 2
     * ).
     */
  double
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
  profileProfile_align_SKL (
    Parameters const& parameters,
    ProfileTypeA const& profile_a,
    ProfileTypeB const& profile_b,
    double const indel_open_cost,
    double const indel_extension_cost,
    vector<uint32_t> & profile_b_advances,
    double const maximum_match_cost // = 0 // by default, no maximum
  ) const
  {
    static const bool be_extra_verbose = false;//true;

    profile_b_advances.clear();
    profile_b_advances.resize( profile_a.length() + 1 );

    // We treat profile b as if it were the 'sequence' (in an analogy to the
    // sequence-profile algorithms), so column j corresponds to position j-1 of
    // profile b (column 0 is before profile b begins), and row i corresponds
    // to position i-1 of profile a (row 0 is before profile a begins).
    uint32_t last_row = profile_a.length();
    uint32_t last_col = profile_b.length();
    register uint32_t row_i, col_i;

    // We make a matrix of BacktraceablePPMatrixCell<double>s to store the costs
    // (we need to keep the whole matrix so we can backtrace).  We store all
    // subcells (match, insertion, deletion) even in the non-affine case,
    // because we need it for backtracing (and because we're not allowing I->D
    // or D->I transitions).
    vector<vector<BacktraceablePPMatrixCell<double> > > matrix( last_row + 1 );

    double cost, skl_distance;
    double self_entropy_a;
    vector<double> self_entropies_b( profile_b.length() );
    for( row_i = 0; row_i <= last_row; row_i++ ) {
      // We have to initialize the matrix row first.
      matrix[ row_i ].resize( last_col + 1 );
      if( row_i > 0 ) {
        // We'll be using this many times, so we compute it here
        self_entropy_a =
          profile_a[ row_i - 1 ].crossEntropy(
            profile_a[ row_i - 1 ]
          );
        // TODO: REMOVE
        assert( self_entropy_a != 0 );
      } else {
        // Shouldn't be used
        self_entropy_a = numeric_limits<double>::quiet_NaN();
      } // row_i > 0

      for( col_i = 0; col_i <= last_col; col_i++ ) {
        if( ( row_i == 0 ) && ( col_i == 0 ) ) {
          matrix[ row_i ][ col_i ][ Match ] = 0;
          matrix[ row_i ][ col_i ][ Insertion ] = numeric_limits<double>::max();
          matrix[ row_i ][ col_i ][ Deletion ] = numeric_limits<double>::max();
          matrix[ row_i ][ col_i ].m_matchFromMatch = true;
          if( be_extra_verbose ) {
            cout << matrix[ row_i ][ col_i ];
          }
          continue;
        }
        if( col_i == 0 ) {
          matrix[ row_i ][ col_i ][ Match ] = numeric_limits<double>::max();
          matrix[ row_i ][ col_i ][ Insertion ] = numeric_limits<double>::max();
        } else {
          // Add the indel costs to the cell from the left
          cost =
            (
              matrix[ row_i ][ col_i - 1 ][ Match ] +
              indel_open_cost
            );
          matrix[ row_i ][ col_i ][ Insertion ] = cost;
          matrix[ row_i ][ col_i ].m_insertionFromMatch = true;
          cost =
            (
              matrix[ row_i ][ col_i - 1 ][ Insertion ] +
              indel_extension_cost
            );
          if( cost <= matrix[ row_i ][ col_i ][ Insertion ] ) {
            matrix[ row_i ][ col_i ].m_insertionFromInsertion = true;
            if( cost < matrix[ row_i ][ col_i ][ Insertion ] ) {
              matrix[ row_i ][ col_i ].m_insertionFromMatch = false;
            }
            matrix[ row_i ][ col_i ][ Insertion ] = cost;
          }
        } // End if col_i > 0
        if( row_i == 0 ) {
          matrix[ row_i ][ col_i ][ Match ] = numeric_limits<double>::max();
          matrix[ row_i ][ col_i ][ Deletion ] = numeric_limits<double>::max();
        } else {
          // Add the indel costs to the cell from above
          cost =
            (
              matrix[ row_i - 1 ][ col_i ][ Match ] +
              indel_open_cost
            );
          matrix[ row_i ][ col_i ][ Deletion ] = cost;
          matrix[ row_i ][ col_i ].m_deletionFromMatch = true;
          cost =
            (
              matrix[ row_i - 1 ][ col_i ][ Deletion ] +
              indel_extension_cost
            );
          if( cost <= matrix[ row_i ][ col_i ][ Deletion ] ) {
            matrix[ row_i ][ col_i ].m_deletionFromDeletion = true;
            if( cost < matrix[ row_i ][ col_i ][ Deletion ] ) {
              matrix[ row_i ][ col_i ].m_deletionFromMatch = false;
            }
            matrix[ row_i ][ col_i ][ Deletion ] = cost;
          }
          if( col_i > 0 ) {
            // Add the symmeterized KL divergence (which is just the sum of the
            // two KL divergences) to the cell from the upper-left.
            if( row_i == 1 ) {
              // Calc and memoize the self-entropy.
              self_entropies_b[ col_i - 1 ] =
                profile_b[ col_i - 1 ].crossEntropy(
                  profile_b[ col_i - 1 ]
                );
            } // End if row_i == 1
            // TODO: REMOVE
            assert( self_entropies_b[ col_i - 1 ] != 0 );

            cost =
              matrix[ row_i - 1 ][ col_i - 1 ][ Match ];
            matrix[ row_i ][ col_i ].m_matchFromMatch = true;
            if( matrix[ row_i - 1 ][ col_i - 1 ][ Insertion ] <= cost ) {
              matrix[ row_i ][ col_i ].m_matchFromInsertion = true;
              if( matrix[ row_i - 1 ][ col_i - 1 ][ Insertion ] < cost ) {
                matrix[ row_i ][ col_i ].m_matchFromMatch = false;
              }
              cost = matrix[ row_i - 1 ][ col_i - 1 ][ Insertion ];
            }
            if( matrix[ row_i - 1 ][ col_i - 1 ][ Deletion ] <= cost ) {
              matrix[ row_i ][ col_i ].m_matchFromDeletion = true;
              if( matrix[ row_i - 1 ][ col_i - 1 ][ Deletion ] < cost ) {
                matrix[ row_i ][ col_i ].m_matchFromMatch = false;
                matrix[ row_i ][ col_i ].m_matchFromInsertion = false;
              }
              cost = matrix[ row_i - 1 ][ col_i - 1 ][ Deletion ];
            }

            // No matter where it's from, the cost involves the symmetric KL
            // distance.
            skl_distance =
              (
                (
                  profile_a[ row_i - 1 ].crossEntropy(
                    profile_b[ col_i - 1 ]
                  ) -
                  self_entropy_a
                ) +
                (
                  profile_b[ col_i - 1 ].crossEntropy(
                    profile_a[ row_i - 1 ]
                  ) -
                  self_entropies_b[ col_i - 1]
                )
              );
            // TODO: REMOVE
            //cout << "skl_distance between " << profile_a[ row_i - 1 ] << " and " << profile_b[ col_i - 1 ] << " is " << skl_distance << endl;
            if(
              ( maximum_match_cost != 0 ) &&
              skl_distance > maximum_match_cost
            ) {
              cost += maximum_match_cost;
            } else {
              cost += skl_distance;
            }
            matrix[ row_i ][ col_i ][ Match ] = cost;
          } // End if col_i > 0
        } // End if row_i > 0
        // That's it.  Nothing else. la da da.
        if( be_extra_verbose ) {
          cout << matrix[ row_i ][ col_i ];
        }
      } // End foreach col_i
      if( be_extra_verbose ) {
        cout << endl;
      }
    } // End foreach row_i

    // Ok now for the backtrace
    row_i = last_row;
    col_i = last_col;
    Subcell lowest_final_subcell;
    // Start in the minimum cost state in the last cell.
    // There could be multiple minimum cost paths.  Prefer the Match, then
    // the Insertion, for lack of any better idea.
    if(
      matrix[ last_row ][ last_col ][ Match ] <=
      matrix[ last_row ][ last_col ][ Insertion ]
    ) {
      if(
        matrix[ last_row ][ last_col ][ Match ] <=
        matrix[ last_row ][ last_col ][ Deletion ]
      ) {
        lowest_final_subcell = Match;
      } else { // Deletion < Match <= Insertion
        lowest_final_subcell = Deletion;
      }
    } else {
      if(
        matrix[ last_row ][ last_col ][ Insertion ] <=
        matrix[ last_row ][ last_col ][ Deletion ]
      ) {
        lowest_final_subcell = Insertion;
      } else { // Deletion < Insertion < Match
        lowest_final_subcell = Deletion;
      }
    }
    Subcell subcell = lowest_final_subcell;
    double last_cost = 0; // For be_extra_verbose
    while( !( ( row_i == 0 ) && ( col_i == 0 ) ) ) {
      if( be_extra_verbose ) {
        cout << ( last_cost - matrix[ row_i ][ col_i ][ subcell ] ) << endl;
        last_cost = matrix[ row_i ][ col_i ][ subcell ];
      }
      if( be_extra_verbose ) {
        cout << "[ " << row_i << ", " << col_i << ", " << subcell << " ] ";// << endl;
      }
      // There could be multiple minimum cost paths.  Prefer the Match, then
      // the Insertion, for lack of any better idea.
      if( subcell == Match ) {
        profile_b_advances[ row_i ] += 1;
        if( matrix[ row_i ][ col_i ].m_matchFromMatch ) {
          // still a Match
        } else if( matrix[ row_i ][ col_i ].m_matchFromInsertion ) {
          subcell = Insertion;
        } else if( matrix[ row_i ][ col_i ].m_matchFromDeletion ) {
          subcell = Deletion;
        } else {
          // TODO: REMOVE?
          cout << "uh-oh: matrix[ " << row_i << " ][ " << col_i << " ].m_matchFrom* false for * in (Match,Insertion,Deletion), yet subcell is Match!" << endl;
          // Uh-oh!
          assert( false );
        }
        row_i -= 1;
        col_i -= 1;
      } else if( subcell == Insertion ) {
        profile_b_advances[ row_i ] += 1;
        if( matrix[ row_i ][ col_i ].m_insertionFromMatch ) {
          subcell = Match;
        } else if( matrix[ row_i ][ col_i ].m_insertionFromInsertion ) {
          // Still an insertion
        //} else if( matrix[ row_i ][ col_i ].m_insertionFromDeletion ) {
        //  // Hrm? We don't allow I->D or D->I transitions..
        //  assert( false );
        } else {
          // Uh-oh!
          assert( false );
        }
        col_i -= 1;
      } else { // subcell == Deletion
        assert( profile_b_advances[ row_i ] == 0 );
        if( matrix[ row_i ][ col_i ].m_deletionFromMatch ) {
          subcell = Match;
        } else if( matrix[ row_i ][ col_i ].m_deletionFromInsertion ) {
          // Hrm? We don't allow I->D or D->I transitions..
          assert( false );
        } else if( matrix[ row_i ][ col_i ].m_deletionFromDeletion ) {
          // Still a deletion
        } else {
          // Uh-oh!
          assert( false );
        }
        row_i -= 1;
      } // End switch current subcell
    } // End backtracing

    return matrix[ last_row ][ last_col ][ lowest_final_subcell ];
  } // profileProfile_align_SKL ( Parameters const &, ProfileTypeA const &, ProfileTypeB const &, double, double, vector<uint32_t> & ) const

  /**
   * Calculate and return the probability that if you drew two paths, one from
   * each of the given profiles, they would coincide exactly.  Note that if the
   * profiles are of different lengths, the probability is 0.
   */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileTypeA, typename ProfileTypeB>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
  calculatePathCooccurrenceProbability (
    Parameters const& parameters,
    ProfileTypeA const& profile_a,
    ProfileTypeB const& profile_b
  ) const
  {
    // I use "row" here to be consistent with the other algorithms, though
    // since we are requiring that both profiles be in the same state at the
    // same time, the row and column is always the same.  row_i==0 corresponds
    // to both profiles being in the pre-align state, and row_i==last_row
    // corresonds to both pofiles being in the post-align state.
    uint32_t last_row = profile_a.length();
    if( last_row != profile_b.length() ) {
      // The paths cannot coincide exactly if the profiles have different
      // lengths.
      return 0;
    }
    uint32_t row_i = 0;

    // The probs of subpaths ending in each subcell at row row_i.
    ScoreType match = 1;
    ScoreType insertion = 1;
    ScoreType deletion = 1;

    ScoreType tmp1, tmp2;

    // Begin in the pre-align state Theoretically there could be an infinite
    // number of pre-aligns.  Each has to occur in both of the profiles
    // simultaneously.  The number of cooccurring pre-aligns follows a
    // Geometric distribution with a probability of continuing ("failing") of
    // the product of the two pre-align->pre-align transition probs times the
    // emission cooccurrence probs of the pre-align emission dists.  The
    // normalizing constant of a Geometric is 1, so this part actually
    // contributes nothing.  But when it "succeeds" to break out, it could be
    // that it was a mutual transition to the Begin state, or that one
    // transitioned out but the other did not, or that they were both
    // pre-aligns but emitted different bases.  We count only the former of
    // these possibilities.  Thus the probability contribution from the
    // pre-align state is the product of the pre-align->begin probs, normalized
    // by 1-(the product of the pre-align->pre-align probs times the
    // cooccurrence prob).
    match =
      profile_a[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ];
    match *=
      profile_b[ Transition::fromPreAlign ][ TransitionFromPreAlign::toBegin ];
    // the normalizing constant is 1-tmp1.
    tmp1 =
      galosh::calculateEmissionCooccurrenceProbability<ScoreType>(
        profile_a,
        profile_b,
        Emission::PreAlignInsertion
      );
    tmp1 *=
      profile_a[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
    tmp1 *=
      profile_b[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ];
    tmp2 = 1.0;
    tmp2 -= tmp1;
    match /= tmp2;
    // That part is shared by deletion.
    deletion = match;

    // But now we differentiate b/n match and deletion
    deletion *=
      profile_a[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
    deletion *=
      profile_b[ Transition::fromBegin ][ TransitionFromBegin::toDeletion ];
    match *=
      profile_a[ Transition::fromBegin ][ TransitionFromBegin::toMatch ];
    match *=
      profile_b[ Transition::fromBegin ][ TransitionFromBegin::toMatch ];
    // We'll handle the match emissions later.

    // And now the row is 1.
    row_i = 1;

    // Now handle all other states..
    for( ; row_i <= last_row; row_i++ ) {
      // Start by calculating the emission from the Match state.
      match *=
        galosh::calculateEmissionCooccurrenceProbability<ScoreType>(
          profile_a[ row_i - 1 ],
          profile_b[ row_i - 1 ],
          Emission::Match
        );
      if( row_i == last_row ) {
        // In the final row there are post-align insertions (only).
        // ERE I AM.  Also todo: transitions to the next state

        // we'll store it in the insertion state here, though it doesn't really
        // matter.
        insertion = match;
        // post-align insertions can follow deletions, too.
        insertion += deletion;

#ifdef USE_END_DISTRIBUTION
        insertion *=
          profile_a[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ];
        insertion *=
          profile_b[ Transition::fromEnd ][ TransitionFromEnd::toPostAlign ];
#endif // USE_END_DISTRIBUTION

        // It doesn't really matter how long the post-align insertion is, so
        // long as everything inserted is the same.  See the discussion above
        // about pre-align insertions.
        insertion *=
          profile_a[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ];
        insertion *=
          profile_b[ Transition::fromPostAlign ][ TransitionFromPostAlign::toTerminal ];
        // the normalizing constant is 1-tmp1.
        tmp1 =
          galosh::calculateEmissionCooccurrenceProbability<ScoreType>(
            profile_a,
            profile_b,
            Emission::PostAlignInsertion
          );
        tmp1 *=
          profile_a[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
        tmp1 *=
          profile_b[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ];
        tmp2 = 1.0;
        tmp2 -= tmp1;
        insertion /= tmp2;

        // And we're done.
        return insertion;
      } else { // if( row_i == last_row ) .. else 
        // Insertion Opens follow matches.
        insertion = match;
        insertion *=
          profile_a[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
        insertion *=
          profile_b[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ];
        insertion *=
          galosh::calculateEmissionCooccurrenceProbability<ScoreType>(
            profile_a[ row_i - 1 ],
            profile_b[ row_i - 1 ],
            Emission::Insertion
          );

        // Again, there could be any number of insertion extensions.  It
        // doesn't really matter, so long as they continue emitting the same
        // thing.  So we need to incorporate the prob that what happened was a
        // joint transition out of the insertion state, conditional on there
        // having been a transition out (or a continuation but a failure to
        // emit the same thing).  As before, the contribution from the
        // insertion extensions themselves ends up being 1, when summed over
        // all possible insertion lengths.
        insertion *=
          profile_a[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
        insertion *=
          profile_b[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
        tmp1 =
          galosh::calculateEmissionCooccurrenceProbability<ScoreType>(
            profile_a[ row_i - 1 ],
            profile_b[ row_i - 1 ],
            Emission::Insertion
          );
        tmp1 *=
          profile_a[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
        tmp1 *=
          profile_b[ Transition::fromInsertion ][ TransitionFromInsertion::toInsertion ];
        tmp2 = 1.0;
        tmp2 -= tmp1;
        insertion /= tmp2;
      } // End if( row_i == last_row ) .. else ..

      // Ok, so now we deal with transitions from this state to the next one.
      tmp2 = match;

      match *= 
        profile_a[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
      match *= 
        profile_b[ Transition::fromMatch ][ TransitionFromMatch::toMatch ];
      tmp1 = deletion;
      tmp1 *=
        profile_a[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ];
      tmp1 *=
        profile_b[ Transition::fromDeletion ][ TransitionFromDeletion::toMatch ];
      match += tmp1;
      tmp1 = insertion;
      tmp1 *=
        profile_a[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
      tmp1 *=
        profile_b[ Transition::fromInsertion ][ TransitionFromInsertion::toMatch ];
      match += tmp1;

      deletion *=
        profile_a[ Transition::fromDeletion ][ TransitionFromDeletion::toDeletion ];
      // tmp2 is the former match value
      tmp2 *=
        profile_a[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ];
      tmp2 *=
        profile_b[ Transition::fromMatch ][ TransitionFromMatch::toDeletion ];
      deletion += tmp2;

      // We will start out the next iteration by setting insertion to match...
      //insertion = 0;
    } // End foreach row_i

    // We should never reach here, since we return at the post-align state.
    assert( false );
    return 0;
  } // calculatePathCooccurrenceProbability( Parameters const&, ProfileTypeA const&, ProfileTypeB const& ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename GlobalCountsType,
            typename PositionCountsType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Given the forward row for row_i and the current SamplerState (with the
     * current column and substate of the sampler for the given sequence,
     * representing the beginning of the partial path from model positions
     * (row_i+1) through (profile.length() - 1), which has already been drawn),
     * draw more of the partial path, corresponding to position row_i of
     * the given profile.  The given target SamplerState reference will be
     * modified to hold the sampler state after drawing the partial path (which
     * is sufficient information, along with the 'current' SamplerState, to
     * reconstruct that partial path).  Note that this method is called in
     * reverse order of the model states, so actually the "target" state will
     * be updated to reflect the first state of the partial path being drawn,
     * and the "current" state is the last state of that partial path.
     *
     * If the GlobalCountsType and/or PositionCountsType argument
     * pointers are non-null, then they should be a Profile type and a
     * ProfilePosition type, storing unsigned integer counts.  These
     * counts will be updated with the number of the observations on the drawn
     * partial path corresponding to each parameter.
     *
     * Note that the added (internal) insertion counts will correspond to the
     * insertions following the Match into profile position (row_i-1), whereas
     * the Match emission counts corresond to the Match into position row_i
     * (so correspond to that profile position, ie use its emission
     * probabilities).
     *
     * If the final argument (the sequence_advances) pointer is non-NULL, then
     * it should be a vector of length ( profile.length() + 4 ).  It will
     * be modified to reflect the drawn partial path.  Note that Matches into
     * row row_i won't be counted until this method is called for the previous
     * row (row_i - 1).
     */
    drawPartialPathForPosition (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& forward_row,
      SamplerState const& current_sampler_state,
      SamplerState & target_sampler_state,
      Random & random,
      GlobalCountsType * global_counts,
      PositionCountsType * position_counts,
      vector<uint32_t> * sequence_advances
    ) const
    {
      const uint32_t last_row = profile.length();
      const uint32_t last_col = sequence.length();

      double unif;
      MatrixValueType tmp_mvt_insertion_to_match, tmp_mvt_deletion_to_match,
        tmp_mvt_match_to_match, tmp_mvt_gap_to_match,
        tmp_mvt_deletion_to_deletion, tmp_mvt_match_to_deletion,
        tmp_mvt_to_deletion, tmp_mvt;

      if( false && ( parameters.debug >= DEBUG_All ) ) {
        cout << "[debug] partialPath: last_row is " << last_row << ", last_col is " << last_col << endl;
      }

      if( row_i == last_row ) { // (No update to position counts in the last row)
        // Note that in the last row, the value in the Match state is the sum
        // of the contribution from the previous Match state in the same row (a
        // postalign insertion) and stuff from the previous row.  We can draw
        // insertions proportional to the previous Match state's contribution
        // to the current Match state's value.
        target_sampler_state.column = last_col;
        ProbabilityType probability_of_insertion;
        bool done = false;
        while( ( target_sampler_state.column > 0 ) && !done ) {
          probability_of_insertion =
            forward_row[ target_sampler_state.column - 1 ][ Match ];
          probability_of_insertion *=
            profile[
              Transition::fromPostAlign
            ][
              TransitionFromPostAlign::toPostAlign
            ];
          probability_of_insertion *=
            getPostAlignEmissionProbability(
              profile,
              sequence,
              ( target_sampler_state.column - 1 )
            );
          probability_of_insertion /=
            forward_row[ target_sampler_state.column ][ Match ];

          if( parameters.debug >= DEBUG_All ) {
            cout << "[debug] partialPath: probability_of_insertion (postalign) is " << probability_of_insertion << endl;
            //cout << "  column is " << target_sampler_state.column << endl;
            //cout << "  forward row is " << forward_row << endl;
            //cout << "  PostAlign->PostAlign is " <<
            //  profile[
            //    Transition::fromPostAlign
            //  ][
            //    TransitionFromPostAlign::toPostAlign
            //  ]
            //     << endl;
            //cout << "   emission probability is " << getPostAlignEmissionProbability(
            //     profile,
            //     sequence,
            //     ( target_sampler_state.column - 1 )
            //   ) << endl;
          }
          if( probability_of_insertion == 0.0 ) {
            unif = 1.0;
          } else if( probability_of_insertion == 1.0 ) {
            unif = 0.0;
          } else {
            unif = random.nextUniform();
            if( parameters.debug >= DEBUG_All ) {
              cout << "[debug] partialPath: unif (postalign) is " << unif << endl;
            }
          }
          if( unif < probability_of_insertion ) {
            if( global_counts != NULL ) {
              ( *global_counts )[
                 Transition::fromPostAlign
               ][
                 TransitionFromPostAlign::toPostAlign
               ] += 1;
              incrementPostAlignEmissionCounts(
                *global_counts,
                sequence,
                target_sampler_state.column - 1
              );
            } // End if global_counts != NULL
            target_sampler_state.column -= 1;
            if( sequence_advances != NULL ) {
              // Post-aligns go in the last sequence_advances index.
              ( *sequence_advances )[ row_i + 3 ] += 1;
            } // End if sequence_advances != NULL
          } else {
            done = true;
          }
        } // End while continuing to do post-align insertions...
        target_sampler_state.subcell = Match; // Always.
        if( global_counts != NULL ) {
          // Also count the observed move to the Terminal state.
          ( *global_counts )[
            Transition::fromPostAlign
          ][
            TransitionFromPostAlign::toTerminal
          ] += 1;
        } // End if global_counts != NULL
        return;
      } // End if this is the last row.

      // If we are in column 0 then we must have deleted the current row.
      if( current_sampler_state.column == 0 ) {
        target_sampler_state.column = 0;
        if( row_i == 0 ) {
          target_sampler_state.subcell = Match;
          if( global_counts != NULL ) {
            ( *global_counts )[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ] += 1;
            ( *global_counts )[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletion
            ] += 1;
          } // End if global_counts != NULL
        } else { // if row_i == 0 .. else ..
          target_sampler_state.subcell = Deletion;
          if( position_counts != NULL ) {
            ( *position_counts )[
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ] += 1;
          } // End if position_counts != NULL
        }
        // Note that we don't increment the sequence_advances, because it
        // should be 0.
        return;
      } // End if current_sampler_state.column == 0

      if( row_i == 0 ) {
        // In the first row, everything is stored in the Match subcell.  In all
        // but column 0 (already handled above), these are actually
        // PreAlignInsertions.  They can be followed by Deletions.
        if( current_sampler_state.subcell == Deletion ) {
          target_sampler_state.subcell = Match; // We store PreAlignInsertions as Matches
          target_sampler_state.column = current_sampler_state.column;
          if( global_counts != NULL ) {
            ( *global_counts )[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ] += 1;
            ( *global_counts )[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletion
            ] += 1;
          } // End if global_counts != NULL

        } else { // current_sampler_state == Match
          target_sampler_state.subcell = Match; // We store PreAlignInsertions as Matches
          target_sampler_state.column = current_sampler_state.column - 1;
          if( global_counts != NULL ) {
            ( *global_counts )[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ] += 1;
            ( *global_counts )[
              Transition::fromBegin
            ][
              TransitionFromBegin::toMatch
            ] += 1;
          } // End if global_counts != NULL

          if( position_counts != NULL ) {
            incrementEmissionCounts(
              *position_counts,
              sequence,
              current_sampler_state.column - 1
            );
          } // End if position_counts != NULL
          if( sequence_advances != NULL ) {
            ( *sequence_advances )[ 2 ] += 1;
          } // End if sequence_advances != NULL

        } // End if current_sampler_state.subcell == Deletion .. else Match ..

        // All the rest in row 0 are PreAlignInsertions.
        while( target_sampler_state.column > 0 ) {
          if( global_counts != NULL ) {
            ( *global_counts )[
               Transition::fromPreAlign
             ][
               TransitionFromPreAlign::toPreAlign
             ] += 1;
            incrementPreAlignEmissionCounts(
              *global_counts,
              sequence,
              target_sampler_state.column - 1
            );
          } // End if global_counts != NULL
          // Increment the sequence advances for the current row, too.
          if( sequence_advances != NULL ) {
            ( *sequence_advances )[ 0 ] += 1;
          } // End if sequence_advances != NULL
        
          // We stay in the Match state, but decrement the column...
          target_sampler_state.column -= 1;
        } // End while( target_sampler_state.column > 0 )

        return;
      } // End if row_i == 0

      if( current_sampler_state.subcell == Deletion ) {
        // Note that we handled the row_i == 0 case already.

        // Then we either continued a deletion or started one.
        tmp_mvt_deletion_to_deletion =
          forward_row[ current_sampler_state.column ][ Deletion ];
        tmp_mvt_deletion_to_deletion *=
          profile[
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toDeletion
          ];
        tmp_mvt_match_to_deletion =
          forward_row[ current_sampler_state.column ][ Match ];
        tmp_mvt_match_to_deletion *=
          profile[
            Transition::fromMatch
          ][
            TransitionFromMatch::toDeletion
          ];

        // Note that we never factored into the tmp_mvt_*_to_deletion the
        // matrixRowScaleFactor or the rabinerInverseScalar, but it doesn't
        // matter since the probability we calculate is a ratio and these
        // factors would affect the numerator and denominator equally.
        tmp_mvt_to_deletion =
          ( tmp_mvt_deletion_to_deletion + tmp_mvt_match_to_deletion );

        ProbabilityType probability_of_deletion =
          tmp_mvt_deletion_to_deletion;
        probability_of_deletion /=
          tmp_mvt_to_deletion;

        if( parameters.debug >= DEBUG_All ) {
          cout << "[debug] partialPath: probability_of_deletion (D->D) is " << probability_of_deletion << endl;
        }
        if( probability_of_deletion == 0 ) {
          unif = 1.0;
        } else if( probability_of_deletion == 1.0 ) {
          unif = 0.0;
        } else {
          unif = random.nextUniform();
          if( parameters.debug >= DEBUG_All ) {
            cout << "[debug] partialPath: unif (D->D) is " << unif << endl;
          }
        }
        if( unif < probability_of_deletion ) {
          target_sampler_state.subcell = Deletion;
          if( position_counts != NULL ) {
            ( *position_counts )[
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ] += 1;
          } // End if position_counts != NULL
        } else { // If ( unif < probability_of_deletion ) .. else ..
          target_sampler_state.subcell = Match;
          if( position_counts != NULL ) {
            ( *position_counts )[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ] += 1;
          } // End if position_counts != NULL
        }
        target_sampler_state.column = current_sampler_state.column;
        // Note that we don't increment the sequence_advances, because it
        // should be 0.
        return;
      } // End if we're in the Deletion state currently

      // If we are using the Plan 7 architecture, the current substate must be
      // either a Deletion or a Match, since Insertions are taken care of on
      // the previous iteration.  Thus by a process of elimination we now know
      // that the subcell is Match.
      assert( current_sampler_state.subcell == Match );

      // These are probabilities of the Match transition (*to* the Match
      // subcell).  Note that we factor out and cancel the Match emission,
      // which appears in all three.
      tmp_mvt_insertion_to_match =
        forward_row[ current_sampler_state.column - 1 ][ Insertion ];
      tmp_mvt_insertion_to_match *= 
        profile[
          Transition::fromInsertion
        ][
          TransitionFromInsertion::toMatch
        ];
      //tmp_mvt_insertion_to_match =
      //  (
      //    forward_row[ current_sampler_state.column - 1 ][ Insertion ] *
      //    (
      //      //profile[ row_i - 1 ][
      //      profile[
      //        Transition::fromInsertion
      //      ][
      //        TransitionFromInsertion::toMatch
      //      ]
      //      //*
      //      //getEmissionProbability(
      //      // profile[ row_i ],
      //      // sequence,
      //      // current_sampler_state.column - 1
      //      //)
      //    )
      //   );
      if( parameters.debug >= DEBUG_All ) {
        if( row_i == ( last_row - 1 ) ) {
          // If we are also needing to test that this isn't a transition to a
          // deletion (which is stored in the Match state in the last row),
          // then we will need to incorporate the emission probability too.
          cout << "[debug] partialPath: tmp_mvt_insertion_to_match (I->M) is " <<
            ( tmp_mvt_insertion_to_match *
              (
                getEmissionProbability(
                  profile[ row_i ],
                  sequence,
                  current_sampler_state.column - 1
                )
              )
            ) << endl;
        } else {
          cout << "[debug] partialPath: tmp_mvt_insertion_to_match (I->M) is " << tmp_mvt_insertion_to_match << endl;
        }
      } // End if DEBUG_ALL
      tmp_mvt_match_to_match =
        forward_row[ current_sampler_state.column - 1 ][ Match ];
      tmp_mvt_match_to_match *=
        profile[
          Transition::fromMatch
        ][
          TransitionFromMatch::toMatch
        ];
      if( parameters.debug >= DEBUG_All ) {
        if( row_i == ( last_row - 1 ) ) {
          // If we are also needing to test that this isn't a transition to a
          // deletion (which is stored in the Match state in the last row),
          // then we will need to incorporate the emission probability too.
          cout << "[debug] partialPath: tmp_mvt_match_to_match (M->M) is " <<
            ( tmp_mvt_match_to_match *
              (
                getEmissionProbability(
                  profile[ row_i ],
                  sequence,
                  current_sampler_state.column - 1
                )
              )
            ) << endl;
        } else {
          cout << "[debug] partialPath: tmp_mvt_match_to_match (M->M) is " << tmp_mvt_match_to_match << endl;
        }
      } // End if DEBUG_ALL
      tmp_mvt_deletion_to_match =
        forward_row[ current_sampler_state.column - 1 ][ Deletion ];
      tmp_mvt_deletion_to_match *=
        profile[
          Transition::fromDeletion
        ][
          TransitionFromDeletion::toMatch
        ];
      if( parameters.debug >= DEBUG_All ) {
        if( row_i == ( last_row - 1 ) ) {
          // If we are also needing to test that this isn't a transition to a
          // deletion (which is stored in the Match state in the last row),
          // then we will need to incorporate the emission probability too.
          cout << "[debug] partialPath: tmp_mvt_deletion_to_match (D->M) is " <<
            ( tmp_mvt_deletion_to_match *
              (
                getEmissionProbability(
                  profile[ row_i ],
                  sequence,
                  current_sampler_state.column - 1
                )
              )
            ) << endl;
        } else {
          cout << "[debug] partialPath: tmp_mvt_deletion_to_match (D->M) is " << tmp_mvt_deletion_to_match << endl;
        }
      } // End if DEBUG_ALL
      tmp_mvt_gap_to_match =
        ( tmp_mvt_deletion_to_match + tmp_mvt_insertion_to_match );

      if( row_i == ( last_row - 1 ) ) {
        // Note, though, that we store Deletions-into-the-final-row in that
        // row's Match state, since postalign insertions are non-affine and can
        // follow deletions.  Thus when row_i == ( last_row - 1 ), we need to
        // consider the possibility of a transition to a deletion as well.

        // Then we either continued a deletion or started one.
        tmp_mvt_deletion_to_deletion =
          forward_row[ current_sampler_state.column ][ Deletion ];
        tmp_mvt_deletion_to_deletion *=
          profile[
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toDeletion
          ];
        if( parameters.debug >= DEBUG_All ) {
          cout << "[debug] partialPath: tmp_mvt_deletion_to_deletion (D->D) is " << tmp_mvt_deletion_to_deletion << endl;
        }
        tmp_mvt_match_to_deletion =
          forward_row[ current_sampler_state.column ][ Match ];
        tmp_mvt_match_to_deletion *=
          profile[
            Transition::fromMatch
          ][
            TransitionFromMatch::toDeletion
          ];

        if( parameters.debug >= DEBUG_All ) {
          cout << "[debug] partialPath: tmp_mvt_match_to_deletion (M->D) is " << tmp_mvt_match_to_deletion << endl;
        }

        tmp_mvt_to_deletion =
          ( tmp_mvt_deletion_to_deletion + tmp_mvt_match_to_deletion );

        // Note that we never factored into the tmp_mvt_*_to_* the
        // matrixRowScaleFactor or the rabinerInverseScalar, but it doesn't
        // matter since the probability we calculate is a ratio and these
        // factors would affect the numerator and denominator equally.
        tmp_mvt =
          tmp_mvt_match_to_match;
        tmp_mvt +=
          tmp_mvt_gap_to_match;
        tmp_mvt *=
          getEmissionProbability(
            profile[ row_i ],
            sequence,
            current_sampler_state.column - 1
          );
        tmp_mvt +=
          tmp_mvt_to_deletion;
        // This is the prob that the transition is to a deletion.  Note that
        // here we do need to multiply in the emission probabilities into the
        // to-Match transitions (see above).
        ProbabilityType probability_of_deletion =
          tmp_mvt_to_deletion;
        probability_of_deletion /=
          tmp_mvt;

        if( parameters.debug >= DEBUG_All ) {
          cout << "[debug] partialPath: probability_of_deletion (?->D) is " << probability_of_deletion << endl;
        }
        if( probability_of_deletion == 0 ) {
          unif = 1.0;
        } else if( probability_of_deletion == 1 ) {
          unif = 0.0;
        } else {
          unif = random.nextUniform();
          if( parameters.debug >= DEBUG_All ) {
            cout << "[debug] partialPath: unif (?->D) is " << unif << endl;
          }
        }
        if( unif < probability_of_deletion ) {
          // Okay so it was a transition to a Deletion.  Now we need to see
          // about the source of the transition.

          // Note that we never factored into the tmp_mvt_*_to_deletion the
          // matrixRowScaleFactor or the rabinerInverseScalar, but it doesn't
          // matter since the probability we calculate is a ratio and these
          // factors would affect the numerator and denominator equally.

          // This is prob that the ?->D transition is a D->D transition
          probability_of_deletion =
            tmp_mvt_deletion_to_deletion;
          probability_of_deletion /=
            tmp_mvt_to_deletion;

          if( parameters.debug >= DEBUG_All ) {
            cout << "[debug] partialPath: probability_of_deletion (D->D) is " << probability_of_deletion << endl;
          }
          if( probability_of_deletion == 0 ) {
            unif = 1.0;
          } else if( probability_of_deletion == 1 ) {
            unif = 0.0;
          } else {
            unif = random.nextUniform();
            if( parameters.debug >= DEBUG_All ) {
              cout << "[debug] partialPath: unif (D->D) is " << unif << endl;
            }
          }
          if( unif < probability_of_deletion ) {
            target_sampler_state.subcell = Deletion;
            if( position_counts != NULL ) {
              ( *position_counts )[
                Transition::fromDeletion
              ][
                TransitionFromDeletion::toDeletion
              ] += 1;
            } // End if position_counts != NULL
          } else {
            target_sampler_state.subcell = Match;
            if( position_counts != NULL ) {
              ( *position_counts )[
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ] += 1;
            } // End if position_counts != NULL
          }
          target_sampler_state.column = current_sampler_state.column;
          // Note that we don't increment the sequence_advances, because this
          // was a deletion.
          return;
        } // End if this is a transition to a Deletion, but stored in the Match
          // state in the last row.
        // .. else it is a transition to a Match, and we can proceed as usual.
      } // End if row_i == ( last_row - 1 )

      // At this point we know that there was a Match to the subsequent row.
      if( sequence_advances != NULL ) {
        ( *sequence_advances )[ row_i + 2 ] += 1;
      } // End if sequence_advances != NULL
      
      ProbabilityType probability_of_match;
      // This is the prob that the *source* is a Match (we know that the target
      // is a Match).
      if( tmp_mvt_gap_to_match == 0 ) {
        probability_of_match = 1;
      } else {
        // Note that we never factored into the tmp_mvt_*_to_match the
        // matrixRowScaleFactor or the rabinerInverseScalar, but it doesn't
        // matter since the probability we calculate is a ratio and these
        // factors would affect the numerator and denominator equally.
        // Nor did we factor in the match emission, for the same reason.
        tmp_mvt =
          tmp_mvt_match_to_match;
        tmp_mvt +=
          tmp_mvt_gap_to_match;
        probability_of_match =
          tmp_mvt_match_to_match;
        probability_of_match /=
          tmp_mvt;

      }
      if( parameters.debug >= DEBUG_All ) {
        cout << "[debug] partialPath: probability_of_match (M->M) is " << probability_of_match << endl;
      }
      if( probability_of_match == 0 ) {
        unif = 1.0;
      } else if( probability_of_match == 1.0 ) {
        unif = 0.0;
      } else {
        unif = random.nextUniform();
        if( parameters.debug >= DEBUG_All ) {
          cout << "[debug] partialPath: unif (M->M) is " << unif << endl;
        }
      }
      if( unif < probability_of_match ) {
        target_sampler_state.subcell = Match;
        target_sampler_state.column =
          current_sampler_state.column - 1;
        if( position_counts != NULL ) {
          ( *position_counts )[
            Transition::fromMatch
          ][
            TransitionFromMatch::toMatch
          ] += 1;
          incrementEmissionCounts(
            *position_counts,
            sequence,
            current_sampler_state.column - 1
          );
        } // End if position_counts != NULL
        return;
      } // End if unif < probability_of_match

      // Okay so if it wasn't a Match, then it was either an Insertion or a
      // Deletion.
      ProbabilityType probability_of_insertion;
      if( current_sampler_state.column == 1 ) {
        // Note that the prob of I->M from col 0 to col 1 is always 0.
        probability_of_insertion = 0;
      } else { // if current_sampler_state.column == 1 .. else .. 
        probability_of_insertion =
          tmp_mvt_insertion_to_match;
        probability_of_insertion /=
          tmp_mvt_gap_to_match;

      } // End if current_sampler_state.column == 1 .. else .. 
      if( parameters.debug >= DEBUG_All ) {
        cout << "[debug] partialPath: probability_of_insertion (I->M) is " << probability_of_insertion << endl;
      }
      if( probability_of_insertion == 0 ) {
        unif = 1.0;
      } else if( probability_of_insertion == 1 ) {
        unif = 0.0;
      } else {
        unif = random.nextUniform();
        if( parameters.debug >= DEBUG_All ) {
          cout << "[debug] partialPath: unif (I->M) is " << unif << endl;
        }
      }
      if( unif < probability_of_insertion ) {
        target_sampler_state.subcell = Insertion;
        target_sampler_state.column =
          current_sampler_state.column - 1;
        if( position_counts != NULL ) {
          ( *position_counts )[
            Transition::fromInsertion
          ][
            TransitionFromInsertion::toMatch
          ] += 1;
        } // End if position_counts != NULL
        if( position_counts != NULL ) {
          incrementEmissionCounts(
            *position_counts,
            sequence,
            current_sampler_state.column - 1
          );
        } // End if position_counts != NULL
      } else { // if unif < probability_of_insertion .. else ..
        target_sampler_state.subcell = Deletion;
        target_sampler_state.column =
          current_sampler_state.column - 1;
        if( position_counts != NULL ) {
          ( *position_counts )[
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toMatch
          ] += 1;
          incrementEmissionCounts(
            *position_counts,
            sequence,
            current_sampler_state.column - 1
          );
        } // End if position_counts != NULL
        return;
      } // End if unif < probability_of_insertion .. else ..

      // Okay so if it was Match then we're done because we know there's no
      // insertions preceding a Match (see above, where we return if it is a
      // Match), but if it was an Insertion, we need to go back to the Match
      // before the insertion began.
      bool done = false;

      while( !done ) {
        if( target_sampler_state.column == 1 ) {
          // Note that the first column has only values in the Deletion subcell
          // (except in row 0, where it only has a value in the Match subcell),
          // so the prob of I->I from col 0 to col 1 is always 0.
          probability_of_insertion = 0;
        } else { // if the probability is 0 .. else .. 
          probability_of_insertion =
            forward_row[ target_sampler_state.column - 1 ][ Insertion ];
          probability_of_insertion *=
            profile[
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toInsertion
            ];
          probability_of_insertion *=
            getInsertionEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              target_sampler_state.column - 1
            );
          probability_of_insertion /=
            forward_row[ target_sampler_state.column ][ Insertion ];
        } // End if the probability is 0 .. else .. 
        if( parameters.debug >= DEBUG_All ) {
          cout << "[debug] partialPath: probability_of_insertion (I->I) is " << probability_of_insertion << endl;
        }
        if( probability_of_insertion == 0 ) {
          unif = 1.0;
        } else if( probability_of_insertion == 1 ) {
          unif = 0.0;
        } else {
          unif = random.nextUniform();
          if( parameters.debug >= DEBUG_All ) {
            cout << "[debug] partialPath: unif (I->I) is " << unif << endl;
          }
        }
        if( unif < probability_of_insertion ) {
          if( position_counts != NULL ) {
            ( *position_counts )[
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toInsertion
            ] += 1;
            incrementInsertionEmissionCounts(
              *position_counts,
              sequence,
              target_sampler_state.column - 1
            );
          } // End if position_counts != NULL
          target_sampler_state.column -= 1;
        } else { // if unif < probability_of_insertion .. else ..
          target_sampler_state.subcell = Match;
          if( position_counts != NULL ) {
            ( *position_counts )[
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ] += 1;
            incrementInsertionEmissionCounts(
              *position_counts,
              sequence,
              target_sampler_state.column - 1
            );
          } // End if position_counts != NULL
          target_sampler_state.column -= 1;
          done = true;
        } // End if unif < probability_of_insertion .. else ..
        // Increment the sequence advances for the current row, too.
        if( sequence_advances != NULL ) {
          ( *sequence_advances )[ row_i + 1 ] += 1;
        } // End if sequence_advances != NULL
      } // End while continuing to do insertions...
      // At this point the subcell should be Match.
      assert( target_sampler_state.subcell == Match );
    } // drawPartialPathForPosition( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, SamplerState const&, SamplerState &, Random &, GlobalCountsType *, PositionCountsType *, vector<uint32_t> * ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename CountsType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Given the forward matrices corresponding to the given profile and
     * sequences, draw paths from the distribution of the paths, given the
     * sequences and profiles; the results are stored in either the given
     * counts arguments or the given multiple alignment argument (one of which
     * should be non-null, or this method is a waste of time).
     *
     * If the CountsType argument pointer is non-null, then it should be a
     * Profile type with the same length as the profile argument, used for
     * storing unsigned integer counts.  These counts will be updated with the
     * number of the observations on the drawn partial path corresponding to
     * each parameter.
     *
     * If the given MultipleAlignment pointer is non-null, then it will be
     * reinitialized and modified to hold the drawn paths.
     */
    drawPaths (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const & sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer const& forward_matrices,
      Random & random,
      CountsType * counts,
      MultipleAlignment<ProfileType, SequenceResidueType> * multiple_alignment
    ) const
    {
      sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );

      if( multiple_alignment != NULL ) {
        multiple_alignment->reinitialize(
          &profile,
          &sequences,
          sequence_count
        );
      } // End if multiple_alignment != NULL, reinitialize it.

      if( counts != NULL ) {
        counts->reinitialize( profile ); // reset the counts, copying length
                                         // from the profile.
      } // End if counts != NULL, reinitialize it.
      

      // Every other time, we must swap which of the SamplerState vectors is
      // current.
      vector<SamplerState> sampler_states_1( sequence_count );
      vector<SamplerState> sampler_states_2( sequence_count );
      bool swap = true;

      uint32_t seq_i;
      uint32_t row_i;
      uint32_t last_row = profile.length();
      typename Matrix::SequentialAccessContainer::const_reverse_iterator forward_matrices_rev_iter;
      for(
        forward_matrices_rev_iter = forward_matrices.rbegin(), row_i = last_row;
        forward_matrices_rev_iter != forward_matrices.rend();
        forward_matrices_rev_iter++, row_i--
      ) {
        // Every other time, we must swap which of the SamplerStates is current.
        swap = !swap;

        for( seq_i = 0; seq_i < sequence_count; seq_i++ ) {
          drawPartialPathForPosition(
            parameters,
            profile,
            sequences[ seq_i ],
            row_i,
            ( *forward_matrices_rev_iter )[ seq_i ],
            ( swap ? sampler_states_2[ seq_i ] : sampler_states_1[ seq_i ] ), // ignored if row_i == last_row
            ( swap ? sampler_states_1[ seq_i ] : sampler_states_2[ seq_i ] ),
            random,
            counts,
            ( ( counts == NULL ) ? NULL : &( ( *counts )[ row_i ] ) ),
            ( ( multiple_alignment == NULL ) ? NULL :
              &( multiple_alignment->m_sequenceAdvances[ seq_i ] ) )
          );
          // TODO: REMOVE
          //cout << "for row_i = " << row_i << ", SamplerState for seq_i = " << seq_i << " is " << ( swap ? sampler_states_1[ seq_i ] : sampler_states_2[ seq_i ] ) << endl;
        } // End foreach seq_i

      } // End foreach row, in reverse

      // TODO: REMOVE
      //if( multiple_alignment != NULL ) {
      //  for( seq_i = 0; seq_i < sequence_count; seq_i++ ) {
      //    cout << "sequence advances for sequence " << seq_i << " is:" << endl;
      //    cout << "( ";
      //    for( row_i = 0; row_i <= ( last_row + 3 ); row_i++ ) {
      //      if( row_i > 0 ) {
      //        cout << ", ";
      //      }
      //      cout << multiple_alignment->m_sequenceAdvances[ seq_i ][ row_i ];
      //    }
      //    cout << ") " << endl;
      //  }
      //}
    } // drawPaths( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t, SequentialAccessContainer const&, Random &, GlobalCountsType *, PositionCountsType *, MultipleAlignment<ProfileType, SequenceResidueType> * ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Given the backward row for row row_i, and the forward row for row_i - 1,
     * calculate the PositionSpecificSequenceScoreCoefficients for position
     * (row_i-1) of the given profile for the given sequence.  Results go in
     * the given PositionSpecificSequenceScoreCoefficients reference.  If
     * parameters.useRabinerScaling is true, the result will be scaled (by the
     * inverse of the product of the rabinerCumulativeInverseScalar values in
     * the given Row references, which is the sequence score).  The inverse of
     * that scalar will be stored in the m_inverseScalar field of the
     * coefficients.  An unscaled version can be retrieved by calling
     * createUnscaledCopy() on the coefficients.
     *
     * NOTE: row_i must be > 0.
     */
    calculatePositionSpecificSequenceScoreCoefficients (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& backward_row,
      PositionSpecificSequenceScoreCoefficients & coefficients
    ) const
    {
      assert( row_i > 0 );

      // TODO: REMOVE
      static const bool do_extra_debugging = false;

      const uint32_t last_row = profile.length();
      const uint32_t last_col = sequence.length();

      SequenceResidueType residue;

      register uint32_t col_i;

      // Either MatrixValueType is big enough, or we are using rabiner scaling,
      // in which case things are scaled so that it is big enough.
      MatrixValueType path_likelihood;
      MatrixValueType tmp;

      coefficients.zero();

      if( false && ( parameters.debug >= DEBUG_All ) ) {
        cout << "[debug] coefficients: last_row is " << last_row << ", last_col is " << last_col << endl;
      }

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      if( row_i >= 2 ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using prev_row_match_distribution == " << prev_row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if row_i >= 2
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      for( col_i = 0; col_i <= last_col; col_i++ ) {

        if( false && ( parameters.debug >= DEBUG_All ) ) {
          cout << "[debug] coefficients: row " << row_i << ", col " << col_i << endl;
        }
        
        // calculate row cell:
        if( col_i > 0 ) {
          // If not in the leftmost column, then there's base emission(s)
          // from this cell.
            residue = sequence[ col_i - 1 ];
//            if( parameters.debug >= DEBUG_All ) {
//              cout << "[debug] coefficients: residue is " << residue << endl;
//            }
        } // End if col_i > 0, get the residue.

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        // to/from DeletionIn
        // from DeletionIn to Match
        if(
          ( row_i > 0 ) && ( col_i >= 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 1 ) 
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toMatch
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp =
            prev_forward_row.m_deletionIn;
#else
          tmp =
            prev_forward_row[ col_i - 1 ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS
          tmp /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the Match( " << residue << " ) value for the del-in end." << endl;
          }
          coefficients[ Emission::Match ][ residue ] +=
            path_likelihood;
        } // End if ( row_i > 0 ) && ( col_i >= 1 ) ( and maybe col_i == 1 )

        // from DeletionIn to DeletionIn (and from Match/Begin to DeletionIn)
        if(
          ( row_i == 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // Actually this is from PreAlign to Begin, from Begin to DeletionIn.
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionIn;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ];
          path_likelihood *=
            profile[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletionIn
            ];
          tmp =
            prev_forward_row[ col_i ][ Match ];
          tmp /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-in open." << endl;
          }
          coefficients.m_constant +=
            path_likelihood;
        } else if(
          ( row_i != 0 ) && ( row_i != last_row )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionIn;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toDeletionIn
            ];

#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp =
            prev_forward_row.m_deletionIn;
#else
          tmp =
            prev_forward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          tmp /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-in extensions." << endl;
          }
          coefficients.m_constant +=
            path_likelihood;
        } // End if ( row_i == 1 ) .. else if( ( row_i != 0 ) && ( row_i != last_row ) ) (and maybe also ( col_i == 0 )).
        // End to/from DeletionIn

        // From DeletionOut
        // from DeletionOut to Match
        // Only happens in the final row, which we've already accounted for in
        // the backward row calculation.

        // from DeletionOut to DeletionOut
        if(
          ( row_i > 2 )
          && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionOut;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toDeletionOut
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp =
            prev_forward_row.m_deletionOut;
#else
          tmp =
            prev_forward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          tmp /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-out extensions." << endl;
          }
          coefficients.m_constant +=
            path_likelihood;
        } // End if ( row_i > 2 ) && ( col_i > 0 ) (and maybe col_i == last_col)
        // End from DeletionOut
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        // from Match
        // from Match to Match
        if( col_i > 0 ) {
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          if( row_i == 1 ) {
            path_likelihood *=
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            path_likelihood *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toMatch
              ];
          } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              prev_row_match_distribution[
                TransitionFromMatch::toMatch
              ];
#else
            path_likelihood *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toMatch
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          }
          tmp =
            prev_forward_row[ col_i - 1 ][ Match ];
          tmp /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          // TODO: REMOVE
          //cout << "[" << col_i << "] M->M: path likelihood is " << path_likelihood << endl;
          coefficients[ Emission::Match ][ residue ] +=
            path_likelihood;
        } // End if col_i == 0, elsif row_i == 0, else ..

        // from Match to Deletion
        if( row_i == 1 ) {
          path_likelihood =
            backward_row[ col_i ][ Deletion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ];
          path_likelihood *=
            profile[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletion
            ];
          tmp = prev_forward_row[ col_i ][ Match ];
          tmp /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          coefficients.m_constant +=
            path_likelihood;
        } else if( row_i == last_row ) {
          // We store the postAlignInsertion stuff in the Match state.
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          path_likelihood *=
            prev_row_match_distribution[
              TransitionFromMatch::toDeletion
            ];
#else
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
#ifdef USE_END_DISTRIBUTION
          path_likelihood *=
            profile[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ];
#endif // USE_END_DISTRIBUTION
          tmp = prev_forward_row[ col_i ][ Match ];
          tmp /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          coefficients.m_constant +=
            path_likelihood;
        } else { // row_i is not 0, 1, or last_row
          path_likelihood =
            backward_row[ col_i ][ Deletion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          path_likelihood *=
            prev_row_match_distribution[
              TransitionFromMatch::toDeletion
            ];
#else
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          tmp = prev_forward_row[ col_i ][ Match ];
          tmp /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          coefficients.m_constant +=
            path_likelihood;
        } // End if it's row 1 .. elsif the last row .. else ..
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        // from Match to DeletionOut
        if(
          ( row_i > 1 ) && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionOut;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;

          // Note that we do not use the scaled distribution, since the
          // backwards row incorporates the rest of the del-out path.  We just
          // want the deletion OPEN probability.
          path_likelihood *=
            profile[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletionOut
            ];
          tmp = prev_forward_row[ col_i ][ Match ];
          tmp /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-out open." << endl;
          }
          coefficients.m_constant +=
            path_likelihood;
        } // End if ( row_i > 1 ) && ( col_i > 0 ) (and maybe col_i == last_col)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        // End from Match

        // from Insertion
        // from Insertion to Match
        if( ( col_i != 0 ) &&
            ( row_i != 1 ) ) { // No ins state in row 0
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toMatch
            ];
          tmp =
            prev_forward_row[ col_i - 1 ][ Insertion ];
          tmp /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
            coefficients[ Emission::Match ][ residue ] +=
              path_likelihood;
        } // End if it's not the first or second row and not the first col ..
        // End from insertion

        // from Deletion
        // from Deletion to Match
        if( row_i == last_row ) {
          // We store post-align stuff in the Match state, even if it comes
          // from a deletion.  So this is really a D->D, stored at M.
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
#ifdef USE_END_DISTRIBUTION
          path_likelihood *=
            profile[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ];
#endif // USE_END_DISTRIBUTION
          tmp = prev_forward_row[ col_i ][ Deletion ];
          tmp /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          coefficients.m_constant +=
            path_likelihood;
        } // End if it's the last row..
        if( ( col_i != 0 ) &&
            ( row_i != 1 ) ) { // no del state in row 0
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toMatch
            ];
          tmp =
            prev_forward_row[ col_i - 1 ][ Deletion ];
          tmp /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          coefficients[ Emission::Match ][ residue ] +=
            path_likelihood;
        } // End if it's not the first or second row and not the first col..
        // from Deletion to Deletion
        if( ( row_i != 1 ) && // no del state in first row
            ( row_i != last_row ) ) {
          path_likelihood =
            backward_row[ col_i ][ Deletion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
          tmp = prev_forward_row[ col_i ][ Deletion ];
          tmp /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp;

          coefficients.m_constant +=
            path_likelihood;
        } // End if it's not the first, second, or last row ..
        // End from Deletion
      } // End for each column, ..

      if( parameters.useRabinerScaling ) {
        coefficients.m_inverseScalar =
          prev_forward_row.m_rabinerCumulativeInverseScalar;
        coefficients.m_inverseScalar *=
          backward_row.m_rabinerCumulativeInverseScalar;
      }
      // else coefficients.m_inverseScalar == 1.0
    } // calculatePositionSpecificSequenceScoreCoefficients( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, Row const&, PositionSpecificSequenceScoreCoefficients & ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType>
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Use the given PositionSpecificSequenceScoreCoefficients reference to
     * update the entente position corresponding to the given profile position.
     * The entente position will be updated by adding the expected number of
     * emissions of each type (divided by the given inverse scalar, if it is
     * non-NULL, or by the calculated sequence score otherwise).  The return
     * value is the (not scaled) sequence score.  Don't forget to zero() the
     * position entente (before calling this for all sequences).  If the
     * inverse_scalar pointer is non-NULL but points to the value 0, the score
     * will be calculated and returned but the position entente will not be
     * updated (it's like having a 0 scalar, not like having a 0 inverse
     * scalar, which would result in a divide-by-0 error).
     */
    updatePositionEntenteForSequence (
      Parameters const& parameters,
      ProfilePosition<ResidueType, ProfileType> const& profile_position,
      PositionSpecificSequenceScoreCoefficients const & coefficients,
      ScoreType const * inverse_scalar,
      PositionEntente & position_entente
    ) const
    {
      PositionEntente position_entente_update;
      ScoreType sequence_score =
        calculatePositionEntenteUpdate(
          parameters,
          profile_position,
          coefficients,
          inverse_scalar,
          position_entente_update
        );

      if( !( ( inverse_scalar != NULL ) && ( *inverse_scalar == 0.0 ) ) ) {
        // Now do the updating, scaling the position_entente appropriately.
        position_entente += position_entente_update;
      }

      return sequence_score;
    } // updatePositionEntenteForSequence( Parameters const&, ProfilePosition const&, PositionSpecificSequenceScoreCoefficients const&, ScoreType const *, PositionEntente & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType>
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Use the given PositionSpecificSequenceScoreCoefficients reference to
     * calculate the update for the entente position corresponding to the given
     * profile position.  The calculated update should be added to the entente
     * position.  This will add the expected number of emissions of each type
     * (divided by the given inverse scalar, if it is non-NULL, or by the
     * calculated sequence score otherwise).  If the inverse_scalar pointer is
     * non-NULL but points to the value 0, the score will be calculated and
     * returned but the position_entente_update will be zero()ed (it's like
     * having a 0 scalar, not like having a 0 inverse scalar, which would
     * result in a divide-by-0 error).  The return value is the (not scaled)
     * sequence score.
     */
    calculatePositionEntenteUpdate (
      Parameters const& parameters,
      ProfilePosition<ResidueType, ProfileType> const& profile_position,
      PositionSpecificSequenceScoreCoefficients const & coefficients,
      ScoreType const * inverse_scalar,
      PositionEntente & position_entente_update
    ) const
    {
      //position_entente_update.zero(); // redundant.  Done in calculateScore(..).

      // TODO: REMOVE
      //cout << "calculatePEU: coefficients are " << coefficients << endl;

      // Calculate the score and the expected emission counts.  Note that if
      // useRabinerScaling is true, the position_entente_update is already
      // scaled by the original sequence score (that is, by
      // 1/coefficients.m_inverseScalar).  Regardless of useRabinerScoring,
      // position_entente_update will be scaled such that is maximum value is
      // 1; dividing by position_entente_update.m_scalar will unscale it
      // completely.
      ScoreType sequence_score =
        coefficients.calculateScore(
          profile_position,
          &position_entente_update
        );
      // TODO: REMOVE
      //cout << "calculatePEU: sequence_score is " << sequence_score << endl;

      // TODO: REMOVE
      //cout << "calculatePEU: position_entente_update (before scaling) is " << position_entente_update << endl;
      //cout << "calculatePEU: unscaled position_entente_update (before scaling) is " << position_entente_update.createUnscaledCopy() << endl;

      if( ( inverse_scalar != NULL ) && ( *inverse_scalar == 0.0 ) ) {
        position_entente_update.zero();
        return sequence_score;
      }

      if( inverse_scalar == NULL ) {
        // TODO: REMOVE
        //cout << "calculatePEU: sequence_score is " << sequence_score << endl;
        //position_entente_update /= sequence_score;
        //( sequence_score * ( position_entente_update.m_scalar / position_entente.m_scalar ) );
        position_entente_update.m_scalar *= sequence_score;
        //if( parameters.useRabinerScaling ) {
        //  // It's already been scaled, but by the former sequence score.
        //  position_entente_update /=
        //    ( sequence_score / ( coefficients.m_inverseScalar * position_entente.m_scalar ) );
        //  //position_entente_update.m_scalar *=
        //  //( sequence_score / coefficients.m_inverseScalar );
        //} else {
        //  position_entente_update /=
        //    ( sequence_score / position_entente.m_scalar );
        //  //position_entente_update.m_scalar *= sequence_score;
        //  // TODO: REMOVE
        //cout << "sequence_score is " << sequence_score << "; position_entente_update.m_scalar is " << position_entente_update.m_scalar << "; ( sequence_score / position_entente_update.m_scalar ) is " << ( sequence_score / position_entente_update.m_scalar ) << endl;
        //} // End if useRabinerScaling .. else ..
        // TODO: REMOVE
        //cout << "calculatePEU: position_entente_update (after scaling by the inverse of the sequence score) is " << position_entente_update << endl;
      } else {
        // TODO: REMOVE
        //cout << "calculatePEU: *inverse_scalar is " << *inverse_scalar << endl;
        //position_entente_update /= *inverse_scalar;
        //( *inverse_scalar * ( position_entente_update.m_scalar / position_entente.m_scalar ) );
        position_entente_update.m_scalar *= *inverse_scalar;
        //if( parameters.useRabinerScaling ) {
        //  // It's already been scaled, but by the former sequence score.
        //  position_entente_update /=
        //    ( *inverse_scalar / ( coefficients.m_inverseScalar * position_entente.m_scalar ) );
        //  //position_entente_update.m_scalar *=
        //  //  ( *inverse_scalar / coefficients.m_inverseScalar );
        //} else {
        //  position_entente_update /=
        //    ( *inverse_scalar / position_entente.m_scalar );
        //  //position_entente_update.m_scalar *= *inverse_scalar;
        //} // End if useRabinerScaling .. else ..
        // TODO: REMOVE
        //cout << "calculatePEU: position_entente_update (after scaling by the inverse of the argument value) is " << position_entente_update << endl;
      }

      return sequence_score;
    } // calculatePositionEntenteUpdate( Parameters const&, ProfilePosition const&, PositionSpecificSequenceScoreCoefficients const&, ScoreType const *, PositionEntente & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Given the forward and backward rows for row row_i, and the forward row
     * for row_i - 1, update the given GlobalEntente for the given sequence,
     * profile pair.  Results go in the given GlobalEntente reference.  Don't
     * forget to zero() the entente (before calling this for all rows and all
     * sequences).
     *
     * NOTE: If useRabinerScaling is true but rabinerScaling_useMaximumValue is
     * also true, then you *MUST* supply a value for inverse_scalar.  We can't
     * calculate it here.
     * 
     * NOTE: If row_i is 0, prev_forward_row is ignored.
     */
    updateGlobalEntenteForSequence (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& forward_row,
      typename Matrix::Row const& backward_row,
      ScoreType const * inverse_scalar,
      GlobalEntente & global_entente
    ) const
    {
      GlobalEntente global_entente_update;
      global_entente_update.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
      // Note that global_entente_update.zero() won't zero the fromPreAlign::toPreAlign when DISALLOW_FLANKING_TRANSITIONS is set, so:
      global_entente_update[ Transition::fromPreAlign ].zero();
      global_entente_update[ Transition::fromPostAlign ].zero();

      // TODO: REMOVE
      assert( profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] == 0 );
      assert( profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] == 0 );
#endif // DISALLOW_FLANKING_TRANSITIONS

      const uint32_t last_row = profile.length();
      const uint32_t last_col = sequence.length();

      // TODO: Add support for 0-length sequences!  We should be learning
      // something about the deletion rate, at least.
      if( last_col == 0 ) {
        cout << "WARNING: 0-length sequences currently are not supported by updateGEFS(..)." << endl;
        assert( false );
        return 0;
      }

      SequenceResidueType residue;

      register uint32_t col_i;

      MatrixValueType path_likelihood, tmp_mvt;
      MatrixValueType matrix_value_type_score = 0;

      // Note that when we are using rabiner scaling, if
      // rabinerScaling_useMaximumValue is false then we do not need to
      // recalculate the score -- we just use the cumulative_inverse_scalars.
      // The score as calculated here ought to be 1 if useRabinerScaling is
      // true and rabinerScaling_useMaximumValue is false, so it can be used as
      // a sanity check.

      // If useRabinerScaling is true but rabinerScaling_useMaximumValue is
      // also true, then you *MUST* supply a value for inverse_scalar.  We
      // can't calculate it here.
      assert( ( parameters.useRabinerScaling && parameters.rabinerScaling_useMaximumValue ) ? ( inverse_scalar != 0 ) : true );
      bool calculate_score = true;
      //!( parameters.useRabinerScaling && !parameters.rabinerScaling_useMaximumValue );

      // NOTE: When useRabinerScaling is true and
      // rabinerScaling_useMaximumValue is false, then (assuming that the
      // matrix rows are up-to-date, which they might not be), the
      // cumulative_inverse_scalar *is* the score.
      ScoreType cumulative_inverse_scalar =
        backward_row.m_rabinerCumulativeInverseScalar;
      if( row_i > 0 ) {
        cumulative_inverse_scalar *=
          prev_forward_row.m_rabinerCumulativeInverseScalar;
      }

      // When we are using rabiner scaling, ->I transitions are within the row,
      // so are scaled by the current row's forward and backward
      // cumulativeInverseScalars (C_t and D_t), instead of the desired product
      // of the current row's backward cumulativeInverseScalar (D_t) with the
      // previous row's forward scalar (C_{t-1}).  To correct for this, we can
      // multiply them by the ratio C_{t-1}/C_t = 1/c_t, where c_t is the
      // rabiner scalar (not cumulative) for the current row.
      // We will call this the correction_factor:
      MatrixValueType rabiner_within_row_correction_factor = 1.0;
      if( parameters.useRabinerScaling ) {
        rabiner_within_row_correction_factor =
          forward_row.m_rabinerInverseScalar;
        assert( rabiner_within_row_correction_factor > 0 );
        // TODO: REMOVE
        //cout << "rabiner_within_row_correction_factor is " << rabiner_within_row_correction_factor;
        //cout << " (should ==) " <<
        //  toDouble( ( forward_row.m_rabinerCumulativeInverseScalar * 
        //              backward_row.m_rabinerCumulativeInverseScalar ) /
        //    cumulative_inverse_scalar ) << endl;

        // TODO: REMOVE?
        //assert(
        //  abs(
        //    toDouble(rabiner_within_row_correction_factor) -
        //    toDouble(
        //      (
        //        forward_row.m_rabinerCumulativeInverseScalar * 
        //        backward_row.m_rabinerCumulativeInverseScalar
        //      ) /
        //      cumulative_inverse_scalar
        //    )
        //  ) < ( 10.0*numeric_limits<double>::epsilon() )
        //);

      } // End if useRabinerScaling

      // TODO: REMOVE
      //if( parameters.useRabinerScaling ) {
      //  cout << "updateGEFS: cumulative_inverse_scalar is " << cumulative_inverse_scalar << endl;
      //  if( row_i == 0 ) {
      //    cout << "\t = backward_row.m_rabinerCumulativeInverseScalar." << endl;
      //  } else {
      //    cout << "\t = ( " <<
      //      prev_forward_row.m_rabinerCumulativeInverseScalar << " * " <<
      //      backward_row.m_rabinerCumulativeInverseScalar << " )" << endl;
      //  }
      //}

      if( false && ( parameters.debug >= DEBUG_All ) ) {
        cout << "[debug] updateGEFS(..) row_i is " << row_i << "; last_row is " << last_row << ", last_col is " << last_col << endl;
      }

      // TODO: REMOVE
      if( false ) {//true ) {
        cout << "updateGEFS: forward_row is " << forward_row << endl;
        cout << "updateGEFS: backward_row is " << backward_row << endl;
        typename Matrix::Row unscaled_row = forward_row;
        if( parameters.useRabinerScaling ) {
          unscaled_row *= forward_row.m_rabinerCumulativeInverseScalar;
          cout << "updateGEFS: rabiner-unscaled forward_row is " << unscaled_row << endl;  
        }
        unscaled_row /= parameters.matrixRowScaleFactor;
        cout << "updateGEFS: unscaled forward_row is " << unscaled_row << endl;  
        unscaled_row = backward_row;
        if( parameters.useRabinerScaling ) {
          unscaled_row *= backward_row.m_rabinerCumulativeInverseScalar;
          cout << "updateGEFS: rabiner-unscaled backward_row is " << unscaled_row << endl;
        }
        unscaled_row /= parameters.matrixRowScaleFactor;
        cout << "updateGEFS: unscaled backward_row is " << unscaled_row << endl;  
      } // End if( false ) (DEBUGGING)

      // TODO: REMOVE
      //cout << "updateGEFS(..): sequence is " << sequence << endl;

      for( col_i = 0; col_i <= last_col; col_i++ ) {

        //if( parameters.debug >= DEBUG_All ) {
        //  cout << "[debug] (entente) row " << row_i << ", col " << col_i << endl;
        //}
        
        // calculate row cell:
        if( col_i > 0 ) {
          // If not in the leftmost column, then there's base emission(s)
          // from this cell.
            residue = sequence[ col_i - 1 ];
//            if( parameters.debug >= DEBUG_All ) {
//              cout << "[debug ] (entente) residue is " << residue << endl;
//            }
        }
    
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        // to/from DeletionIn
        // from DeletionIn to Match
        if(
          ( row_i > 0 ) && ( col_i >= 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS // & USE_DEL_IN_DEL_OUT
          && ( col_i == 1 ) 
#endif // DISALLOW_FLANKING_TRANSITIONS ( & USE_DEL_IN_DEL_OUT )
        ) {
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toMatch
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS // & USE_DEL_IN_DEL_OUT
          tmp_mvt =
            prev_forward_row.m_deletionIn;
#else
          tmp_mvt =
            prev_forward_row[ col_i - 1 ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS ( & USE_DEL_IN_DEL_OUT )
          tmp_mvt /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
#ifdef USE_END_DISTRIBUTION // & USE_DEL_IN_DEL_OUT
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION ( & USE_DEL_IN_DEL_OUT )

          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Match ][ residue ];

          global_entente_update[ 
            Transition::fromDeletionIn
          ][
            TransitionFromDeletionIn::toMatch
          ] += path_likelihood;

          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            global_entente_update[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else {
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
            } // End if calculate_score
          } // End if row_i == last_row .. else ..
        } // End if ( row_i > 0 ) && ( col_i >= 1 ) ( and maybe col_i == 1 )

        // from DeletionIn to DeletionIn (and from Match/Begin to DeletionIn)
        if(
          ( row_i == 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // Actually this is from PreAlign to Begin, from Begin to DeletionIn.
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionIn;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ];
          path_likelihood *=
            profile[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletionIn
            ];
          tmp_mvt =
            prev_forward_row[ col_i ][ Match ];
          tmp_mvt /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          global_entente_update[
            Transition::fromPreAlign
          ][
            TransitionFromPreAlign::toBegin
          ] += path_likelihood;
          global_entente_update[
            Transition::fromBegin
          ][
            TransitionFromBegin::toDeletionIn
          ] += path_likelihood;

          if( calculate_score ) {
            matrix_value_type_score += path_likelihood;
          } // End if calculate_score
        } else if(
          ( row_i != 0 ) && ( row_i != last_row )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionIn;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toDeletionIn
            ];

#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_mvt =
            prev_forward_row.m_deletionIn;
#else
          tmp_mvt =
            prev_forward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          tmp_mvt /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          global_entente_update[
            Transition::fromDeletionIn
          ][
            TransitionFromDeletionIn::toDeletionIn
          ] += path_likelihood;

          if( calculate_score ) {
            matrix_value_type_score += path_likelihood;
          }
        } // End if ( row_i == 1 ) .. else if( ( row_i != 0 ) && ( row_i != last_row ) ) (and maybe also ( col_i == 0 )).
        // End to/from DeletionIn

        // From DeletionOut
        // from DeletionOut to DeletionOut
        // and from DeletionOut to End
        if(
          ( row_i > 2 )
          && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionOut;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toDeletionOut
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_mvt =
            prev_forward_row.m_deletionOut;
#else
          tmp_mvt =
            prev_forward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          // Note we don't need to multiply in the DeletionOut::toEnd or the
          // End::toPostAlign probabilities, because these are already
          // incorporated into the backward_row.  But we do need to count them.
          if( row_i == last_row ) {
            // TODO: REMOVE
            //cout << "For W->W INTO last_row, adding " << path_likelihood << " to the W->E GEU value." << endl;
            global_entente_update[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toEnd
            ] += path_likelihood;
#ifdef USE_END_DISTRIBUTION
            // If this is the last row, that transition includes from End to
            // PostAlign.
            global_entente_update[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } // End if row_i == last_row

          global_entente_update[
            Transition::fromDeletionOut
          ][
            TransitionFromDeletionOut::toDeletionOut
          ] += path_likelihood;

          if( calculate_score && ( row_i != last_row ) ) {
            matrix_value_type_score += path_likelihood;
          }
        } // End if ( row_i > 2 ) && ( col_i > 0 ) (and maybe col_i == last_col)
        // End from DeletionOut
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        // from Match
        // from Match to Match
        if( ( col_i == 0 ) && ( last_col != 0 ) ) {
          if( row_i == 0 ) {
            if( calculate_score ) {
              // Note in row 0 the score is completely calculated using col 0.
              matrix_value_type_score = // assuming forward row has a 1 here
                backward_row[ 0 ][ Match ];
              // Note that even if USE_DEL_IN_DEL_OUT, we need not add anything
              // else, since those values are already rolled into the Match
              // cell here.
              matrix_value_type_score *=
                parameters.matrixRowScaleFactor;
            } // End if calculate_score
          }
          // Do nothing else.
        } else if( row_i == 0 ) { // col_i > 0 && row_i == 0
          // Note in row 0 the score is completely calculated using col 0.

          // In the top row, we use preAlignInsertions.
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toPreAlign
            ];
          path_likelihood *=
            getPreAlignEmissionProbability(
              profile,
              sequence,
              col_i - 1
            );
          // In the top row we use preAlignInsertions, which are not
          // affine, so we store the values in the Match cell instead of
          // the Insertion cell.
          tmp_mvt = forward_row[ col_i - 1 ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          path_likelihood *= rabiner_within_row_correction_factor;

          global_entente_update[
            Transition::fromPreAlign
          ][
            TransitionFromPreAlign::toPreAlign
          ] +=
            path_likelihood;

            global_entente_update[ Emission::PreAlignInsertion ][ residue ] +=
              path_likelihood;
        } else { // col_i > 0 (or last_col == col_i == 0), row_i > 0
          if( row_i == last_row ) {
            if( col_i == last_col ) {
              // Note in the last row the sequence score is calculated entirely
              // using the last column.
              // Include in the entente the C->T transitions, which are in
              // every path so it's just the total sequence score

              //// This is how we used to do it, but now that we have the
              //// rabinerScaling, we need the *scaled* C->T, which is just the
              //// backward value there.
              //matrix_value_type_score =
                //(
                // profile[
                //   Transition::fromPostAlign
                // ][
                //   TransitionFromPostAlign::toTerminal
                // ] *
                // forward_row[ col_i ][ Match ]
                //);
              // So here is the new way:
              path_likelihood =
                backward_row[ col_i ][ Match ];
              path_likelihood /=
                parameters.matrixRowScaleFactor;
              tmp_mvt = forward_row[ col_i ][ Match ];
              tmp_mvt /= parameters.matrixRowScaleFactor;
              path_likelihood *= tmp_mvt;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS // & USE_DEL_IN_DEL_OUT
              tmp_mvt = forward_row.m_deletionOut;
#else // !DISALLOW_FLANKING_TRANSITIONS
              tmp_mvt = forward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS ( & USE_DEL_IN_DEL_OUT )
              tmp_mvt /=
                parameters.matrixRowScaleFactor;
              // Note we don't use the backward_row deletion out value here
              // because it includes the delout end and end->postalign values,
              // which are also included in the forward_row value.  All that we
              // need is the postalign->terminal part.
              tmp_mvt *=
                profile[
                  Transition::fromPostAlign
                ][
                  TransitionFromPostAlign::toTerminal
                ];

              path_likelihood += tmp_mvt;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

              path_likelihood *= rabiner_within_row_correction_factor;

              if( calculate_score ) {
                matrix_value_type_score = path_likelihood;
              } // End if calculate_score

              global_entente_update[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toTerminal
              ] +=
                path_likelihood;
              //cout << "\tafter: " << global_entente_update[
              //  Transition::fromPostAlign
              //][
              //  TransitionFromPostAlign::toTerminal
              //] << endl;
            } // End if it's the bottom-right cell.
            // In the bottom row, we use postAlignInsertions.  So this is
            // really I->I (postalign, so technically C->C), stored in M
            // subcell.
            path_likelihood =
              backward_row[ col_i ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ];
            path_likelihood *=
              getPostAlignEmissionProbability(
                profile,
                sequence,
                col_i - 1
              );
            // In the bottom row we use postAlignInsertions, which are not
            // affine, so we store the values in the Match cell instead of
            // the Insertion cell.
            tmp_mvt = forward_row[ col_i - 1 ][ Match ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
            tmp_mvt +=
              forward_row[ col_i - 1 ][ DeletionOut ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT && !DISALLOW_FLANKING_TRANSITIONS
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;
            path_likelihood *= rabiner_within_row_correction_factor;

            global_entente_update[
              Transition::fromPostAlign
            ][
              TransitionFromPostAlign::toPostAlign
            ] += path_likelihood;
            // TODO: REMOVE
            //cout << "\tafter: " << global_entente_update[
            //  Transition::fromPostAlign
            //][
            //  TransitionFromPostAlign::toPostAlign
            //] << endl;
            //cout << "\tafter (more): " << global_entente_update[
            //  Transition::fromPostAlign
            //] << endl;
            //cout << "\tafter (whole shebang): " << global_entente_update << endl;
            // TODO: REMOVE!
            //exit( 0 );

              global_entente_update[ Emission::PostAlignInsertion ][ residue ] +=
                path_likelihood;

            // Note in the last row the sequence score is calculated entirely
            // using the last column.
          } // End if it's the last row..
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;

          if( row_i == 1 ) {
            // TODO: REMOVE
            //cout << "N->B prob is " <<
            //  profile[
            //    Transition::fromPreAlign
            //  ][
            //    TransitionFromPreAlign::toBegin
            //  ] << endl;
            path_likelihood *=
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            path_likelihood *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toMatch
              ];
          } else { // if row_i == 1 .. else ..
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              prev_row_match_distribution[
                TransitionFromMatch::toMatch
              ];
#else
            path_likelihood *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toMatch
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          } // if row_i == 1 .. else ..
          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Match ][ residue ];
#ifdef USE_END_DISTRIBUTION
            if( row_i == last_row ) {
              path_likelihood *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i - 1 ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          if( row_i == 1 ) {
            // TODO: REMOVE
            //cout << "adding " << path_likelihood << " to the N->B value of the GEU." << endl;
            global_entente_update[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ] += path_likelihood;
            global_entente_update[
              Transition::fromBegin
            ][
              TransitionFromBegin::toMatch
            ] += path_likelihood;
          } else {
            global_entente_update[
              Transition::fromMatch
            ][
              TransitionFromMatch::toMatch
            ] += path_likelihood;
          } // end if row_i == 1 .. else ..
          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            global_entente_update[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else {
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *2*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score

          } // End if row_i == last_row .. else ..
        } // End if col_i == 0, elsif row_i == 0, else ..
        // from Match to Insertion
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != last_row ) ) {
          path_likelihood =
            backward_row[ col_i ][ Insertion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          path_likelihood *=
            row_match_distribution[
              TransitionFromMatch::toInsertion
            ];
#else
          path_likelihood *=
            profile[ row_i - 1 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Insertion ][ residue ];
          tmp_mvt = forward_row[ col_i - 1 ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          path_likelihood *= rabiner_within_row_correction_factor;

          global_entente_update[
            Transition::fromMatch
          ][
            TransitionFromMatch::toInsertion
          ] += path_likelihood;
            // TODO: REMOVE
            //cout << "updateGEFS(..): *1* about to add " << path_likelihood << " to the insertion emission distribution for residue " << residue << endl;
            global_entente_update[ Emission::Insertion ][ residue ] +=
              path_likelihood;
          // Don't include this path_likelihood in the matrix_value_type_score,
          // since this is a transition to an Insertion state, so the path will
          // be counted elsewhere.
        } // End if this is not the first column and not the first or last row ..
        // from Match to Deletion
        if( row_i != 0 ) {
          if( row_i == 1 ) {
            path_likelihood =
              backward_row[ col_i ][ Deletion ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            path_likelihood *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletion
              ];
            tmp_mvt = prev_forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            global_entente_update[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ] += path_likelihood;
            global_entente_update[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletion
            ] += path_likelihood;
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *3*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score

          } else if( row_i == last_row ) {
            // We store the postAlignInsertion stuff in the Match state.
            path_likelihood =
              backward_row[ col_i ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              prev_row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            path_likelihood *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
#ifdef USE_END_DISTRIBUTION
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            tmp_mvt = prev_forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            global_entente_update[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ] += path_likelihood;
#ifdef USE_END_DISTRIBUTION
            global_entente_update[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else { // row_i is not 0, 1, or last_row
            path_likelihood =
              backward_row[ col_i ][ Deletion ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              prev_row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            path_likelihood *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            tmp_mvt = prev_forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            global_entente_update[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ] += path_likelihood;
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *4*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score

          } // End if it's row 1 .. elsif the last row .. else ..
        } // End if it's not the first row ..
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        // from Match to DeletionOut
        if(
          ( row_i > 1 ) && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS // & USE_DEL_IN_DEL_OUT
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS ( & USE_DEL_IN_DEL_OUT )
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS // & USE_DEL_IN_DEL_OUT
         path_likelihood =
           backward_row.m_deletionOut;
#else
         path_likelihood =
           backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS ( & USE_DEL_IN_DEL_OUT ) .. else ..
         assert( path_likelihood != 0 );
         path_likelihood /=
           parameters.matrixRowScaleFactor;

         // Note that we do not use the scaled distribution, since the
         // backwards row incorporates the rest of the del-out path.  We just
         // want the deletion OPEN probability.
         path_likelihood *=
           profile[
             Transition::fromMatch
           ][
             TransitionFromMatch::toDeletionOut
           ];
         // TODO: REMOVE (Need to account for machine epsilon to use this)
         //assert( path_likelihood == prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] );
//#ifndef NDEBUG
//         if( path_likelihood != prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] ) {
//           bool path_likelihood_gt_prmd = ( static_cast<ScaledMatchDistributionProbabilityType>( path_likelihood ) > prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] );
//           ScaledMatchDistributionProbabilityType diff = ( path_likelihood_gt_prmd ? ( static_cast<ScaledMatchDistributionProbabilityType>( path_likelihood ) - prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] ) : ( prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ]  - static_cast<ScaledMatchDistributionProbabilityType>( path_likelihood ) ) );
//           if( ( diff / ( path_likelihood_gt_prmd ? static_cast<ScaledMatchDistributionProbabilityType>( path_likelihood ) : prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] ) ) > 1E-5 ) { // TODO: DEHACKIFY MAGIC # MACHINE EPSILON!
//             cout << "UH OH: path_likelihood is " << path_likelihood << " but prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] is " << prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] << ".  Their diff is ( path_likelihood - prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] ) = " << ( path_likelihood_gt_prmd ? ' ' : '-' ) << diff << ", and abs(diff) / (the greater value) is " << ( diff / ( path_likelihood_gt_prmd ? static_cast<ScaledMatchDistributionProbabilityType>( path_likelihood ) : prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] ) ) << "." << endl;
//             assert(  toDouble( ( diff / ( path_likelihood_gt_prmd ? static_cast<ScaledMatchDistributionProbabilityType>( path_likelihood ) : prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ] ) ) < 1E-5 ) ); // TODO: DEHACIFY MAGIC # MACHINE EPSILON!
//           }
//         }
//#endif // !NDEBUG
          // TODO: REMOVE!!! TESTING!!! ERE I AM!!!
         //path_likelihood = prev_row_match_distribution[ TransitionFromMatch::toDeletionOut ];

          tmp_mvt = prev_forward_row[ col_i ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          // TODO: REMOVE
          if( false && ( row_i > 2 ) && ( path_likelihood == static_cast<MatrixValueType>( 0 ) ) ) {
            cout << "path_likelihood is now 0.  row_i is " << row_i << ", col_i is " << col_i << "; backward_row.m_deletionOut is " << backward_row.m_deletionOut << "; prev-fwd-row-match at this column is " << prev_forward_row[ col_i ][ Match ] << endl;
            //assert( ( path_likelihood != 0 ) );
            assert( ( path_likelihood != static_cast<MatrixValueType>( 0 ) ) );
          }

          // Note we don't need to multiply in the DeletionOut::toEnd or the
          // End::toPostAlign probabilities, because these are already
          // incorporated into the backward_row.  But we do need to count them.
          if( row_i == last_row ) {
            // TODO: REMOVE
            //cout << "For M->W INTO last_row, adding " << path_likelihood << " to the W->E GEU value." << endl;
            global_entente_update[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toEnd
            ] += path_likelihood;
#ifdef USE_END_DISTRIBUTION // & USE_DEL_IN_DEL_OUT
            // If this is the last row, that transition includes from End to
            // PostAlign.
            global_entente_update[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION ( & USE_DEL_IN_DEL_OUT )
          } // End if row_i == last_row

          if( calculate_score && ( row_i != last_row ) ) {
            matrix_value_type_score += path_likelihood;
          } // End if calculate_score
  
          // TODO: PUT BACK.. TESTING.
          global_entente_update[
            Transition::fromMatch
          ][
            TransitionFromMatch::toDeletionOut
          ] += path_likelihood;

          // TODO: REMOVE
          //cout << "for row_i == " << row_i << ", contribution from match to deletionOut is " << path_likelihood << endl;

          // TODO: REMOVE.. TESTING
          // TODO: Especially remove this, which makes NO sense!
          //path_likelihood /=
          //  profile[
          //    Transition::fromMatch
          //  ][
          //    TransitionFromMatch::toDeletionOut
          //  ];

          // Divide out the extensions and the close.
          //( ( ( row_i - 2 ) + 2 ) < profile.length() ) {
          //// then there's del-out extensions to consider (from the prev pos)
          //tmp_mvt =
          //  pow( profile[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ], ( ProbabilityType )( profile.length() - ( ( row_i - 2 ) + 2 ) ) );
          //path_likelihood /=
          //  tmp_mvt;
          //// End if there's extensions
          //th_likelihood /=
          //profile[
          //  Transition::fromDeletionOut
          //][
          //  TransitionFromDeletionOut::toEnd
          //];
#ifdef USE_END_DISTRIBUTION // & USE_DEL_IN_DEL_OUT
          //path_likelihood /=
          //  profile[
          //    Transition::fromEnd
          //  ][
          //    TransitionFromEnd::toPostAlign
          //  ];
#endif // USE_END_DISTRIBUTION ( & USE_DEL_IN_DEL_OUT )
          //path_likelihood /=
          //  profile[
          //    Transition::fromPostAlign
          //  ][
          //    TransitionFromPostAlign::toTerminal
          //  ];
          // TODO: REMOVE
          //cout << "for row_i == " << row_i << ", contribution from match to deletionOut is " << path_likelihood << endl;
          //global_entente_update[
          //  Transition::fromMatch
          //][
          //  TransitionFromMatch::toDeletionOut
          //] += path_likelihood;
      } // End if ( row_i > 1 ) && ( col_i > 0 ) (and maybe col_i == last_col)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        // End from Match

        // from Insertion
        // from Insertion to Match
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != 1 ) ) { // No ins state in row 0
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toMatch
            ];

            path_likelihood *=
              profile[ row_i - 1 ][ Emission::Match ][ residue ];
#ifdef USE_END_DISTRIBUTION
            if( row_i == last_row ) {
              path_likelihood *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i - 1 ][ Insertion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          global_entente_update[
            Transition::fromInsertion
          ][
            TransitionFromInsertion::toMatch
          ] += path_likelihood;
          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            global_entente_update[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else {
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *5*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score

          }
        } // End if it's not the first or second row and not the first col ..
        // from Insertion to Insertion
        if( ( row_i != 0 ) &&
            ( row_i != last_row ) &&
            ( col_i != 0 ) ) {
          path_likelihood =
            backward_row[ col_i ][ Insertion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 1 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toInsertion
            ];
            path_likelihood *=
              profile[ row_i - 1 ][ Emission::Insertion ][ residue ];
          tmp_mvt = forward_row[ col_i - 1 ][ Insertion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          path_likelihood *= rabiner_within_row_correction_factor;

          global_entente_update[
            Transition::fromInsertion
          ][
            TransitionFromInsertion::toInsertion
          ] += path_likelihood;
            // TODO: REMOVE
            //cout << "updateGEFS(..): *2* about to add " << path_likelihood << " to the insertion emission distribution for residue " << residue << endl;
            global_entente_update[ Emission::Insertion ][ residue ] +=
              path_likelihood;
          // Don't include this path in the matrix_value_type_score, since this
          // is a transition to an Insertion..
        } // End if this is not the first or last row, and not the first col..
        // End from insertion

        // from Deletion
        // from Deletion to Match
        if( row_i == last_row ) {
          // We store post-align stuff in the Match state, even if it comes
          // from a deletion.  So this is really a D->D, stored at M.
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
#ifdef USE_END_DISTRIBUTION
          path_likelihood *=
            profile[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ];
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i ][ Deletion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          global_entente_update[
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toDeletion
          ] += path_likelihood;
#ifdef USE_END_DISTRIBUTION
          global_entente_update[
            Transition::fromEnd
          ][
            TransitionFromEnd::toPostAlign
          ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
        }
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != 1 ) ) { // no del state in row 0
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toMatch
            ];

            path_likelihood *=
              profile[ row_i - 1 ][ Emission::Match ][ residue ];
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i - 1 ][ Deletion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          global_entente_update[
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toMatch
          ] += path_likelihood;
          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            global_entente_update[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else { // row_i < last_row
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *6*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score
          } // End if row_i == last_row .. else ..

          // D'oh!  Why was this here?!
          // TODO: REMOVE!
          //global_entente_update[
          //  Transition::fromDeletion
          //][
          //  TransitionFromDeletion::toDeletion
          //] += path_likelihood;

        } // End if it's not the first or second row and not the first col..
        // from Deletion to Deletion
        if( ( row_i != 0 ) &&
            ( row_i != 1 ) && // no del state in first row
            ( row_i != last_row ) ) {
          path_likelihood =
            backward_row[ col_i ][ Deletion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
          tmp_mvt = prev_forward_row[ col_i ][ Deletion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          global_entente_update[
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toDeletion
          ] += path_likelihood;

          if( calculate_score ) {
            matrix_value_type_score += path_likelihood;
          }
        } // End if it's not the first, second, or last row ..
        // End from Deletion
      } // End for each column, ..

#ifdef DISALLOW_FLANKING_TRANSITIONS
      // TODO: ERE I AM.  Why would this make sense?  First of all, wouldn't these be 0 anyway?  Secondly, aren't we enforcing the DISALLOW_FLANKING_TRANSITIONS by constraining the profile?  Shouldn't these be asserts, changed to row_i == 0 ?  See below.
//     if( row_i != 0 ) {
//       global_entente_update[
//         Transition::fromPreAlign
//       ].zero();
//     }
//     if( row_i != last_row ) {
//       global_entente_update[
//         Transition::fromPostAlign
//       ].zero();
//     }
#endif // DISALLOW_FLANKING_TRANSITIONS

      // TODO: REMOVE
      //cout << "updateGEFS: global_entente_update (before scaling) is " << global_entente_update << endl;

      // TODO: REMOVE
      //cout << "updateGEFS: matrix_value_type_score is:" << matrix_value_type_score << endl;

      if( inverse_scalar == NULL ) {
        if( parameters.useRabinerScaling ) {
          // It's already scaled by the cumulative_inverse_scalar, but that
          // might not be the same thing as the score..
          global_entente_update.m_scalar *= matrix_value_type_score;
          global_entente_update.m_scalar /= cumulative_inverse_scalar;
          //cout << "updateGEFS: NOT scaling by score calculated from rabiner scalars:" << cumulative_inverse_scalar << endl;
        } else {
          //cout << "updateGEFS: scaling by score:" << matrix_value_type_score << endl;
          global_entente_update.m_scalar *= matrix_value_type_score;
        } // End if useRabinerScaling .. else ..
        // TODO: REMOVE
        //cout << "updateGEFS: global_entente_update (after scaling by the inverse of the sequence score) is " << global_entente_update << endl;
      } else if( *inverse_scalar != 0.0 ) {
        //cout << "updateGEFS: scaling by *inverse_scalar:" << *inverse_scalar << endl;
        if( parameters.useRabinerScaling ) {
          // It's already been scaled, but by the cumulative_inverse_scalar.
          global_entente_update.m_scalar *= *inverse_scalar;
          global_entente_update.m_scalar /= cumulative_inverse_scalar;
        } else {
          global_entente_update.m_scalar *= *inverse_scalar;
        } // End if useRabinerScaling .. else ..
        // TODO: REMOVE
        //cout << "updateGEFS: global_entente_update (after scaling by the inverse of the argument value) is " << global_entente_update << endl;
      }
      // TODO: REMOVE
      //if( isnan( global_entente_update[
      //        Transition::fromMatch
      //      ][
      //        TransitionFromMatch::toMatch
      //        ] ) ) {
      //  cout << "Uh-oh, updateGEFS(..): got nan!" << endl;
      //  cout << "\t row_i is " << row_i << endl;
      //  cout << "\t last_row is " << last_row << endl;
      //  cout << "\t last_col is " << last_col << endl;
      //  if( inverse_scalar == NULL ) {
      //    cout << "\t inverse_scalar is NULL" << endl;
      //    cout << "\t matrix_value_type_score is " << matrix_value_type_score << endl;
      //  } else {
      //    cout << "\t inverse_scalar is " << *inverse_scalar << endl;
      //    if( *inverse_scalar == 0.0 ) {
      //      cout << "\t cumulative_inverse_scalar = " << cumulative_inverse_scalar << endl;
      //    }
      //  }
      //  cout << "\t global_entente.m_scalar = " << global_entente.m_scalar << endl;
      //  cout << "\t global_entente_update is " << global_entente_update << endl;
      //} // end if isnan..

      // Now do the updating, scaling the global_entente appropriately.
      // NEW:
      //  cout << "\t global_entente_update is " << global_entente_update << endl;
#ifdef DISALLOW_FLANKING_TRANSITIONS
      if( row_i == 0 ) {
#ifndef NDEBUG
        if( global_entente_update.maximumValue() != 0 ) {
          cout << endl << "Unexpectedly > 0: global_entente_update for row 0 = " << global_entente_update << endl;
          cout << "Profile globals are ";
          profile.writeExceptPositions( cout );
          cout << endl;
        }
        assert( global_entente_update.maximumValue() == 0 );
#endif // !NDEBUG
      } else {
#endif // DISALLOW_FLANKING_TRANSITIONS
#ifdef DISALLOW_FLANKING_TRANSITIONS
        // NOTE: See above; the row_i isn't 0.
        //if( row_i == 0 ) {
        //  assert( max( global_entente_update.maximumValue(), global_entente_update[ Transition::fromPreAlign ].maximumValue() ) != 0 );
        //} else
        if( row_i == last_row ) {
          assert( max( global_entente_update.maximumValue(), global_entente_update[ Transition::fromPostAlign ].maximumValue() ) != 0 );
        } else {
          assert( global_entente_update.maximumValue() != 0 );
        }
#else
        assert( global_entente_update.maximumValue() != 0 );
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..

        // TODO: REMOVE
        //cout << "global_entente_update is " << global_entente_update.createUnscaledCopy() << endl;

        global_entente += global_entente_update;
        
        // NOTE: See alternative below
#ifdef DISALLOW_FLANKING_TRANSITIONS
        global_entente[ Transition::fromPreAlign ] +=
         global_entente_update[ Transition::fromPreAlign ];
        global_entente[ Transition::fromPostAlign ] +=
          global_entente_update[ Transition::fromPostAlign ];
        // NOTE: I've commented this out because it was helping me debug what I thought was a problem: the C->T shows up as the score, not as 1 as would happen without DISALLOW_FLANKING_TRANSITIONS -- but in the end, the C->T is irrelevant when DISALLOW_FLANKING_TRANSITIONS, so other than reassuring me that everything is kosher, the following is not necessary, and it can be confusing because the unscale() doesn't apply to C->T either, so after unscaling, things look funky no matter what.
//        // NOTE That += doesn't have an effect on the pre-align and post-align transitions when DISALLOW_FLANKING_TRANSITIONS is #defined.
//        // NOTE HACK:  these are not scaled additions, so we manually do the scaling !!!!
//        // TODO: Redesign so this isn't necessary.  eg. make += add everything...
//        ScoreType scale_ratio_st = global_entente.m_scalar;
//        scale_ratio_st /= global_entente_update.m_scalar;
//        global_entente_update[ Transition::fromPreAlign ] *= scale_ratio_st;
//        global_entente[ Transition::fromPreAlign ] +=
//          global_entente_update[ Transition::fromPreAlign ];
//        global_entente_update[ Transition::fromPostAlign ] *= scale_ratio_st;
//        global_entente[ Transition::fromPostAlign ] +=
//          global_entente_update[ Transition::fromPostAlign ];
#endif // DISALLOW_FLANKING_TRANSITIONS
#ifdef DISALLOW_FLANKING_TRANSITIONS
      }
#endif // DISALLOW_FLANKING_TRANSITIONS



      // OLD:
      //if( !( ( inverse_scalar != NULL ) && ( *inverse_scalar == 0.0 ) ) ) {
      //
      //  // Incorporate the updated entente, accounting for scaling.
      //
      //  // TODO: REMOVE
      //  //cout << "updateGEFS(..): global_entente before is " << global_entente << endl;
      //  //cout << "updateGEFS(..): global_entente_update is " << global_entente_update << endl;
      //
      //  // First account for the different scales
      //  ScoreType scale_ratio_st = global_entente.m_scalar;
      //  scale_ratio_st /= global_entente_update.m_scalar;
      //  MatrixValueType scale_ratio = scale_ratio_st;
      //
      //  // TODO: REMOVE
      //  //cout << "updateGEFS(..): scale_ratio is " << scale_ratio << endl;
      //  if( scale_ratio != 0 ) {
      //    global_entente_update *= scale_ratio;
      //    global_entente += global_entente_update;
      //    // Make max value 1.0, updating m_scalar so that it will unscale to
      //    // the correct values.
      //    global_entente.rescale();
      //  }
      //  // TODO: REMOVE
      //  //cout << "updateGEFS(..): global_entente after is " << global_entente << "; m_scalar is " << global_entente.m_scalar << "." << endl;
      //  //global_entente_update = global_entente;
      //  //global_entente_update /= global_entente.m_scalar;
      //  //cout << "\t unscaled: " << global_entente_update << endl;
      //}

      return matrix_value_type_score;
    } // updateGlobalEntenteForSequence( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, Row const&, Row const&, GlobalEntente & ) const


    /**
     * Given the forward matrices for a set of Sequences, calculate their
     * AlignmentProfiles, and place the result in the given vector ref.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
    template <typename ProfileType, typename SequenceType>
  GALOSH_INLINE_TRIVIAL
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    calculateAlignmentProfiles (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<SequenceType> const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer const& forward_matrices,
      vector<AlignmentProfile> & alignment_profiles
    ) const
    {
      sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );

      // Allocate some temporary backward matrices...
      typename Matrix::RowVector backward_rows_1 =
        typename Matrix::RowVector(
          sequences,
          sequence_count
        );
      typename Matrix::RowVector backward_rows_2 =
        typename Matrix::RowVector(
          sequences,
          sequence_count
        );
      calculateAlignmentProfiles(
        parameters,
        profile,
        sequences,
        sequence_count,
        NULL,
        forward_matrices,
        backward_rows_1,
        backward_rows_2,
        alignment_profiles
      );
    } // calculateAlignmentProfiles ( Parameters const&, ProfileType const&, vector<SequenceType> const&, uint32_t const&, SequentialAccessContainer const&, vector<AlignmentProfile> & ) const

    /**
     * Given the forward matrices for a set of Sequences, calculate their
     * AlignmentProfiles, and place the result in the given vector ref.  The
     * given RowVectors are for temporary storage -- they will be clobbered.
     * If the scores are not provided, and if useRabinerScaling is true and
     * rabinerScaling_useMaximumValue is also true, then they will be
     * calculated (each is the cumulativeInverseScalar of the last row of the
     * forward matrix) -- otherwise they are not needed.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
    template <typename ProfileType, typename SequenceType>
  GALOSH_INLINE_ALGORITHM
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    calculateAlignmentProfiles (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<SequenceType> const& sequences,
      uint32_t const& sequence_count,
      vector<ScoreType> const * sequence_scores,
      typename Matrix::SequentialAccessContainer const& forward_matrices,
      typename Matrix::RowVector & backward_rows_1,
      typename Matrix::RowVector & backward_rows_2,
      vector<AlignmentProfile> & alignment_profiles
    ) const
    {
      typename Matrix::RowVector * backward_rows_ptr = &backward_rows_1;
      typename Matrix::RowVector * next_backward_rows_ptr = &backward_rows_2;
      typename Matrix::RowVector * temp_backward_rows_ptr;

      uint32_t last_seq =
        ( ( sequence_count == 0 ) ? ( sequences.size() - 1 ) : ( min( static_cast<size_t>( sequence_count ), sequences.size() ) - 1 ) );
      // TODO: REMOVE
      //cout << "sequence_count is " << sequence_count << "; last_seq is " << last_seq << endl;

      bool we_are_calculating_scores = false;
      vector<ScoreType> local_sequence_scores;
      if( parameters.useRabinerScaling && parameters.rabinerScaling_useMaximumValue ) {
        // We need to pass scores to the calculateAlignmentProfilePosition(..)
        // method.
        if( sequence_scores == 0 ) {
          we_are_calculating_scores = true;
          local_sequence_scores.resize( sequences.size() );
          sequence_scores = &local_sequence_scores;
        }
      } // End if we need the scores

      uint32_t seq_i;
      typename Matrix::SequentialAccessContainer::const_reverse_iterator forward_matrices_rev_iter = forward_matrices.rbegin();
      typename Matrix::SequentialAccessContainer::const_reverse_iterator prev_forward_matrices_rev_iter = forward_matrices_rev_iter;
      ++prev_forward_matrices_rev_iter;

      uint32_t last_row = profile.length();
      uint32_t row_i = last_row;
      do { // While row_i >= 0
        // Rotate the backward_rows
        temp_backward_rows_ptr = backward_rows_ptr;
        backward_rows_ptr = next_backward_rows_ptr;
        next_backward_rows_ptr = temp_backward_rows_ptr;

        // TODO: REMOVE
        //cout << "row_i is " << row_i << endl;

        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {

          // TODO: REMOVE
          //cout << "seq_i is " << seq_i << endl;

          if( ( row_i == last_row ) && we_are_calculating_scores ) {
            local_sequence_scores[ seq_i ] =
              ( *forward_matrices_rev_iter )[ seq_i ][ sequences[ seq_i ].length() ][ Match ];
            local_sequence_scores[ seq_i ] *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toTerminal
              ];
            local_sequence_scores[ seq_i ] /=
              parameters.matrixRowScaleFactor;
            if( parameters.useRabinerScaling ) {
              local_sequence_scores[ seq_i ] *=
                ( *forward_matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
            }
            // TODO: REMOVE
            cout << "in calculateAlignmentProfiles(..), sequence " << seq_i << " score is " << ( *sequence_scores )[ seq_i ] << endl;
          } // End if we need to calc the score

          // Calculate the backward row for seq_i
          if( parameters.useRabinerScaling ) {
            ( *backward_rows_ptr )[ seq_i ].m_rabinerInverseScalar =
              ( *forward_matrices_rev_iter )[ seq_i ].m_rabinerInverseScalar;
          } // End if useRabinerScaling
          // TODO: REMOVE
          //cout << "Calling backward_calculateRow" << endl;
          backward_calculateRow(
             parameters,
             profile,
             sequences[ seq_i ],
             row_i,
             ( *next_backward_rows_ptr )[ seq_i ],
             ( *backward_rows_ptr )[ seq_i ]
           );
          if( parameters.useRabinerScaling ) {
            if( row_i == last_row ) {
              ( *backward_rows_ptr )[ seq_i ].m_rabinerCumulativeInverseScalar =
                ( *backward_rows_ptr )[ seq_i ].m_rabinerInverseScalar;
            } else {
              ( *backward_rows_ptr )[ seq_i ].m_rabinerCumulativeInverseScalar =
                (
                 ( *next_backward_rows_ptr )[ seq_i ].m_rabinerCumulativeInverseScalar *
                 ( *backward_rows_ptr )[ seq_i ].m_rabinerInverseScalar
                );
            } // End if row_i == last_row .. else ..
          } // End if useRabinerScaling

          // TODO: REMOVE
          //cout << "Calling calculateAlignmentProfilePosition" << endl;

          calculateAlignmentProfilePosition(
            parameters,
            profile,
            sequences[ seq_i ],
            row_i,
            ( ( row_i == 0 ) ?
              ( *prev_forward_matrices_rev_iter )[ seq_i ] /*ignored*/ :
              ( *prev_forward_matrices_rev_iter )[ seq_i ] ),
            ( *forward_matrices_rev_iter )[ seq_i ],
            backward_rows_ptr->operator[](  seq_i ),
            next_backward_rows_ptr->operator[](  seq_i ),
            ( ( sequence_scores == 0 ) ? NULL : &( *sequence_scores )[ seq_i ] ),
            alignment_profiles[ seq_i ][ row_i ],
            NULL
          );
        } // End foreach seq_i

        //prev_forward_matrices_rev_iter = forward_matrices_rev_iter;
        forward_matrices_rev_iter = prev_forward_matrices_rev_iter;
        // Increment matrices iterator, which is going in reverse through the
        // RowVectors
        ++prev_forward_matrices_rev_iter;
      } while( row_i-- > 0 );
    } // calculateAlignmentProfiles ( Parameters const&, ProfileType const&, vector<SequenceType> const&, uint32_t const&, vector<ScoreType> const *, SequentialAccessContainer const&, Matrix<RowVector> &, Matrix<RowVector> &, vector<AlignmentProfile> & ) const

    /**
     * Given the forward and backward rows for row row_i, and the forward row
     * for ( row_i - 1 ) and the backward row for ( row_i + 1 ), calculate the
     * expected uses of each profile parameter at profile position ( row_i - 1
     * ) for the given sequence, profile pair.  Results will be divided by the
     * sequence score or by the given inverse_scalar (if it is
     * non-null). Results go in the given AlignmentProfilePosition reference.
     * If the coefficients pointer is non-NULL, the coefficients will be
     * calculated too (as if you'd called
     * calculatePositionSpecificSequenceScoreCoefficients(..).
     *
     * NOTE: If useRabinerScaling is true but rabinerScaling_useMaximumValue is
     * also true, then you *MUST* supply a value for inverse_scalar.  We can't
     * calculate it here.  Likewise when the backward_row is not up-to-date, as
     * happens in ProfileTrainer.hpp.
     * 
     * NOTE: If row_i is 0, prev_forward_row is ignored, and coefficients will
     * be all 0.  If row_i == last_row, next_backward_row is ignored.
     */
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    calculateAlignmentProfilePosition (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row const& forward_row,
      typename Matrix::Row const& backward_row,
      typename Matrix::Row const& next_backward_row,
      ScoreType const * inverse_scalar,
      AlignmentProfilePosition & pos,
      PositionSpecificSequenceScoreCoefficients * coefficients
    ) const
    {
      pos.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
      // Note that global_entente_update.zero() won't zero the fromPreAlign::toPreAlign when DISALLOW_FLANKING_TRANSITIONS is set, so:
      pos[ Transition::fromPreAlign ].zero();
      pos[ Transition::fromPostAlign ].zero();

      // TODO: REMOVE
      assert( profile[ Transition::fromPreAlign ][ TransitionFromPreAlign::toPreAlign ] == 0 );
      assert( profile[ Transition::fromPostAlign ][ TransitionFromPostAlign::toPostAlign ] == 0 );
#endif // DISALLOW_FLANKING_TRANSITIONS

      // NOTE: This is rather complicated because Match emissions for prof pos
      // (row_i-1) are calculated using transitions INTO row_i, while the
      // transitions out of prof pos (row_i-1) are calculated using transitions
      // OUT OF row_i.

      static const bool do_extra_debugging = false;

      // TODO: REMOVE
      //cout << "hi 0" << endl;

      // TODO: REMOVE
      //cout << "hi 1" << endl;

      const uint32_t last_row = profile.length();
      //cout << "hi 1.2" << endl;
      const uint32_t last_col = sequence.length();
      //cout << "hi 1.7" << endl;
      // TODO: Add support for 0-length sequences!  We should be learning
      // something about the deletion rate, at least.
      assert( last_col > 0 );

      SequenceResidueType residue;
      SequenceResidueType next_residue;

      register uint32_t col_i;

      MatrixValueType path_likelihood;
      MatrixValueType matrix_value_type_score = 0;
      MatrixValueType tmp_mvt, tmp_mvt2;

      // Note that when we are using rabiner scaling, if
      // rabinerScaling_useMaximumValue is false then we do not need to
      // recalculate the score -- we just use the cumulative_inverse_scalars.
      // The score as calculated here ought to be 1 if useRabinerScaling is
      // true and rabinerScaling_useMaximumValue is false, so it can be used as
      // a sanity check.

      // If useRabinerScaling is true but rabinerScaling_useMaximumValue is
      // also true, then you *MUST* supply a value for inverse_scalar.  We
      // can't calculate it here.  Likewise if the backward_row is not
      // up-to-date (as happens in ProfileTrainer.hpp).

      // TODO: Couldn't we calculate the score using a method similar to the
      // forward_score(..) method that takes the forward and backward rows?

      assert( ( parameters.useRabinerScaling && parameters.rabinerScaling_useMaximumValue ) ? ( inverse_scalar != 0 ) : true );
      // TODO: PUT BACK
      bool calculate_score = true;
      //!( parameters.useRabinerScaling && !parameters.rabinerScaling_useMaximumValue );

      // NOTE: When useRabinerScaling is true and
      // rabinerScaling_useMaximumValue is false, then (assuming that the
      // matrix rows are up-to-date, which they might not be), the
      // cumulative_inverse_scalar *is* the score.
      ScoreType cumulative_inverse_scalar =
        backward_row.m_rabinerCumulativeInverseScalar;
      if( row_i > 0 ) {
        cumulative_inverse_scalar *=
          prev_forward_row.m_rabinerCumulativeInverseScalar;
      }

      // When we are using rabiner scaling, ->I transitions are within the row,
      // so are scaled by the current row's forward and backward
      // cumulativeInverseScalars (C_t and D_t), instead of the desired product
      // of the current row's backward cumulativeInverseScalar (D_t) with the
      // previous row's forward scalar (C_{t-1}).  To correct for this, we can
      // multiply them by the ratio C_{t-1}/C_t = 1/c_t, where c_t is the
      // rabinerInverseScalar (not cumulative) for the current row.
      // We will call this the correction_factor:
      MatrixValueType rabiner_within_row_correction_factor = 1.0;
      if( parameters.useRabinerScaling ) {
        rabiner_within_row_correction_factor =
          forward_row.m_rabinerInverseScalar;
        assert( rabiner_within_row_correction_factor > 0 );
      }

#ifndef NDEBUG
      if( rabiner_within_row_correction_factor > 1E10 ) {
        cout << "rabiner_within_row_correction_factor is " << rabiner_within_row_correction_factor << endl;
        cout << "This should be the same as: " <<
        toDouble( ( forward_row.m_rabinerCumulativeInverseScalar * 
                    backward_row.m_rabinerCumulativeInverseScalar ) /
          cumulative_inverse_scalar ) << endl;
        assert( false );
      }
#endif // !NDEBUG

      // TODO: REMOVE
      //if( parameters.useRabinerScaling ) {
      //  cout << "calculateAlignmentProfilePosition: cumulative_inverse_scalar is " << cumulative_inverse_scalar << endl;
      //  if( row_i == 0 ) {
      //    cout << "\t = backward_row.m_rabinerCumulativeInverseScalar." << endl;
      //  } else {
      //    cout << "\t = ( " <<
      //      prev_forward_row.m_rabinerCumulativeInverseScalar << " * " <<
      //      backward_row.m_rabinerCumulativeInverseScalar << " )" << endl;
      //  }
      //}

      // TODO: REMOVE
      //cout << "hi 2" << endl;
//      if( false || ( parameters.debug >= DEBUG_All ) ) {
//        cout << "[debug] calculateAlignmentProfilePosition(..) row_i is " << row_i << "; last_row is " << last_row << ", last_col is " << last_col << endl;
//      }

      // TODO: REMOVE
//      if( false ) {
//        cout << "calculateAlignmentProfilePosition: prev_forward_row is " << prev_forward_row << endl;
//        cout << "calculateAlignmentProfilePosition: forward_row is " << forward_row << endl;
//        cout << "calculateAlignmentProfilePosition: backward_row is " << backward_row << endl;
//        cout << "calculateAlignmentProfilePosition: next_backward_row is " << next_backward_row << endl;
//        //typename Matrix::Row unscaled_row = forward_row;
//        //if( parameters.useRabinerScaling ) {
//        //  unscaled_row *= forward_row.m_rabinerCumulativeInverseScalar;
//        //  cout << "calculateAlignmentProfilePosition: rabiner-unscaled forward_row is " << unscaled_row << endl;  
//        //}
//        //unscaled_row /= parameters.matrixRowScaleFactor;
//        //cout << "calculateAlignmentProfilePosition: unscaled forward_row is " << unscaled_row << endl;  
//        //unscaled_row = backward_row;
//        //if( parameters.useRabinerScaling ) {
//        //  unscaled_row *= backward_row.m_rabinerCumulativeInverseScalar;
//        //  cout << "calculateAlignmentProfilePosition: rabiner-unscaled backward_row is " << unscaled_row << endl;
//        //}
//        //unscaled_row /= parameters.matrixRowScaleFactor;
//        //cout << "calculateAlignmentProfilePosition: unscaled forward_row is " << unscaled_row << endl;  
//      } // End if( false ) (DEBUGGING)

      // TODO: REMOVE
      //cout << "calculateAlignmentProfilePosition(..): sequence is " << sequence << endl;

      // Zero the coefficients, if non-NULL.
      if( coefficients != 0 ) {
        coefficients->zero();
      }

      // TODO: REMOVE
      //if( ( row_i + 1 ) == last_row ) {
      //  cout << "hi from the second-to-last-row." << endl;
      //  cout << "prev_forward_row is " << prev_forward_row << endl;
      //  cout << "forward_row is " << forward_row << endl;
      //  cout << "backward_row is " << backward_row << endl;
      //  cout << "next_backward_row is " << next_backward_row << endl;
      //  cout << "*inverse_scalar is " << ( ( inverse_scalar == 0 ) ? 0 : *inverse_scalar ) << endl;
      //}

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      if( row_i >= 2 ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using prev_row_match_distribution == " << prev_row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if row_i >= 2
      if( row_i >= 1 ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using row_match_distribution == " << row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if row_i >= 1
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      for( col_i = 0; col_i <= last_col; col_i++ ) {

        if( false && ( parameters.debug >= DEBUG_All ) ) {
          cout << "[debug] calculateAlignmentProfilePosition(..): row " << row_i << ", col " << col_i << endl;
        }
        
        // calculate row cell:
        // If not in the leftmost column, then there's base emission(s)
        // from this cell.
          if( col_i > 0 ) {
            residue = sequence[ col_i - 1 ];
//            if( parameters.debug >= DEBUG_All ) {
//              cout << "[debug ] (entente) residue is " << residue << endl;
//            }
          }
          if( col_i < last_col ) {
            next_residue = sequence[ col_i ];
//            if( parameters.debug >= DEBUG_All ) {
//              cout << "[debug ] (entente) next_residue is " << next_residue << endl;
//            }
          }

        ////// INTO row_i //////
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        // to/from DeletionIn
        // from DeletionIn to Match
        if(
          ( row_i > 0 ) && ( col_i >= 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 1 ) 
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toMatch
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_mvt =
            prev_forward_row.m_deletionIn;
#else
          tmp_mvt =
            prev_forward_row[ col_i - 1 ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS
          tmp_mvt /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the Match( " << residue << " ) value for the del-in end." << endl;
          }
          if( coefficients != 0 ) {
            ( *coefficients )[ Emission::Match ][ residue ] +=
              path_likelihood;
          } // End if also computing coefficients

          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Match ][ residue ];

          pos[ Emission::Match ][ residue ] +=
            path_likelihood;
          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            pos[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else {
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
            } // End if calculate_score
          } // End if row_i == last_row .. else ..
        } // End if ( row_i > 0 ) && ( col_i >= 1 ) ( and maybe col_i == 1 )

        // from DeletionIn to DeletionIn (and from Match/Begin to DeletionIn)
        if(
          ( row_i == 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // Actually this is from PreAlign to Begin, from Begin to DeletionIn.
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionIn;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ];
          path_likelihood *=
            profile[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletionIn
            ];
          tmp_mvt =
            prev_forward_row[ col_i ][ Match ];
          tmp_mvt /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-in open." << endl;
          }
          if( coefficients != 0 ) {
            coefficients->m_constant +=
              path_likelihood;
          } // End if we are calculating coefficients

          if( calculate_score ) {
            matrix_value_type_score += path_likelihood;
          } // End if calculate_score
        } else if(
          ( row_i != 0 ) && ( row_i != last_row )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionIn;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toDeletionIn
            ];

#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_mvt =
            prev_forward_row.m_deletionIn;
#else
          tmp_mvt =
            prev_forward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          tmp_mvt /=
            parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-in extensions." << endl;
          }
          if( coefficients != 0 ) {
            coefficients->m_constant +=
              path_likelihood;
          } // End if we are calculating coefficients

          if( calculate_score ) {
            matrix_value_type_score += path_likelihood;
          }
        } // End if ( row_i == 1 ) .. else if( ( row_i != 0 ) && ( row_i != last_row ) ) (and maybe also ( col_i == 0 )).
        // End to/from DeletionIn

        // From DeletionOut
        // from DeletionOut to DeletionOut
        // and from DeletionOut to End
        if(
          ( row_i > 2 )
          && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionOut;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toDeletionOut
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_mvt =
            prev_forward_row.m_deletionOut;
#else
          tmp_mvt =
            prev_forward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          // Note we don't need to multiply in the DeletionOut::toEnd or the
          // End::toPostAlign probabilities, because these are already
          // incorporated into the backward_row.  But we do need to count them.
          if( row_i == last_row ) {
            // TODO: REMOVE
            //cout << "For W->W INTO last_row, adding " << path_likelihood << " to the W->E AlignmentProfile value." << endl;
            pos[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toEnd
            ] += path_likelihood;
#ifdef USE_END_DISTRIBUTION
            // If this is the last row, that transition includes from End to
            // PostAlign.
            pos[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } // End if row_i == last_row

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-out extensions." << endl;
          }

          if( coefficients != 0 ) {
            coefficients->m_constant +=
              path_likelihood;
          } // End if we are calculating coefficients

          if( calculate_score && ( row_i != last_row ) ) {
            matrix_value_type_score += path_likelihood;
          }
        } // End if ( row_i > 2 ) && ( col_i > 0 ) (and maybe col_i == last_col)
        // End from DeletionOut
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        // from Match
        // from Match to Match
        if( ( col_i == 0 ) && ( last_col != 0 ) ) {
          if( row_i == 0 ) {
            if( calculate_score ) {
              // Note in row 0 the score is completely calculated using col 0.
              matrix_value_type_score = // assuming forward row has a 1 here
                backward_row[ 0 ][ Match ];
              matrix_value_type_score /=
                parameters.matrixRowScaleFactor;

              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *0*, matrix_value_type_score is still 0 after setting it to " << backward_row[ 0 ][ Match ] << endl;
              //}
            } // End if calculate_score
          }
          // Do nothing else.
        } else if( row_i == 0 ) { // col_i > 0 && row_i == 0
          // Note in row 0 the score is completely calculated using col 0.

          // In the top row, we use preAlignInsertions.
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toPreAlign
            ];
          path_likelihood *=
            getPreAlignEmissionProbability(
              profile,
              sequence,
              col_i - 1
            );
          // In the top row we use preAlignInsertions, which are not
          // affine, so we store the values in the Match cell instead of
          // the Insertion cell.
          tmp_mvt = forward_row[ col_i - 1 ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          path_likelihood *= rabiner_within_row_correction_factor;

          // TODO: REMOVE
          assert( path_likelihood < 1E10 );

          // In the alignment profile we store "pre-align Opens" in the Match
          // distribution (as Match->Insertion), and "pre-align Extensions" in
          // the pre-align distribution, so we can differentiate.
          if( col_i == 1 ) {
            pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ] += path_likelihood;
              // TODO: REMOVE
              //cout << "calculateAlignmentProfilePosition(..): *1* about to add " << path_likelihood << " to the insertion emission distribution for residue " << residue << endl;
            pos[ Emission::Insertion ][ residue ] +=
              path_likelihood;
          } else {
            pos[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toPreAlign
            ] +=
              path_likelihood;
            pos[ Emission::PreAlignInsertion ][ residue ] +=
              path_likelihood;
          } // End if col_i == 1 .. else ..
        } else { // col_i > 0 (or last_col == col_i == 0), row_i > 0
          if( row_i == last_row ) {
            if( col_i == last_col ) {
              // Note in the last row the sequence score is calculated entirely
              // using the last column.
              path_likelihood =
                backward_row[ col_i ][ Match ];
              path_likelihood /=
                parameters.matrixRowScaleFactor;
              tmp_mvt = forward_row[ col_i ][ Match ];
              tmp_mvt /= parameters.matrixRowScaleFactor;
              path_likelihood *= tmp_mvt;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
              tmp_mvt = forward_row.m_deletionOut;
#else // !DISALLOW_FLANKING_TRANSITIONS
              tmp_mvt = forward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
              tmp_mvt /=
                parameters.matrixRowScaleFactor;
              // Note we don't use the backward_row deletion out value here
              // because it includes the delout end and end->postalign values,
              // which are also included in the forward_row value.  All that we
              // need is the postalign->terminal part.
              tmp_mvt *=
                profile[
                  Transition::fromPostAlign
                ][
                  TransitionFromPostAlign::toTerminal
                ];

              path_likelihood += tmp_mvt;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

              path_likelihood *= rabiner_within_row_correction_factor;
              pos[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toTerminal
              ] +=
                path_likelihood;
              if( calculate_score ) {
                matrix_value_type_score = path_likelihood;
              } // End if calculate_score
            } // End if it's the bottom-right cell.
            // In the bottom row, we use postAlignInsertions.  So this is
            // really I->I (postalign, so technically C->C), stored in M
            // subcell.

            // For alignment profiles, are being tricksy: we need to separately
            // store the "post-align Opens" in the Match->Insertion probability
            // (and the Insertion emission distribution) and the "post-align
            // Extensions" in the post-align distributions.  This is rather
            // tricky since we've lost the distinguishing information by
            // storing everything in the Match state.  So here we need to
            // consider the contributions to the forward matrix row separately.

            // This part is the same, regardless of the source of the Match.
            path_likelihood =
              backward_row[ col_i ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ];
            path_likelihood *=
              getPostAlignEmissionProbability(
                profile,
                sequence,
                col_i - 1
              );
            // And we do *not* need to correct for the post-align transition being within-row..
            // path_likelihood *= rabiner_within_row_correction_factor;

            // In the bottom row we use postAlignInsertions, which are not
            // affine, so we store the values in the Match cell instead of
            // the Insertion cell.
            // See the note above about separating the sources.
            //tmp_mvt = forward_row[ col_i - 1 ][ Match ];
            //tmp_mvt /= parameters.matrixRowScaleFactor;

            // Note that these are all recreating the values contributing to
            // forward_row[ col_i - 1 ][ Match ] -- this is the previous column
            // in the current row.
            tmp_mvt = 0;
            if( col_i > 1 ) {
              // From prev-row Match->Match
              tmp_mvt2 = prev_forward_row[ col_i - 2 ][ Match ];
              tmp_mvt2 /= parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              tmp_mvt2 *=
                prev_row_match_distribution[
                  TransitionFromMatch::toMatch
                ];
#else
              tmp_mvt2 *=
                profile[ row_i - 2 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toMatch
                ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
              tmp_mvt = tmp_mvt2;
              // From prev-row Insertion->Match
              tmp_mvt2 = prev_forward_row[ col_i - 2 ][ Insertion ];
              tmp_mvt2 /= parameters.matrixRowScaleFactor;
              tmp_mvt2 *=
                profile[ row_i - 2 ][
                  Transition::fromInsertion
                ][
                  TransitionFromInsertion::toMatch
              ];
              tmp_mvt += tmp_mvt2;
              // From prev-row Deletion->Match
              tmp_mvt2 = prev_forward_row[ col_i - 2 ][ Deletion ];
              tmp_mvt2 /= parameters.matrixRowScaleFactor;
              tmp_mvt2 *=
                profile[ row_i - 2 ][
                  Transition::fromDeletion
                ][
                  TransitionFromDeletion::toMatch
              ];
              tmp_mvt += tmp_mvt2;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              // From prev-row DeletionIn->Match
              if(
                true
#ifdef DISALLOW_FLANKING_TRANSITIONS
                && ( col_i == 2 )
#endif // DISALLOW_FLANKING_TRANSITIONS
              ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
                tmp_mvt2 = prev_forward_row.m_deletionIn;
#else
                tmp_mvt2 = prev_forward_row[ col_i - 2 ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
                tmp_mvt2 /= parameters.matrixRowScaleFactor;
                tmp_mvt2 *=
                  profile[
                    Transition::fromDeletionIn
                  ][
                    TransitionFromDeletionIn::toMatch
                ];
                tmp_mvt += tmp_mvt2;
              } // End if( true ) (and maybe col_i == 1)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
              // All of those emitted a Match.
              tmp_mvt *=
                getEmissionProbability(
                  profile[ row_i - 1 ],
                  sequence,
                  col_i - 2
                );
            } // End if col_i > 1
            // Note: we can assume col_i > 0

            // From prev-row Match->Deletion
            tmp_mvt2 = prev_forward_row[ col_i - 1 ][ Match ];
            tmp_mvt2 /= parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            tmp_mvt2 *=
              prev_row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            tmp_mvt2 *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            tmp_mvt += tmp_mvt2;
            // From prev-row Deletion->Deletion
            tmp_mvt2 = prev_forward_row[ col_i - 1 ][ Deletion ];
            tmp_mvt2 /= parameters.matrixRowScaleFactor;
            tmp_mvt2 *=
              profile[ row_i - 2 ][
                Transition::fromDeletion
              ][
                TransitionFromDeletion::toDeletion
              ];
            tmp_mvt += tmp_mvt2;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
            // From DeletionOut to Post-Align
            // Note that DeletionOut is the only thing we do store separately
            // from post-aligns in the bottom row.
            tmp_mvt2 = forward_row[ col_i - 1 ][ DeletionOut ];
            tmp_mvt2 /= parameters.matrixRowScaleFactor;
            // In this case we *do* need to multiply in the correction factor,
            // since it doesn't originate in the previous row.
            tmp_mvt2 *= rabiner_within_row_correction_factor;
            tmp_mvt += tmp_mvt2;
#endif //defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )

            path_likelihood *= tmp_mvt;

            // TODO: REMOVE
            assert( path_likelihood < 1E10 );

            // TODO: REMOVE
            //cout << "[row " << row_i << ", col " << col_i << "] C->C OPEN: " << path_likelihood << endl;
            // Store that at the M->I location, and in the insertion emission slot.
            pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ] += path_likelihood;
            pos[ Emission::Insertion ][ residue ] +=
              path_likelihood;

            // Okay, now for the "post-align Extension" part
            // This part is the same, regardless of the source of the Match.
            tmp_mvt2 =
              backward_row[ col_i ][ Match ];
            tmp_mvt2 /=
              parameters.matrixRowScaleFactor;
            tmp_mvt2 *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ];
            tmp_mvt2 *=
              getPostAlignEmissionProbability(
                profile,
                sequence,
                col_i - 1
              );
            // And we need to correct for the post-align transition being within-row..
            tmp_mvt2 *= rabiner_within_row_correction_factor;

            // From current row, prev column Match spot (this is an extension).
            if( col_i > 1 ) {
              tmp_mvt = forward_row[ col_i - 2 ][ Match ];
              tmp_mvt /= parameters.matrixRowScaleFactor;
              tmp_mvt *=
                profile[
                  Transition::fromPostAlign
                ][
                  TransitionFromPostAlign::toPostAlign
                ];
              tmp_mvt *=
                getPostAlignEmissionProbability(
                  profile,
                  sequence,
                  col_i - 2
                );
  
              tmp_mvt2 *= tmp_mvt;

              // TODO: REMOVE
              //cout << "[row " << row_i << ", col " << col_i << "] C->C Extension: " << tmp_mvt2 << endl;

              pos[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ] += tmp_mvt2;
  
              pos[ Emission::PostAlignInsertion ][ residue ] +=
                tmp_mvt2;
            } else { // if col_i > 1 .. else ..
              tmp_mvt2 = 0;
            } // End if col_i > 1 .. else ..
            // Note in the last row the sequence score is calculated entirely
            // using the last column.

            // TODO: REMOVE
            if( false ) {
              // Okay, so it should be the case that path_likelihood + tmp_mvt2
              // is what we whould have calculated.  Let's see.
              tmp_mvt2 += path_likelihood;
              //cout << "total C->C: " << tmp_mvt2 << endl;
              
              // Old:
              path_likelihood =
                backward_row[ col_i ][ Match ];
              path_likelihood /=
                parameters.matrixRowScaleFactor;
              path_likelihood *=
                profile[
                  Transition::fromPostAlign
                ][
                  TransitionFromPostAlign::toPostAlign
                ];
              path_likelihood *=
                getPostAlignEmissionProbability(
                  profile,
                  sequence,
                  col_i - 1
                );
              // And we need to correct for the post-align transition being within-row..
              path_likelihood *= rabiner_within_row_correction_factor;
              tmp_mvt = forward_row[ col_i - 1 ][ Match ];
              tmp_mvt /= parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
              // From DeletionOut to Post-Align
              // Note that DeletionOut is the only thing we do store separately
              // from post-aligns in the bottom row.
              tmp_mvt2 = forward_row[ col_i - 1 ][ DeletionOut ];
              tmp_mvt2 /= parameters.matrixRowScaleFactor;
              tmp_mvt += tmp_mvt2;
#endif //defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
              path_likelihood *= tmp_mvt;
              //cout << "\t Should be the same as: " << path_likelihood << endl;
              cout << path_likelihood << " should equal " << tmp_mvt2 << endl;
              cout << "\tthe difference is " << ( toDouble(path_likelihood) - toDouble(tmp_mvt2) ) << endl;
              //assert( path_likelihood == tmp_mvt2 );
              assert( abs( toDouble(path_likelihood) - toDouble(tmp_mvt2) ) < (10*numeric_limits<double>::epsilon()) );
            } // End double-checking that we get the same result either way

          } // End if it's the last row..
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          if( row_i == 1 ) {
            path_likelihood *=
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            path_likelihood *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toMatch
              ];
          } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              prev_row_match_distribution[
                TransitionFromMatch::toMatch
              ];
#else
            path_likelihood *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toMatch
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          }
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i - 1 ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          if( coefficients != 0 ) {
            ( *coefficients )[ Emission::Match ][ residue ] +=
              path_likelihood;
          } // End if also computing coefficients
          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Match ][ residue ];

          pos[ Emission::Match ][ residue ] +=
            path_likelihood;
          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            pos[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else {
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *2*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score
          } // End if row_i != last_row
        } // End if col_i == 0, elsif row_i == 0, else ..
        // from Match to Insertion
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != last_row ) ) {
          path_likelihood =
            backward_row[ col_i ][ Insertion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          path_likelihood *=
            row_match_distribution[
              TransitionFromMatch::toInsertion
            ];
#else
          path_likelihood *=
            profile[ row_i - 1 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Insertion ][ residue ];
          tmp_mvt = forward_row[ col_i - 1 ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          path_likelihood *= rabiner_within_row_correction_factor;

          // TODO: REMOVE
          if( path_likelihood >= 1E10 ) {
            cout << "backward_row[ " << col_i << " ][ Insertion ]: " << backward_row[ col_i ][ Insertion ] << endl;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            cout << "row_match_distribution[ Transition::fromMatch ][ TransitionFromMatch::toInsertion ]: " << 
            row_match_distribution[
              TransitionFromMatch::toInsertion
            ] << endl;
#else
            cout << "profile[ " << row_i << " - 1 ][ Transition::fromMatch ][ TransitionFromMatch::toInsertion ]: " << 
            profile[ row_i - 1 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ] << endl;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            cout << "profile[ row_i - 1 ][ Emission::Insertion ][ residue ]: " << profile[ row_i - 1 ][ Emission::Insertion ][ residue ] << endl;
            cout << "forward_row[ col_i - 1 ][ Match ]: " << forward_row[ col_i - 1 ][ Match ] << endl;
            cout << "rabiner_within_row_correction_factor: " << rabiner_within_row_correction_factor << endl;
            assert( path_likelihood < 1E10 );
          }
          pos[
            Transition::fromMatch
          ][
            TransitionFromMatch::toInsertion
          ] += path_likelihood;
            // TODO: REMOVE
            //cout << "calculateAlignmentProfilePosition(..): *1* about to add " << path_likelihood << " to the insertion emission distribution for residue " << residue << endl;
          pos[ Emission::Insertion ][ residue ] +=
            path_likelihood;
          // Don't include this path_likelihood in the matrix_value_type_score,
          // since this is a transition to an Insertion state, so the path will
          // be counted elsewhere.
        } // End if this is not the first column and not the first or last row ..
        // from Match to Deletion
        if( row_i != 0 ) {
          if( row_i == 1 ) {
            path_likelihood =
              backward_row[ col_i ][ Deletion ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            path_likelihood *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletion
              ];
            tmp_mvt = prev_forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            if( coefficients != 0 ) {
              coefficients->m_constant +=
                path_likelihood;
            } // End if we are calculating coefficients

            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *3*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score

          } else if( row_i == last_row ) {
            // We store the postAlignInsertion stuff in the Match state.
            path_likelihood =
              backward_row[ col_i ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              prev_row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            path_likelihood *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
#ifdef USE_END_DISTRIBUTION
            path_likelihood *=
              // And we transition from the end to the postalign.
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            tmp_mvt = prev_forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            if( coefficients != 0 ) {
              coefficients->m_constant +=
                path_likelihood;
            } // End if we are calculating coefficients

#ifdef USE_END_DISTRIBUTION
            pos[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else { // row_i is not 0, 1, or last_row
            path_likelihood =
              backward_row[ col_i ][ Deletion ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              prev_row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            path_likelihood *=
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            tmp_mvt = prev_forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            if( coefficients != 0 ) {
              coefficients->m_constant +=
                path_likelihood;
            } // End if we are calculating coefficients

            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *4*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score

          } // End if it's row 1 .. elsif the last row .. else ..
        } // End if it's not the first row ..
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
        // from Match to DeletionOut
        if(
          ( row_i > 1 ) && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          path_likelihood =
            backward_row.m_deletionOut;
#else
          path_likelihood =
            backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          path_likelihood /=
            parameters.matrixRowScaleFactor;

          // Note that we do not use the scaled distribution, since the
          // backwards row incorporates the rest of the del-out path.  We just
          // want the deletion OPEN probability.
          path_likelihood *=
            profile[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletionOut
            ];
          tmp_mvt = prev_forward_row[ col_i ][ Match ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          // Note we don't need to multiply in the DeletionOut::toEnd or the
          // End::toPostAlign probabilities, because these are already
          // incorporated into the backward_row.  But we do need to count them.
          if( row_i == last_row ) {
            // TODO: REMOVE
            //cout << "For M->W INTO last_row, adding " << path_likelihood << " to the W->E AlignmentProfile value." << endl;
            pos[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toEnd
            ] += path_likelihood;
#ifdef USE_END_DISTRIBUTION
            // If this is the last row, that transition includes from End to
            // PostAlign.
            pos[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } // End if row_i == last_row

          if( do_extra_debugging ) {
            cout << "col " << col_i << ": adding " << path_likelihood << " to the constant for the del-out open." << endl;
          }
          if( coefficients != 0 ) {
            coefficients->m_constant +=
              path_likelihood;
          } // End if we are calculating coefficients
  
          if( calculate_score && ( row_i != last_row ) ) {
            matrix_value_type_score += path_likelihood;
          } // End if calculate_score
      } // End if ( row_i > 1 ) && ( col_i > 0 ) (and maybe col_i == last_col)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        // End from Match

        // from Insertion
        // from Insertion to Match
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != 1 ) ) { // No ins state in row 0
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toMatch
            ];
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i - 1 ][ Insertion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          if( coefficients != 0 ) {
            ( *coefficients )[ Emission::Match ][ residue ] +=
              path_likelihood;
          } // End if also computing coefficients
            path_likelihood *=
              profile[ row_i - 1 ][ Emission::Match ][ residue ];

            pos[ Emission::Match ][ residue ] +=
              path_likelihood;
          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            pos[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else {
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *5*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score
          } // End if row_i != last_row
        } // End if it's not the first or second row and not the first col ..
        // from Insertion to Insertion
        if( ( row_i != 0 ) &&
            ( row_i != last_row ) &&
            ( col_i != 0 ) ) {
          path_likelihood =
            backward_row[ col_i ][ Insertion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 1 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toInsertion
            ];
          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Insertion ][ residue ];
          tmp_mvt = forward_row[ col_i - 1 ][ Insertion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          path_likelihood *= rabiner_within_row_correction_factor;

          // TODO: REMOVE
          assert( path_likelihood < 1E10 );
          pos[
            Transition::fromInsertion
          ][
            TransitionFromInsertion::toInsertion
          ] += path_likelihood;
            // TODO: REMOVE
            //cout << "calculateAlignmentProfilePosition(..): *2* about to add " << path_likelihood << " to the insertion emission distribution for residue " << residue << endl;
          pos[ Emission::Insertion ][ residue ] +=
            path_likelihood;
          // Don't include this path in the matrix_value_type_score, since this
          // is a transition to an Insertion..
        } // End if this is not the first or last row, and not the first col..
        // End from insertion

        // from Deletion
        // from Deletion to Match
        if( row_i == last_row ) {
          // We store post-align stuff in the Match state, even if it comes
          // from a deletion.  So this is really a D->D, stored at M.
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
#ifdef USE_END_DISTRIBUTION
          path_likelihood *=
            // And we transition from the end to the postalign.
            profile[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ];
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i ][ Deletion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          if( coefficients != 0 ) {
            coefficients->m_constant +=
              path_likelihood;
          } // End if we are calculating coefficients

#ifdef USE_END_DISTRIBUTION
          pos[
            Transition::fromEnd
          ][
            TransitionFromEnd::toPostAlign
          ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
        } // End if row_i == last_row
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != 1 ) ) { // no del state in row 0
          path_likelihood =
            backward_row[ col_i ][ Match ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toMatch
            ];

#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            path_likelihood *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_mvt = prev_forward_row[ col_i - 1 ][ Deletion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;
          if( coefficients != 0 ) {
            ( *coefficients )[ Emission::Match ][ residue ] +=
              path_likelihood;
          } // End if also computing coefficients

          path_likelihood *=
            profile[ row_i - 1 ][ Emission::Match ][ residue ];

          pos[ Emission::Match ][ residue ] +=
            path_likelihood;
          if( row_i == last_row ) {
#ifdef USE_END_DISTRIBUTION
            pos[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ] += path_likelihood;
#endif // USE_END_DISTRIBUTION
          } else {
            if( calculate_score ) {
              matrix_value_type_score += path_likelihood;
              // TODO: REMOVE
              //if( matrix_value_type_score == 0 ) {
              //  cout << "uh-oh *6*, matrix_value_type_score is still 0 after adding in " << path_likelihood << endl;
              //}
            } // End if calculate_score
          } // End if row_i == last_row .. else ..
        } // End if it's not the first or second row and not the first col..
        // from Deletion to Deletion
        if( ( row_i != 0 ) &&
            ( row_i != 1 ) && // no del state in first row
            ( row_i != last_row ) ) {
          path_likelihood =
            backward_row[ col_i ][ Deletion ];
          path_likelihood /=
            parameters.matrixRowScaleFactor;
          path_likelihood *=
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
          tmp_mvt = prev_forward_row[ col_i ][ Deletion ];
          tmp_mvt /= parameters.matrixRowScaleFactor;
          path_likelihood *= tmp_mvt;

          if( coefficients != 0 ) {
            coefficients->m_constant +=
              path_likelihood;
          } // End if we are calculating coefficients

          if( calculate_score ) {
            matrix_value_type_score += path_likelihood;
          }
        } // End if it's not the first, second, or last row ..
        // End from Deletion

        ////// OUT OF row_i //////
        if( row_i == last_row ) {
          // Then do nothing.  No "OUT OF" last_row transitions.
        } else { // if row_i == last_row .. else .. 
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          // to/from DeletionIn
          // from DeletionIn to Match
          if(
            ( ( row_i + 1 ) > 0 ) && ( ( col_i + 1 ) >= 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            && ( ( col_i + 1 ) == 1 ) 
#endif // DISALLOW_FLANKING_TRANSITIONS
          ) {
            path_likelihood =
              next_backward_row[ ( col_i + 1 ) ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromDeletionIn
              ][
                TransitionFromDeletionIn::toMatch
              ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
            tmp_mvt =
              forward_row.m_deletionIn;
#else
            tmp_mvt =
              forward_row[ ( col_i + 1 ) - 1 ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS
            tmp_mvt /=
              parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;
#ifdef USE_END_DISTRIBUTION
            if( ( row_i + 1 ) == last_row ) {
              // Note this is counted in the "OUT OF" for row_i == last_row;
              // see above.
              path_likelihood *=
                profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
            }
#endif // USE_END_DISTRIBUTION
            path_likelihood *=
              profile[ ( row_i + 1 ) - 1 ][ Emission::Match ][ next_residue ];

            pos[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toMatch
            ] += path_likelihood;
          } // End if ( ( row_i + 1 ) > 0 ) && ( ( col_i + 1 ) >= 1 ) ( and maybe ( col_i + 1 ) == 1 )

          // from DeletionIn to DeletionIn (and from Match/Begin to DeletionIn)
          if(
            ( ( row_i + 1 ) == 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
          ) {
            // Actually this is from PreAlign to Begin, from Begin to DeletionIn.
#ifdef DISALLOW_FLANKING_TRANSITIONS
            path_likelihood =
              next_backward_row.m_deletionIn;
#else
            path_likelihood =
              next_backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            path_likelihood *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletionIn
              ];
            tmp_mvt =
              forward_row[ col_i ][ Match ];
            tmp_mvt /=
              parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ] += path_likelihood;
            pos[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletionIn
            ] += path_likelihood;
          } else if(
            ( ( row_i + 1 ) != 0 ) && ( ( row_i + 1 ) != last_row )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
          ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
            path_likelihood =
              next_backward_row.m_deletionIn;
#else
            path_likelihood =
              next_backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromDeletionIn
              ][
                TransitionFromDeletionIn::toDeletionIn
              ];

#ifdef DISALLOW_FLANKING_TRANSITIONS
            tmp_mvt =
              forward_row.m_deletionIn;
#else
            tmp_mvt =
              forward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
            tmp_mvt /=
              parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toDeletionIn
            ] += path_likelihood;
          } // End if ( ( row_i + 1 ) == 1 ) .. else if( ( ( row_i + 1 ) != 0 ) && ( ( row_i + 1 ) != last_row ) ) (and maybe also ( col_i == 0 )).
          // End to/from DeletionIn

          // From DeletionOut
          // from DeletionOut to End
          // Not an OUT_OF transition.  It is counted above..

          // from DeletionOut to DeletionOut
          if(
            ( ( row_i + 1 ) > 2 )
            && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
          ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
            path_likelihood =
              next_backward_row.m_deletionOut;
#else
            path_likelihood =
              next_backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromDeletionOut
              ][
                TransitionFromDeletionOut::toDeletionOut
              ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
            tmp_mvt =
              forward_row.m_deletionOut;
#else
            tmp_mvt =
              forward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            // TODO: REMOVE
            //cout << "For W->W OUT OF row " << row_i << ", adding " << path_likelihood << " to the W->W AlignmentProfile value." << endl;
            pos[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toDeletionOut
            ] += path_likelihood;
          } // End if ( ( row_i + 1 ) > 2 ) && ( col_i > 0 ) (and maybe col_i == last_col)
          // End from DeletionOut
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

          // from Match
          // from Match to Match
          if( col_i != last_col ) {
            path_likelihood =
              next_backward_row[ ( col_i + 1 ) ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;

            if( ( row_i + 1 ) == 1 ) {
              path_likelihood *=
                profile[
                  Transition::fromPreAlign
                ][
                  TransitionFromPreAlign::toBegin
                ];
              path_likelihood *=
                profile[
                  Transition::fromBegin
                ][
                  TransitionFromBegin::toMatch
                ];
            } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              path_likelihood *=
                row_match_distribution[
                  TransitionFromMatch::toMatch
                ];
#else
              path_likelihood *=
                profile[ ( row_i + 1 ) - 2 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toMatch
                ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            }
#ifdef USE_END_DISTRIBUTION
            if( ( row_i + 1 ) == last_row ) {
              path_likelihood *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION
            tmp_mvt = forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;
            path_likelihood *=
              profile[ ( row_i + 1 ) - 1 ][ Emission::Match ][ next_residue ];

            // TODO: REMOVE
            //assert( !isnan( path_likelihood ) );
            if( ( row_i + 1 ) == 1 ) {
              pos[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ] += path_likelihood;
              pos[
                Transition::fromBegin
              ][
                TransitionFromBegin::toMatch
              ] += path_likelihood;
              // In the alignment profile we are tricksy: we store "pre-align
              // Opens" in the Match distribution (as Match->Insertion), and
              // "pre-align Extensions" in the pre-align distribution, so we
              // can differentiate.  We also keep track of the M->M and M->D
              // transitions in the first column, to know how many
              // non-pre-align-opens there were.
              if( col_i == 1 ) {
                pos[
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toMatch
                ] += path_likelihood;
              }
            } else {
              pos[
                Transition::fromMatch
              ][
                TransitionFromMatch::toMatch
              ] += path_likelihood;
              // TODO: REMOVE
              //assert( !isnan( pos[
              //  Transition::fromMatch
              //][
              //  TransitionFromMatch::toMatch
              //] ) );
            } // end if ( row_i + 1 ) == 1 .. else ..
            
            // from Match to Insertion
            // Not an "OUT OF" transition
          } // End if col_i != last_col
          
          // from Match to Deletion
          if( ( row_i + 1 ) == 1 ) {
            path_likelihood =
              next_backward_row[ col_i ][ Deletion ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            path_likelihood *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletion
              ];
            tmp_mvt = forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ] += path_likelihood;
            pos[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletion
            ] += path_likelihood;
            // In the alignment profile we store "pre-align Opens" in the
            // Match distribution (as Match->Insertion), and "pre-align
            // Extensions" in the pre-align distribution, so we can
            // differentiate.  We also keep track of the M->M and M->D
            // transitions in the first column, to know how many
            // non-pre-align-opens there were.
            if( col_i == 0 ) {
              pos[
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ] += path_likelihood;
            }
          } else if( ( row_i + 1 ) == last_row ) {
            // We store the postAlignInsertion stuff in the Match state.
            path_likelihood =
              next_backward_row[ col_i ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            path_likelihood *=
              profile[ ( row_i + 1 ) - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
#ifdef USE_END_DISTRIBUTION
            path_likelihood *=
              // And we transition from the end to the postalign.
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            tmp_mvt = forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ] += path_likelihood;
          } else { // ( row_i + 1 ) is not 0, 1, or last_row
            path_likelihood =
              next_backward_row[ col_i ][ Deletion ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            path_likelihood *=
              row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            path_likelihood *=
              profile[ ( row_i + 1 ) - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            tmp_mvt = forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ] += path_likelihood;
          } // End if it's row 1 .. elsif the last row .. else ..
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          // from Match to DeletionOut
          if(
            ( ( row_i + 1 ) > 1 ) && ( col_i > 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
          ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
            path_likelihood =
              next_backward_row.m_deletionOut;
#else
            path_likelihood =
              next_backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
            path_likelihood /=
              parameters.matrixRowScaleFactor;

            // Note that we do not use the scaled distribution, since the
            // backwards row incorporates the rest of the del-out path.  We just
            // want the deletion OPEN probability.
            path_likelihood *=
              profile[
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletionOut
              ];
            tmp_mvt = forward_row[ col_i ][ Match ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            // TODO: REMOVE
            //cout << "For M->W OUT OF row " << row_i << ", adding " << path_likelihood << " to the M->W AlignmentProfile value." << endl;
            pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletionOut
            ] += path_likelihood;
          } // End if ( ( row_i + 1 ) > 1 ) && ( col_i > 0 ) (and maybe col_i == last_col)
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          // End from Match

          // from Insertion
          // from Insertion to Match
          if( ( col_i != last_col ) &&
              ( ( row_i + 1 ) != 1 ) ) {
            path_likelihood =
              next_backward_row[ ( col_i + 1 ) ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[ ( row_i + 1 ) - 2 ][
                Transition::fromInsertion
              ][
                TransitionFromInsertion::toMatch
              ];

              path_likelihood *=
                profile[ ( row_i + 1 ) - 1 ][ Emission::Match ][ next_residue ];
#ifdef USE_END_DISTRIBUTION
            if( ( row_i + 1 ) == last_row ) {
              path_likelihood *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION
            tmp_mvt = forward_row[ ( col_i + 1 ) - 1 ][ Insertion ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toMatch
            ] += path_likelihood;
          } // End if it's not the first row and not the last col ..
          
          // from Insertion to Insertion: Not an OUT OF transition
          // End from insertion
          
          // from Deletion
          // from Deletion to Match
          if( ( row_i + 1 ) == last_row ) {
            // We store post-align stuff in the Match state, even if it comes
            // from a deletion.  So this is really a D->D, stored at M.
            path_likelihood =
              next_backward_row[ col_i ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[ ( row_i + 1 ) - 2 ][
                Transition::fromDeletion
              ][
                TransitionFromDeletion::toDeletion
              ];
#ifdef USE_END_DISTRIBUTION
            path_likelihood *=
              // And we transition from the end to the postalign.
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            tmp_mvt = forward_row[ col_i ][ Deletion ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ] += path_likelihood;
          } // End if ( row_i + 1 == last_row )
          if( ( col_i != last_col ) &&
              ( ( row_i + 1 ) != 1 ) ) {
            path_likelihood =
              next_backward_row[ ( col_i + 1 ) ][ Match ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[ ( row_i + 1 ) - 2 ][
                Transition::fromDeletion
              ][
                TransitionFromDeletion::toMatch
              ];

              path_likelihood *=
                profile[ ( row_i + 1 ) - 1 ][ Emission::Match ][ next_residue ];
#ifdef USE_END_DISTRIBUTION
            if( ( row_i + 1 ) == last_row ) {
              path_likelihood *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION
            tmp_mvt = forward_row[ ( col_i + 1 ) - 1 ][ Deletion ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;

            pos[
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toMatch
            ] += path_likelihood;
          } // End if it's not the first row and not the last col..
          // from Deletion to Deletion
          if( ( ( row_i + 1 ) != 1 ) && // no del state in first row
              ( ( row_i + 1 ) != last_row ) ) {
            path_likelihood =
              next_backward_row[ col_i ][ Deletion ];
            path_likelihood /=
              parameters.matrixRowScaleFactor;
            path_likelihood *=
              profile[ ( row_i + 1 ) - 2 ][
                Transition::fromDeletion
              ][
                TransitionFromDeletion::toDeletion
              ];
            tmp_mvt = forward_row[ col_i ][ Deletion ];
            tmp_mvt /= parameters.matrixRowScaleFactor;
            path_likelihood *= tmp_mvt;
          
            pos[
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ] += path_likelihood;

          } // End if it's not the first, second, or last row ..
          // End from Deletion
        } // End if row_i != last_row, compute OUT OF transitions

      } // End for each column, ..

//#ifdef DISALLOW_FLANKING_TRANSITIONS
//      if( row_i != 0 ) {
//        pos[
//          Transition::fromPreAlign
//        ].zero();
//      }
//      if( row_i != last_row ) {
//        pos[
//          Transition::fromPostAlign
//        ].zero();
//      }
//#endif // DISALLOW_FLANKING_TRANSITIONS

      // TODO: REMOVE
      //if( row_i == last_row ) {
      //cout << "calculateAlignmentProfilePosition: pos (before scaling) is " << pos << endl;
      //}

      // TODO: REMOVE
      if( do_extra_debugging ) {
        if( calculate_score ) {
          cout << "calculateAlignmentProfilePosition: matrix_value_type_score is:" << matrix_value_type_score << endl;
        } else {
          if( parameters.useRabinerScaling ) {
            if( !parameters.rabinerScaling_useMaximumValue ) {
              cout << "in calculateAlignmentProfilePosition(..): score (calculated using the cumulative_inverse_scalars) is " << cumulative_inverse_scalar << endl;
            }
          } // End if useRabinerScaling
          cout << "in calculateAlignmentProfilePosition(..): inverse_scalar is " << ( ( inverse_scalar == 0 ) ? ScoreType( 0 ) : *inverse_scalar ) << endl;
        }
      } // End if do_extra_debugging

      // TODO: REMOVE.
      // ASSUMING inverse_scalar is the score.
      // Ultimately we can turn off calculate_score when inverse_scalar == NULL.
      if( false && !parameters.useRabinerScaling && calculate_score && ( inverse_scalar != NULL ) && ( *inverse_scalar != 0 ) ) {
        // TODO: REMOVE!
        ScoreType alt_score =
          forward_score(
            parameters,
            profile,
            row_i,
            forward_row,
            backward_row
          );

        /*
        ScoreType mvts_is_diff = ( ( matrix_value_type_score >= *inverse_scalar ) ? ( static_cast<ScoreType>( matrix_value_type_score ) - *inverse_scalar ) : ( *inverse_scalar - static_cast<ScoreType>( matrix_value_type_score ) ) );
        ScoreType mvts_is_pct = mvts_is_diff;
        mvts_is_pct /= ( ( matrix_value_type_score >= *inverse_scalar ) ? *inverse_scalar : static_cast<ScoreType>( matrix_value_type_score ) );
        ScoreType as_is_diff = ( ( alt_score >= *inverse_scalar ) ? ( alt_score - *inverse_scalar ) : ( *inverse_scalar - alt_score ) );
        ScoreType as_is_pct = as_is_diff;
        as_is_pct /= ( ( alt_score >= *inverse_scalar ) ? *inverse_scalar : alt_score );
        // TODO: REMOVE
        if( false && ( row_i == last_row ) ) {
          if( matrix_value_type_score >= *inverse_scalar ) {
            cout << "matrix_value_type_score - *inverse_scalar is " << ( matrix_value_type_score - *inverse_scalar ) << endl;
          } else {
            cout << "*inverse_scalar - matrix_value_type_score is " << ( *inverse_scalar - matrix_value_type_score ) << endl;
          }
          cout << "\t mvts_is_pct is " << mvts_is_pct << endl;
          if( alt_score >= *inverse_scalar ) {
            cout << "alt_score - *inverse_scalar is " << ( alt_score - *inverse_scalar ) << endl;
          } else {
            cout << "*inverse_scalar - alt_score is " << ( *inverse_scalar - alt_score ) << endl;
          }
          cout << "\t as_is_pct is " << as_is_pct << endl;
          cout << "\t row_i is " << row_i << endl;
          cout << "\t last_row is " << last_row << endl;
          cout << "\t last_col is " << last_col << endl;
          //cout << "\t sequence is " << sequence << endl;
          cout << "\t score, calculated here, is " << matrix_value_type_score << endl;
          cout << "\t *inverse_scalar is " << *inverse_scalar << endl;
          cout << "\t score, recalculated using forward_score( Parameters, bool, Row, Row ), is " << alt_score << endl;
          //cout << "\t forward_row.m_rabinerInverseScalar is " << forward_row.m_rabinerInverseScalar << endl;
          //cout << "\t forward_row.m_rabinerCumulativeInverseScalar is " << forward_row.m_rabinerCumulativeInverseScalar << endl;
          //cout << "\t backward_row.m_rabinerInverseScalar is " << backward_row.m_rabinerInverseScalar << endl;
          //cout << "\t backward_row.m_rabinerCumulativeInverseScalar is " << backward_row.m_rabinerCumulativeInverseScalar << endl;
        } // End if row_i == 0
        if(
          ( mvts_is_pct > 1E-4 ) ||
          ( as_is_pct > 1E-4 )
        ) {
          if( mvts_is_pct > 1E-4 ) {
            cout << "ERROR: score calculated in calculateAlignmentProfilePosition(..) is different than the score given in *inverse_scalar: " << matrix_value_type_score << " != " << *inverse_scalar << endl;
          }
          if( as_is_pct > 1E-4 ) {
            cout << "ERROR: score given in *inverse_scalar is different than score calculated using forward_score( Parameters, bool, Row, Row ): " << *inverse_scalar << " != " << alt_score << endl;
          }
          cout << "\t row_i is " << row_i << endl;
          cout << "\t last_row is " << last_row << endl;
          cout << "\t last_col is " << last_col << endl;
          cout << "\t sequence is " << sequence << endl;
          //cout << "forward_row.m_rabinerInverseScalar is " << forward_row.m_rabinerInverseScalar << endl;
          //cout << "forward_row.m_rabinerCumulativeInverseScalar is " << forward_row.m_rabinerCumulativeInverseScalar << endl;
          //cout << "backward_row.m_rabinerInverseScalar is " << backward_row.m_rabinerInverseScalar << endl;
          //cout << "backward_row.m_rabinerCumulativeInverseScalar is " << backward_row.m_rabinerCumulativeInverseScalar << endl;
          assert( mvts_is_pct <= 1E-4 );
          assert( as_is_pct <= 1E-4 );
        }
        */
      } // End if calculate_score && ( inverse_scalar != NULL )

      if( inverse_scalar == NULL ) {
        if( parameters.useRabinerScaling ) {
          // We need to scale it, since it is presently scaled by
          // the cumulative_inverse_scalar, not by the sequence score.
          //cout << "calculateAlignmentProfilePosition: scaling by score / cumulative_inverse_scalar = " << matrix_value_type_score << " / " << cumulative_inverse_scalar << endl;
          pos.m_scalar *= matrix_value_type_score;
          pos.m_scalar /= cumulative_inverse_scalar;
          //cout << "calculateAlignmentProfilePosition: NOT scaling by score calculated from rabiner scalars:" << cumulative_inverse_scalar << endl;
          //pos.m_scalar *= cumulative_inverse_scalar;
            //( matrix_value_type_score ) / cumulative_inverse_scalar;
            //matrix_value_type_score;
          //( ( forward_row.m_rabinerCumulativeInverseScalar *
          //      backward_row.m_rabinerCumulativeInverseScalar ) /
          //    cumulative_inverse_scalar );
        } else {
          //cout << "calculateAlignmentProfilePosition: scaling by score = " << matrix_value_type_score << endl;
          pos.m_scalar *= matrix_value_type_score;
        } // End if useRabinerScaling .. else ..
        // TODO: REMOVE
        //cout << "calculateAlignmentProfilePosition: pos (after scaling by the inverse of the sequence score) is " << pos << endl;
      } else if( *inverse_scalar != 0.0 ) {
        // TODO: REMOVE.  TESTING
#ifndef NDEBUG
        if( calculate_score ) {
          if( matrix_value_type_score != *inverse_scalar ) {
            bool mvt_score_gt_inv_scalar = ( static_cast<ScoreType>( matrix_value_type_score ) > *inverse_scalar );
            ScoreType diff = ( mvt_score_gt_inv_scalar ? ( static_cast<ScoreType>( matrix_value_type_score ) - *inverse_scalar ) : ( *inverse_scalar  - static_cast<ScoreType>( matrix_value_type_score ) ) );
            if( ( diff / ( mvt_score_gt_inv_scalar ? static_cast<ScoreType>( matrix_value_type_score ) : *inverse_scalar ) ) > 1E-5 ) { // TODO: DEHACKIFY MAGIC # MACHINE EPSILON!
              cout << "UH OH: matrix_value_type_score is " << matrix_value_type_score << " but *inverse_scalar is " << *inverse_scalar << ".  Their diff is ( matrix_value_type_score - *inverse_scalar ) = " << ( mvt_score_gt_inv_scalar ? ' ' : '-' ) << diff << ", and abs(diff) / (the greater value) is " << ( diff / ( mvt_score_gt_inv_scalar ? static_cast<ScoreType>( matrix_value_type_score ) : *inverse_scalar ) ) << "." << endl;
              // NOTE: When using anchor rows and columns, the score can drift slightly, and that's ok.  It should suffice to check that the scores are on the same scale (eg. < 1 rather than < 1E-5 after diving out the scale)
              assert(  toDouble( ( diff / ( mvt_score_gt_inv_scalar ? static_cast<ScoreType>( matrix_value_type_score ) : *inverse_scalar ) ) < 1 ) ); // TODO: DEHACIFY MAGIC # 1 !!!
              //assert(  toDouble( ( diff / ( mvt_score_gt_inv_scalar ? static_cast<ScoreType>( matrix_value_type_score ) : *inverse_scalar ) ) < 1E-5 ) ); // TODO: DEHACIFY MAGIC # MACHINE EPSILON!
            }
          }
        }
#endif // !NDEBUG
        // TODO: REMOVE
        //cout << "calculateAlignmentProfilePosition: scaling by *inverse_scalar:" << *inverse_scalar << endl;
        if( parameters.useRabinerScaling ) {
          // It's already been scaled, but by the cumulative_inverse_scalar
          pos.m_scalar *= *inverse_scalar;
          pos.m_scalar /= cumulative_inverse_scalar;
        } else {
          pos.m_scalar *= *inverse_scalar;
        } // End if useRabinerScaling .. else ..
        // TODO: REMOVE
        //cout << "calculateAlignmentProfilePosition: pos (after scaling by the inverse of the argument value) is " << pos << endl;
      } // End if we don't have an inverse scalar .. else if we do ..

      // TODO: REMOVE
      assert( !isnan( pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toMatch
            ] ) );
      // TODO: REMOVE
      assert( !isinf( pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toMatch
            ] ) );
      // TODO: REMOVE
      assert( !isinf( pos.m_scalar ) );

      // TODO: REMOVE!
      //cout << "before rescale(), pos is " << pos << endl;
      //cout << "before rescale(), pos.createUnscaledCopy() is " << pos.createUnscaledCopy() << endl;
      // Rescale to avoid underflow later.
      // TODO: PUT BACK!
      pos.rescale();
      // TODO: REMOVE!
      //cout << "after rescale(), pos is " << pos << endl;
      //cout << "after rescale(), pos.createUnscaledCopy() is " << pos.createUnscaledCopy() << endl;

      if( coefficients != 0 ) {
        coefficients->m_inverseScalar = cumulative_inverse_scalar;
      }

      // TODO: REMOVE
      assert( !isinf( pos.m_scalar ) );
      // TODO: REMOVE
      if( isnan( pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toMatch
              ] ) ) {
        cout << "Uh-oh, calculateAlignmentProfilePosition(..): got nan!" << endl;
        cout << "\t row_i is " << row_i << endl;
        cout << "\t last_row is " << last_row << endl;
        cout << "\t last_col is " << last_col << endl;
        if( inverse_scalar == NULL ) {
          cout << "\t inverse_scalar is NULL" << endl;
          cout << "\t matrix_value_type_score is " << matrix_value_type_score << endl;
        } else {
          cout << "\t inverse_scalar is " << *inverse_scalar << endl;
          if( *inverse_scalar == 0.0 ) {
            cout << "\t cumulative_inverse_scalar = " << cumulative_inverse_scalar << endl;
          }
        }
        if( row_i < last_row ) {
          cout << "\t profile[ ( row_i + 1 ) - 1 ] is " << profile[ ( row_i + 1 ) - 1 ] << endl;
        }
        cout << "\t pos is " << pos << endl;
        assert( !isnan( pos[
              Transition::fromMatch
            ][
              TransitionFromMatch::toMatch
              ] ) );
        exit( 1 );
      } // end if isnan..

      // TODO: REMOVE
      /*
      if(
        false &&
        ( row_i != 0 ) && ( row_i != last_row )
      ) {
        MatrixValueType transition_total =
          pos[ Transition::fromMatch ].total();
        MatrixValueType emission_total = pos[ Emission::Match ].total();
        //cout << "after rescale, transition total is " << transition_total << ", emission_total is " << emission_total << ", and the difference is " << ( ( transition_total > emission_total ) ? ( transition_total - emission_total ) : ( emission_total - transition_total ) ) << endl;
        if(
          (
            ( transition_total > emission_total ) &&
            ( ( transition_total - emission_total ) > 1E-4 )
          ) ||
          (
            ( transition_total < emission_total ) &&
            ( ( emission_total - transition_total ) > 1E-4 )
          )
        ) {
          cout << "ERROR: The total emissions from the Match state is " << emission_total << ", but the total transitions from the Match state out of this row is " << transition_total << ", and the difference is " << ( ( transition_total > emission_total ) ? ( transition_total - emission_total ) : ( emission_total - transition_total ) ) << endl;
          cout << "row_i is " << row_i << "; last_row is " << last_row << endl;
          cout << "\t last_col is " << last_col << endl;
          cout << "pos is " << pos << endl;
          cout << "sequence is " << sequence << endl;
          //cout << "prev_forward_row is " << prev_forward_row << endl;
          //cout << "forward_row is " << forward_row << endl;
          //typename Matrix::Row new_forward_row( forward_row );
          //forward_calculateRow(
          //  parameters,
          //  false, // not viterbi
          //  profile,
          //  sequence,
          //  row_i,
          //  prev_forward_row,
          //  new_forward_row
          //);
          //cout << "After recalculating forward row, it is " << new_forward_row << endl;
          //cout << "next_backward_row is " << next_backward_row << endl;
          //cout << "backward_row is " << backward_row << endl;
          //typename Matrix::Row new_backward_row( backward_row );
          //backward_calculateRow(
          //  parameters,
          //  false, // not viterbi
          //  profile,
          //  sequence,
          //  row_i,
          //  next_backward_row,
          //  new_backward_row
          //);
          //cout << "After recalculating backward row, it is " << new_backward_row << endl;

          assert(
            pos[ Emission::Match ].total() ==
            pos[ Transition::fromMatch ].total()
          );
        } // End if the total emitted from the match state is not the same as
          // the total transitions out of the match state...
      } // End if row_i is not the first nor the last row...
      */

      return matrix_value_type_score;
    } // calculateAlignmentProfilePosition ( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, Row const&, Row const&, AlignmentProfilePosition & ) const

#ifdef ALLOW_BOLTZMANN_GIBBS
  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType>
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Use the given ProfilePosition reference
     * and Coefficients reference to update the
     * PositionBoltzmannGibbs corresponding to the given profile position
     * with Baldi's gradient ascent step (equation 3 from Pierre Baldi, Yves
     * Chauvin).  If the given inverse scalar is non-NULL, it will be used
     * instead of the calculated sequence score for scaling the change.  The
     * return value is the (not scaled) sequence score.  Note that the learning
     * rate will *not* be incorporated; the caller is responsible for
     * multiplying the resulting PositionBoltzmannGibbs by the learning
     * rate afterwards (or equivalently dividing its m_scalar by the learning rate). 
     *
     * Don't forget to zero() the PositionBoltzmannGibbs (before calling
     * this for all sequences).  If the inverse_scalar pointer is non-NULL but
     * points to the value 0, the score will be calculated and returned but the
     * position_boltzmann_gibbs_change will not be updated (it's like having a 0 scalar, not
     * like having a 0 inverse scalar, which would result in a divide-by-0
     * error).
     */
    updatePositionBoltzmannGibbsChangeForSequence (
      Parameters const& parameters,
      ProfilePosition<ResidueType, ProfileType> const& profile_position,
      PositionSpecificSequenceScoreCoefficients const& coefficients,
      ScoreType const * inverse_scalar,
      PositionBoltzmannGibbs<ResidueType, ScoreType> & position_boltzmann_gibbs_change
    ) const
    {
      // WARNING: Not thread safe.  I'm thinking this is faster, but maybe it
      // would be optimized better without the static?
      static PositionEntente sequence_position_entente =
        PositionEntente();
      // WARNING: Not thread safe.  I'm thinking this is faster, but maybe it
      // would be optimized better without the static?
      static ProfilePosition<ResidueType, ProfileType> tmp_position =
        ProfilePosition<ResidueType, ProfileType>();

      sequence_position_entente.zero();
      ScoreType sequence_score =
        calculatePositionEntenteUpdate(
          parameters,
          profile_position,
          coefficients,
          inverse_scalar,
          sequence_position_entente
        );

      if( ( inverse_scalar != NULL ) && ( *inverse_scalar == 0.0 ) ) {
        return sequence_score;
      }
      // Now do the updating:
      position_boltzmann_gibbs_change += sequence_position_entente;

      // TODO: REMOVE
      //cout << "updatePositionBoltzmannGibbsChangeForSequence: just added " << sequence_position_entente << " (which unscaled is " << sequence_position_entente.createUnscaledCopy() << ") to position_boltzmann_gibbs_change (which is now " << position_boltzmann_gibbs_change << "; unscaled: " << position_boltzmann_gibbs_change.createUnscaledCopy() << ")" << endl;

      tmp_position = profile_position;

      // First of all, we need to make sure that the ententes are on the same
      // scale.  We just normalized the position entente, so it is already
      // unscaled.
      sequence_position_entente.unscale();

      // TODO: REMOVE
      //cout << "updatePositionBoltzmannGibbsChangeForSequence: sequence_position_entente[ Emission::Match ].total() is " << sequence_position_entente[ Emission::Match ].total() << endl;

      // Now weigh it...
      tmp_position[ Emission::Match ] *=
        sequence_position_entente[ Emission::Match ].total();
      // TODO: If PositionSpecificParameters ever includes other types, do them too (ie see AlignmentProfilePosition, which includes "global" params too).

      position_boltzmann_gibbs_change -= tmp_position;

      // TODO: REMOVE
      //cout << "updatePositionBoltzmannGibbsChangeForSequence: just subtracted " << tmp_position << " from position_boltzmann_gibbs_change, which is now " << position_boltzmann_gibbs_change << " (unscaled: " << position_boltzmann_gibbs_change.createUnscaledCopy() << ")" << endl;

      return sequence_score;
    } // updatePositionBoltzmannGibbsChangeForSequence( Parameters const&, ProfilePosition const&, PositionEntente const&, PositionSpecificSequenceScoreCoefficients const&, ScoreType const *, PositionBoltzmannGibbs & ) const
#endif // ALLOW_BOLTZMANN_GIBBS


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  GALOSH_INLINE_ALGORITHM
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate the distances among the given sequence ententes, and replace
     * the values of the given PositionEntenteDistances object with these
     * distances.
     */
    updatePositionEntenteDistances (
      vector<PositionEntente> const& position_ententes,
      PositionEntenteDistances & position_entente_distances
    ) const
    {
      bool distance_metric_is_symmetric = true;
      uint32_t num_sequences = position_ententes.size();
      uint32_t from_seq_i, to_seq_i;
      double distance;
      for( from_seq_i = 0; from_seq_i < num_sequences; from_seq_i++ ) {
        
        // TODO: REMOVE
        //cout << "Position entente for seq " << from_seq_i << " is " << position_ententes[ from_seq_i ] << endl;

        for( to_seq_i = 0; to_seq_i < num_sequences; to_seq_i++ ) {
          if( to_seq_i == from_seq_i ) {
            position_entente_distances[ from_seq_i ][ to_seq_i ] = 0;
            if( distance_metric_is_symmetric ) {
              break;
            }
          }
          distance =
            position_ententes[ from_seq_i ].euclideanDistance(
              position_ententes[ to_seq_i ]
            );
          position_entente_distances[ from_seq_i ][ to_seq_i ] = distance;
          if( distance_metric_is_symmetric ) {
            position_entente_distances[ to_seq_i ][ from_seq_i ] = distance;
          }
        }
      }
    } // updatePositionEntenteDistances ( vector<PositionEntente> const&, PositionEntenteDistances & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileTreeType, typename SequenceType>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile tree, using the forward or viterbi algorithm.  The
     * sequences vector must be of length profile_tree.nodeCount().  Use the
     * given matrices, which must be large enough to accommodate the largest
     * set of sequences.
     */
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileTreeType const& profile_tree,
      vector<SequenceType> const& sequences,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      // TODO: REMOVE
      //if( use_viterbi ) {
      //  cout << "root node 0 is " << endl;
      //  cout << *profile_tree.getProfileTreeRoot() << endl;
      //}
      ScoreType score =
        forward_score(
          parameters,
          use_viterbi,
          ( *profile_tree.getProfileTreeRoot() ),
          sequences[ 0 ],
          sequences[ 0 ].size(),
          matrices
        );
      uint32_t tree_size = profile_tree.nodeCount();
      for( uint32_t internal_node_i = 1; internal_node_i < tree_size; internal_node_i++ ) {
        // TODO: REMOVE
        //if( use_viterbi ) {
        //  cout << "internal node " << internal_node_i << " is " << endl;
        //  cout << profile_tree.getProfileTreeInternalNode( internal_node_i ) << endl;
        //}
        score *=
          forward_score(
            parameters,
            use_viterbi,
            profile_tree.getProfileTreeInternalNode( internal_node_i ),
            sequences[ internal_node_i ],
            sequences[ internal_node_i ].size(),
            matrices
          );          
        // TODO: REMOVE
        //if( use_viterbi ) {
        //  cout << "again, internal node " << internal_node_i << " is " << endl;
        //  cout << profile_tree.getProfileTreeInternalNode( internal_node_i ) << endl;
        //}
      } // End foreach internal_node .. 
      return score;
    } // forward_score( Parameters const&, bool, ProfileTree const&, vector<SequenceType> const&, SequentialAccessContainer & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate the likelihood of the sequences whose forward and backward
     * rows are given.  They are assumed to have been pre-calculated (see
     * forward_calculateRow(..) and backward_calculateRow(..)), for the same
     * profile, sequences, and row number.  As a side effect, the given
     * ScoreType vector (if non-NULL) becomes filled with the scores of the
     * sequences, and the largest_score pointer (if non-null) will be filled
     * with the largest of the scores of the sequences.  We need to know if it
     * is the last row because we store things differently there..
     */
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      uint32_t const& sequence_count,
      uint32_t const row_i,
      typename Matrix::RowVector const & forward_rows,
      typename Matrix::RowVector const & backward_rows,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score
    ) const
  {
    assert( sequence_scores->size() == sequence_count );
    assert( forward_rows.size() == sequence_count );
    assert( backward_rows.size() == sequence_count );

    ScoreType score = 1.0;

    if( largest_score != 0 ) {
      *largest_score = 0;
    }
    ScoreType sequence_score;
    for( uint32_t seq_i = 0; seq_i < sequence_count; seq_i++ ) {
      sequence_score =
        forward_score(
          parameters,
          profile,
          row_i,
          forward_rows[ seq_i ],
          backward_rows[ seq_i ]
        );
      if( sequence_scores != 0 ) {
        ( *sequence_scores )[ seq_i ] = sequence_score;
      }
      if( ( largest_score != 0 ) && ( sequence_score > *largest_score ) ) {
        *largest_score = sequence_score;
      }
      score *= sequence_score;
    } // End foreach seq_i

    return score;
  } // forward_score( Parameters const&, bool, uint32_t const&, RowVector &, RowVector &, vector<ScoreType> *, ScoreType * ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate the likelihood of the sequence whose forward and backward rows
     * are given.  They are assumed to have been pre-calculated (see
     * forward_calculateRow(..) and backward_calculateRow(..)), for the same
     * profile, sequence, and row number.  We need to know if it is the last
     * row because we store things differently there..
     */
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      uint32_t const row_i,
      typename Matrix::Row const & forward_row,
      typename Matrix::Row const & backward_row
    ) const
  {
    MatrixValueType sequence_score_mvt = 0.0;

    // TODO: REMOVE
    //cout << "forward_score( parameters, forward_row, backward_row ): forward_row is " << forward_row << endl;
    //cout << "\t backward_row is " << backward_row << endl;

    uint32_t last_col = ( forward_row.size() - 1 );
    assert( backward_row.size() == ( last_col + 1 ) );
    uint32_t last_row = profile.length();

    if( row_i == last_row ) {
      sequence_score_mvt =
        forward_row[ last_col ][ Match ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
      sequence_score_mvt +=
        forward_row.m_deletionOut;
#else
      sequence_score_mvt +=
        forward_row[ last_col ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      sequence_score_mvt *=
        profile[
          Transition::fromPostAlign
        ][
          TransitionFromPostAlign::toTerminal
        ];
      sequence_score_mvt /=
        parameters.matrixRowScaleFactor;
    } else if( row_i == 0 ) {
      sequence_score_mvt = backward_row[ 0 ][ Match ];
      sequence_score_mvt /= parameters.matrixRowScaleFactor;
    } else { // if row_i == last_row .. else if row_i == 0 .. else ..
      MatrixValueType tmp_mvt;
#if defined( USE_DEL_IN_DEL_OUT) && defined( DISALLOW_FLANKING_TRANSITIONS )
      tmp_mvt =
        forward_row.m_deletionIn;
      tmp_mvt /= parameters.matrixRowScaleFactor;
      tmp_mvt *=
        backward_row.m_deletionIn;
      tmp_mvt /= parameters.matrixRowScaleFactor;
      sequence_score_mvt += tmp_mvt;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT && DISALLOW_FLANKING_TRANSITIONS
      for( uint32_t col_i = 0; col_i <= last_col; col_i++ ) {
        // All paths go through exactly one Deletion or exactly one Match
        // subcell.  Or maybe a DelOut subcell.
#if defined( USE_DEL_IN_DEL_OUT) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        tmp_mvt =
          forward_row[ col_i ][ DeletionIn ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        tmp_mvt *=
          backward_row[ col_i ][ DeletionIn ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        sequence_score_mvt += tmp_mvt;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT && !DISALLOW_FLANKING_TRANSITIONS
        tmp_mvt =
          forward_row[ col_i ][ Match ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        tmp_mvt *=
          backward_row[ col_i ][ Match ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        sequence_score_mvt += tmp_mvt;
        tmp_mvt =
          forward_row[ col_i ][ Deletion ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        tmp_mvt *=
          backward_row[ col_i ][ Deletion ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        sequence_score_mvt += tmp_mvt;
#if defined( USE_DEL_IN_DEL_OUT) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        tmp_mvt =
          forward_row[ col_i ][ DeletionOut ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        tmp_mvt *=
          backward_row[ col_i ][ DeletionOut ];
        tmp_mvt /= parameters.matrixRowScaleFactor;
        sequence_score_mvt += tmp_mvt;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT && !DISALLOW_FLANKING_TRANSITIONS
      } // End foreach col_i
#if defined( USE_DEL_IN_DEL_OUT) && defined( DISALLOW_FLANKING_TRANSITIONS )
      tmp_mvt =
        forward_row.m_deletionOut;
      tmp_mvt /= parameters.matrixRowScaleFactor;
      tmp_mvt *=
        backward_row.m_deletionOut;
      tmp_mvt /= parameters.matrixRowScaleFactor;
      sequence_score_mvt += tmp_mvt;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT && DISALLOW_FLANKING_TRANSITIONS
    } // End if row_i == last_row .. else if row_i == 0 .. else ..last_row

    // TODO: REMOVE
    //cout << "sequence_score_mvt is " << sequence_score_mvt << endl;
    // TODO: REMOVE
    //cout << "forward_row rabinerCumulativeInverseScalar is " << forward_row.m_rabinerCumulativeInverseScalar << endl;
    //cout << "backward_row rabinerCumulativeInverseScalar is " << backward_row.m_rabinerCumulativeInverseScalar << endl;

    // Account for scaling.
    ScoreType sequence_score = sequence_score_mvt;
    if( parameters.useRabinerScaling ) {
      if( row_i > 0 ) {
        sequence_score *= forward_row.m_rabinerCumulativeInverseScalar;
      }
      if( row_i < last_row ) {
        sequence_score *= backward_row.m_rabinerCumulativeInverseScalar;
      }
    }
    // TODO: REMOVE
    //cout << "sequence_score is " << sequence_score << endl;

    return sequence_score;
  } // forward_score( Parameters const&, bool, uint32_t const&, Row &, Row & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  bool
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate the likelihood of the given sequence, given the given profile
     * model, using the given Rows as temporary storage (these will be modified
     * such that afterwards, if the return value is false then forward_row will
     * hold the forward matrix row values for the last row (at index
     * profile.length()), and prev_forward_row will hold the values for the row
     * just before that (at index profile.length()-1)), if the return value is
     * true then these rows will be swapped.  The side-effect, besides changing
     * the forward rows, is to change the value in the score reference
     * argument.  The return value indicates whether the rows are swapped.
     * too.
     */
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      typename Matrix::Row & prev_forward_row,
      typename Matrix::Row & forward_row,
      ScoreType & score
    ) const
    {
      // NOTE: The return value will be true iff profile.length() is odd.
      uint32_t last_row = profile.length();
      uint32_t row_i;

      bool swap = true;
      for( row_i = 0; row_i <= last_row; row_i++ ) {
        // Every other time, we must swap which forward row is current.
        swap = !swap;

        forward_calculateRow(
          parameters,
          profile,
          sequence,
          row_i,
          ( swap ? forward_row : prev_forward_row ),
          ( swap ? prev_forward_row : forward_row )
        );
        if( parameters.useRabinerScaling ) {
          if( row_i == 0 ) {
            ( swap ? prev_forward_row : forward_row ).m_rabinerCumulativeInverseScalar =
              ( swap ? prev_forward_row : forward_row ).m_rabinerInverseScalar;
          } else {
            ( swap ? prev_forward_row : forward_row ).m_rabinerCumulativeInverseScalar =
              (
               ( swap ? forward_row : prev_forward_row ).m_rabinerCumulativeInverseScalar *
               ( swap ? prev_forward_row : forward_row ).m_rabinerInverseScalar
              );
          }
        }
      } // End foreach row.

      if( swap ) {
        score =
          prev_forward_row[ sequence.length() ][ Match ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
        score +=
          prev_forward_row.m_deletionOut;
#else
        score +=
          prev_forward_row[ sequence.length() ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        score *=
          profile[
            Transition::fromPostAlign
          ][
            TransitionFromPostAlign::toTerminal
          ];
        score /=
          parameters.matrixRowScaleFactor;
      } else { // not swap
        score =
          forward_row[ sequence.length() ][ Match ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
        score +=
          forward_row.m_deletionOut;
#else
        score +=
          forward_row[ sequence.length() ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        score *=
          profile[
            Transition::fromPostAlign
          ][
            TransitionFromPostAlign::toTerminal
          ];
        score /=
          parameters.matrixRowScaleFactor;
      } // End if swap .. else ..
      if( parameters.useRabinerScaling ) {
        score *=
          ( swap ? prev_forward_row : forward_row ).m_rabinerCumulativeInverseScalar;
      }
      return swap;
    } // forward_score( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, Row&, Row&, ScoreType & ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.  As a side effect,
     * the given matrices references becomes filled with the forward matrix
     * values, and the given ScoreType vector (if non-NULL) becomes filled
     * with the scores of the sequences, and the largest_score pointer (if
     * non-null) will be filled with the largest of the scores of the
     * sequences.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     */
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer & matrices,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score
    ) const
    {
      sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );

      // TODO: REMOVE
      //cout << "forward_score(..): profile is " << profile << endl;
      static const bool do_extra_debugging = false;//true;

      ScoreType sequence_score;
      ScoreType score;

      typename Matrix::SequentialAccessContainer::iterator prev_matrices_iter =
        matrices.end(); // will be set before it is used.
      typename Matrix::SequentialAccessContainer::iterator matrices_iter;
      uint32_t last_row = profile.length();
      register uint32_t row_i;

      uint32_t last_seq = sequence_count - 1;
      register uint32_t seq_i;

      for(
        row_i = 0, matrices_iter = matrices.begin();
        row_i <= last_row;
        row_i++, prev_matrices_iter = matrices_iter++
      ) {
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          forward_calculateRow(
            parameters,
            use_viterbi,
            profile,
            sequences[ seq_i ],
            row_i,
            ( ( row_i == 0 ) ?
            ( *prev_matrices_iter )[ seq_i ] : // ignored
              ( *prev_matrices_iter )[ seq_i ] ),
            ( *matrices_iter )[ seq_i ]
          );
          // TODO: REMOVE
          //if( seq_i == 0 ) {
          //  cout << "seq 0 forward row " << row_i << ":" << endl;
          //  cout << ( *matrices_iter )[ 0 ] << endl;
          //}
          if( parameters.useRabinerScaling ) {
            // TODO: REMOVE
            //if( seq_i == 0 ) {
            //  cout << "inverse scalar is " << ( *matrices_iter )[ seq_i ].m_rabinerInverseScalar << "! row_i = " << row_i << ", seq_i = " << seq_i << endl;
            //}
            if( do_extra_debugging && isnan( ( *matrices_iter )[ seq_i ].m_rabinerInverseScalar ) ) {
              cout << "inverse scalar is nan! row_i = " << row_i << ", seq_i = " << seq_i << endl;
            }
            if( row_i == 0 ) {
              ( *matrices_iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                ( *matrices_iter )[ seq_i ].m_rabinerInverseScalar;
            } else {
              ( *matrices_iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                ( *prev_matrices_iter )[ seq_i ].m_rabinerCumulativeInverseScalar;
              ( *matrices_iter )[ seq_i ].m_rabinerCumulativeInverseScalar *=
                ( *matrices_iter )[ seq_i ].m_rabinerInverseScalar;
            }
          } // End if useRabinerScaling
        } // End foreach sequence.
      } // End foreach row.

      score = 1.0;
      if( largest_score != NULL ) {
        *largest_score = numeric_limits<double>::min();
      }
      for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
        sequence_score =
          ( *( matrices.rbegin() ) )[ seq_i ][ sequences[ seq_i ].length() ][ Match ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          if( use_viterbi ) {
            if( ( *( matrices.rbegin() ) )[ seq_i ].m_deletionOut > sequence_score ) {
              sequence_score =
                ( *( matrices.rbegin() ) )[ seq_i ].m_deletionOut;
            }
          } else {
            sequence_score +=
              ( *( matrices.rbegin() ) )[ seq_i ].m_deletionOut;
          }

#else
          if( use_viterbi ) {
            if( ( *( matrices.rbegin() ) )[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ] > sequence_score ) {
              sequence_score =
                ( *( matrices.rbegin() ) )[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ];
            }
          } else {
            sequence_score +=
              ( *( matrices.rbegin() ) )[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ];
          }
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        if( parameters.useRabinerScaling ) {
          sequence_score *=
            ( *( matrices.rbegin() ) )[ seq_i ].m_rabinerCumulativeInverseScalar;
        } // End if useRabinerScaling
        sequence_score *=
          profile[
            Transition::fromPostAlign
          ][
            TransitionFromPostAlign::toTerminal
          ];
        sequence_score /=
          parameters.matrixRowScaleFactor;

        // TODO: REMOVE
        //cout << "SEQUENCE " << seq_i << " SCORE IS " << sequence_score << endl;
        //cout << "\t as double:" << toDouble( sequence_score ) << endl;

        // TODO: Debug mode only...
        if( do_extra_debugging && ( isinf( sequence_score ) || isnan( sequence_score ) ) ) {
          cout << "SEQUENCE " << seq_i << " SCORE IS " << sequence_score << endl;
          cout << "The sequence is:" << sequences[ seq_i ] << endl;
          cout << "The last row of the forward matrix for this sequence is " << endl << ( *( matrices.rbegin() ) )[ seq_i ] << endl;
          cout << "The m_rabinerCumulativeInverseScalar value is " <<
            ( *( matrices.rbegin() ) )[ seq_i ].m_rabinerCumulativeInverseScalar << endl;
          //cout << "PROFILE IS:" << endl << profile << endl;
          cout << "( *( matrices.rbegin() ) )[ seq_i ][ sequences[ seq_i ].length() ][ Match ] returns " << ( *( matrices.rbegin() ) )[ seq_i ][ sequences[ seq_i ].length() ][ Match ] << endl;
          cout << "(\
                profile[\
                  Transition::fromPostAlign\
                ][\
                  TransitionFromPostAlign::toTerminal\
                ]\
              ) returns " << (
                profile[
                  Transition::fromPostAlign
                ][
                  TransitionFromPostAlign::toTerminal
                ]
                                                              ) << endl;
          exit( 0 );
        }
        score *= sequence_score;
        if( sequence_scores != NULL ) {
          ( *sequence_scores )[ seq_i ] = sequence_score;
        } // End if sequence_scores != NULL
        if( ( largest_score != NULL ) && ( sequence_score > *largest_score ) ) {
          *largest_score = sequence_score;
        }
      } // End foreach sequence
      return score;
    } // forward_score( Parameters const&, bool, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer &, vector<ScoreType> *, ScoreType * ) const



  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the forward algorithm.
     */
    forward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count
    ) const
    {
      sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );

      typename Matrix::RowVector prev_forward_rows =
        typename Matrix::RowVector( sequences, sequence_count );
      typename Matrix::RowVector forward_rows =
        typename Matrix::RowVector( sequences, sequence_count );

      // Initialize the rows.
      uint32_t last_seq = sequence_count - 1;
      for( uint32_t seq_i = 0; seq_i <= last_seq; seq_i++ ) {
        prev_forward_rows[ seq_i ] =
          typename Matrix::Row( sequences[ seq_i ].length() + 1 );
        forward_rows[ seq_i ] =
          typename Matrix::Row( sequences[ seq_i ].length() + 1 );
      }

      ScoreType score;
      forward_score(
        parameters,
        false, // don't use viterbi
        profile,
        sequences,
        sequence_count,
        prev_forward_rows,
        forward_rows,
        score // changes as a side-effect
      );
      return score;
    } // forward_score( Parameters const&, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const& ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    template <typename ProfileType,
              typename SequenceResidueType,
              typename ScaledMatchDistributionProbabilityType>
#else
    template <typename ProfileType,
              typename SequenceResidueType>
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
  GALOSH_INLINE_ALGORITHM
  bool
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate the total likelihood of the given sequences, given the given
     * profile model, using the given Rows as temporary storage (these will be
     * modified such that afterwards, if the return value is false then
     * forward_rows will hold the forward matrix row values for the last rows
     * (at index profile.length()), and prev_forward_rows will hold the values
     * for the rows just before that (at index profile.length()-1)).  If the
     * return value is true then these rows will be swapped.  One side-effect,
     * besides changing the forward rows, is to change the value in the score
     * reference argument.  The return value indicates whether the rows are
     * swapped. If the given ScoreType vector is non-NULL, it will become
     * filled with the scores of the sequences, and the largest_score pointer
     * (if non-null) will be filled with the largest of the scores of the
     * sequences.  If parameters.useRabinerScaling is true, the
     * rabinerInverseScaling and rabinerCumulativeInverseScalar values of the
     * forward rows will be set, too.  If the anchor_columns pointer is
     * non-null, then anchor columns (to be used with
     * forward_reverseCalculateRow(..) will be stored also (using the member
     * var anchor_columns->m_storeEveryNthColumn).  If the anchor_rows pointer
     * is non-null, then anchor rows will be stored also (using the member var
     * anchor_rows->m_storeEveryNthRow).  If the use_viterbi argument is true,
     * calculate most-likely-path likelihoods rather than full (all path)
     * likelihoods.
     */
    forward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      vector<MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> > const & scaled_match_distributions,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer * anchor_columns,
      typename Matrix::SequentialAccessContainer * anchor_rows,
      typename Matrix::RowVector & prev_forward_rows,
      typename Matrix::RowVector & forward_rows,
      vector<ScoreType> * sequence_scores,
      ScoreType * largest_score,
      ScoreType & score
    ) const
    {
      sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );
      // NOTE: The return value will be true iff profile.length() is odd.

      uint32_t last_row = profile.length();
      uint32_t row_i;
      uint32_t last_seq = sequence_count - 1;
      uint32_t seq_i;

      typename Matrix::SequentialAccessContainer::iterator anchor_rows_iter;
      typename Matrix::SequentialAccessContainer::iterator anchor_columns_iter;
      if(
        ( anchor_rows != 0 ) &&
        ( anchor_rows->m_storeEveryNthRow != 0 )
      ) {
        anchor_rows_iter = anchor_rows->begin();
      }
      if(
        ( anchor_columns != 0 ) &&
        ( anchor_columns->m_storeEveryNthColumn != 0 )
      ) {
        anchor_columns_iter = anchor_columns->begin();
      }

      if( largest_score != 0 ) {
        *largest_score = numeric_limits<double>::min();
      }

      bool swap = true;
      for( row_i = 0; row_i <= last_row; row_i++ ) {
        if(
          ( anchor_rows != 0 ) &&
          ( anchor_rows->m_storeEveryNthRow != 0 ) &&
          ( row_i != 0 ) &&
          ( ( row_i % anchor_rows->m_storeEveryNthRow ) == 0 )
        ) {
          anchor_rows_iter++;
          assert( anchor_rows_iter != anchor_rows->end() );
        } // End if storing anchor rows too and this is an anchor row
        if(
          ( anchor_columns != 0 ) &&
          ( row_i != 0 ) &&
          ( anchor_columns->m_storeEveryNthColumn != 0 )
        ) {
          anchor_columns_iter++;
        }

        // Every other time, we must swap which forward row is current.
        swap = !swap;

        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          forward_calculateRow(
            parameters,
            use_viterbi,
            profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            (
              ( ( row_i <= 1 ) || ( ( row_i - 2 ) >= ( last_row - 1 ) ) ) ?
              scaled_match_distributions[ row_i ] : // ignored
              scaled_match_distributions[ row_i - 2 ]
            ),
            (
              ( ( row_i == 0 ) || ( ( row_i - 1 ) >= ( last_row - 1 ) ) ) ?
              scaled_match_distributions[ row_i ] : // ignored
              scaled_match_distributions[ row_i - 1 ]
            ),
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
            sequences[ seq_i ],
            row_i,
            ( swap ? forward_rows[ seq_i ] : prev_forward_rows[ seq_i ] ),
            ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] )
          );
          if(
            ( anchor_columns != 0 ) &&
            ( anchor_columns->m_storeEveryNthColumn != 0 )
          ) {
            // This won't copy the rabiner scalars -- see below for that.
            ( *anchor_columns_iter )[ seq_i ].storeColumns(
              ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ),
              anchor_columns->m_storeEveryNthColumn
            );
          } // End if storing anchor columns too.
          if(
            ( anchor_rows != 0 ) &&
            ( anchor_rows->m_storeEveryNthRow != 0 ) &&
            ( ( row_i % anchor_rows->m_storeEveryNthRow ) == 0 )
          ) {
            // This won't copy rabiner scalars.  See below for that.
            ( *anchor_rows_iter )[ seq_i ] =
              ( swap ? prev_forward_rows : forward_rows )[ seq_i ];
          } // End if storing anchor rows too and this is an anchor row

          if( parameters.useRabinerScaling ) {
            if( row_i == 0 ) {
              ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ).m_rabinerCumulativeInverseScalar =
                ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ).m_rabinerInverseScalar;
            } else {
              ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ).m_rabinerCumulativeInverseScalar =
                (
                 ( swap ? forward_rows[ seq_i ] : prev_forward_rows[ seq_i ] ).m_rabinerCumulativeInverseScalar *
                  ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ).m_rabinerInverseScalar
                );
            }
            if( ( anchor_columns != 0 ) && ( anchor_columns->m_storeEveryNthColumn != 0 ) ) {
              ( *anchor_columns_iter )[ seq_i ].m_rabinerInverseScalar =
                ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ).m_rabinerInverseScalar;
              ( *anchor_columns_iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ).m_rabinerCumulativeInverseScalar;
            } // End if anchor_columns != 0;
            if(
              ( anchor_rows != 0 ) &&
              ( anchor_rows->m_storeEveryNthRow != 0 ) &&
              ( ( row_i % anchor_rows->m_storeEveryNthRow ) == 0 )
            ) {
              ( *anchor_rows_iter )[ seq_i ].m_rabinerInverseScalar =
                ( swap ? prev_forward_rows : forward_rows )[ seq_i ].m_rabinerInverseScalar;
              ( *anchor_rows_iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                ( swap ? prev_forward_rows : forward_rows )[ seq_i ].m_rabinerCumulativeInverseScalar;
            } // End if storing anchor rows too and this is an anchor row
          } // End if useRabinerScaling

          // TODO: REMOVE
          //if( seq_i == 1 ) {
          //  if( row_i == 1 ) {
          //    cout << "forward row 1 for seq 1 is " << ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ) << endl;
          //  }
          //  if( row_i == 2 ) {
          //    cout << "forward row 2 for seq 1 is " << ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ) << endl;
          //  }
          //  if( row_i == 3 ) {
          //    cout << "forward row 3 for seq 1 is " << ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ) << endl;
          //  }
          //  if( row_i == 4 ) {
          //    cout << "forward row 4 for seq 1 is " << ( swap ? prev_forward_rows[ seq_i ] : forward_rows[ seq_i ] ) << endl;
          //  }
          //}

        } // End foreach sequence.

      } // End foreach row.

      score = 1.0;
      ScoreType sequence_score;
      for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
        if( swap ) {
          sequence_score =
            prev_forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ Match ];
          // TODO: REMOVE
          //cout << "[swap] seq " << seq_i << " Match: " << sequence_score << endl;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          if( use_viterbi ) {
            if( prev_forward_rows[ seq_i ].m_deletionOut > sequence_score ) {
              sequence_score =
                prev_forward_rows[ seq_i ].m_deletionOut;
            }
          } else {
            sequence_score +=
              prev_forward_rows[ seq_i ].m_deletionOut;
          }

#else
          if( use_viterbi ) {
            if( prev_forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ] > sequence_score ) {
              sequence_score =
                prev_forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ];
            }
          } else {
            sequence_score +=
              prev_forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ];
          }
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          // TODO: REMOVE
          //cout << "[swap] seq " << seq_i << " Match+DelOut: " << sequence_score << endl;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          sequence_score *=
            profile[
              Transition::fromPostAlign
            ][
              TransitionFromPostAlign::toTerminal
            ];
          // TODO: REMOVE
          //cout << "[swap] seq " << seq_i << " after C->T: " << sequence_score << endl;
          sequence_score /=
            parameters.matrixRowScaleFactor;
          // TODO: REMOVE
          //cout << "[swap] seq " << seq_i << " after /= matrixRowScaleFactor: " << sequence_score << endl;
          if( parameters.useRabinerScaling ) {
            sequence_score *=
              prev_forward_rows[ seq_i ].m_rabinerCumulativeInverseScalar;
          }
          // TODO: REMOVE
          //cout << "[swap] seq " << seq_i << " score: " << sequence_score << endl;
        } else { // not swap
          sequence_score =
            forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ Match ];
          // TODO: REMOVE
          //cout << "seq " << seq_i << " Match: " << sequence_score << endl;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          if( use_viterbi ) {
            if( forward_rows[ seq_i ].m_deletionOut > sequence_score ) {
              sequence_score =
                forward_rows[ seq_i ].m_deletionOut;
            }
          } else {
            sequence_score +=
              forward_rows[ seq_i ].m_deletionOut;
          }

#else
          if( use_viterbi ) {
            if( forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ] > sequence_score ) {
              sequence_score =
                forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ];
            }
          } else {
            sequence_score +=
              forward_rows[ seq_i ][ sequences[ seq_i ].length() ][ DeletionOut ];
          }
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          // TODO: REMOVE
          //cout << "seq " << seq_i << " Match+DelOut: " << sequence_score << endl;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          sequence_score *=
            profile[
              Transition::fromPostAlign
            ][
              TransitionFromPostAlign::toTerminal
            ];
          // TODO: REMOVE
          //cout << "seq " << seq_i << " after C->T: " << sequence_score << endl;
          sequence_score /=
            parameters.matrixRowScaleFactor;
          // TODO: REMOVE
          //cout << "seq " << seq_i << " after /= matrixRowScaleFactor: " << sequence_score << endl;
          if( parameters.useRabinerScaling ) {
            sequence_score *=
              forward_rows[ seq_i ].m_rabinerCumulativeInverseScalar;
          }
          // TODO: REMOVE
          //cout << "seq " << seq_i << " score: " << sequence_score << endl;
        } // End if swap .. else ..
        if( sequence_scores != 0 ) {
          ( *sequence_scores )[ seq_i ] = sequence_score;
        }
        if( ( largest_score != 0 ) && ( sequence_score > *largest_score ) ) {
          *largest_score = sequence_score;
        }
        score *= sequence_score;
        // TODO: REMOVE.
        //cout << "forward_calculateScore(..): seq " << seq_i << "'s score is " << sequence_score << endl;
        //cout << "forward_calculateScore(..): after seq " << seq_i << ", the total score is " << score << endl;
      } // End foreach sequence...

      // TODO: REMOVE!
      //exit( 0 );

      return swap;
    } // forward_score( Parameters const&, bool, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer *, uint32_t const &, RowVector &, RowVector &, ScoreType & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
     /**
     * Using the forward row (of the Forward algorithm matrix) for (row_i - 1),
     * calculate the forward row for row_i.  Return the given forward_row
     * reference.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     * If parameters.useRabinerScaling is true, then
     * the calculated forward row will be scaled by the sum of all of its
     * values (or by the max value if parameters.rabinerScaling_useMaximumValue
     * is true).
     *
     * NOTE: when row_i is last_row, forward_row.m_rabinerInverseScalar is the
     * last Match value times the postalign-to-terminal probability.
     * NOTE: when row_i is 0, prev_forward_row is ignored.
     * NOTE: For now this requires that the profile length be > 1.
     */
    forward_calculateRow (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & prev_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& prev_forward_row,
      typename Matrix::Row & forward_row
    ) const
    {
      const bool do_extra_debugging = false;//( row_i >= 3 );//false;
      const bool be_extra_verbose = false;//( row_i >= 3 );//false;

      const uint32_t last_row = profile.length();
      const uint32_t last_col = sequence.length();

      assert( last_row > 1 ); // TODO: Implement short-profile stuff.

      // TODO: REMOVE?  Testing.
      if( row_i > 0 ) {
        assert( profile[ row_i - 1 ].getProfileTreeRoot() == &profile );
      }

      const bool use_rabiner_scaling = parameters.useRabinerScaling;

      if( be_extra_verbose || ( parameters.debug == DEBUG_Special ) ) {
        cout << "row_i: " << row_i << endl;
        cout << "SEQUENCE: " << endl;
        cout << sequence << endl;
      }

      register uint32_t col_i;

      MatrixValueType col_i_match_new, col_i_insertion_new, col_i_deletion_new;
      MatrixValueType col_i_deletion_in_new, col_i_deletion_out_new;
      MatrixValueType tmp_new;

      if( do_extra_debugging ) {
        cout << "profile globals are ";
        profile.writeExceptPositions( cout );
        cout << endl;
        if( ( row_i >= 2 ) && ( ( row_i - 2 ) < ( profile.length() - 1 ) ) ) {
          // TODO: REMOVE
            cout << "profile[ row_i - 2 ] is " << profile[ row_i - 2 ] << endl;
        } // End if row_i >= 2
        if( ( row_i >= 1 ) && ( ( row_i - 1 ) < ( profile.length() - 1 ) ) ) {
          // TODO: REMOVE
            cout << "profile[ row_i - 1 ] is " << profile[ row_i - 1 ] << endl;
        } // End if row_i >= 1
      } // End if do_extra_debugging
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      // We use scaled versions of the TransitionFromMatch distributions.
      if( ( row_i >= 2 ) && ( ( row_i - 2 ) < ( profile.length() - 1 ) ) ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using prev_row_match_distribution == " << prev_row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if row_i >= 2
      if( ( row_i >= 1 ) && ( ( row_i - 1 ) < ( profile.length() - 1 ) ) ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using row_match_distribution == " << row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if row_i >= 1
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      MatrixValueType rabiner_inverse_scalar_mvt;
      if( use_rabiner_scaling ) {
        rabiner_inverse_scalar_mvt = 1.0;
      } // End if use_rabiner_scaling

      if( be_extra_verbose || ( parameters.debug == DEBUG_Special ) ) {
        cout << "[forward] last_row is " << last_row << ", last_col is " << last_col << endl;
      }

#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
        forward_row.m_deletionIn = 0;
        forward_row.m_deletionOut = 0;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )
      for( col_i = 0; col_i <= last_col; col_i++ ) {

        col_i_match_new = 0;
        col_i_insertion_new = 0;
        col_i_deletion_new = 0;
        col_i_deletion_in_new = 0;
        col_i_deletion_out_new = 0;
        tmp_new = 0;

        // TODO: REMOVE
        //if( parameters.debug == DEBUG_Special ) {
        //  cout << "[forward] row " << row_i << ", col " << col_i << endl;
        //}
        
        // calculate row cell:
        forward_row[ col_i ][ Match ] = 0;
        forward_row[ col_i ][ Insertion ] = 0;
        forward_row[ col_i ][ Deletion ] = 0;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        forward_row[ col_i ][ DeletionIn ] = 0;
        forward_row[ col_i ][ DeletionOut ] = 0;
#endif // ( USE_DEL_IN_DEL_OUT && !DISALLOW_FLANKING_TRANSITIONS )

#if defined( USE_DEL_IN_DEL_OUT )
        // from DeletionIn
        // from DeletionIn to Match
        if( row_i > 1 ) {
          if(
            ( col_i >= 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            && ( col_i == 1 )
#endif // DISALLOW_FLANKING_TRANSITIONS
          ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_SWENTRY_SWEXIT
            tmp_new =
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            tmp_new *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletionIn
              ];
            tmp_new *= ( 1.0 / ( profile.length() - 1 ) );
#else // USE_SWENTRY_SWEXIT ( & DISALLOW_FLANKING_TRANSITIONS ) .. else ..
            tmp_new = prev_forward_row.m_deletionIn;
#endif // USE_SWENTRY_SWEXIT ( & DISALLOW_FLANKING_TRANSITIONS ) .. else ..
#else
            // TODO: Implement support for USE_SWENTRY_SWEXIT when !DISALLOW_FLANKING_TRANSITIONS
            tmp_new = prev_forward_row[ col_i - 1 ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else .. 
#ifndef USE_SWENTRY_SWEXIT
            tmp_new *=
              profile[
                Transition::fromDeletionIn
              ][
                TransitionFromDeletionIn::toMatch
              ];
#endif // !USE_SWENTRY_SWEXIT
            tmp_new *=
              getEmissionProbability(
                profile[ row_i - 1 ],
                sequence,
                col_i - 1
              );
#ifdef USE_END_DISTRIBUTION
            if( row_i == last_row ) {
              tmp_new *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION

            // TODO: REMOVE
            if( do_extra_debugging ) {
              if( row_i == last_row ) {
                cout << "col " << col_i << " Z->M->C Match contribution: " << tmp_new << endl;
              } else {
                cout << "col " << col_i << " Z->M Match contribution: " << tmp_new << endl;
              }
            } // End if do_extra_debugging
            if( use_viterbi ) {
              if( tmp_new > col_i_match_new ) {
                col_i_match_new = tmp_new;
              }
            } else {
              col_i_match_new += tmp_new;
            }
          } // End if col_i > 0 (and maybe col_i == 1)
        } // End if row_i > 1

#ifndef USE_SWENTRY_SWEXIT
        // from DeletionIn to DeletionIn
        if( ( row_i == 0 ) || ( row_i == last_row ) ) {
          // Note that deletion-ins must end with a Match into some position,
          // so in the last row it's always 0.

          // Then deletion in is just 0.
          tmp_new = 0;
        } else if(
          ( row_i == 1 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // In other columns we start with the pre-aligns.
          tmp_new =
            prev_forward_row[ col_i ][ Match ];
          tmp_new *=
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ];
          tmp_new *=
            profile[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletionIn
            ];
        } else // row_i > 1
#ifdef DISALLOW_FLANKING_TRANSITIONS
          if( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        {
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_new =
            prev_forward_row.m_deletionIn;
#else
          tmp_new =
            prev_forward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else .. 
          tmp_new *=
            profile[
              Transition::fromDeletionIn
            ][
              TransitionFromDeletionIn::toDeletionIn
            ];
        } // End if row_i == 0 .. else 1 .. else ..
        // Don't worry about viterbi, since there's only one path to del-in in
        // this cell.
#ifdef DISALLOW_FLANKING_TRANSITIONS
        if( col_i == 0 ) {
          forward_row.m_deletionIn = tmp_new;
        }
#else
        col_i_deletion_in_new = tmp_new;
#endif // DISALLOW_FLANKING_TRANSITIONS .. else .. 
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

        // from Match
        // from Match to Match
        if( col_i == 0 ) {
          if( row_i == 0 ) {
            col_i_match_new =
              parameters.matrixRowScaleFactor;
          }        
        } else if( row_i == 0 ) { // col_i > 0 && row_i == 0
          // In the top row, we use preAlignInsertions
          tmp_new =
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toPreAlign
            ];
          // TODO: REMOVE
          //cout << "about to call getPreAlignEmissionProbability(" << endl;
          //cout << profile << ", " << endl;
          //cout << sequence << ", " << endl;
          //cout << "( " << col_i << " - 1 ), " << endl;
          //cout << ")" << endl;
          //cout << "sequence[ col_i - 1 ] is " << sequence[ col_i - 1 ] << endl;
          tmp_new *=
            getPreAlignEmissionProbability(
              profile,
              sequence,
              col_i - 1
            );
          tmp_new *=
              // In the top row we use preAlignInsertions, which are not
              // affine, so we store the values in the Match cell instead of
              // the Insertion cell.
            forward_row[ col_i - 1 ][ Match ];

          // TODO: REMOVE.
          if( do_extra_debugging && isnan( tmp_new ) ) {
            cout << "tmp_new just became nan = ( "
                 << profile[
                      Transition::fromPreAlign
                    ][
                      TransitionFromPreAlign::toPreAlign
                    ] << " * "
                 << getPreAlignEmissionProbability(
                      profile,
                      sequence,
                      col_i - 1
                    ) << " * "
                 << forward_row[ col_i - 1 ][ Match ] << " )" << endl;
          } // End if isnan( tmp_new )
          // TODO: REMOVE
          if( do_extra_debugging ) {
            cout << "col " << col_i << " N->N Match contribution: " << tmp_new << endl;
          } // End if do_extra_debugging
          if( use_viterbi ) {
            if( tmp_new > col_i_match_new ) {
              col_i_match_new = tmp_new;
            }
          } else {
            col_i_match_new += tmp_new;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
              cout << "forward_calculateRow: *1* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
            }
            
          }
        } else { // col_i > 0, row_i > 0
          if( row_i == last_row ) { // col_i > 0 && row_i == last_row
            // In the bottom row, we use postAlignInsertions.  So this is
            // really I->I (postalign, so technically C->C), stored in M
            // subcell.

            // In the bottom row we use postAlignInsertions, which are not
            // affine, so we store the values in the Match cell instead of
            // the Insertion cell.
            tmp_new =
              forward_row[ col_i - 1 ][ Match ];

            tmp_new *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ];
            tmp_new *=
              getPostAlignEmissionProbability(
                profile,
                sequence,
                col_i - 1
              );

            // TODO: REMOVE
            if( do_extra_debugging ) {
              cout << "col " << col_i << " C->C Match contribution, Match state: " << tmp_new << endl;
            } // End if do_extra_debugging
            if( use_viterbi ) {
              if( tmp_new > col_i_match_new ) {
                col_i_match_new = tmp_new;
              }
            } else {
              col_i_match_new += tmp_new;
              // TODO: REMOVE
              if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
                cout << "forward_calculateRow: *2* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              }
            }
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            // Then we need to include paths that deleted out early.  Note that
            // if there's no flanking transitions, then this is all 0 anyway.
#ifndef DISALLOW_FLANKING_TRANSITIONS
            tmp_new =
              forward_row[ col_i - 1 ][ DeletionOut ];
            tmp_new *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ];
            tmp_new *=
              getPostAlignEmissionProbability(
                profile,
                sequence,
                col_i - 1
              );

            // TODO: REMOVE
            if( do_extra_debugging ) {
              cout << "col " << col_i << " C->C Match contribution, DeletionOut state: " << tmp_new << endl;
            } // End if do_extra_debugging
            if( use_viterbi ) {
              if( tmp_new > col_i_match_new ) {
                col_i_match_new = tmp_new;
              }
            } else {
              col_i_match_new += tmp_new;
              // TODO: REMOVE
              if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
                cout << "forward_calculateRow: *2b* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              }
            }
#endif // !DISALLOW_FLANKING_TRANSITIONS
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          } // End if it's the last row..

          if( row_i == 1 ) {
            tmp_new =
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            tmp_new *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toMatch
              ];
          } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            // Use the scaled match distribution
            tmp_new =
              prev_row_match_distribution[
                TransitionFromMatch::toMatch
              ];
#else
            tmp_new =
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toMatch
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          }
          tmp_new *=
            getEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              col_i - 1
            );
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            tmp_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_new *=
            prev_forward_row[ col_i - 1 ][ Match ];

          // TODO: REMOVE
          if( do_extra_debugging ) {
            if( row_i == 1 ) {
              cout << "col " << col_i << " N->M Match contribution: " << tmp_new << endl;
            } else if( row_i == last_row ) {
              cout << "col " << col_i << " M->M->C Match contribution: " << tmp_new << endl;
            } else {
              cout << "col " << col_i << " M->M Match contribution: " << tmp_new << endl;
            }
          } // End if do_extra_debugging
          if( use_viterbi ) {
            if( tmp_new > col_i_match_new ) {
              col_i_match_new = tmp_new;
            }
          } else {
            col_i_match_new += tmp_new;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
              cout << "forward_calculateRow: *3* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              cout << "tmp_new = ( "
                << ( ( row_i == 1 ) ?
                (
                  profile[
                    Transition::fromPreAlign
                  ][
                    TransitionFromPreAlign::toBegin
                  ] *
                  profile[
                    Transition::fromBegin
                  ][
                    TransitionFromBegin::toMatch
                  ]
                ) :
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                // Use the scaled match distribution
                static_cast<ProbabilityType>( prev_row_match_distribution[
                                                               TransitionFromMatch::toMatch
                                                            ] )
#else
                profile[ row_i - 2 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toMatch
                ]
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                     ) << " * " << 
             getEmissionProbability(
               profile[ row_i - 1 ],
               sequence,
               col_i - 1
                                    ) << " * " <<
             ( ( row_i == last_row ) ?
               (
#ifdef USE_END_DISTRIBUTION
                 ( MatrixValueType )profile[
                   Transition::fromEnd
                 ][
                   TransitionFromEnd::toPostAlign
                 ] *
#endif // USE_END_DISTRIBUTION
                 prev_forward_row[ col_i - 1 ][ Match ]
               ) :
               prev_forward_row[ col_i - 1 ][ Match ] ) << " );" << endl;
              cout << "profile[ " << ( row_i - 1 ) << " ] is: " <<
                profile[ row_i - 1 ] << endl;
            }
            
          }
        } // End if col_i == 0, elsif row_i == 0, else ..
        // from Match to Insertion
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != last_row ) ) {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          // Use the scaled match distribution
          tmp_new =
            row_match_distribution[
              TransitionFromMatch::toInsertion
            ];
#else
          tmp_new =
            profile[ row_i - 1 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          tmp_new *=
            getInsertionEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              col_i - 1
            );
          tmp_new *= forward_row[ col_i - 1 ][ Match ];

          // TODO: REMOVE
          //cout << "M->I tmp_new is " << tmp_new << "." << endl;

          if( use_viterbi ) {
            if( tmp_new > col_i_insertion_new ) {
              col_i_insertion_new = tmp_new;
            }
          } else {
            col_i_insertion_new += tmp_new;
            // TODO: REMOVE
            //cout << "col_i_insertion_new is now " << col_i_insertion_new << "." << endl;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_insertion_new ) || isinf( col_i_insertion_new ) || ( col_i_insertion_new < 0 ) ) ) {
              cout << "forward_calculateRow: *iM* col_i_insertion_new has just become " << col_i_insertion_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              cout << "tmp_new = ( " <<
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                row_match_distribution[
                  TransitionFromMatch::toInsertion
                ] <<
#else
                profile[ row_i - 1 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toInsertion
                ] <<
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                " * " << 
                getInsertionEmissionProbability(
                  profile[ row_i - 1 ],
                  sequence,
                  col_i - 1
                ) << " = " <<
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                ( row_match_distribution[
                  TransitionFromMatch::toInsertion
                ] *
#else
                ( profile[ row_i - 1 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toInsertion
                ] *
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                getInsertionEmissionProbability(
                  profile[ row_i - 1 ],
                  sequence,
                  col_i - 1
                ) ) << " ) * " <<
                forward_row[ col_i - 1 ][ Match ] << endl;
            }
          }
        } // End if this is not the first column and not the first or last row ..
        // from Match to Deletion
        if( row_i != 0 ) {
          if( row_i == 1 ) {
            tmp_new =
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            tmp_new *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletion
              ];
            tmp_new *=
              prev_forward_row[ col_i ][ Match ];

            if( use_viterbi ) {
              if( tmp_new > col_i_deletion_new ) {
                col_i_deletion_new = tmp_new;
              }
            } else {
              col_i_deletion_new += tmp_new;
            }

          } else if( row_i == last_row ) {
            // We store the postAlignInsertion stuff in the Match state.
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            // Use the scaled match distribution
            tmp_new =
              prev_row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            tmp_new =
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
#ifdef USE_END_DISTRIBUTION
            tmp_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            tmp_new *=
              prev_forward_row[ col_i ][ Match ];

            // TODO: REMOVE
            if( do_extra_debugging ) {
              cout << "col " << col_i << " M->D->C Match contribution: " << tmp_new << endl;
            } // End if do_extra_debugging
            if( use_viterbi ) {
              if( tmp_new > col_i_match_new ) {
                col_i_match_new = tmp_new;
              }
            } else {
              col_i_match_new += tmp_new;
              // TODO: REMOVE
              if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
                cout << "forward_calculateRow: *4* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              }
            
            }
          } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            // Use the scaled match distribution
            tmp_new =
              prev_row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            tmp_new =
              profile[ row_i - 2 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            tmp_new *=
              prev_forward_row[ col_i ][ Match ];

            if( use_viterbi ) {
              if( tmp_new > col_i_deletion_new ) {
                col_i_deletion_new = tmp_new;
              }
            } else {
              col_i_deletion_new += tmp_new;
              // TODO: REMOVE
              if( do_extra_debugging && ( isnan( col_i_deletion_new ) || isinf( col_i_deletion_new ) || ( col_i_deletion_new < 0 ) ) ) {
                cout << "forward_calculateRow: *dM* col_i_deletion_new has just become " << col_i_deletion_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
                cout << "tmp_new = ( " <<
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
                  prev_row_match_distribution[
                    TransitionFromMatch::toDeletion
                  ] <<
#else
                  profile[ row_i - 2 ][
                    Transition::fromMatch
                  ][
                    TransitionFromMatch::toDeletion
                  ] <<
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
                  " * " << 
                  prev_forward_row[ col_i ][ Match ] << " ) " << endl;
              }
            }
          } // End if it's row 1 .. elsif the last row .. else ..
        } // End if it's not the first row ..
#if defined( USE_DEL_IN_DEL_OUT )
        // From Match to DeletionOut
        if(
          ( row_i > 1 ) && ( col_i != 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // Note that we don't use the scaled disribution, since we just want
          // to include the part of the path that ends here, and the scaled
          // distribution includes the entire path.

          tmp_new =
            profile[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletionOut
            ];
          tmp_new *=
            prev_forward_row[ col_i ][ Match ];

          if( row_i == last_row ) {
#ifndef USE_SWENTRY_SWEXIT
            tmp_new *=
              profile[
                Transition::fromDeletionOut
              ][
                TransitionFromDeletionOut::toEnd
              ];
#endif // USE_SWENTRY_SWEXIT
#ifdef USE_END_DISTRIBUTION
            tmp_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
          } // End if( row_i == last_row )

#ifndef USE_SWENTRY_SWEXIT
#ifdef DISALLOW_FLANKING_TRANSITIONS
          // TODO: REMOVE
          assert( col_i == last_col );
          forward_row.m_deletionOut = tmp_new;
#else
          col_i_deletion_out_new = tmp_new;
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..

          // from DeletionOut to DeletionOut
          if( row_i > 2 ) {
            tmp_new =
              profile[
                Transition::fromDeletionOut
              ][
                TransitionFromDeletionOut::toDeletionOut
              ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
            tmp_new *= prev_forward_row.m_deletionOut;
#else
            tmp_new *= prev_forward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
            if( row_i == last_row ) {
              tmp_new *=
                profile[
                  Transition::fromDeletionOut
                ][
                  TransitionFromDeletionOut::toEnd
                ];
#ifdef USE_END_DISTRIBUTION
              tmp_new *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
#endif // USE_END_DISTRIBUTION
            } // End if( row_i == last_row )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            if( use_viterbi ) {
              if( tmp_new > forward_row.m_deletionOut ) {
                forward_row.m_deletionOut = tmp_new;
              }
            } else {
              forward_row.m_deletionOut += tmp_new;
            }
#else
            if( use_viterbi ) {
              if( tmp_new > col_i_deletion_out_new ) {
                col_i_deletion_out_new = tmp_new;
              }
            } else {
              col_i_deletion_out_new += tmp_new;
            }
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          } // End if ( row_i > 2 )
#endif // !USE_SWENTRY_SWEXIT
        } // End if ( row_i > 1 ) && ( col_i != 0 ) ( and maybe col_i == last_col )
#endif // USE_DEL_IN_DEL_OUT
        // End from Match

        // from Insertion
        // from Insertion to Match
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != 1 ) ) { // No ins state in row 0
          tmp_new =
            profile[ row_i - 2 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toMatch
            ];
          tmp_new *=
            getEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              col_i - 1
            );
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            tmp_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_new *=
            prev_forward_row[ col_i - 1 ][ Insertion ];

          // TODO: REMOVE
          if( do_extra_debugging ) {
            if( row_i == last_row ) {
              cout << "col " << col_i << " I->M->C Match contribution: " << tmp_new << endl;
            } else {
              cout << "col " << col_i << " I->M Match contribution: " << tmp_new << endl;
            }
          } // End if do_extra_debugging
          if( use_viterbi ) {
            if( tmp_new > col_i_match_new ) {
              col_i_match_new = tmp_new;
            }
          } else {
            col_i_match_new += tmp_new;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
              cout << "forward_calculateRow: *5* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              cout << "tmp_new = ( " <<
                profile[ row_i - 2 ][
                  Transition::fromInsertion
                  ][
                  TransitionFromInsertion::toMatch
                  ] << " * " << 
                getEmissionProbability(
                  profile[ row_i - 1 ],
                  sequence,
                  col_i - 1
                ) << " ) * ( ( ( row_i == last_row ) = " <<
                ( row_i == last_row ) << " ) ? ( " << 
#ifdef USE_END_DISTRIBUTION
                   profile[
                      Transition::fromEnd
                    ][
                      TransitionFromEnd::toPostAlign
                    ]
#else
                1.0
#endif // USE_END_DISTRIBUTION
                << " ) * " <<
                  prev_forward_row[ col_i - 1 ][ Insertion ]
                   << " ) : " <<
                prev_forward_row[ col_i - 1 ][ Insertion ] << " ) )" << endl;
            }
            
          }
        } // End if it's not the first or second row and not the first col ..
        // from Insertion to Insertion
        if( ( row_i != 0 ) &&
            ( row_i != last_row ) &&
            ( col_i != 0 ) ) {
          tmp_new =
            profile[ row_i - 1 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toInsertion
            ];
          tmp_new *=
            getInsertionEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              col_i - 1
            );
          tmp_new *=
            forward_row[ col_i - 1 ][ Insertion ];

          // TODO: REMOVE
          //cout << "I->I tmp_new is " << tmp_new << "." << endl;

          if( use_viterbi ) {
            if( tmp_new > col_i_insertion_new ) {
              col_i_insertion_new = tmp_new;
            }
          } else {
            col_i_insertion_new += tmp_new;
            // TODO: REMOVE
            //cout << "col_i_insertion_new is now " << col_i_insertion_new << "." << endl;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_insertion_new ) || isinf( col_i_insertion_new ) || ( col_i_insertion_new < 0 ) ) ) {
              cout << "forward_calculateRow: *iI* col_i_insertion_new has just become " << col_i_insertion_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              tmp_new =
                (
               profile[ row_i - 1 ][
                  Transition::fromInsertion
                ][
                  TransitionFromInsertion::toInsertion
                ] *
                getInsertionEmissionProbability(
                  profile[ row_i - 1 ],
                  sequence,
                  col_i - 1
                )
                );
              cout << "tmp_new before mutliplying by forward_row[ col_i - 1 ][ Insertion ] is " << tmp_new << endl;
              tmp_new *=
                (
                  forward_row[ col_i - 1 ][ Insertion ]
                );
              cout << "tmp_new after mutliplying by forward_row[ col_i - 1 ][ Insertion ] is " << tmp_new << endl;
            }
          }
        } // End if this is not the first or last row, and not the first col..
        // End from insertion

        // from Deletion
        // from Deletion to Match
        if( row_i == last_row ) {
          // We store post-align stuff in the Match state, even if it comes
          // from a deletion.  So this is really a D->D, stored at M.
          tmp_new =
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
#ifdef USE_END_DISTRIBUTION
          tmp_new *=
            profile[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ];
#endif // USE_END_DISTRIBUTION
          tmp_new *=
            prev_forward_row[ col_i ][ Deletion ];

          // TODO: REMOVE
          if( do_extra_debugging ) {
            cout << "col " << col_i << " D->D->C Match contribution: " << tmp_new << endl;
          } // End if do_extra_debugging
          if( use_viterbi ) {
            if( tmp_new > col_i_match_new ) {
              col_i_match_new = tmp_new;
            }
          } else {
            col_i_match_new += tmp_new;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
              cout << "forward_calculateRow: *6* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
            }
            
          }
        } // End if it's the last row..
        if( ( col_i != 0 ) &&
            ( row_i != 0 ) &&
            ( row_i != 1 ) ) { // no del state in row 0
          tmp_new =
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toMatch
            ];
          tmp_new *=
            getEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              col_i - 1
            );
#ifdef USE_END_DISTRIBUTION
          if( row_i == last_row ) {
            tmp_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_new *=
            prev_forward_row[ col_i - 1 ][ Deletion ];

          // TODO: REMOVE
          if( do_extra_debugging ) {
            if( row_i == last_row ) {
              cout << "col " << col_i << " D->M->C Match contribution: " << tmp_new << endl;
            } else {
              cout << "col " << col_i << " D->M Match contribution: " << tmp_new << endl;
            }
          } // End if do_extra_debugging
          if( use_viterbi ) {
            if( tmp_new > col_i_match_new ) {
              col_i_match_new = tmp_new;
            }
          } else {
            col_i_match_new += tmp_new;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
              cout << "forward_calculateRow: *7* col_i_match_new has just become " << col_i_match_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              cout << "tmp_new = ( ( " <<
               profile[ row_i - 2 ][
                 Transition::fromDeletion
               ][
                 TransitionFromDeletion::toMatch
               ] << " * " <<
               getEmissionProbability(
                 profile[ row_i - 1 ],
                 sequence,
                 col_i - 1
               ) << " ) * " <<
                " ( ( row_i == last_row ) = " << ( row_i == last_row ) << " ? ( ( "
                   << 
#ifdef USE_END_DISTRIBUTION
                profile[
                     Transition::fromEnd
                   ][
                     TransitionFromEnd::toPostAlign
                   ]
#else
                1.0
#endif // USE_END_DISTRIBUTION
                   << " ) * " <<
                prev_forward_row[ col_i - 1 ][ Deletion ] <<
                " ) : " <<
                prev_forward_row[ col_i - 1 ][ Deletion ] << " ) )" << endl;
            }
            
          }
        } // End if it's not the first or second row and not the first col..
        // from Deletion to Deletion
        if( ( row_i != 0 ) &&
            ( row_i != 1 ) && // no del state in first row
            ( row_i != last_row ) ) {
          tmp_new =
            profile[ row_i - 2 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
          tmp_new *=
            prev_forward_row[ col_i ][ Deletion ];

          if( use_viterbi ) {
            if( tmp_new > col_i_deletion_new ) {
              col_i_deletion_new = tmp_new;
            }
          } else {
            col_i_deletion_new += tmp_new;
            // TODO: REMOVE
            if( do_extra_debugging && ( isnan( col_i_deletion_new ) || isinf( col_i_deletion_new ) || ( col_i_deletion_new < 0 ) ) ) {
              cout << "forward_calculateRow: *dD* col_i_deletion_new has just become " << col_i_deletion_new << "!  tmp_new is " << tmp_new << ".  row_i = " << row_i << ", col_i = " << col_i << endl;
              cout << "Profile is " << &profile;
              cout << "profile[ row_i - 2 ].m_root is " << profile[ row_i - 2 ].m_root << endl;
              // TODO: REMOVE!
              exit( 0 );
            }
          }
        } // End if it's not the first, second, or last row ..
        // End from Deletion

        forward_row[ col_i ][ Match ] = col_i_match_new;
        forward_row[ col_i ][ Insertion ] = col_i_insertion_new;
        forward_row[ col_i ][ Deletion ] = col_i_deletion_new;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
        forward_row[ col_i ][ DeletionIn ] = col_i_deletion_in_new;
        forward_row[ col_i ][ DeletionOut ] = col_i_deletion_out_new;
#endif // ( USE_DEL_IN_DEL_OUT && !DISALLOW_FLANKING_TRANSITIONS )

        // TODO: REMOVE  Testing.
        if( do_extra_debugging && ( isnan( col_i_match_new ) || isinf( col_i_match_new ) || ( col_i_match_new < 0 ) ) ) {
          cout << "forward_calculateRow: col_i_match_new is " << col_i_match_new << "!  row_i = " << row_i << ", col_i = " << col_i << endl;
          cout << "forward_row is " << forward_row << endl;
          //cout << "The profile is " << profile << endl;
          //cout << "The sequence is " << sequence << endl;
          // TODO: REMOVE!
          exit( 0 );
        }
        if( do_extra_debugging && ( isnan( col_i_insertion_new ) || isinf( col_i_insertion_new ) || ( col_i_insertion_new < 0 ) ) ) {
          cout << "forward_calculateRow: col_i_insertion_new is " << col_i_insertion_new << "!  row_i = " << row_i << ", col_i = " << col_i << endl;
          cout << "The profile is " << profile << endl;
          cout << "The sequence is " << sequence << endl;
        }
        if( do_extra_debugging && ( isnan( col_i_deletion_new ) || isinf( col_i_deletion_new ) || ( col_i_deletion_new < 0 ) ) ) {
          cout << "forward_calculateRow: col_i_deletion_new is " << col_i_deletion_new << "!  row_i = " << row_i << ", col_i = " << col_i << endl;
          cout << "The profile is " << profile << endl;
          cout << "The sequence is " << sequence << endl;
        }

        if( use_rabiner_scaling && ( row_i != last_row ) ) {
          // It turns out that rabinerScaling can use any scale (except in the
          // last row, in which the cumulativeInverseScalar must be the score).
          // Rabiner suggested scaling each time separately, by the total of
          // the matrix values.  We rotate, so we scale each state separately
          // (not each time) -- but actually we group together the states in
          // each position, so we scale the positions.  We can either scale
          // them by the total of the values there (or by the total of just the
          // Match and Deletion values) or by the maximum value, which lets us
          // use more of the range of the MatrixValueType.
          if( parameters.rabinerScaling_useMaximumValue ) {
            if( col_i_match_new > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = col_i_match_new;
            }
            if( col_i_insertion_new > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = col_i_insertion_new;
            }
            if( col_i_deletion_new > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = col_i_deletion_new;
            }
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            if( col_i == last_col ) {
              if( forward_row.m_deletionIn > rabiner_inverse_scalar_mvt ) {
                rabiner_inverse_scalar_mvt = forward_row.m_deletionIn;
              }
              if( forward_row.m_deletionOut > rabiner_inverse_scalar_mvt ) {
                rabiner_inverse_scalar_mvt = forward_row.m_deletionOut;
              }
            }
#else
            if( col_i_deletion_in_new > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = col_i_deletion_in_new;
            }
            if( col_i_deletion_out_new > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = col_i_deletion_out_new;
            }
#endif // DISALLOW_FLANKING_TRANSITIONS
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
          } else { // if rabinerScaling_useMaximumValue .. else ..
            rabiner_inverse_scalar_mvt +=
              col_i_match_new;
            rabiner_inverse_scalar_mvt +=
              col_i_insertion_new;
            rabiner_inverse_scalar_mvt +=
              col_i_deletion_new;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifdef DISALLOW_FLANKING_TRANSITIONS
            if( col_i == last_col ) {
              rabiner_inverse_scalar_mvt +=
                forward_row.m_deletionIn;
              rabiner_inverse_scalar_mvt +=
                forward_row.m_deletionOut;
            }
#else
            rabiner_inverse_scalar_mvt +=
              col_i_deletion_in_new;
            rabiner_inverse_scalar_mvt +=
              col_i_deletion_out_new;
#endif // DISALLOW_FLANKING_TRANSITIONS
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
            // TODO: REMOVE
            if( do_extra_debugging && isnan( rabiner_inverse_scalar_mvt ) ) {
              cout << "in forward_calculateRow (row_i = " << row_i << ", col_i = " << col_i << "): m_rabinerInverseScalar is nan after adding in: " << endl;
              cout << "\t+= " << col_i_match_new << endl;
              cout << "\t+= " << col_i_insertion_new << endl;
              cout << "\t+= " << col_i_deletion_new << endl;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( DISALLOW_FLANKING_TRANSITIONS )
              cout << "\t+= " << col_i_deletion_in_new << endl;
              cout << "\t+= " << col_i_deletion_out_new << endl;
#endif // ( USE_DEL_IN_DEL_OUT && !DISALLOW_FLANKING_TRANSITIONS )

            }
          } // End if rabinerScaling_useMaximumValue .. else ..
        } // End if use_rabiner_scaling && ( row_i != last_row )
      } // End for each column, ..

      if( use_rabiner_scaling ) {
        if( row_i == last_row ) {
          rabiner_inverse_scalar_mvt =
            forward_row[ last_col ][ Match ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          // Note that we don't count del-ins here, because del-ins have to end
          // by the last position -- they have to delete into some position
          // (with a Match).
#ifdef DISALLOW_FLANKING_TRANSITIONS
          rabiner_inverse_scalar_mvt +=
            forward_row.m_deletionOut;
#else
          rabiner_inverse_scalar_mvt +=
            forward_row[ last_col ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

          rabiner_inverse_scalar_mvt *=
            profile[
              Transition::fromPostAlign
            ][
              TransitionFromPostAlign::toTerminal
            ];

          // TODO: REMOVE
          if( do_extra_debugging && isnan( rabiner_inverse_scalar_mvt ) ) {
            cout << "in forward_calculateRow: rabiner_inverse_scalar_mvt is nan after setting it to ( " << 
              forward_row[ last_col ][ Match ] << " * " << 
             profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toTerminal
              ]
                 << " )." << endl;
          }
          // TODO: REMOVE
          if( do_extra_debugging && ( rabiner_inverse_scalar_mvt == 0 ) ) {
            cout << "in forward_calculateRow: rabiner_inverse_scalar_mvt is 0 after setting it to ( " << 
              forward_row[ last_col ][ Match ] << " * " << 
              (
                profile[
                  Transition::fromPostAlign
                ][
                  TransitionFromPostAlign::toTerminal
                 ]
               )
                 << " )." << endl;
            if( forward_row[ last_col ][ Match ] == 0 ) {
              cout << "forward row is " << forward_row << endl;
            }
          }
        } else { // if( row_i == last_row ) .. else ..
          // TODO: REMOVE
           if( do_extra_debugging && isnan( rabiner_inverse_scalar_mvt ) ) {
             cout << "in forward_calculateRow: m_rabinerInverseScalar is nan." << endl;
           }
        } // End if( row_i == last_row ) .. else ..
        if( be_extra_verbose ) {
          cout << "Before scaling by " << parameters.matrixRowScaleFactor << ", rabiner_inverse_scalar_mvt is " << rabiner_inverse_scalar_mvt << endl;
        }
        rabiner_inverse_scalar_mvt /=
          parameters.matrixRowScaleFactor;
        if( be_extra_verbose ) {
          cout << "After scaling by " << parameters.matrixRowScaleFactor << ", rabiner_inverse_scalar_mvt is " << rabiner_inverse_scalar_mvt << endl;
        }
        // Note this can still cause crazy problems if the mvt is within some
        // epsilon of 0.
        // TODO: DEHACKIFY!  MAGIC #s
        if( ( rabiner_inverse_scalar_mvt == 0.0 ) || ( rabiner_inverse_scalar_mvt < 1E-250 ) ) {
          rabiner_inverse_scalar_mvt = 1.0;
          rabiner_inverse_scalar_mvt /=
            parameters.matrixRowScaleFactor;
        }
        forward_row /= rabiner_inverse_scalar_mvt;
        forward_row.m_rabinerInverseScalar =
          rabiner_inverse_scalar_mvt;
      } // End if use_rabiner_scaling

      // TODO: REMOVE  Testing.
      //for( col_i = 0; col_i <= last_col; col_i++ ) {
      //
      //  if( do_extra_debugging && ( isnan( forward_row[ col_i ][ Match ] ) || isinf( forward_row[ col_i ][ Match ] ) || ( forward_row[ col_i ][ Match ] < 0 ) ) ) {
      //    cout << "forward_calculateRow: forward_row[ col_i ][ Match ] is " << forward_row[ col_i ][ Match ] << "!  row_i = " << row_i << ", col_i = " << col_i << endl;
      //    cout << "forward_row is " << forward_row << endl;
      //    cout << "rabiner_inverse_scalar_mvt is " << rabiner_inverse_scalar_mvt << endl;
      //    //cout << "The profile is " << profile << endl;
      //    //cout << "The sequence is " << sequence << endl;
      //    // TODO: REMOVE!
      //    exit( 0 );
      //  }
      //  if( do_extra_debugging && ( isnan( forward_row[ col_i ][ Insertion ] ) || isinf( forward_row[ col_i ][ Insertion ] ) || ( forward_row[ col_i ][ Insertion ] < 0 ) ) ) {
      //    cout << "forward_calculateRow: forward_row[ col_i ][ Insertion ] is " << forward_row[ col_i ][ Insertion ] << "!  row_i = " << row_i << ", col_i = " << col_i << endl;
      //    cout << "forward_row is " << forward_row << endl;
      //    cout << "rabiner_inverse_scalar_mvt is " << rabiner_inverse_scalar_mvt << endl;
      //    //cout << "The profile is " << profile << endl;
      //    //cout << "The sequence is " << sequence << endl;
      //    // TODO: REMOVE!
      //    exit( 0 );
      //  }
      //  if( do_extra_debugging && ( isnan( forward_row[ col_i ][ Deletion ] ) || isinf( forward_row[ col_i ][ Deletion ] ) || ( forward_row[ col_i ][ Deletion ] < 0 ) ) ) {
      //    cout << "forward_calculateRow: forward_row[ col_i ][ Deletion ] is " << forward_row[ col_i ][ Deletion ] << "!  row_i = " << row_i << ", col_i = " << col_i << endl;
      //    cout << "forward_row is " << forward_row << endl;
      //    cout << "rabiner_inverse_scalar_mvt is " << rabiner_inverse_scalar_mvt << endl;
      //    //cout << "The profile is " << profile << endl;
      //    //cout << "The sequence is " << sequence << endl;
      //    // TODO: REMOVE!
      //    exit( 0 );
      //  }
      //} // End foreach col_i

      if( be_extra_verbose ) {
        cout << "forward row is " << forward_row << endl;
      }
      return forward_row;
    } // forward_calculateRow( Parameters const&, bool, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, Row & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
     /**
     * Using the forward row (of the Forward algorithm matrix) for (row_i + 1),
     * calculate the forward row for row_i.  Return the given forward_row
     * reference.  row_i must be not be 0, last_row, or last_row - 1.  The
     * values for the last column must be set ahead of time, as must the values
     * at colums evenly divisible by the store_every_Nth_column argument if it
     * is not 0.  The rabinerInverseScalars should be correctly set already (if
     * parameters.useRabinerScaling), as should the rows' deletionIn and
     * deletionOut values -- note that we do not support USE_DEL_IN_DEL_OUT
     * unless DISALLOW_FLANKING_TRANSITIONS is true.  Note that this is *not*
     * backward_calculateRow, which computes "backward" values of the
     * "forward-backward" algorithm.  This method computes "forward" values,
     * but computes earlier forward values from later ones (ie. in "reverse").
     */
    forward_reverseCalculateRow (
      Parameters const& parameters,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & next_row_match_distribution,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      uint32_t const & store_every_Nth_column,
      typename Matrix::Row const& next_forward_row,
      typename Matrix::Row & forward_row
    ) const
    {
#if defined( USE_DEL_IN_DEL_OUT) && !defined( DISALLOW_FLANKING_TRANSITIONS )
      // Not implemented!
      assert( false && "forward_reverseCalculateRow(..) is not implemented for teh case of USE_DEL_IN_DEL_OUT but not DISALLOW_FLANKING_TRANSITIONS" );
#endif // USE_DEL_IN_DEL_OUT && !DISALLOW_FLANKING_TRANSITIONS

      static const bool do_extra_debugging = false;
      const bool be_extra_verbose = false;

      const uint32_t last_row = profile.length();
      const uint32_t last_col = sequence.length();

      // We assume that the last col values have already been
      // set.  Also assume row_i != last_row and row_i != ( last_row - 1 ).
      // Also assume row_i != 0.  Also assume that the rabiner scalars are
      // already set.
      assert( row_i < ( last_row - 1 ) );
      assert( row_i > 0 );

      // Note: so far it doesn't seem to affect training! TODO: Experiment later.

      if( be_extra_verbose || ( parameters.debug == DEBUG_Special ) ) {
        cout << "row_i: " << row_i << endl;
        cout << "SEQUENCE: " << endl;
        cout << sequence << endl;
      }

      register uint32_t col_i;

      MatrixValueType col_i_match_new, col_i_insertion_new, col_i_deletion_new;
      MatrixValueType tmp_ItoM, tmp_DtoM, tmp_ItoI, tmp_MtoI, tmp_DtoD, tmp_MtoD;
      MatrixValueType tmp_new, tmp_new_2;

      MatrixValueType next_rabiner_inverse_scalar_mvt;
      MatrixValueType rabiner_inverse_scalar_mvt;
      if( parameters.useRabinerScaling ) {
        next_rabiner_inverse_scalar_mvt =
          next_forward_row.m_rabinerInverseScalar;
        rabiner_inverse_scalar_mvt =
          forward_row.m_rabinerInverseScalar;
      } // End if parameters.useRabinerScaling

      if( be_extra_verbose || ( parameters.debug == DEBUG_Special ) ) {
        cout << "[forward] last_row is " << last_row << ", last_col is " << last_col << endl;
      }
      if( be_extra_verbose ) {
        cout << "next forward row is " << next_forward_row << endl;
        cout << "last column of target forward row is " << forward_row[ last_col ] << endl;
      }

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      // We need to use scaled versions of the TransitionFromMatch distributions.
      if( row_i >= 1 ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using row_match_distribution == " << row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if row_i >= 1
      if( ( row_i + 1 ) >= 1 ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using next_row_match_distribution == " << next_row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if ( row_i + 1 ) >= 1
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      // row_i and col_i are *target* row and col.
      col_i = ( last_col - 1 );
      do { // count cols down to 0
        if(
          ( store_every_Nth_column != 0 ) &&
          ( ( col_i % store_every_Nth_column ) == 0 )
        ) {
          if( be_extra_verbose ) {
            cout << "Anchor column: " << col_i << endl;
          }
          continue;
        } // End if this is an anchor column

        // calculate row cell:
        col_i_match_new = 0;
        col_i_insertion_new = 0;
        col_i_deletion_new = 0;

        // TODO: REMOVE
        //if( parameters.debug == DEBUG_Special ) {
        //  cout << "[forward] row " << row_i << ", col " << col_i << endl;
        //}
        
        // Match
        // wrt Match cell in row_i + 1, col_i + 1
        if( col_i == 0 ) {
          // The Match value in the first column is always 0 (except in row 0,
          // but we assume row_i > 0)
          col_i_match_new = 0;
        } else {
          col_i_match_new =
            next_forward_row[ col_i + 1 ][ Match ];
          if( parameters.useRabinerScaling ) {
            col_i_match_new *=
              next_rabiner_inverse_scalar_mvt;
          } // End if parameters.useRabinerScaling
        } // End if col_i == 0 .. else ..

        // TODO: REMOVE
        //if( col_i == 497 ) {
        //  cout << "col_i_match_new is " << col_i_match_new << endl;
        //}

        // wrt Insertion cell in row_i, col_i + 1
        if( col_i <= 1 ) {
          // Insertion cell values are 0 in the first two columns.
          tmp_ItoM = 0;
          tmp_ItoI = 0;
        } else {
          tmp_ItoM =
            profile[ row_i - 1 ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toMatch
            ];
          tmp_ItoM *=
            getEmissionProbability(
              profile[ ( row_i + 1 ) - 1 ],
              sequence,
              ( col_i + 1 ) - 1
            );
          tmp_ItoI =
            profile[ row_i ][
              Transition::fromInsertion
            ][
              TransitionFromInsertion::toInsertion
            ];
          tmp_ItoI *=
            getInsertionEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              ( col_i + 1 ) - 1
            );
          // Mathematically, tmp_ItoI will not be 0, but numerically it might
          // be:
          if( tmp_ItoI != 0 ) {
            tmp_new =
              forward_row[ col_i + 1 ][ Insertion ];
            // TODO: Comment out again
            //if( parameters.useRabinerScaling ) {
            //  tmp_new *=
            //    rabiner_inverse_scalar_mvt;
            //} // End if parameters.useRabinerScaling
            tmp_new *=
              tmp_ItoM;
            tmp_new /=
              tmp_ItoI;
            if( tmp_new < col_i_match_new ) {
              col_i_match_new -=
                tmp_new;
            } else {
              col_i_match_new = 0;
            }
            // TODO: REMOVE
            //if( col_i == 497 ) {
            //  cout << "forward_row[ col_i + 1 ][ Insertion ] is " << forward_row[ col_i + 1 ][ Insertion ] << endl;
            //  cout << "tmp_ItoM is " << tmp_ItoM << endl;
            //  cout << "tmp_ItoI is " << tmp_ItoI << endl;
            //  cout << "tmp_new is " << tmp_new << endl;
            //  cout << "col_i_match_new is " << col_i_match_new << endl;
            //}

          } else {
            // TODO: REMOVE
            cerr << "Warning in reverseCalculateRow at "  << __LINE__ << ": tmp_ItoI is " << tmp_ItoI << endl;
          } // End if ( tmp_ItoI != 0 )
        } // End if col_i <= 1 .. else ..

        // wrt Deletion cell in row_i + 1, col_i
        tmp_DtoM =
          profile[ row_i - 1 ][
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toMatch
          ];
        tmp_DtoM *=
          getEmissionProbability(
            profile[ ( row_i + 1 ) - 1 ],
            sequence,
            ( col_i + 1 ) - 1
          );
        tmp_DtoD =
          profile[ row_i - 1 ][
            Transition::fromDeletion
          ][
            TransitionFromDeletion::toDeletion
          ];
        // Mathematically, tmp_DtoD will not be 0, but numerically it might
        // be:
        if( col_i > 0 ) {
          if( tmp_DtoD != 0 ) {
            tmp_new =
              next_forward_row[ col_i ][ Deletion ];
            if( parameters.useRabinerScaling ) {
              tmp_new *=
                next_rabiner_inverse_scalar_mvt;
            } // End if parameters.useRabinerScaling
            tmp_new *=
              tmp_DtoM;
            tmp_new /=
              tmp_DtoD;
            if( tmp_new < col_i_match_new ) {
              col_i_match_new -=
                tmp_new;
            } else {
              col_i_match_new = 0;
            }
            // TODO: REMOVE
            //if( col_i == 497 ) {
            //  cout << "forward_row[ col_i ][ Deletion ] is " << forward_row[ col_i ][ Deletion ] << endl;
            //  cout << "tmp_DtoM is " << tmp_DtoM << endl;
            //  cout << "tmp_DtoD is " << tmp_DtoD << endl;
            //  cout << "tmp_new is " << tmp_new << endl;
            //  cout << "col_i_match_new is " << col_i_match_new << endl;
            //}
          } else {
            // TODO: REMOVE
            cerr << "Warning in reverseCalculateRow at "  << __LINE__ << ": tmp_DtoD is " << tmp_DtoD << endl;
          }
        } // End if col_i > 0

        // And now for the denomenator
        if( col_i > 0 ) {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          tmp_new = // MtoM
            row_match_distribution[
              TransitionFromMatch::toMatch
            ];
#else
          tmp_new = // MtoM
            profile[ row_i - 1 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toMatch
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          tmp_new *= // MtoM
            getEmissionProbability(
              profile[ ( row_i + 1 ) - 1 ],
              sequence,
              ( col_i + 1 ) - 1
            );

          if( col_i > 1 ) {
            tmp_new_2 =
              tmp_ItoM;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            tmp_MtoI =
              next_row_match_distribution[
                TransitionFromMatch::toInsertion
              ];
#else
            tmp_MtoI =
              profile[ row_i ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toInsertion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            tmp_MtoI *=
              getInsertionEmissionProbability(
                profile[ row_i - 1 ],
                sequence,
                ( col_i + 1 ) - 1
              );
            tmp_new_2 *=
              tmp_MtoI;
            // Mathematically, tmp_ItoI will not be 0, but numerically it might
            // be:
            if( tmp_ItoI != 0 ) {
              tmp_new_2 /=
                tmp_ItoI;
              if( tmp_new_2 < tmp_new ) {
                tmp_new -=
                  tmp_new_2;
              } else {
                tmp_new = 0;
              }
            } else {
              // TODO: REMOVE
              cerr << "Warning in reverseCalculateRow at "  << __LINE__ << ": tmp_ItoI is " << tmp_ItoI << endl;
            }
          } // End if col_i > 1
  
          tmp_new_2 =
            tmp_DtoM;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          tmp_MtoD =
            row_match_distribution[
              TransitionFromMatch::toDeletion
            ];
#else
          tmp_MtoD =
            profile[ row_i - 1 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletion
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          tmp_new_2 *=
            tmp_MtoD;
          // Mathematically, tmp_DtoD will not be 0, but numerically it might
          // be:
          if( tmp_DtoD != 0 ) {
            tmp_new_2 /=
              tmp_DtoD;
            if( tmp_new_2 < tmp_new ) {
              tmp_new -=
                tmp_new_2;
            } else {
              tmp_new = 0;
            }
          } else {
            // TODO: REMOVE
            cerr << "Warning in reverseCalculateRow at "  << __LINE__ << ": tmp_DtoD is " << tmp_DtoD << endl;
          }

          if( tmp_new == 0 ) {
            col_i_match_new = 0;
          } else {
            col_i_match_new /=
              tmp_new;
          }

          // TODO: REMOVE
          //if( col_i == 497 ) {
          //  cout << "tmp_new is " << tmp_new << endl;
          //  cout << "col_i_match_new is " << col_i_match_new << endl;
          //}

        } // End if col_i > 0

        // Insertion
        if( ( col_i <= 1 ) || ( tmp_ItoI == 0 ) ) {
          col_i_insertion_new = 0;
        } else {
          col_i_insertion_new =
            forward_row[ col_i + 1 ][ Insertion ];
          // TODO: Comment out again
          //if( parameters.useRabinerScaling ) {
          //  col_i_insertion_new *=
          //    rabiner_inverse_scalar_mvt;
          //} // End if parameters.useRabinerScaling
          tmp_new =
            col_i_match_new;
          tmp_new *=
            tmp_MtoI;
          if( tmp_new >= col_i_insertion_new ) {
            //assert( tmp_new <= col_i_insertion_new );
            col_i_insertion_new = 0;
          } else {
            col_i_insertion_new -=
              tmp_new;
            col_i_insertion_new /=
              tmp_ItoI;
          }
        } // End if col_i <= 1 .. else ..

        // Deletion
        if( tmp_DtoD == 0 ) {
          col_i_deletion_new = 0;
        } else {
          col_i_deletion_new =
            next_forward_row[ col_i ][ Deletion ];
          if( parameters.useRabinerScaling ) {
            col_i_deletion_new *=
              next_rabiner_inverse_scalar_mvt;
          } // End if parameters.useRabinerScaling
          if( col_i > 0 ) {
            tmp_new =
              col_i_match_new;
            tmp_new *=
              tmp_MtoD;
            // TODO: REMOVE
            if( tmp_new >= col_i_deletion_new ) {
              //cout << "col_i_match_new is " << col_i_match_new << endl;
              //cout << "col_i_insertion_new is " << col_i_insertion_new << endl;
              //cout << "The rest of the forward row is ";
              //for( uint32_t col_j = ( col_i + 1 ); col_j <= last_col; col_j++ ) {
              //  cout << forward_row[ col_j ] << endl;
              //}
              //cout << "The rest of the next forward row is ";
              //for( uint32_t col_j = col_i; col_j <= last_col; col_j++ ) {
              //  cout << next_forward_row[ col_j ] << endl;
              //}
              //assert( tmp_new <= col_i_deletion_new );
              col_i_deletion_new = 0;
            } else {
              col_i_deletion_new -=
                tmp_new;
            }
          } // End if col_i > 0
          if( col_i_deletion_new != 0 ) {
            col_i_deletion_new /=
              tmp_DtoD;
          }
        } // End if tmp_DtoD == 0 .. else ..

        forward_row[ col_i ][ Match ] = col_i_match_new;
        forward_row[ col_i ][ Insertion ] = col_i_insertion_new;
        forward_row[ col_i ][ Deletion ] = col_i_deletion_new;
        //if( parameters.useRabinerScaling ) {
        //  forward_row[ col_i ] /=
        //    rabiner_inverse_scalar_mvt;
        //} // End if parameters.useRabinerScaling
      } while( col_i-- > 0 ); // End for each column, downto 0 (including 0)

      if( be_extra_verbose ) {
        cout << "forward row is now " << forward_row << endl;
      }
      return forward_row;
    } // forward_reverseCalculateRow( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, uint32_t const &, Row const&, Row & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate and return the likelihood of the given sequence, given the
     * given profile model, using the backward algorithm.  If the given
     * rabiner_inverse_scalars vector pointer is non-null, it should have
     * length one greater than the length of the profile, and should contain
     * the m_rabinerInverseScalars of the forward rows, calculated using
     * forward_calculateRow(..).
     * NOTE: If parameters.useRabinerScaling is true, you should provide a
     * non-null rabiner_inverse_scalars vector pointer, and vice-versa.  We
     * just assume that if it is non-null, then you mean us to use it.
     */
    backward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const vector<MatrixValueType> * const rabiner_inverse_scalars
    ) const
    {
      typename Matrix::Row prev_backward_row =
        typename Matrix::Row( sequence.length() + 1 );
      typename Matrix::Row backward_row =
        typename Matrix::Row( sequence.length() + 1 );
      ScoreType score;

      backward_score(
        parameters,
        profile,
        sequence,
        rabiner_inverse_scalars,
        prev_backward_row,
        backward_row,
        score
      );
      return score;
    } // backward_score( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, const vector<MatrixValueType> * const ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  bool
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate the likelihood of the given sequence, given the given profile
     * model, using the given Rows as temporary storage (these will be modified
     * such that afterwards, if the return value is false then backward_row
     * will hold the backward matrix row values for the first row (at index 0),
     * and prev_backward_row will hold the values for the row just after that
     * (at index 1), if the return value is true then these rows will be
     * swapped.  The side-effect, besides changing the backward rows, is to
     * change the value in the score reference argument.  The return value
     * indicates whether the rows are swapped.  If the given
     * rabiner_inverse_scalars vector pointer is non-null, it should have
     * length one greater than the length of the profile, and should contain
     * the m_rabinerInverseScalars of the forward rows, calculated using
     * forward_calculateRow(..).
     * NOTE: If parameters.useRabinerScaling is true, you should provide a
     * non-null rabiner_inverse_scalars vector pointer, and vice-versa.  We
     * just assume that if it is non-null, then you mean us to use it.
     */
    backward_score (
      Parameters const& parameters,
      ProfileType const& profile,
      Sequence<SequenceResidueType> const& sequence,
      const vector<MatrixValueType> * const rabiner_inverse_scalars,
      typename Matrix::Row & prev_backward_row,
      typename Matrix::Row & backward_row,
      ScoreType & score
    ) const
    {
      // NOTE: The return value will be true iff profile.length() is odd.

      uint32_t last_row = profile.length();
      uint32_t row_i = last_row;

      bool done = false;
      bool swap = true;
      while( !done ) { // Unsigned index going down, so no for loop.
        // Every other time, we must swap which backward row is current.
        swap = !swap;

        if( rabiner_inverse_scalars != NULL ) {
          ( swap ? prev_backward_row : backward_row ).m_rabinerInverseScalar =
            ( *rabiner_inverse_scalars )[ row_i ];
        } // End if useRabinerScaling

        backward_calculateRow(
          parameters,
          profile,
          sequence,
          row_i,
          ( swap ? backward_row : prev_backward_row ),
          ( swap ? prev_backward_row : backward_row )
        );

        if( rabiner_inverse_scalars != NULL ) {
          if( row_i == last_row ) {
            ( swap ? backward_row : prev_backward_row ).m_rabinerCumulativeInverseScalar =
              ( swap ? prev_backward_row : backward_row ).m_rabinerInverseScalar;
          } else {
            ( swap ? backward_row : prev_backward_row ).m_rabinerCumulativeInverseScalar =
              (
               ( swap ? prev_backward_row : backward_row ).m_rabinerCumulativeInverseScalar *
               ( swap ? prev_backward_row : backward_row ).m_rabinerInverseScalar
              );
          }
        } // End if useRabinerScaling

        if( row_i == 0 ) {
          break;
        }
        row_i--;
      } // End foreach row, from the end up to 0

      if( swap ) {
        score =
          prev_backward_row[ 0 ][ Match ];
        score /=
          parameters.matrixRowScaleFactor;
      } else {
        score =
          backward_row[ 0 ][ Match ];
        score /=
          parameters.matrixRowScaleFactor;
      }
      if( rabiner_inverse_scalars != NULL ) { // parameters.useRabinerScaling
        score *=
          ( swap ? backward_row : prev_backward_row ).m_rabinerCumulativeInverseScalar;
      }
      return swap;
    } // backward_score( Parameters const&, Profile const&, Sequence<SequenceResidueType> const&, const vector<MatrixValueType> * const, Row &, Row &, ScoreType & ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  ScoreType
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Calculate and return the total likelihood of the given sequences, given
     * the given profile model, using the backward algorithm.  As a side
     * effect, the given matrices references becomes filled with the backward
     * matrix values.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     * If useRabinerScaling is true, the m_rabinerInverseScalars of the
     * matrices will be used, so you must first compute the forward matrices
     * (using, eg. forward_score(..)), then call
     * matrices.copyRabinerInverseScalarsFrom( forward_matrices ) BEFORE
     * calling this.
     */
    backward_score (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
      vector<Sequence<SequenceResidueType> > const& sequences,
      uint32_t sequence_count,
      typename Matrix::SequentialAccessContainer & matrices
    ) const
    {
      sequence_count = ( ( sequence_count == 0 ) ? sequences.size() : min( static_cast<size_t>( sequence_count ), sequences.size() ) );

      // NOTE: The return value will be true iff profile.length() is odd.
      ScoreType score;

      typename Matrix::SequentialAccessContainer::reverse_iterator matrices_rev_iter =
        matrices.rbegin();
      typename Matrix::SequentialAccessContainer::reverse_iterator prev_matrices_rev_iter =
        matrices.rbegin();
      uint32_t last_row = profile.length();
      uint32_t row_i = last_row;
      uint32_t last_seq = sequence_count - 1;
      uint32_t seq_i;
      bool done = false;

      // TODO: REMOVE
      //cout << "Before Matrices: " << matrices << endl;

      while( !done ) {
        for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
          //if( row_i == last_row ) {
          //  // TODO: REMOVE
          //  cout << "about to call backward_calculateRow with prev row: (ignored)" << endl;
          //} else {
          //  // TODO: REMOVE
          //  cout << "about to call backward_calculateRow with prev row: " <<
          //    ( *( prev_matrices_rev_iter ) )[ seq_i ] << endl;
          //}
          backward_calculateRow(
            parameters,
            use_viterbi,
            profile,
            sequences[ seq_i ],
            row_i,
            ( ( row_i == last_row ) ?
              ( *matrices_rev_iter )[ seq_i ] : // ignored
              ( *( prev_matrices_rev_iter ) )[ seq_i ] ),
            ( *matrices_rev_iter )[ seq_i ]
          );

          // Note that the backard inverse scalars can actually be different
          // from the forward ones (there's a hack switch in
          // backward_calculateRow to calculate and use different rabiner
          // scalars there).
          if( parameters.useRabinerScaling ) {
            if( row_i == last_row ) {
              ( *matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                ( *matrices_rev_iter )[ seq_i ].m_rabinerInverseScalar;
            } else {
              ( *matrices_rev_iter )[ seq_i ].m_rabinerCumulativeInverseScalar =
                (
                  ( *( prev_matrices_rev_iter ) )[ seq_i ].m_rabinerCumulativeInverseScalar *
                  ( *matrices_rev_iter )[ seq_i ].m_rabinerInverseScalar
                );
            }
          } // End if useRabinerScaling
        } // End foreach sequence.
        if( row_i == 0 ) {
          done = true;
        } else {
          if( row_i != last_row ) {
            prev_matrices_rev_iter++;
          }
          matrices_rev_iter++;
          row_i--;
        }
      } // End foreach row.

      // TODO: REMOVE
      //cout << "After Matrices: " << matrices << endl;

      score = 1;
      ScoreType sequence_score;
      for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
        sequence_score =
          ( *( matrices.begin() ) )[ seq_i ][ 0 ][ Match ];
        sequence_score /=
          parameters.matrixRowScaleFactor;

        if( parameters.useRabinerScaling ) {
          sequence_score *=
            ( *( matrices.begin() ) )[ seq_i ].m_rabinerCumulativeInverseScalar;
        } // End if useRabinerScaling
        // TODO: REMOVE
        cout << "in backward_score(..): seq " << seq_i << " score: " <<
          sequence_score << endl;
        score *= sequence_score;
      } // End foreach sequence
      return score;
    } // backward_score( Parameters const&, bool, Profile const&, vector<Sequence<SequenceResidueType> > const&, uint32_t const&, SequentialAccessContainer& ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
  template <typename ProfileType,
            typename SequenceResidueType,
            typename ScaledMatchDistributionProbabilityType>
#else
  template <typename ProfileType,
            typename SequenceResidueType>
#endif
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  typename DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Matrix::Row &
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Using the backward row (of the Backward algorithm matrix) for (row_i +
     * 1), calculate the backward row for row_i.  Return the given backward_row
     * reference.  If the use_viterbi argument is true, calculate
     * most-likely-path likelihoods rather than full (all path) likelihoods.
     * If parameters.useRabinerScaling is true, then
     * the calculated forward row will be scaled by
     * backward_row.m_rabinerInverseScalar.
     * backward_row.m_rabinerInverseScalar is the sum of the *forward* row
     * values for row_i, calculated using forward_calculateRow(..).
     * NOTE: when row_i is last_row, backward_row.m_rabinerInverseScalar is the
     * last Match value times the postalign-to-terminal probability.
     */
    backward_calculateRow (
      Parameters const& parameters,
      bool use_viterbi,
      ProfileType const& profile,
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> const & row_match_distribution,
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      Sequence<SequenceResidueType> const& sequence,
      const uint32_t row_i,
      typename Matrix::Row const& next_backward_row,
      typename Matrix::Row & backward_row
    ) const
    {
      static const bool do_extra_debugging = false;
      static const bool be_extra_verbose = false;
      // Note we would normally use the rabinerInverseScalar from the
      // corresponding forward row, set ahead of time, but since the actual
      // scalar doesn't matter, I'm trying to avoid underflow/overflow issues
      // by using a separate scalar for each backward row.  It works!
      // NOTE: With this true, it is actually *not* necessary to first set the
      // rabiner scalars to the forward scalar values (by calling
      // copyRabinerInverseScalarsFrom(..)), nor is it necessary to pass the
      // forward scalars to the backward_score(..) function.  But of course,
      // with this true, it becomes necessary to update the
      // rabinerCumulativeInverseScalar values after calling this.
      // TODO: Make this a parameter?
      static const bool calculate_rabiner_inverse_scalar_here = false;//true;

      const uint32_t last_row = profile.length();
      const uint32_t last_col = sequence.length();

      const bool use_rabiner_scaling = parameters.useRabinerScaling;

      // We would use a for loop here, but since col_i is
      // an unsigned type, and since we're decrementing as we go, we can't
      // check that it's hit zero *after* it's gone below zero, so we use
      // a while loop instead.
      register bool done = false;

      // TODO: REMOVE
      /*
      if( backward_row.length() != ( last_col + 1 ) ) {
        cout << "Uh oh: sequence is " << sequence << "; its length is " << last_col << ", which should be one less than the length of the backward row, which is " << backward_row << ", and has length " << backward_row.length() << "." << endl;
        //assert( 0 );
      }
      */
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      // We use scaled versions of the TransitionFromMatch distributions.
      if( ( row_i >= 1 ) && ( row_i != last_row ) ) {
        // TODO: REMOVE
        if( do_extra_debugging ) {
          cout << "using row_match_distribution == " << row_match_distribution << endl;
        } // End if do_extra_debugging
      } // End if ( row_i >= 1 ) && ( row_i != last_row )
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      MatrixValueType rabiner_inverse_scalar_mvt;
      if( use_rabiner_scaling && calculate_rabiner_inverse_scalar_here ) {
        rabiner_inverse_scalar_mvt = 0;
      } // End if use_rabiner_scaling && calculate_rabiner_inverse_scalar_here

#if defined( USE_DEL_IN_DEL_OUT ) && defined( DISALLOW_FLANKING_TRANSITIONS )
      backward_row.m_deletionIn = 0;
      backward_row.m_deletionOut = 0;
#endif // ( USE_DEL_IN_DEL_OUT && DISALLOW_FLANKING_TRANSITIONS )

      register uint32_t col_i = last_col;
      MatrixValueType col_i_new, tmp_new;
      while( !done ) {
    
        if( false && ( parameters.debug >= DEBUG_All ) ) {
          cout << "[backward] row " << row_i << ", col " << col_i << endl;
        }

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
#ifndef DISALLOW_FLANKING_TRANSITIONS
        backward_row[ col_i ][ DeletionIn ] = 0;
        backward_row[ col_i ][ DeletionOut ] = 0;
#endif // !DISALLOW_FLANKING_TRANSITIONS
        // From DeletionIn
        col_i_new = 0;
        if(
          ( col_i != last_col )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // from DeletionIn to Match
          if( !( ( row_i == 0 ) || ( row_i == last_row ) ) ) {
            tmp_new =
              profile[
                Transition::fromDeletionIn
              ][
                TransitionFromDeletionIn::toMatch
              ];
            tmp_new *=
              getEmissionProbability(
                profile[ row_i ],
                sequence,
                col_i
              );
            tmp_new *=
              next_backward_row[ col_i + 1 ][ Match ];
            if( row_i == ( last_row - 1 ) ) {
              // Then this is actually a transition into the last row, so
              // there's also post-align stuff to deal with.
#ifdef USE_END_DISTRIBUTION
              tmp_new *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
#endif // USE_END_DISTRIBUTION
            } // End if ( row_i == ( last_row - 1 ) )
            col_i_new = tmp_new;
          } // End if this is not the first row nor the last row ..
          // from DeletionIn to DeletionIn
          if( ( row_i != 0 ) && ( row_i < ( last_row - 1 ) ) ) {
            // There's no deletionIns of the last position, since deletionIns
            // have to end in a Match.
            tmp_new =
              profile[
                Transition::fromDeletionIn
              ][
                TransitionFromDeletionIn::toDeletionIn
              ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
            tmp_new *= next_backward_row.m_deletionIn;
#else
            tmp_new *= next_backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..

            if( use_viterbi ) {
              if( tmp_new > col_i_new ) {
                col_i_new = tmp_new;
              }
            } else {
              col_i_new += tmp_new;
            }
          } // End if ( ( row_i > 0 ) && ( row_i < ( last_row - 1 ) ) )
        } // End if ( col_i != last_col ) ( and maybe col_i == 0 )
        // End from DeletionIn
#ifdef DISALLOW_FLANKING_TRANSITIONS
        if( col_i == 0 ) {
          backward_row.m_deletionIn = col_i_new;
        }
#else
        backward_row[ col_i ][ DeletionIn ] = col_i_new;
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
        
        // From DeletionOut
        col_i_new = 0;
        if( ( col_i == 0 ) || ( row_i == 0 ) ) {
          // No delOuts here..
        } else if(
          ( row_i == last_row )
        ) {
          if( col_i == last_col ) {
            // from DeletionOut to End to PostAlign to Terminal
            col_i_new =
              profile[
                Transition::fromDeletionOut
              ][
                TransitionFromDeletionOut::toEnd
              ];
#ifdef USE_END_DISTRIBUTION
            col_i_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            // deletionOuts directly into the final cell don't go through the
            // post-align (Match) cell, so we need to incorporate the terminal
            // part ourselves, and multiply in the scale factor.
            col_i_new *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toTerminal
              ];
            col_i_new *=
              parameters.matrixRowScaleFactor;
          } // End if col_i == last_col
#ifndef DISALLOW_FLANKING_TRANSITIONS
          else { // col_i < last_col
            // from DeletionOut to Match
            // Actually from DeletionOut to End, End to PostAlign, PostAlign to
            // PostAlign.
            // In the last row the del-out paths go to / through post-aligns.
            col_i_new =
              profile[
                Transition::fromDeletionOut
              ][
                TransitionFromDeletionOut::toEnd
              ];
#ifdef USE_END_DISTRIBUTION
            col_i_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            col_i_new *=
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ];
            col_i_new *=
              getPostAlignEmissionProbability(
                profile,
                sequence,
                col_i
              );
            col_i_new *=
               // In the bottom row we use postAlignInsertions, which are not
               // affine, so we store the values in the Match cell instead of
               // the Insertion cell.
              backward_row[ col_i + 1 ][ Match ];
          } // End if col_i < last_col
#endif // !DISALLOW_FLANKING_TRANSITIONS
        } else
#ifdef DISALLOW_FLANKING_TRANSITIONS
          if( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        {
          // from DeletionOut to DeletionOut
#ifdef DISALLOW_FLANKING_TRANSITIONS
          col_i_new =
            next_backward_row.m_deletionOut;
#else
          col_i_new =
            next_backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          col_i_new *=
            profile[
              Transition::fromDeletionOut
            ][
              TransitionFromDeletionOut::toDeletionOut
            ];
        }
        // Don't worry about viterbi, since there's only one path from
        // deletionOut.
#ifdef DISALLOW_FLANKING_TRANSITIONS
        if( col_i == last_col ) {
          backward_row.m_deletionOut = col_i_new;
        }
#else
        backward_row[ col_i ][ DeletionOut ] = col_i_new;
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
        // End from DeletionOut
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

        if( col_i < last_col ) {
          // If not in the rightmost column, then there's base emission(s)
          // from the subsequent cell (col_i - 1 would be from *this* cell).
        } else if( row_i == last_row ) {
          // This is the very bottom-rightmost cell.
#ifdef USE_END_DISTRIBUTION
          backward_row[ col_i ][ Match ] =
            profile[
              Transition::fromEnd
            ][
              TransitionFromEnd::toPostAlign
            ];
#else
          backward_row[ col_i ][ Match ] = 1;
#endif // USE_END_DISTRIBUTION .. else ..
          backward_row[ col_i ][ Match ] *=
            profile[
              Transition::fromPostAlign
            ][
              TransitionFromPostAlign::toTerminal
            ];
          backward_row[ col_i ][ Match ] *=
            parameters.matrixRowScaleFactor;

          backward_row[ col_i ][ Insertion ] =
            backward_row[ col_i ][ Deletion ] =
            0;

          // Move on to the next col_i, descending...
          if( col_i == 0 ) {
            done = true;
          } else {
            col_i--;
          }
          continue;
        } // end if col_i > 0 .. else if row_i == last_row ..
    
        // from Match
        col_i_new = 0;
        tmp_new = 0;
        // from Match to Match
        if( col_i != last_col ) {
          if( row_i == last_row ) {
            // In the bottom row, we use postAlignInsertions.
            tmp_new =
              profile[
                Transition::fromPostAlign
              ][
                TransitionFromPostAlign::toPostAlign
              ];
            tmp_new *=
              getPostAlignEmissionProbability(
                profile,
                sequence,
                col_i
              );
            tmp_new *=
               // In the bottom row we use postAlignInsertions, which are not
               // affine, so we store the values in the Match cell instead of
               // the Insertion cell.
              backward_row[ col_i + 1 ][ Match ];

            if( use_viterbi ) {
              if( tmp_new > col_i_new ) {
                col_i_new = tmp_new;
              }
            } else {
              col_i_new += tmp_new;
              // TODO: REMOVE.  TESTING.
              //if(
              //  ( col_i_new == 0 ) ||
              //  ( col_i_new.val.f < 1E-18 ) ||
              //  ( tmp_new.val.f < 1E-18 )
              //) {
              //  cout << " col_i is " << col_i << endl;
              //  cout << " col_i_new is " << col_i_new << endl;
              //  cout << " col_i_new.val.f is " << col_i_new.val.f << endl;
              //  cout << " col_i_new.val.e is " << col_i_new.val.e << endl;
              //  cout << " tmp_new is " << tmp_new;
              //  cout << " tmp_new.val.f is " << tmp_new.val.f << endl;
              //  cout << " tmp_new.val.e is " << tmp_new.val.e << endl;
              //  cout << "   post-align-to-post-align prob is " <<
              //    profile[
              //      Transition::fromPostAlign
              //    ][
              //      TransitionFromPostAlign::toPostAlign
              //    ] << endl;
              //  tmp_new = 
              //    getPostAlignEmissionProbability(
              //      profile,
              //      sequence,
              //      col_i
              //    );
              //  cout << "   post-align ins prob is " << tmp_new << endl;
              //  cout << "                     f is " << tmp_new.val.f << endl;
              //  cout << "                     e is " << tmp_new.val.e << endl;
              //  tmp_new =
              //    backward_row[ col_i + 1 ][ Match ];
              //  cout << "    backward_row[ col_i + 1 ][ Match ] is " << tmp_new << endl;
              //  cout << "                                     f is " << tmp_new.val.f << endl;
              //  cout << "                                     e is " << tmp_new.val.e << endl;
              //  tmp_new =
              //    profile[
              //      Transition::fromPostAlign
              //    ][
              //      TransitionFromPostAlign::toPostAlign
              //    ];
              //  assert( tmp_new != 0 );
              //  tmp_new *=
              //    getPostAlignEmissionProbability(
              //      profile,
              //      sequence,
              //      col_i
              //    );
              //  assert( tmp_new != 0 );
              //  cout << " tmp_new.val.f is " << tmp_new.val.f << endl;
              //  cout << " tmp_new.val.e is " << tmp_new.val.e << endl;
              //  tmp_new *=
              //    // In the bottom row we use postAlignInsertions, which are not
              //    // affine, so we store the values in the Match cell instead of
              //    // the Insertion cell.
              //    backward_row[ col_i + 1 ][ Match ];
              //  cout << " tmp_new.val.f is " << tmp_new.val.f << endl;
              //  cout << " tmp_new.val.e is " << tmp_new.val.e << endl;
              //  if( tmp_new == 0 ) {
              //    assert( tmp_new != 0 );
              //  }
              //  assert( col_i_new != 0 );
              //}
            }
          } else {
            if( row_i == 0 ) {
              // In the top row, we use preAlignInsertions
              tmp_new =
                profile[
                  Transition::fromPreAlign
                ][
                  TransitionFromPreAlign::toPreAlign
                ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
              // TODO: REMOVE
              assert( tmp_new == 0 );
#endif // DISALLOW_FLANKING_TRANSITIONS
              tmp_new *=
                getPreAlignEmissionProbability(
                  profile,
                  sequence,
                  col_i
                );
              tmp_new *=
                  // In the top row we use preAlignInsertions, which are not
                  // affine, so we store the values in the Match cell instead of
                  // the Insertion cell.
                backward_row[ col_i + 1 ][ Match ];

              if( use_viterbi ) {
                if( tmp_new > col_i_new ) {
                  col_i_new = tmp_new;
                }
              } else {
                col_i_new += tmp_new;
              }
            } // End if it's the first row, include preAlign insertions
            //// TODO: There is an assumption that row 1 is not the last row.
            if( row_i == 0 ) {
              tmp_new =
                profile[
                  Transition::fromPreAlign
                ][
                  TransitionFromPreAlign::toBegin
                ];
              tmp_new *=
                profile[
                  Transition::fromBegin
                ][
                  TransitionFromBegin::toMatch
                ];
            } else { // row_i > 0
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
              tmp_new =
                row_match_distribution[
                  TransitionFromMatch::toMatch
                ];
#else
              tmp_new =
                profile[ row_i - 1 ][
                  Transition::fromMatch
                ][
                  TransitionFromMatch::toMatch
                ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            }
            tmp_new *=
              getEmissionProbability(
                profile[ row_i ],
                sequence,
                col_i
              );
#ifdef USE_END_DISTRIBUTION
            if( row_i == ( last_row - 1 ) ) {
              tmp_new *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION
            tmp_new *=
              next_backward_row[ col_i + 1 ][ Match ];

            if( use_viterbi ) {
              if( tmp_new > col_i_new ) {
                col_i_new = tmp_new;
              }
            } else {
              col_i_new += tmp_new;
            }
          } // End if this is the first row .. elsif the last row .. else ..
        } // End if this is not the last column .. 
        // from Match to Insertion
        if( ( col_i != last_col ) &&
            ( row_i != 0 ) &&
            ( row_i != last_row ) ) {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
          tmp_new =
            row_match_distribution[
              TransitionFromMatch::toInsertion
            ];
#else
          tmp_new =
            profile[ row_i - 1 ][
              Transition::fromMatch
            ][
              TransitionFromMatch::toInsertion
            ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
          tmp_new *=
            getInsertionEmissionProbability(
              profile[ row_i - 1 ],
              sequence,
              col_i
            );
          tmp_new *=
            backward_row[ col_i + 1 ][ Insertion ];

          if( use_viterbi ) {
            if( tmp_new > col_i_new ) {
              col_i_new = tmp_new;
            }
          } else {
            col_i_new += tmp_new;
          }
        } // End if this is not the last column and not the first or last row ..
        // from Match to Deletion
        if( row_i != last_row ) {
          if( row_i == 0 ) {
            tmp_new =
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            tmp_new *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletion
              ];
            tmp_new *=
              next_backward_row[ col_i ][ Deletion ];

            if( use_viterbi ) {
              if( tmp_new > col_i_new ) {
                col_i_new = tmp_new;
              }
            } else {
              col_i_new += tmp_new;
            }
          } else {
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            tmp_new =
              row_match_distribution[
                TransitionFromMatch::toDeletion
              ];
#else
            tmp_new =
              profile[ row_i - 1 ][
                Transition::fromMatch
              ][
                TransitionFromMatch::toDeletion
              ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            
            // We store the postAlignInsertion stuff in the Match state.
            if( row_i == ( last_row - 1 ) ) {
#ifdef USE_END_DISTRIBUTION
              tmp_new *=
                profile[
                      Transition::fromEnd
                    ][
                      TransitionFromEnd::toPostAlign
                    ];
#endif // USE_END_DISTRIBUTION
              tmp_new *=
                next_backward_row[ col_i ][ Match ];
            } else {
              tmp_new *=
                next_backward_row[ col_i ][ Deletion ];
            }

            if( use_viterbi ) {
              if( tmp_new > col_i_new ) {
                col_i_new = tmp_new;
              }
            } else {
              col_i_new += tmp_new;
            }
          } // End if it's the first row .. else ..
        } // End if it's not the last row ..
#if defined( USE_DEL_IN_DEL_OUT )
        // from Match to DeletionIn
        // Really this is from Begin to DeletionIn, and incorporates from
        // PreAlign to Begin..
        if(
          ( row_i == 0 )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == 0 )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // Implicitly we're leaving the PreAlign state too.
          tmp_new =
            profile[
              Transition::fromPreAlign
            ][
              TransitionFromPreAlign::toBegin
            ];
#ifdef USE_SWENTRY_SWEXIT
          ///// TODDO: ERE I AM.  The problem is that I think we will need to store a m_deletionIn field for backwardRows if we insist on being able to calculate the score, so we can keep track of the mass that passes through the B->Z transition.  For forward calculation, the missing mass is the remaining fraction: B->Z * ( ( profile.length() - 1 ) - ( pos.i + 1 ) ) / ( profile.length() - 1 ).  We don't need to keep track of it.  Nor do we ever need to keep track of the delOut mass, because backwards it just goes straight to the PostAlign state, and forwards ... hrm maybe we do need to keep track of it for the same reason we keep track of the delIn here.  Yes. /// I could just keep all of them and use the redundant one anyway, sort-of-simplifying things.  Then rather than B->Z*1/(profile_length-1), we'd use prev_row.m_deletionIn * ( 1 / ( ( profile.length() - 1 ) - ( pos.i + 1 ) ) ), and each m_deletionIn would be the remaining fraction as described above: that's the total prob of all paths that haven't yet deleted-in.  That maintains its interpretation from before.

          // mark
            tmp_new =
              profile[
                Transition::fromPreAlign
              ][
                TransitionFromPreAlign::toBegin
              ];
            tmp_new *=
              profile[
                Transition::fromBegin
              ][
                TransitionFromBegin::toDeletionIn
              ];
            tmp_new *= ( 1.0 / ( profile.length() - 1 ) );
#else // USE_SWENTRY_SWEXIT .. else ..
          tmp_new *=
            profile[
              Transition::fromBegin
            ][
              TransitionFromBegin::toDeletionIn
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_new *= next_backward_row.m_deletionIn;
#else
          tmp_new *= next_backward_row[ col_i ][ DeletionIn ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
#endif // USE_SWENTRY_SWEXIT .. else ..
          if( use_viterbi ) {
            if( tmp_new > col_i_new ) {
              col_i_new = tmp_new;
            }
          } else {
            col_i_new += tmp_new;
          }
        } // End if ( row_i == 0 )
        // from Match to DeletionOut
        if(
          ( row_i != 0 ) && ( row_i < last_row )
#ifdef DISALLOW_FLANKING_TRANSITIONS
          && ( col_i == last_col )
#endif // DISALLOW_FLANKING_TRANSITIONS
        ) {
          // Note that we don't use the scaled distribution -- we *want* just
          // the DeletionOut OPEN probability here, since the rest is already
          // multiplied into the next backward row's value.
          tmp_new =
            profile[
              Transition::fromMatch
            ][
              TransitionFromMatch::toDeletionOut
            ];
#ifdef DISALLOW_FLANKING_TRANSITIONS
          tmp_new *= next_backward_row.m_deletionOut;
#else
          tmp_new *= next_backward_row[ col_i ][ DeletionOut ];
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..
          
          if( use_viterbi ) {
            if( tmp_new > col_i_new ) {
              col_i_new = tmp_new;
            }
          } else {
            col_i_new += tmp_new;
          }
        } // End if ( row_i != 0 ) && ( row_i != last_row ) (and maybe col_i == last_col )
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
        // End from Match
        backward_row[ col_i ][ Match ] = col_i_new;
    
        // from Insertion
        col_i_new = 0;
        tmp_new = 0;
        if( ( row_i == 0 ) || ( row_i == last_row ) ) {
          // In the top row we use preAlignInsertions, which are not affine, so we
          // store the values in the Match cell instead of the Insertion cell.
          // Same with the bottom row, where we use postAlignInsertions.
    
          // (So do nothing else.)
        } else {
          // from Insertion to Match
          if( ( col_i != 0 ) && ( col_i != last_col ) ) {
            tmp_new =
              profile[ row_i - 1 ][
                Transition::fromInsertion
              ][
                TransitionFromInsertion::toMatch
              ];
            tmp_new *=
              getEmissionProbability(
                profile[ row_i ],
                sequence,
                col_i
              );
#ifdef USE_END_DISTRIBUTION
            if( row_i == ( last_row - 1 ) ) {
              tmp_new *=
                profile[
                  Transition::fromEnd
                ][
                  TransitionFromEnd::toPostAlign
                ];
            }
#endif // USE_END_DISTRIBUTION
            tmp_new *=
              next_backward_row[ col_i + 1 ][ Match ];

            if( use_viterbi ) {
              if( tmp_new > col_i_new ) {
                col_i_new = tmp_new;
              }
            } else {
              col_i_new += tmp_new;
            }
          } // End if it's not the first nor the last col.
          // from Insertion to Insertion
          if( ( col_i != 0 ) && ( col_i != last_col ) ) {
            tmp_new =
              profile[ row_i - 1 ][
                Transition::fromInsertion
              ][
                TransitionFromInsertion::toInsertion
              ];
            tmp_new *=
              getInsertionEmissionProbability(
                profile[ row_i - 1 ],
                sequence,
                col_i
              );
            tmp_new *=
              backward_row[ col_i + 1 ][ Insertion ];

            if( use_viterbi ) {
              if( tmp_new > col_i_new ) {
                col_i_new = tmp_new;
              }
            } else {
              col_i_new += tmp_new;
            }
          } // End if it's not the first nor the last col.
        } // End if this is the topmost or bottommost row, else ..
        // End from insertion
        backward_row[ col_i ][ Insertion ] = col_i_new;
    
        // from Deletion
        col_i_new = 0;
        tmp_new = 0;
        // from Deletion to Match
        if( ( col_i != last_col ) &&
            ( row_i != 0 ) &&
            ( row_i != last_row ) ) {
          tmp_new =
            profile[ row_i - 1 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toMatch
            ];
          tmp_new *=
            getEmissionProbability(
              profile[ row_i ],
              sequence,
              col_i
            );
#ifdef USE_END_DISTRIBUTION
          if( row_i == ( last_row - 1 ) ) {
            tmp_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
          }
#endif // USE_END_DISTRIBUTION
          tmp_new *=
            next_backward_row[ col_i + 1 ][ Match ];

          if( use_viterbi ) {
            if( tmp_new > col_i_new ) {
              col_i_new = tmp_new;
            }
          } else {
            col_i_new += tmp_new;
          }
        } // End if it's not the last col, and neither the first nor last row..
        // from Deletion to Deletion
        if( ( row_i != last_row ) &&
            ( row_i != 0 ) ) {
          tmp_new =
            profile[ row_i - 1 ][
              Transition::fromDeletion
            ][
              TransitionFromDeletion::toDeletion
            ];
             // We store the postAlignInsertion stuff in the Match state.
          if( row_i == ( last_row - 1 ) ) {
#ifdef USE_END_DISTRIBUTION
            tmp_new *=
              profile[
                Transition::fromEnd
              ][
                TransitionFromEnd::toPostAlign
              ];
#endif // USE_END_DISTRIBUTION
            tmp_new *=
              next_backward_row[ col_i ][ Match ];
          } else {
            tmp_new *=
              next_backward_row[ col_i ][ Deletion ];
          }

          if( use_viterbi ) {
            if( tmp_new > col_i_new ) {
              col_i_new = tmp_new;
            }
          } else {
            col_i_new += tmp_new;
          }
        } // End if it's neither the first nor the last row
        // End from Deletion
        backward_row[ col_i ][ Deletion ] = col_i_new;
    
        if(
          use_rabiner_scaling &&
          calculate_rabiner_inverse_scalar_here
        ) {
          // It turns out that rabinerScaling can use any scale.
          // Rabiner suggested scaling each time separately, by the total of
          // the matrix values.  We rotate, so we scale each state separately
          // (not each time) -- but actually we group together the states in
          // each position, so we scale the positions.  We can either scale
          // them by the total of the values there (or by the total of just the
          // Match and Deletion values) or by the maximum value, which lets us
          // use more of the range of the MatrixValueType.
          if( parameters.rabinerScaling_useMaximumValue ) {
            if( backward_row[ col_i ][ Match ] > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = backward_row[ col_i ][ Match ];
            }
            if( backward_row[ col_i ][ Deletion ] > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = backward_row[ col_i ][ Deletion ];
            }
            if( backward_row[ col_i ][ Insertion ] > rabiner_inverse_scalar_mvt ) {
              rabiner_inverse_scalar_mvt = backward_row[ col_i ][ Insertion ];
            }
          } else { // if rabinerScaling_useMaximumValue .. else ..
            rabiner_inverse_scalar_mvt +=
              backward_row[ col_i ][ Match ];
            rabiner_inverse_scalar_mvt +=
              backward_row[ col_i ][ Insertion ];
            rabiner_inverse_scalar_mvt +=
              backward_row[ col_i ][ Deletion ];
          } // End if rabinerScaling_useMaximumValue .. else ..

          // TODO: REMOVE
           if( do_extra_debugging && isnan( rabiner_inverse_scalar_mvt ) ) {
             cout << "in backward_calculateRow (row_i = " << row_i << ", col_i = " << col_i << "): m_rabinerInverseScalar is nan after adding in: " << endl;
             cout << "\t+= " << backward_row[ col_i ][ Match ] << endl;
             cout << "\t+= " << backward_row[ col_i ][ Insertion ] << endl;
             cout << "\t+= " << backward_row[ col_i ][ Deletion ] << endl;
           }
        } // End if use_rabiner_scaling && calculate_rabiner_inverse_scalar_here

        if( col_i == 0 ) {
          done = true;
        } else {
          col_i--;
        }
      } // For each col_i, right to left.

      if( use_rabiner_scaling && calculate_rabiner_inverse_scalar_here ) {
        // TODO: REMOVE
        if( do_extra_debugging && isnan( rabiner_inverse_scalar_mvt ) ) {
          cout << "in backward_calculateRow: m_rabinerInverseScalar is nan." << endl;
        }
        rabiner_inverse_scalar_mvt /=
          parameters.matrixRowScaleFactor;
        if( be_extra_verbose ) {
          cout << "After scaling by " << parameters.matrixRowScaleFactor << ", rabiner_inverse_scalar_mvt is " << rabiner_inverse_scalar_mvt << endl;
        }
        if( be_extra_verbose ) {
          cout << "Before scaling by " << parameters.matrixRowScaleFactor << ", rabiner_inverse_scalar_mvt is " << rabiner_inverse_scalar_mvt << endl;
        }
        // Note this can still cause crazy problems if the mvt is within some
        // epsilon of 0.
        // TODO: DEHACKIFY!  MAGIC #s
        if( ( rabiner_inverse_scalar_mvt == 0.0 ) || ( rabiner_inverse_scalar_mvt < 1E-250 ) ) {
          rabiner_inverse_scalar_mvt = 1.0;
          rabiner_inverse_scalar_mvt /=
            parameters.matrixRowScaleFactor;
        }
        backward_row /= rabiner_inverse_scalar_mvt;
        backward_row.m_rabinerInverseScalar =
          rabiner_inverse_scalar_mvt;
      } else if( use_rabiner_scaling && ( backward_row.m_rabinerInverseScalar > 0 ) ) {
        // Scale it.
        backward_row /= backward_row.m_rabinerInverseScalar;
      }
      return backward_row;
    } // backward_calculateRow( Parameters const&, bool, Profile const&, Sequence<SequenceResidueType> const&, uint32_t, Row const&, Row & ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM
  uint32_t
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Generate (randomly) sequences from the given Profile model.  Actually,
     * paths are drawn, and in addition to modifying the given Fasta reference
     * to hold the sequences, the given MultipleAlignment reference is modified
     * to hold the paths.  The Fasta type includes sequence description
     * strings.  These descriptions will be set to be ( description_prefix +
     * seq_i ) for each sequence seq_i (so the first sequence would have name
     * "Randomly generated 0" if the description_prefix == "Randomly generated
     * ").
     *
     * NOTE: For now we will only include sequences of length > 0.  We will
     * keep drawing until we get a sequence of non-zero length, if we have to.
     * The return value is the number of 0-length sequences drawn.
     */
    drawSequences (
      Parameters const& parameters,
      ProfileType const & profile,
      uint32_t const sequence_count,
      string const & description_prefix,
      Random & random,
      Fasta<SequenceResidueType> & fasta,
      MultipleAlignment<ProfileType, SequenceResidueType> & multiple_alignment
    ) const
    {
      fasta.reinitialize(
        sequence_count
      );
      multiple_alignment.reinitialize(
        &profile,
        &fasta,
        sequence_count
      );

      bool restart = false;
      uint32_t zero_length_seq_count = 0;
      for( uint32_t seq_i = 0; seq_i < sequence_count; seq_i++ ) {
        if( restart ) {
          seq_i = 0;
          restart = false;
        }
        // TODO: REMOVE.  Debugging.
        assert( multiple_alignment.m_sequenceAdvances.size() > seq_i );
        assert( fasta.size() > seq_i );
        drawSequence(
          parameters,
          profile,
          random,
          multiple_alignment.m_sequenceAdvances[ seq_i ],
          fasta[ seq_i ]
        );
        if( fasta[ seq_i ].length() == 0 ) {
          // TODO: REMOVE?
          cout << "WARNING: drawSequences(..) is redrawing a sequence because it has length 0, which we do not allow." << endl;
          zero_length_seq_count += 1;
          // Try again.
          if( seq_i == 0 ) {
            restart = true;
          } else {
            seq_i--;
          }
        } else {
          assert( fasta.m_descriptions.size() > seq_i );
          fasta.m_descriptions[ seq_i ] =
            ( description_prefix + boost::lexical_cast<string>( seq_i ) );
        }
      } // End foreach sequence

      return zero_length_seq_count;
    } // drawSequences( Parameters const &, ProfileType const &, uint32_t const ) const


  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  uint32_t
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    /**
     * Generate (randomly) a sequence from the given Profile model.  Both the
     * sequence and the particular path drawn are returned to the caller via
     * modifying the provided references.  The sequence_advances int-vector
     * will have length two greater than that of the profile, and the integer
     * at index i will indicate the number of sequence positions advanced at
     * profile position (i-1) by the path.  So if a profile position's integer
     * is 0, that profile position was a deletion in the path generating the
     * sequence.  If it is 1, then it was a match.  Any number n greater than 1
     * indicates that there were n-1 insertions following a match at that
     * position.  Index 0 indicates the number of pre-align insertions.  The
     * last index indicates the number of post-align insertions.  Together, the
     * sequence and the sequence_advances convey all of the information in
     * the path.  The return value is the length of the generated sequence
     * (sequence.length()).
     */
    drawSequence (
      Parameters const & parameters,
      ProfileType const & profile,
      Random & random,
      vector<uint32_t> & sequence_advances, // the path to be created
      Sequence<SequenceResidueType> & sequence // the sequence to be created
    ) const
    {
      // TODO: PUT BACK.
      DebugLevel debug = parameters.debug;
      // TODO: REMOVE.
      //DebugLevel debug = DEBUG_All;

      uint32_t profile_length = profile.length();
      if( sequence_advances.size() != ( profile_length + 4 ) ) {
        sequence_advances.resize( profile_length + 4 );
      }
      sequence_advances.assign( ( profile_length + 4 ), 0 );

      SequenceResidueType residue;

      uint32_t seq_adv_i;
      uint32_t last_seq_adv_i = profile_length + 3;
      uint32_t pos_i;
      char current_transition;
      Subcell current_substate;

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, MatrixValueType> match_distribution;
      ProbabilityType del_out_extensions_prob;
      ProbabilityType old_scale, new_scale;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      // Start in PreAlign state.
      current_transition =
        profile[ Transition::fromPreAlign ].draw( random );
      while( current_transition == TransitionFromPreAlign::toPreAlign ) {
        if( debug >= DEBUG_All ) {
          cout << "[drawSequence] current_transition == " << TransitionFromPreAlign::toPreAlign << endl;
        }
        // Count one more preAlign.
        sequence_advances[ 0 ] += 1;
        // Emit from the preAlignInsertion distribution
        residue = profile[ Emission::PreAlignInsertion ].draw( random );
        seqan::fill( sequence, ( sequence.length() + 1 ), residue );
        current_transition =
          profile[ Transition::fromPreAlign ].draw( random );
      } // End while we continue to insert from the preAlign distribution
      if( debug >= DEBUG_All ) {
        cout << "[drawSequence] current_transition == " << TransitionFromPreAlign::toBegin << endl;
      } // End if debug >= DEBUG_All

      // Now we are in the Begin state, since this comes after preAlign.
      current_transition =
        profile[ Transition::fromBegin ].draw( random );
      if( current_transition == TransitionFromBegin::toMatch ) {
        current_substate.setToMatch();
        if( debug >= DEBUG_All ) {
          cout << "[drawSequence] current_transition == " << TransitionFromBegin::toMatch << endl;
        } // End if debug >= DEBUG_All
      } else if( current_transition == TransitionFromBegin::toDeletion ) {
        current_substate.setToDeletion();
        if( debug >= DEBUG_All ) {
          cout << "[drawSequence] current_transition == " << TransitionFromBegin::toDeletion << endl;
        } // End if debug >= DEBUG_All
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      } else if( current_transition == TransitionFromBegin::toDeletionIn ) {
        current_substate.setToDeletionIn();
        if( debug >= DEBUG_All ) {
          cout << "[drawSequence] current_transition == " << TransitionFromBegin::toDeletionIn << endl;
        } // End if debug >= DEBUG_All
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
      } else {
        // ?!
        cerr << "ERROR! current_transition is " << current_transition << ", which is unrecognized!" << endl;
        assert( false );
      }

      pos_i = 0;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      while( current_substate.isDeletionIn() ) {
        // Count one more deletionIn
        sequence_advances[ 1 ] += 1;
        // And advance the profile position.
        pos_i += 1;
        current_transition =
          profile[ Transition::fromDeletionIn ].draw( random );
        if( current_transition == TransitionFromDeletionIn::toDeletionIn ) {
          if( debug >= DEBUG_All ) {
            cout << "[drawSequence] current_transition == " << TransitionFromDeletionIn::toDeletionIn << endl;
          } // End if debug >= DEBUG_All
        } else { // TransitionFromDeletionIn::toMatch
          current_substate.setToMatch();
          if( debug >= DEBUG_All ) {
            cout << "[drawSequence] current_transition == " << TransitionFromDeletionIn::toMatch << endl;
          } // End if debug >= DEBUG_All
        }
      } // End while( current_substate.isDeletionIn() )
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      for( ; pos_i < profile_length; pos_i++ ) {
        seq_adv_i = ( pos_i + 2 );
        if( debug >= DEBUG_All ) {
          cout << "[drawSequence] pos_i == " << pos_i << endl;
        }

        if( current_substate.isMatch() ) {
          // First emit something from the Match distribution.
          // Count one for the pos.
          sequence_advances[ seq_adv_i ] = 1;
          residue = profile[ pos_i ][ Emission::Match ].draw( random );
          seqan::fill( sequence, ( sequence.length() + 1 ), residue );
          if( pos_i < ( profile_length - 1 ) ) {
            // Now see about transitioning..
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            // In this case we need to custom-make a distribution with the
            // right Del-Out prob (that incorporates all of the necessary
            // del-out extensions).
            profile.createScaledMatchDistributionForPosition(
              pos_i,
              match_distribution
            );

            if( debug >= DEBUG_All ) {
              cout << "[drawSequence] match distribution == " << match_distribution << endl;
            } // End if debug >= DEBUG_All
            current_transition =
              match_distribution.draw( random );
#else
            current_transition =
              profile[ pos_i ][ Transition::fromMatch ].draw( random );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT .. else ..
            if( current_transition == TransitionFromMatch::toMatch ) {
              if( debug >= DEBUG_All ) {
                cout << "[drawSequence] current_transition == " << TransitionFromMatch::toMatch << endl;
              } // End if debug >= DEBUG_All
              current_substate.setToMatch();
            } else if( current_transition == TransitionFromMatch::toInsertion ) {
              if( debug >= DEBUG_All ) {
                cout << "[drawSequence] current_transition == " << TransitionFromMatch::toInsertion << endl;
              } // End if debug >= DEBUG_All
              current_substate.setToInsertion();
            } else if( current_transition == TransitionFromMatch::toDeletion ) {
              if( debug >= DEBUG_All ) {
                cout << "[drawSequence] current_transition == " << TransitionFromMatch::toDeletion << endl;
              } // End if debug >= DEBUG_All
              current_substate.setToDeletion();
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            } else if( current_transition == TransitionFromMatch::toDeletionOut ) {
              if( debug >= DEBUG_All ) {
                cout << "[drawSequence] current_transition == " << TransitionFromMatch::toDeletionOut << endl;
              } // End if debug >= DEBUG_All
              current_substate.setToDeletionOut();
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
            } else {
              // ?!
              cerr << "ERROR! current_transition is " << current_transition << ", which is unrecognized!" << endl;
              assert( false );
            }
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
            if( current_substate.isDeletionOut() ) {
              // Okay then all the rest are deletions.  Count the del-outs.
              sequence_advances[ last_seq_adv_i - 1 ] =
                ( profile_length - ( pos_i + 1 ) );
              pos_i = ( profile_length - 1 );
            } // End if isDeletionOut()
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
            while( current_substate.isInsertion() ) {
              // Count one more for the row.
              sequence_advances[ seq_adv_i ]++;
              // Emit from the insertion distribution
              residue = profile[ pos_i ][ Emission::Insertion ].draw( random );
              seqan::fill( sequence, ( sequence.length() + 1 ), residue );
              current_transition =
                profile[ pos_i ][ Transition::fromInsertion ].draw( random );
              if( current_transition == TransitionFromInsertion::toMatch ) {
                current_substate.setToMatch();
              }
              if( debug >= DEBUG_All ) {
                if( current_substate.isMatch() ) {
                  cout << "[drawSequence] current_transition == " << TransitionFromInsertion::toMatch << endl;
                } else {
                  cout << "[drawSequence] current_transition == " << TransitionFromInsertion::toInsertion << endl;
                }
              } // End if debug >= DEBUG_All
            } // End while we continue to insert
            // At this point, current_transition is either toMatch or
            // toDeletion..

            // Advance pos_i in next iteration...
          } // End if( pos_i != ( profile_length - 1 ) )
          // else the transition must be to End...
        } else if( pos_i != ( profile_length - 1 ) ) { // If current_substate.isMatch() .. else ..
          // Current substate is deletion.
          assert( current_substate.isDeletion() );
          current_transition =
            profile[ pos_i ][ Transition::fromDeletion ].draw( random );
          if( current_transition == TransitionFromDeletion::toMatch ) {
            current_substate.setToMatch();
          }
          if( debug >= DEBUG_All ) {
            if( current_substate.isMatch() ) {
              cout << "[drawSequence] current_transition == " << TransitionFromDeletion::toMatch << endl;
            } else {
              cout << "[drawSequence] current_transition == " << TransitionFromDeletion::toDeletion << endl;
            }
          } // End if debug >= DEBUG_All
          // Advance pos_i in next iteration...
        } // End if current_substate.isMatch() .. else if( not last pos ) ..
      } // End foreach pos_i

      assert( pos_i == profile_length );
      // The transition from the last position to End has probability 1.

      // Okay, now we're in the End state.
#ifdef USE_END_DISTRIBUTION
      // For now this is NEVER toLoop.  There is no point in even drawing
      // it.
      //current_transition =
      //  profile[ Transition::fromEnd ].draw( random );
      
      if( debug >= DEBUG_All ) {
        cout << "[drawSequence] current_transition == " << TransitionFromEnd::toPostAlign << endl;
      } // End if debug >= DEBUG_All
#endif // USE_END_DISTRIBUTION

      // In PostAlign state.
      current_transition =
        profile[ Transition::fromPostAlign ].draw( random );
      while( current_transition == TransitionFromPostAlign::toPostAlign ) {
        if( debug >= DEBUG_All ) {
          cout << "[drawSequence] current_transition == " << TransitionFromPostAlign::toPostAlign << endl;
        } // End if debug >= DEBUG_All
        // Count one more for postAligns
        sequence_advances[ last_seq_adv_i ]++;
        // Emit from the postAlignInsertion distribution
        residue = profile[ Emission::PostAlignInsertion ].draw( random );
        seqan::fill( sequence, ( sequence.length() + 1 ), residue );
        current_transition =
          profile[ Transition::fromPostAlign ].draw( random );
      } // End while we continue to insert from the postAlign distribution

      // And now we're done! (in Terminal state)
      if( debug >= DEBUG_All ) {
        cout << "[drawSequence] current_transition == " << TransitionFromPostAlign::toTerminal << endl;
      } // End if debug >= DEBUG_All

      if( debug >= DEBUG_All ) {
        cout << "[drawSequence] Drew sequence of length " << sequence.length() << endl;
        cout << sequence << endl;
      } // End if debug >= DEBUG_All

      return sequence.length();
    } // drawSequence( Parameters const&, ProfileType const&, vector<int> &, Sequence<SequenceResidueType> & ) const

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <typename ProfileTreeType,
            typename SequenceResidueType>
  GALOSH_INLINE_ALGORITHM_INNERLOOP
  void
  DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::
    alignLeaves (
      Parameters const& parameters,
      TreeMultipleAlignment<ProfileTreeType, SequenceResidueType> & multiple_alignment
    ) const
    {
      /// Note that this disregards del-in and del-out counts -- it treats all
      /// deletions the same.

      static const bool be_extra_verbose = false;//true;
      /// \todo Make this a parameter
      static const bool use_DFS_sequence_order = true;

      const uint32_t root_length =
        multiple_alignment.m_profileTree->getProfileTreeRoot()->length();
      //uint32_t num_sequences = multiple_alignment.m_sequences->size();

      const uint32_t node_count = multiple_alignment.m_sequenceIndices->size();
      //cout << "NODE COUNT is " << node_count << endl;
      /// \todo make a multiple_alignment.m_profileTree->nodeCount() method, and assert we get the same result.

      // Calculate all of the profile-profile alignments.
      vector<vector<uint32_t> > profile_profile_alignments( node_count - 1 );

      // If use_sibling_alignments is true, use the m_profileProfileAlignments
      // vector (which is the alignments between each internal node's
      // children).  Otherwise calculate the alignments between children and
      // their parents here.
      const bool use_sibling_alignments = true;

      uint32_t seq_i;
      if( use_sibling_alignments ) {
        uint32_t child_1_vertex;
        uint32_t child_2_vertex;
        uint32_t parent_length;
        uint32_t profile_profile_alignment_index;
        // For each internal node (not leaves), take the alignment between its
        // two children and turn it into alignments between it and its
        // children.
        for( uint32_t parent_node_i = 0; parent_node_i < node_count; parent_node_i++ ) {
          if(
            ( multiple_alignment.m_profileTree->childCount( parent_node_i ) == 0 )
          ) {
            // Leaf.
            continue;
          }
          // Assuming two children.
          assert( multiple_alignment.m_profileTree->childCount( parent_node_i ) == 2 );

          child_1_vertex =
            multiple_alignment.m_profileTree->getChildVertex(
              parent_node_i,
              1
            );
          child_2_vertex =
            multiple_alignment.m_profileTree->getChildVertex(
              parent_node_i,
              2
            );

          parent_length =
            multiple_alignment.m_profileTree->getProfileTreeInternalNode( parent_node_i ).length();
          profile_profile_alignments[ child_1_vertex - 1 ].resize(
            parent_length + 1
          );
          profile_profile_alignments[ child_2_vertex - 1 ].resize(
            parent_length + 1
          );

          // By definition there are no pre-align insertions.
          profile_profile_alignments[ child_1_vertex - 1 ][ 0 ] = 0;
          profile_profile_alignments[ child_2_vertex - 1 ][ 0 ] = 0;
          profile_profile_alignment_index = 1;
          for(
            uint32_t i = 0;
            i < ( *multiple_alignment.m_profileProfileAlignments )[ parent_node_i ][ 0 ];
            i++
          ) {
            // These are all in child_2, but not in child_1.
            profile_profile_alignments[ child_1_vertex - 1 ][ profile_profile_alignment_index ] = 0;
            profile_profile_alignments[ child_2_vertex - 1 ][ profile_profile_alignment_index ] = 1;
            profile_profile_alignment_index += 1;
          } // End foreach 'prealign' relative to child_1.
          for(
            uint32_t child_1_pos_plus1 = 1;
            child_1_pos_plus1 < ( *multiple_alignment.m_profileProfileAlignments )[ parent_node_i ].size();
            child_1_pos_plus1++
          ) {
            if( ( *multiple_alignment.m_profileProfileAlignments )[ parent_node_i ][ child_1_pos_plus1 ] == 0 ) {
              // In child_1 but not child_2.
              profile_profile_alignments[ child_1_vertex - 1 ][ profile_profile_alignment_index ] = 1;
              profile_profile_alignments[ child_2_vertex - 1 ][ profile_profile_alignment_index ] = 0;
              profile_profile_alignment_index += 1;
            } else {
              // It's at least 1, which means it's in both children.
              profile_profile_alignments[ child_1_vertex - 1 ][ profile_profile_alignment_index ] = 1;
              profile_profile_alignments[ child_2_vertex - 1 ][ profile_profile_alignment_index ] = 1;
              profile_profile_alignment_index += 1;
              for(
                uint32_t i = 2;
                i <= ( *multiple_alignment.m_profileProfileAlignments )[ parent_node_i ][ child_1_pos_plus1 ];
                i++
              ) {
                // These are all in child_2, but not in child_1.
                profile_profile_alignments[ child_1_vertex - 1 ][ profile_profile_alignment_index ] = 0;
                profile_profile_alignments[ child_2_vertex - 1 ][ profile_profile_alignment_index ] = 1;
                profile_profile_alignment_index += 1;
              } // End foreach child_2 advance at child_1_pos_plus1 after the 1st..
            } // End if there's 0 advances .. else ..
          } // End foreach child_1_pos_plus1
          // Should be pointing to the end of the newly-created alignments.
          assert( profile_profile_alignment_index == ( parent_length + 1 ) );
        } // End foreach parent_node_i
      } else { // if use_sibling_alignments .. else ..
        /// \todo Make these parameters
        const double indel_open_cost = 1.0;//20.0;//10.0;//.25;
        const double indel_extension_cost = 1.0;//20.0;//10.0;//.25;

        for( uint32_t node_i = 1; node_i < node_count; node_i++ ) {
          if(
            be_extra_verbose &&
            ( multiple_alignment.m_profileTree->childCount( node_i ) == 0 )
          ) {
            // It's a leaf.  Do some debug checking...
            for( seq_i = 0; seq_i < ( *multiple_alignment.m_sequenceIndices )[ node_i ].size(); seq_i++ ) {
              // TODO: REMOVE?  We are assuming that the profile associated with
              // a sequence has the same length as that sequence.
              if( !( multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i ).length() == ( *multiple_alignment.m_sequences )[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ].length() ) ) {
                cout << "Length mismatch in node " << node_i << ", which has length " << multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i ).length() << ", but a sequence (index " << seq_i << " of " << ( *multiple_alignment.m_sequenceIndices )[ node_i ].size() << " in it: sequence " << ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] << ") has length " << ( *multiple_alignment.m_sequences )[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ].length() << endl;
                assert( multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i ).length() == ( *multiple_alignment.m_sequences )[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ].length() );
                cout << "length mismatch in DynamicProgramming::alignLeaves()" << endl;
                exit( 0 );
              } // End if there's a problem.. (testing/debugging)
            }
          } // End if it's a leaf. (debugging)
        
          // Do profile-profile alignment between the two profiles
          profileProfile_align_SKL(
            parameters,
            multiple_alignment.m_profileTree->getParent( node_i ),
            multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i ),
            indel_open_cost,
            indel_extension_cost,
            profile_profile_alignments[ node_i - 1 ],
            ( indel_open_cost * 2 ) // do allow "gap matches"
          );
          if( be_extra_verbose ) {
            cout << "The alignment between node " << node_i << " and its parent, node " << multiple_alignment.m_profileTree->getParentVertex( node_i ) << ", is ( " << profile_profile_alignments[ node_i - 1 ][ 0 ];
            for( uint32_t i = 1; i < profile_profile_alignments[ node_i - 1 ].size(); i++ ) {
              cout << ", " << profile_profile_alignments[ node_i - 1 ][ i ];
            }  // End foreach alignment position
            cout << " )" << endl;
          } // End if be_extra_verbose
        } // End foreach node_i, calculate the profile-profile alignment between that node and its parent.
      } // End if use_sibling_alignments .. else ..

      // Now do a depth-first-traversal, generating root index map information
      // along the way (at non-leaf nodes), and alignments at leaves.
      vector<vector<uint32_t> > root_index_maps;
      root_index_maps.resize( node_count );
      root_index_maps[ 0 ].resize( root_length + 1, 0 );
      for( uint32_t node_i = 1; node_i < node_count; node_i++ ) {
        uint32_t node_length =
          multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i ).length();
        // TODO: REMOVE
        //cout << "node " << node_i << " length is " << node_length << endl;
        root_index_maps[ node_i ].resize( node_length + 1, 0 );
      }

      // For some reason I can't get this to work using bool type!
      vector<vector<uint8_t> > root_index_map_isInsertions;
      //vector<vector<bool> > root_index_map_isInsertions;
      root_index_map_isInsertions.resize( node_count );
      root_index_map_isInsertions[ 0 ].resize( root_length, false );
      for( uint32_t node_i = 1; node_i < node_count; node_i++ ) {
        uint32_t node_length =
          multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i ).length();
        root_index_map_isInsertions[ node_i ].resize( node_length, false );
      }

      // If node_i is not the root, which child of its parent is it? (min 1)
      uint32_t node_i_which_child = 0;
      uint32_t node_i_parent = 0;
      uint32_t node_length = 0, parent_length = 0;
      uint32_t node_pos_i = 0, parent_pos_i = 0;
      int tmp_counter = 0;
      // Do a depth-first traversal
      uint32_t node_i = 0;
      uint32_t leaves_visited = 0;
      do {
        if( be_extra_verbose ) {
          cout << "For node_i " << node_i << ":" << endl;
        }
        // Is it an internal node (or the root)?
        if( multiple_alignment.m_profileTree->childCount( node_i ) > 0 ) {
          if( be_extra_verbose ) {
            cout << "It is an internal node." << endl;
          }
          if( node_i == 0 ) {
            node_length = root_length;
            //root_index_maps[ node_i ].resize( root_length + 1, 0 );
            //root_index_map_isInsertions[ node_i ].resize( root_length, false );
            for( node_pos_i = 0; node_pos_i <= root_length; node_pos_i++ ) {
              root_index_maps[ node_i ][ node_pos_i ] = node_pos_i;
              if( node_pos_i > 0 ) {
                root_index_map_isInsertions[ node_i ][ node_pos_i - 1 ] = false;
              }
            }
          } else { // if node_i is root .. else ..
            node_length =
              multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i ).length();

            // first pos is always pre-align
            root_index_maps[ node_i ][ 0 ] = 0;
            // Assuming node_parent_length is already set.
            node_pos_i = 1;
            for( parent_pos_i = 0;
                 (
                   ( parent_pos_i <= parent_length ) &&
                   ( node_pos_i <= node_length )
                 );
                 parent_pos_i++
            ) {
              if( false && be_extra_verbose ) {
                cout << "[node " << node_i << ", node_pos " << node_pos_i << ", parent_pos " << parent_pos_i << "]: pp alignment value is " << profile_profile_alignments[ node_i - 1 ][ parent_pos_i ] << endl;
              } // End if be_extra_verbose
              if( profile_profile_alignments[ node_i - 1 ][ parent_pos_i ] > 0 ) {
                tmp_counter =
                  profile_profile_alignments[ node_i - 1 ][ parent_pos_i ];
                root_index_maps[ node_i ][ node_pos_i ] =
                  root_index_maps[ node_i_parent ][ parent_pos_i ];
                // First count is an insertion if the parent pos was an
                // insertion.
                root_index_map_isInsertions[ node_i ][ node_pos_i - 1 ] =
                  root_index_map_isInsertions[ node_i_parent ][ parent_pos_i - 1 ];
                node_pos_i += 1;
                while( --tmp_counter > 0 ) {
                  root_index_maps[ node_i ][ node_pos_i ] =
                    root_index_maps[ node_i_parent ][ parent_pos_i ];
                  // Subsequent counts are always insertions
                  root_index_map_isInsertions[ node_i ][ node_pos_i - 1 ] =
                    true;
                  node_pos_i += 1;
                } // End while assigning more of node_i's positions to the same
                  // root position as the parent is assigned to at
                  // parent_pos_i.
              } // End if parent's p-p alignment value is non-zero at parent_pos_i.
            } // End foreach parent_pos_i
          } // End if node_i is root .. else ..

          if( be_extra_verbose ) {
            cout << "[node " << node_i << "] The root index map is (";
            cout << root_index_maps[ node_i ][ 0 ];
            for( node_pos_i = 1; node_pos_i <= node_length; node_pos_i++ ) {
              cout << "," << root_index_maps[ node_i ][ node_pos_i ];
            }
            cout << ")" << endl;
            cout << "[node " << node_i << "] The root index map (isInsertion) is (";
            cout << ( ( root_index_map_isInsertions[ node_i ][ 0 ] == true ) ? "T" : "F" );
            for( node_pos_i = 2; node_pos_i <= node_length; node_pos_i++ ) {
              cout << "," << ( ( root_index_map_isInsertions[ node_i ][ node_pos_i - 1 ] == true ) ? "T" : "F" );
            }
            cout << ")" << endl;
          } // End if be_extra_verbose
          // Walk down.
          node_i_which_child = 1;
          node_i_parent = node_i;
          parent_length = node_length;
          node_i =
            multiple_alignment.m_profileTree->getChildVertex( node_i, node_i_which_child );
          if( be_extra_verbose ) {
            cout << "After decending from " << node_i_parent << ", node_i is now " << node_i << endl;
          } // End if be_extra_verbose
        } else { // if this is an internal node .. else ..
          if( be_extra_verbose ) {
            cout << "It is a leaf." << endl;
          }
          // A leaf.  Make the alignment.
          // Note node_i > 0 (since 0 is the root).
          assert( node_i > 0 );

          if( use_DFS_sequence_order ) {
            for( seq_i = 0; seq_i < ( *multiple_alignment.m_sequenceIndices )[ node_i ].size(); seq_i++ ) {
              multiple_alignment.m_sequenceOrder[ leaves_visited ] =
                ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ];
              if( be_extra_verbose ) {
                cout << "Sequence " << ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] << " at node " << node_i << " gets to be in position " << leaves_visited << " of the sequenceOrder." << endl;
              }

              leaves_visited += 1;
            } // End foreach seq_i
          } // End if use_DFS_sequence_order

          for( parent_pos_i = 0; parent_pos_i <= parent_length; parent_pos_i++ ) {
            if( false && be_extra_verbose ) {
              cout << "[node " << node_i << ", parent_pos " << parent_pos_i << "]: pp alignment value is " << profile_profile_alignments[ node_i - 1 ][ parent_pos_i ] << endl;
            } // End if be_extra_verbose
            if( profile_profile_alignments[ node_i - 1 ][ parent_pos_i ] == 0 ) {
              // If it's not an insertion cell, note that we have a deletion.
              if(
                !(
                  ( parent_pos_i == 0 ) ||
                  ( root_index_maps[ node_i_parent ][ parent_pos_i ] == 0 ) ||
                  ( root_index_map_isInsertions[ node_i_parent ][ parent_pos_i - 1 ] == true )
                )
              ) {
                for( seq_i = 0; seq_i < ( *multiple_alignment.m_sequenceIndices )[ node_i ].size(); seq_i++ ) {
                  multiple_alignment.m_insertionsAfterPosition[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ root_index_maps[ node_i_parent ][ parent_pos_i ] ] = 0;
                  multiple_alignment.m_matchIndicators[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ root_index_maps[ node_i_parent ][ parent_pos_i ] - 1 ] = false;
                } // End foreach sequence assigned to the leaf node_i
              } // End if it's not an insertion cell, it's a deletion.
              continue;
            } // End if there's nothing aligned here, move along.
            // It's either a count that's predetermined to be an insertion (eg
            // it's pre-align, or it's been designated by the
            // root_index_map_isInsertions flag) ... or it's not ..
            if(
              ( parent_pos_i == 0 ) ||
              ( root_index_maps[ node_i_parent ][ parent_pos_i ] == 0 ) ||
              ( root_index_map_isInsertions[ node_i_parent ][ parent_pos_i - 1 ] == true )
            ) {
              for( seq_i = 0; seq_i < ( *multiple_alignment.m_sequenceIndices )[ node_i ].size(); seq_i++ ) {
                multiple_alignment.m_insertionsAfterPosition[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ root_index_maps[ node_i_parent ][ parent_pos_i ] ] +=
                  profile_profile_alignments[ node_i - 1 ][ parent_pos_i ];
              } // End foreach sequence assigned to the leaf node_i
            } else { // if this is an insertion count .. else ..
              for( seq_i = 0; seq_i < ( *multiple_alignment.m_sequenceIndices )[ node_i ].size(); seq_i++ ) {
                if( profile_profile_alignments[ node_i - 1 ][ parent_pos_i ] == 0 ) {
                  multiple_alignment.m_matchIndicators[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ root_index_maps[ node_i_parent ][ parent_pos_i ] - 1 ] = false;
                  multiple_alignment.m_insertionsAfterPosition[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ root_index_maps[ node_i_parent ][ parent_pos_i ] ] = 0;
                } else {
                  multiple_alignment.m_matchIndicators[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ root_index_maps[ node_i_parent ][ parent_pos_i ] - 1 ] = true;
  
                  // Ok, but maybe there's also insertions.
                  if( profile_profile_alignments[ node_i - 1 ][ parent_pos_i ] > 1 ) {
                    for( seq_i = 0; seq_i < ( *multiple_alignment.m_sequenceIndices )[ node_i ].size(); seq_i++ ) {
                      multiple_alignment.m_insertionsAfterPosition[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ root_index_maps[ node_i_parent ][ parent_pos_i ] ] +=
                        ( profile_profile_alignments[ node_i - 1 ][ parent_pos_i ] - 1 );
                    } // End foreach sequence assigned to the leaf node_i
                  } // End if there's any insertions following the Match.
                } // End if this is a Match to the root, maybe followed by insertions..
              } // End foreach sequence assigned to this leaf node
            } // End if this is an insertion count .. else ..
          } // End foreach parent_pos_i

          if( be_extra_verbose ) {
            for( seq_i = 0; seq_i < ( *multiple_alignment.m_sequenceIndices )[ node_i ].size(); seq_i++ ) {
              cout << "Aligned sequence " << ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] << " at node " << node_i << ":" << endl;
              cout << ( *multiple_alignment.m_sequences )[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ] << endl;
              cout << " ";
              for( parent_pos_i = 1; // actually root pos.
                   parent_pos_i <= root_length;
                   parent_pos_i++
              ) {
                cout << ( ( multiple_alignment.m_matchIndicators[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ parent_pos_i - 1 ] == true ) ? "T" : "F" );
              }
              cout << endl;
              for( parent_pos_i = 0; // actually root pos with pre-align offset
                   parent_pos_i <= root_length;
                   parent_pos_i++
              ) {
                cout << multiple_alignment.m_insertionsAfterPosition[ ( *multiple_alignment.m_sequenceIndices )[ node_i ][ seq_i ] ][ parent_pos_i ];
              }
              cout << endl;
            } // End foreach seq at this node
          } // End if be_extra_verbose

          // Move along
          // Don't descend further here.
          while( ( node_i != 0 ) && ( node_i_which_child == multiple_alignment.m_profileTree->childCount( node_i_parent ) ) ) {
            if( be_extra_verbose ) {
              cout << "This is the last child of its parent.  Ascending tree." << endl;
            }
            if( be_extra_verbose ) {
              cout << "node_i is currently " << node_i << endl;
              cout << "node_i_parent is currently " << node_i_parent << endl;
              cout << "node_i_which_child is currently " << node_i_which_child << endl;
            } // End if be_extra_verbose

            // We're at the last child of its parent.  Move up.
            node_i = node_i_parent;
            node_i_parent = multiple_alignment.m_profileTree->getParentVertex( node_i );
            if( node_i_parent == 0 ) {
              parent_length = root_length;
            } else {
              parent_length =
                multiple_alignment.m_profileTree->getProfileTreeInternalNode( node_i_parent ).length();
            }
            node_i_which_child =
              multiple_alignment.m_profileTree->getChildIndexInParent( node_i_parent, node_i );
          } // While we're at the last child of our parent, move up.
          if( node_i == 0 ) {
            // Can't move up.  Done.
            break;
          }
          if( be_extra_verbose ) {
            cout << "Moving to next sibling." << endl;
          }
          // Move along to the next child.
          node_i_which_child += 1;
          node_i =
            multiple_alignment.m_profileTree->getChildVertex( node_i_parent, node_i_which_child );
          if( be_extra_verbose ) {
            cout << "node_i is now " << node_i << endl;
          } // End if be_extra_verbose
        } // End if this is an internal node .. else ..
      } while( true ); // We break when done.
 
      return;
    } // alignLeaves( Parameters const &, ProfileTreeType const &, MultipleAlignment & )


    /**
     * \class AlignmentProfileAccessor
     * \author Ted Holzman
     * \date 2/29/2012
     *
     * This is a class to access the inner class DynamicProgramming::AlignmentProfile.
     * It contains input and accessor functions.
     */

    namespace io = boost::iostreams;
    using namespace std;
    using namespace galosh;
    namespace bx = boost::xpressive;
    using namespace bx;

    template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
    class
    AlignmentProfileAccessor :
       public DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::AlignmentProfile
    {
       private:
          int m_orig_nseq;

       public:
          AlignmentProfileAccessor () {
             m_orig_nseq = 0;
          }

          AlignmentProfileAccessor(int nseq) : m_orig_nseq(nseq)
          {}     

          typedef ResidueType APAResidueType;
          /**
           * \var std::map<std::string,std::string> > kvPairs;
           *
           * When parsing an Alignment Profile, there may be lines whose first non-blank character is a pound-sign.
       	   * Everything after the pound-sign is ignored <em> in so far as the Alignment Profile proper</em> is concerned.
           * However, the comment lines may contain key-value pairs of the form <key>=<value>.  These can be accessed
           * later for any purpose.  Our first use is to supply the n-of-sequences which went into the generation of
           * the alignment profile, in a format like this:
           *
           *    nseq=10
           */

          std::map<std::string,std::string> kvPairs;
          /**
           * Stream reader/parsing routines
           *
           * \operator >>
           * Clear current AlignmentProfile.  Read in a new one.
           */
           friend std::istream &
           operator>> (
              std::istream & normalStream,
              AlignmentProfileAccessor<ResidueType, ProbabilityType, ScoreType, MatrixValueType> & prof
           )
           {
              static io::filtering_istream is;
              galosh::input_comment_diversion_filter * icdf;
              is.push(galosh::input_comment_diversion_filter());
              icdf = is.component<input_comment_diversion_filter>((int)(is.size()-1));
              is.push(normalStream);

              /**
               * Don't even begin if the input stream is at EOF or in a failed state.
               * Otherwise, clear the current vector and iterate through the stream.  Each
               * line is an <AlignmentProfilePosition> in square brackets.  It consists of
               * MatchEmissionParameters followed by GlobalParameters
               */
               assert( !is.fail() && !is.eof());
               prof.clear();
               int lcount = 0;  //debug and/or error reporting
               char c;
               while( !is.eof()) {
            	  typename DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::AlignmentProfilePosition curPosition;
            	      //is >> "[ ";
            	  c = is.get();
            	  if(c == '\n' || c == -1) continue;
                  assert(c == '[');
                  assert(is.get() == ' ');
                  MatchEmissionParameters<ResidueType,ProbabilityType> tmpMEP;
                  tmpMEP.readMatchEmissionParameters(is);
                  curPosition.MatchEmissionParameters<ResidueType,ProbabilityType>::copyFrom(tmpMEP);
                  assert(is.get() ==',');
                  assert(is.get() ==' ');
                  GlobalParameters<ResidueType,ProbabilityType> tmpGP;
                  tmpGP.readGlobalParameters(is);
                  curPosition.GlobalParameters<ResidueType,ProbabilityType>::copyFrom(tmpGP);
                  c = is.get();
                  if(c != ' ') {
                     std::cerr << "Problem parsing alignment profile on line " << lcount+1 << std::endl;
                     std::cerr << "Expected ' ]', saw: " << std::endl;
                     do {
                        std::cerr << c << std::endl;
                     } while ((c = is.get()) != '\n');
                  }
                  prof.push_back(curPosition);
                  ++lcount;
                  is.ignore( 100000, '\n' );
              }

              ///TAH 3/12  debug code:
              //std::cerr << "Read " << ++lcount << " alignment profile lines.\n"; std::cerr.flush();
              /**
               * having read in and created an AlignmentProfile, make sure it the same when we
               * output it again
               */
              /*
                std::ofstream testout("test.profile.out");
                testout << static_cast<typename DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::AlignmentProfile>(prof);
                testout.close();
              */
              ///TAH 3/12 end debug code

              /**
               * For detailed information on the "placeholder" thing, see the documentation for xpressive.  Basically it is a way
               * of allowing the semantic portion of a xpressive regular expression (the part in braces) to access a variable that
               * has no lvalue at compile time.
               */

              placeholder<std::map<std::string,std::string> > dummy;
              // match a quoted string
              sregex quoteStr = as_xpr('"') >> *~(as_xpr('"')) >> '"';
              // match a kv pair, where the v can be a quoted string
              sregex pair = -*_ >> ( (s1= +_w) >> "=" >> (s2= (+_w | quoteStr)) )
                 [ dummy[s1] = bx::as<std::string>(s2) ];
              // Match one or more token=token pairs, separated
              // by anything
              sregex pairs = pair >> *(-+_ >> pair);
              bx::smatch what;
              what.let(dummy = prof.kvPairs);
              regex_search(icdf->get_raw_comments(),what,pairs);
              /// TAH 4/12 debug code:
              /*
              std::cerr << "Contents of raw comment string\n";
              std::cerr << icdf->get_raw_comments() << std::endl; cerr.flush();
              std::cerr << "Contents of key/value map\n";
              for(map<string,string>::iterator  it = prof.kvPairs.begin(); it != prof.kvPairs.end(); it++ ) {
                 std::cerr << "key: " << it->first << ", val: " << it->second << endl;
              }
              std::cerr.flush();
              */
              //end debug code
              /// Side effect. The idea here:  if somebody pre-set the number of original sequences (before the alignment
              /// profile is read in, then that number will take precedence over any comment in the sequence.
              /// Otherwise, nseq=10 will set the number of aligned sequences to 10.
              if(prof.kvPairs.count("nseq")) {
            	  if(prof.m_orig_nseq == 0) prof.m_orig_nseq = atoi(prof.kvPairs["nseq"].c_str());
              }

              /**
               * \note - it is an interesting question whether we should return \a is (the filtering iStream) or \a normalStream
               * the non-filtering one. I think that, for the time being, it is fine to return either.  It may be useful, in the
               * future to capture the returned filtering_istream, get its input_comment_diversion_filter, and get the raw comments
               * from it.  But currently we're only interested in the key=value pairs, which are made available in kvPairs.
               **/
               return is;
           }; // operator>>(std::istream,alignmentprofile)

    	   bool
           fromFile (
    	      std::istream & is,
              AlignmentProfileAccessor<ResidueType, ProbabilityType, ScoreType, MatrixValueType> & prof
           ) 
           {
              is >> prof;
              return true;
           }  // fromFile(stream, AlignmentProfile) wrapper around >>

           bool
           fromFile (
              char *fn,
              AlignmentProfileAccessor<ResidueType, ProbabilityType, ScoreType, MatrixValueType> & prof
           )
           {
               std::ifstream profileFile(const_cast<char *>(fn));
        	   return fromFile(profileFile,prof);
    	   } // fromFile(stream, AlignmentProfile) wrapper around >>

           /// \todo check that the next two functions work when they're
           /// actually compiled.

           bool
           fromFile (const char *fn)
           {
              return fromFile(fn,this);
           } // fromFile(char *) Read this AlignmentProfile

           bool
           fromFile (std::istream & is) {
              return fromFile(is,this);
           } // fromFile(std::istream) Read this AlignmentProfile

           int
           orig_nseq () { // access original number of sequences in alignment
              return m_orig_nseq;
           } // orig_nseq
   }; // AlignmentProfileAccessor

} // End namespace galosh

#endif // __GALOSH_DYNAMICPROGRAMMING_HPP__
