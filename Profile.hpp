/**
 * \file Profile.hpp
 * \author  D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 * \brief
 *      Class definition for the Galosh Profile HMM class.  
 *
 *      A Profile is a data
 *      structure for the model parameters.  Conceptually, for every position
 *      of the profile there are a bunch of parameters for each kind of Plan 7
 *      profile transition/emission.  We make this a vector of maps from
 *      ProfileKeys (eg. Transition) to maps from ProfileKey-specific params
 *      (eg. Transition_M_to_D) to the values of those params.
 *
 *      Also defines dirichlet priors for each set of profile parameters.
 *
 *      More about Plan7 can be found at
 *      http://www.csb.yale.edu/userguides/seq/hmmer/docs/node11.html
 *      (or any other HMMER docs mirror, node 11)
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

#ifndef __GALOSH_PROFILE_HPP__
#define __GALOSH_PROFILE_HPP__

#include "Prolific.hpp"

#include "ProfileHMM.hpp"
#include "Parameters.hpp"
#include "MultinomialDistribution.hpp"
#include "Sequence.hpp"

#include <string>
using std::string;
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;
#include <sstream>
#include <fstream>
#include <map>
using std::map;
#include <vector>
using std::vector;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <boost/lexical_cast.hpp>
// TODO: REMOVE
//#include <cassert>

#include <seqan/basic.h>

/// Note that in the ProfileTreeRoot code we differentiate between size(), which returns the 
/// actual underlying size of the vector that holds the ProfilePositions, and length(), which
/// might be overridden to return something other than size() (as it is in the case of the 
/// PositionEntente in DynamicProgramming.hpp).  Also for this reason we differentiate between 
/// operator[] (which might be overridden) and vector<ProfilePosition<ResidueType, ProbabilityType> >\::operator[]
/// (which always returns the corresponding index into the underlying array).

namespace galosh {

  // TODO: I want to support global parameter sets as well as missing values.

  struct profile_distribution_tag
  {
    // Identifying tag for profile distributions.
  };
  struct profile_transition_distribution_tag : public profile_distribution_tag
  {
    // Identifying tag for transition distributions.
  };
  struct profile_Match_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the Match state.
  };
  struct profile_Insertion_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the Insertion state.
  };
  struct profile_Deletion_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the Deletion state.
  };
  struct profile_PreAlign_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the PreAlign (N) state.
  };
  struct profile_Begin_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the Begin (B) state.
  };
#ifdef USE_END_DISTRIBUTION
  struct profile_End_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the End (E) state.
  };
#endif // USE_END_DISTRIBUTION
  struct profile_PostAlign_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the PostAlign (C) state.
  };
  struct profile_Loop_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the Loop (J) state.
  };
#ifdef USE_DEL_IN_DEL_OUT
  struct profile_DeletionIn_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the DeletionIn (Z) state.
  };
  struct profile_DeletionOut_distribution_tag :
    public profile_transition_distribution_tag
  {
    // Identifying tag for the transition distribution from the DeletionOut (W) state.
  };
#endif // USE_DEL_IN_DEL_OUT

  struct profile_emission_distribution_tag : public profile_distribution_tag
  {
    // Identifying tag for emission distributions.
  };
  struct profile_Match_emission_distribution_tag : public profile_emission_distribution_tag
  {
    // Identifying tag for match emission distributions.
  };
  struct profile_Insertion_emission_distribution_tag : public profile_emission_distribution_tag
  {
    // Identifying tag for Insertion emission distributions.
  };
  struct profile_PreAlign_emission_distribution_tag :
    public profile_emission_distribution_tag
  {
    // Identifying tag for PreAlign Insertion emission distributions.
  };
  struct profile_PostAlign_emission_distribution_tag :
    public profile_emission_distribution_tag
  {
    // Identifying tag for PostAlign Insertion emission distributions.
  };

  class Emission : public profile_emission_distribution_tag
  {
  public:
    static profile_Match_emission_distribution_tag const Match;
    static profile_Insertion_emission_distribution_tag const Insertion;
    static profile_PreAlign_emission_distribution_tag const PreAlignInsertion;
    static profile_PostAlign_emission_distribution_tag const PostAlignInsertion;
  
    Emission ()
    {
      // Do nothing.
    } // <init>()
  }; // End class galosh::Emission

const profile_Match_emission_distribution_tag Emission::Match =
  profile_Match_emission_distribution_tag();
const profile_Insertion_emission_distribution_tag Emission::Insertion =
  profile_Insertion_emission_distribution_tag();
const profile_PreAlign_emission_distribution_tag Emission::PreAlignInsertion =
  profile_PreAlign_emission_distribution_tag();
const profile_PostAlign_emission_distribution_tag Emission::PostAlignInsertion =
  profile_PostAlign_emission_distribution_tag();

  //////////////
  // These are holdovers to keep the old DynamicProgramming.hpp stuff happy.
  class TransitionFromMatch
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type toMatch, toInsertion, toDeletion;
#ifdef USE_DEL_IN_DEL_OUT
    static const galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type toDeletionOut;
#endif // USE_DEL_IN_DEL_OUT
  }; // End class galosh::TransitionFromMatch
  
  galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type const
    TransitionFromMatch::toMatch( 'M' ),
    TransitionFromMatch::toInsertion( 'I' ),
    TransitionFromMatch::toDeletion( 'D' );
#ifdef USE_DEL_IN_DEL_OUT
  galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type const
    TransitionFromMatch::toDeletionOut( 'W' );
#endif // USE_DEL_IN_DEL_OUT

  class TransitionFromInsertion
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::InsertionStateLabel, galosh::Plan7>::Type toMatch, toInsertion;
  }; // End class galosh::TransitionFromInsertion
  
  galosh::StateLabelTransitionTargets<galosh::InsertionStateLabel, galosh::Plan7>::Type const
    TransitionFromInsertion::toMatch( 'M' ),
    TransitionFromInsertion::toInsertion( 'I' );

  class TransitionFromDeletion
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::DeletionStateLabel, galosh::Plan7>::Type toMatch, toDeletion;
  }; // End class galosh::TransitionFromDeletion
  
  galosh::StateLabelTransitionTargets<galosh::DeletionStateLabel, galosh::Plan7>::Type const
    TransitionFromDeletion::toMatch( 'M' ),
    TransitionFromDeletion::toDeletion( 'D' );

  class TransitionFromPreAlign
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::PreAlignStateLabel, galosh::Plan7>::Type toPreAlign, toBegin;
  }; // End class galosh::TransitionFromPreAlign
  
  galosh::StateLabelTransitionTargets<galosh::PreAlignStateLabel, galosh::Plan7>::Type const
    TransitionFromPreAlign::toPreAlign( 'N' ),
    TransitionFromPreAlign::toBegin( 'B' );

  class TransitionFromBegin
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::BeginStateLabel, galosh::Plan7>::Type toMatch, toDeletion;
#ifdef USE_DEL_IN_DEL_OUT
    static const galosh::StateLabelTransitionTargets<galosh::BeginStateLabel, galosh::Plan7>::Type toDeletionIn;
#endif // USE_DEL_IN_DEL_OUT
  }; // End class galosh::TransitionFromBegin
  
  galosh::StateLabelTransitionTargets<galosh::BeginStateLabel, galosh::Plan7>::Type const
    TransitionFromBegin::toMatch( 'M' ),
    TransitionFromBegin::toDeletion( 'D' );
#ifdef USE_DEL_IN_DEL_OUT
  galosh::StateLabelTransitionTargets<galosh::BeginStateLabel, galosh::Plan7>::Type const
    TransitionFromBegin::toDeletionIn( 'Z' );
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_END_DISTRIBUTION
  class TransitionFromEnd
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::EndStateLabel, galosh::Plan7>::Type toPostAlign, toLoop;
  }; // End class galosh::TransitionFromEnd

  galosh::StateLabelTransitionTargets<galosh::EndStateLabel, galosh::Plan7>::Type const
    TransitionFromEnd::toPostAlign( 'C' ),
    TransitionFromEnd::toLoop( 'J' );
#endif // USE_END_DISTRIBUTION

  class TransitionFromPostAlign
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::PostAlignStateLabel, galosh::Plan7>::Type toPostAlign, toTerminal;
  }; // End class galosh::TransitionFromPostAlign
  
  galosh::StateLabelTransitionTargets<galosh::PostAlignStateLabel, galosh::Plan7>::Type const
    TransitionFromPostAlign::toPostAlign( 'C' ),
    TransitionFromPostAlign::toTerminal( 'T' );

#ifdef USE_DEL_IN_DEL_OUT
  class TransitionFromDeletionIn
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::DeletionInStateLabel, galosh::Plan7>::Type toDeletionIn, toMatch;
  }; // End class galosh::TransitionFromDeletionIn
  
  galosh::StateLabelTransitionTargets<galosh::DeletionInStateLabel, galosh::Plan7>::Type const
    TransitionFromDeletionIn::toDeletionIn( 'Z' ),
    TransitionFromDeletionIn::toMatch( 'M' );

  class TransitionFromDeletionOut
  {
    public:
    static const galosh::StateLabelTransitionTargets<galosh::DeletionOutStateLabel, galosh::Plan7>::Type toDeletionOut, toEnd;
  }; // End class galosh::TransitionFromDeletionOut
  
  galosh::StateLabelTransitionTargets<galosh::DeletionOutStateLabel, galosh::Plan7>::Type const
    TransitionFromDeletionOut::toDeletionOut( 'W' ),
    TransitionFromDeletionOut::toEnd( 'E' );
#endif // USE_DEL_IN_DEL_OUT

  // End TransitionFrom ...
  //////////////

  class Transition : public profile_transition_distribution_tag
  {
#define TRANSITION_CONSTS( state ) \
    class From ## state : public profile_ ## state ## _distribution_tag \
    { \
    public: \
      From ## state () \
      { \
      } /* <init>() */ \
    }; /* End inner class Transition::From ## state */ \
    static From ## state const from ## state;
  
  public:
    TRANSITION_CONSTS( Match )
    TRANSITION_CONSTS( Insertion )
    TRANSITION_CONSTS( Deletion )
    TRANSITION_CONSTS( PreAlign )
    TRANSITION_CONSTS( Begin )
#ifdef USE_END_DISTRIBUTION
    TRANSITION_CONSTS( End )
#endif // USE_END_DISTRIBUTION
    TRANSITION_CONSTS( PostAlign )
    TRANSITION_CONSTS( Loop )
#ifdef USE_DEL_IN_DEL_OUT
    TRANSITION_CONSTS( DeletionIn )
    TRANSITION_CONSTS( DeletionOut )
#endif // USE_DEL_IN_DEL_OUT
  
    Transition ()
    {
      // Do nothing.
    } // <init>()
  }; // End class galosh::Transition

  /**
   * This static template fn applies to anything that supports the operator[](
   * tag ) method.  Calculate and return the probability that if
   * you drew one emission from each object, the emission from the given first
   * object would be identical to the emission from the given second object,
   * using the emission values of the two objects as the relevant
   * probabilities.  Note that this is symmetric, so the order of the arguments
   * does not matter.
   */
  template <typename ProbabilityType,
            typename EmissionContainerType1,
            typename EmissionContainerType2,
            typename EmissionType>
  static inline
  ProbabilityType
  calculateEmissionCooccurrenceProbability (
    EmissionContainerType1 const & obj1,
    EmissionContainerType2 const & obj2,
    EmissionType const & tag
  )
  {
    return
      obj1[ tag ].calculateExpectedValue(
        obj2[ tag ]
      );
  } // static template ProbabilityType calculateEmissionCooccurrenceProbability ( EmissionContainerType1 const&, EmissionContainerType2 const&, EmissionType const & )

  // Forward declarations
  //template <typename InternalNodeOrRoot> class ProfileTreeNode;
  template <typename ResidueType, typename ProbabilityType> class ProfileTreeRoot;
  template <typename ResidueType, typename ProbabilityType> class ProfileTreeInternalNode;
  template <typename ProfileType> struct profile_traits;
  template <typename ScoreType, typename ParameterCollectionType> class ScalableParameterCollection;

// Declarations and accessors for the transition parameter collections
#define GALOSH_OPERATOR_DECLARE_TRANSITION( probability_type, state ) \
    /*protected:                                                        */ \
    public: \
    galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets< state ## StateLabel, galosh::Plan7>::Type, probability_type> m_ ## state ## _Distribution; \
    public: \
    galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets< state ## StateLabel, galosh::Plan7>::Type, probability_type> &        \
    operator[] ( profile_ ## state ## _distribution_tag const ) \
    { \
      return m_ ## state ## _Distribution; \
    } /* operator[] ( profile_ ## state ## _distribution_tag const ) */ \
    galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets< state ## StateLabel, galosh::Plan7>::Type, probability_type> const&   \
    operator[] ( profile_ ## state ## _distribution_tag const ) const\
    { \
      return m_ ## state ## _Distribution; \
    } /* operator[] ( profile_ ## state ## _distribution_tag const ) */ \

// Accessors for derived classes of the transition parameter collections
#define GALOSH_OPERATOR_TRANSITION( collection_type, probability_type, state ) \
    public: \
    galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets< state ## StateLabel, galosh::Plan7>::Type, probability_type> & \
    operator[] ( profile_ ## state ## _distribution_tag const tag ) \
    { \
      return this->collection_type<ResidueType, probability_type>::operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _distribution_tag const ) */ \
    galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets< state ## StateLabel, galosh::Plan7>::Type, probability_type> const& \
    operator[] ( profile_ ## state ## _distribution_tag const tag ) const\
    { \
      return this->collection_type<ResidueType, probability_type>::operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _distribution_tag const ) */ \

// Accessors for classes delegating to a transition parameter collection
#define GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, state ) \
    public: \
    galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets< state ## StateLabel, galosh::Plan7>::Type, probability_type> & \
    operator[] ( profile_ ## state ## _distribution_tag const tag ) \
    { \
      return obj .operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _distribution_tag const ) */ \
    galosh::MultinomialDistribution<galosh::StateLabelTransitionTargets< state ## StateLabel, galosh::Plan7>::Type, probability_type> const& \
    operator[] ( profile_ ## state ## _distribution_tag const tag ) const\
    { \
      return obj .operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _distribution_tag const ) */ \

// Declarations and accessors for the emission parameter collections
#define GALOSH_OPERATOR_DECLARE_EMISSION( state, probability_type ) \
  /* protected: */                                                      \
    public: \
    MultinomialDistribution<ResidueType, probability_type> m_ ## state ## _Emission_Distribution; \
  GALOSH_OPERATOR_DECLARE_EMISSION_OPERATORS( state, state, probability_type )

#define GALOSH_OPERATOR_DECLARE_EMISSION_OPERATORS( tag_state, dist_state, probability_type ) \
    public: \
    galosh::MultinomialDistribution<ResidueType, probability_type> & \
    operator[] ( profile_ ## tag_state ## _emission_distribution_tag const ) \
    { \
      return m_ ## dist_state ## _Emission_Distribution; \
    } /* operator[] ( profile_ ## tag_state ## _emission_distribution_tag const ) */ \
    galosh::MultinomialDistribution<ResidueType, probability_type> const& \
    operator[] ( profile_ ## tag_state ## _emission_distribution_tag const ) const\
    { \
      return m_ ## dist_state ## _Emission_Distribution; \
    } /* operator[] ( profile_ ## tag_state ## _emission_distribution_tag const ) */ \

// TODO: REMOVE?  TESTING.
#define GALOSH_OPERATOR_DECLARE_EMISSION_USEINSERTION( state, probability_type ) \
  /* protected: */                                                      \
    public: \
    MultinomialDistribution<ResidueType, probability_type> m_ ## state ## _Emission_Distribution;

// Accessors for derived classes of the emission parameter collections
#define GALOSH_OPERATOR_EMISSION( collection_type, state, probability_type ) \
  public: \
    galosh::MultinomialDistribution<ResidueType, probability_type> & \
    operator[] ( profile_ ## state ## _emission_distribution_tag const tag ) \
    { \
      return this->collection_type<ResidueType, probability_type>::operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _emission_distribution_tag const ) */ \
    galosh::MultinomialDistribution<ResidueType, probability_type> const& \
    operator[] ( profile_ ## state ## _emission_distribution_tag const tag ) const\
    { \
      return this->collection_type<ResidueType, probability_type>::operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _emission_distribution_tag const ) */ \

// Accessors for classes delegating to an emission parameter collection
#define GALOSH_OPERATOR_DELEGATING_EMISSION( obj, state, probability_type ) \
  public: \
    galosh::MultinomialDistribution<ResidueType, probability_type> & \
    operator[] ( profile_ ## state ## _emission_distribution_tag const tag ) \
    { \
      return obj .operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _emission_distribution_tag const ) */ \
    galosh::MultinomialDistribution<ResidueType, probability_type> const& \
    operator[] ( profile_ ## state ## _emission_distribution_tag const tag ) const\
    { \
      return obj .operator[]( tag ); \
    } /* operator[] ( profile_ ## state ## _emission_distribution_tag const ) */ \

  /**
   * These are the emission parameters for the Match distribution.
   */
  template <typename ResidueType, typename ProbabilityType>
  class MatchEmissionParameters
  {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_NVP( m_Match_Emission_Distribution );
    } // serialize( Archive &, const unsigned int )

  protected:

// Accessors and declarations for emissions
    GALOSH_OPERATOR_DECLARE_EMISSION( Match, ProbabilityType )

// Accessors for derived classes
#define GALOSH_OPERATORS_MatchEmissionParameters( probability_type ) \
    GALOSH_OPERATOR_EMISSION( MatchEmissionParameters, Match, probability_type ) \

// Accessors for delegating classes
#define GALOSH_OPERATORS_DELEGATING_MatchEmissionParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, Match, probability_type ) \


  public:

    MatchEmissionParameters ();

    template <typename AnyProbabilityType>
    MatchEmissionParameters (
      MatchEmissionParameters<ResidueType, AnyProbabilityType> const & copy_from
    );

    MatchEmissionParameters (
      MatchEmissionParameters const & copy_from
    );

    void
    reinitialize ();

    template <typename AnyProbabilityType>
    MatchEmissionParameters<ResidueType, ProbabilityType> &
    operator= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    MatchEmissionParameters<ResidueType, ProbabilityType> &
    operator= ( MatchEmissionParameters<ResidueType, ProbabilityType> const& other_pos );

    template <typename AnyProbabilityType>
    void
    copyFrom ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Divide each contained distribution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    MatchEmissionParameters &
    operator/= ( AnyProbabilityType const& denominator );

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    MatchEmissionParameters<ResidueType, ProbabilityType> &
    operator*= ( AnyProbabilityType const& scalar );

    /**
     * Add to each contained distribution the values in the given other
     * MatchEmissionParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    MatchEmissionParameters<ResidueType, ProbabilityType> &
    operator+= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Subtract from each contained distribution the values in the given other
     * MatchEmissionParameters.  Note that this may violate the rule that the
     * probabilities are greater than 0.
     */
    template <typename AnyProbabilityType>
    MatchEmissionParameters<ResidueType, ProbabilityType> &
    operator-= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).  Note that
     * the cross entropy is non-symmetric (calling other_pos.crossEntropy(
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
    template <typename AnyProbabilityType>
    double
    crossEntropy (
      MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos
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
    template <typename AnyProbabilityType,
              typename AnyMatchEmissionParametersType>
    double
    crossEntropy (
      MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyMatchEmissionParametersType const * const weights
    ) const;

    /**
     * Calculate and return the Euclidean distance between this set of emission
     * parameters and another set (treating every probability as an orthogonal
     * dimension).
     */
    double
    euclideanDistance (
      MatchEmissionParameters const& other_pos
    ) const;

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of emission parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    double
    euclideanDistanceSquared (
      MatchEmissionParameters const& other_pos
    ) const;

    /**
     * How many free parameters are there?  This is the sum of the free
     * parameters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const;

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ();

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
    void
    even ();

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const;

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
    normalize ( ProbabilityType const & min );

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random );

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
    template <typename AnyMatchEmissionParametersType>
    void
    dirichlet (
      AnyMatchEmissionParametersType const & counts,
      Random & random
    );

    /**
     * Set all values such that each distrubution is distributed according to
     * the given DirichletMixture distribution (which should be a
     * DynamicProgramming::DirichletMixtureMatchEmissionPrior class).
     */
    template <typename DirichletMixtureType>
    void
    dirichletMixture (
      DirichletMixtureType const & me_prior,
      Random & random
    );

    /**
     * Return the largest value in the contained distributions.
     */
    ProbabilityType
    maximumValue () const;

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream & is,
      MatchEmissionParameters & pos
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      pos.readMatchEmissionParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // operator>> ( istream &, MatchEmissionParameters & )

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readMatchEmissionParameters (
      std::istream& is
    );

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readParameterCollection (
      std::istream & is
    )
    {
      readMatchEmissionParameters( is );
    } // readParameterCollection( istream & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      MatchEmissionParameters const& pos
    )
    {
      os << "[ ";
      pos.writeMatchEmissionParameters( os );
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, MatchEmissionParameters const& )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
    template<class CharT, class Traits>
    void writeMatchEmissionParameters (
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
      writeMatchEmissionParameters( os );
    } // writeParameterCollection( basic_ostream & ) const

  }; // End class MatchEmissionParameters

  /**
   * These are the emission parameters for the Insertion distribution.
   */
template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType = seqan::True>
  class InsertionEmissionParameters
  {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_NVP( m_Insertion_Emission_Distribution );
    } // serialize( Archive &, const unsigned int )

  protected:

// Accessors and declarations for emissions
    GALOSH_OPERATOR_DECLARE_EMISSION( Insertion, ProbabilityType )

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
// Accessors for derived classes
#define GALOSH_OPERATORS_InsertionEmissionParameters( probability_type ) \
    GALOSH_OPERATOR_EMISSION( InsertionEmissionParameters, Insertion, probability_type ) \

// Accessors for delegating classes
#define GALOSH_OPERATORS_DELEGATING_InsertionEmissionParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, Insertion, probability_type ) \

#else // !USE_FLANKING_EMISSION_DISTRIBUTIONS
    GALOSH_OPERATOR_DECLARE_EMISSION_OPERATORS( PreAlign, Insertion, ProbabilityType )
    GALOSH_OPERATOR_DECLARE_EMISSION_OPERATORS( PostAlign, Insertion, ProbabilityType )

// Accessors for derived classes
#define GALOSH_OPERATORS_InsertionEmissionParameters( probability_type ) \
    GALOSH_OPERATOR_EMISSION( InsertionEmissionParameters, Insertion, probability_type ) \
    GALOSH_OPERATOR_EMISSION( InsertionEmissionParameters, PreAlign, probability_type ) \
    GALOSH_OPERATOR_EMISSION( InsertionEmissionParameters, PostAlign, probability_type ) \


// Accessors for delegating classes
#define GALOSH_OPERATORS_DELEGATING_InsertionEmissionParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, Insertion, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PreAlign, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PostAlign, probability_type ) \

#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS .. else ..

  public:

    InsertionEmissionParameters ();

    InsertionEmissionParameters (
      InsertionEmissionParameters const & copy_from
    );

    void
    reinitialize ();

    template <typename AnyProbabilityType>
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
    operator= ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos );

    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
    operator= ( InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos );

    template <typename AnyProbabilityType>
    void
    copyFrom ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos );

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    InsertionEmissionParameters &
    operator/= ( AnyProbabilityType const& denominator );

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
    operator*= ( AnyProbabilityType const& scalar );

    /**
     * Add to each contained distribution the values in the given other
     * InsertionEmissionParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
    operator+= ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos );

    /**
     * Subtract from each contained distribution the values in the given other
     * InsertionEmissionParameters.  Note that this may violate the rule that the
     * probabilities are greater than 0.
     */
    template <typename AnyProbabilityType>
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
    operator-= ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos );

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).  Note that
     * the cross entropy is non-symmetric (calling other_pos.crossEntropy(
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
    template <typename AnyProbabilityType>
    double
    crossEntropy (
      InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos
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
    template <typename AnyProbabilityType,
              typename AnyInsertionEmissionParametersType>
    double
    crossEntropy (
      InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos,
      AnyInsertionEmissionParametersType const * const weights
    ) const;

    /**
     * Calculate and return the Euclidean distance between this set of emission
     * parameters and another set (treating every probability as an orthogonal
     * dimension).
     */
    double
    euclideanDistance (
      InsertionEmissionParameters const& other_pos
    ) const;

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of emission parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    double
    euclideanDistanceSquared (
      InsertionEmissionParameters const& other_pos
    ) const;

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const;

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ();

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
    void
    even ();

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const;

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
    void normalize ( ProbabilityType const & min );

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random );

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
    template <typename AnyInsertionEmissionParametersType>
    void
    dirichlet (
      AnyInsertionEmissionParametersType const & counts,
      Random & random
    );

    /**
     * Return the largest value in the contained distributions.
     */
    ProbabilityType
    maximumValue () const;

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream & is,
      InsertionEmissionParameters & pos
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      pos.readInsertionEmissionParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // operator>> ( istream &, InsertionEmissionParameters & )

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readInsertionEmissionParameters (
      std::istream& is
    );

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readParameterCollection (
      std::istream & is
    )
    {
      readInsertionEmissionParameters( is );
    } // readParameterCollection( istream & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      InsertionEmissionParameters const& pos
    )
    {
      os << "[ ";
      pos.writeInsertionEmissionParameters( os );
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, InsertionEmissionParameters const& )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
    template<class CharT, class Traits>
    void
    writeInsertionEmissionParameters (
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
      writeInsertionEmissionParameters( os );
    } // writeParameterCollection( basic_ostream & ) const

  }; // End class InsertionEmissionParameters

  /**
   * Each profile position stores certain position-specific transition
   * parameters.  These are they.
   */
  template <typename ResidueType, typename ProbabilityType>
  class PositionTransitionParameters
  {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_NVP( m_Match_Distribution )
         & BOOST_SERIALIZATION_NVP( m_Insertion_Distribution )
         & BOOST_SERIALIZATION_NVP( m_Deletion_Distribution );
    } // serialize( Archive &, const unsigned int )

// Accessors and declarations for transitions
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, Match )
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, Insertion )
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, Deletion )

// Accessors for derived classes
#define GALOSH_OPERATORS_PositionTransitionParameters( probability_type ) \
  GALOSH_OPERATOR_TRANSITION( PositionTransitionParameters, probability_type, Match ) \
  GALOSH_OPERATOR_TRANSITION( PositionTransitionParameters, probability_type, Deletion ) \
  GALOSH_OPERATOR_TRANSITION( PositionTransitionParameters, probability_type, Insertion ) \

// Accessors for delegating classes
#define GALOSH_OPERATORS_DELEGATING_PositionTransitionParameters( obj, probability_type ) \
  GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, Match ) \
  GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, Deletion ) \
  GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, Insertion ) \

  public:

    PositionTransitionParameters ();

    PositionTransitionParameters (
      PositionTransitionParameters const & copy_from
    );

    void
    reinitialize ();

    template <typename AnyProbabilityType>
    PositionTransitionParameters<ResidueType, ProbabilityType> &
    operator= ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    PositionTransitionParameters<ResidueType, ProbabilityType> &
    operator= ( PositionTransitionParameters<ResidueType, ProbabilityType> const& other_pos );

    template <typename AnyProbabilityType>
    void
    copyFrom ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PositionTransitionParameters &
    operator/= ( AnyProbabilityType const& denominator );

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PositionTransitionParameters<ResidueType, ProbabilityType> &
    operator*= ( AnyProbabilityType const& scalar );

    /**
     * Add to each contained distribution the values in the given other
     * PositionTransitionParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PositionTransitionParameters<ResidueType, ProbabilityType> &
    operator+= ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Subtract from each contained distribution the values in the given other
     * PositionTransitionParameters.  Note that this may violate the
     * rule that the probabilities are greater than 0.
     */
    template <typename AnyProbabilityType>
    PositionTransitionParameters<ResidueType, ProbabilityType> &
    operator-= ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).  Note that
     * the cross entropy is non-symmetric (calling other_pos.crossEntropy(
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
    template <typename AnyProbabilityType>
    double
    crossEntropy (
      PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos
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
    template <typename AnyProbabilityType,
              typename AnyPositionTransitionParametersType>
    double
    crossEntropy (
      PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPositionTransitionParametersType const * const weights
    ) const;

    /**
     * Calculate and return the Euclidean distance between this set of
     * transition parameters and another set (treating every probability as an
     * orthogonal dimension).
     */
    double
    euclideanDistance (
      PositionTransitionParameters const& other_pos
    ) const;

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of transition parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    // Templated this because ProfileTreeInternalNode is not a
    // PositionTransitionParameters but supports its interface.
    template <typename TransitionParametersType>
    double
    euclideanDistanceSquared (
      TransitionParametersType const& other_pos
    ) const;

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const;

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ();

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
    void
    even ();

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const;

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
    template <typename other_type>
    void normalize ( other_type const & min );

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
    void
    normalize ( ProbabilityType const & min );

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random );

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
    template <typename AnyPositionTransitionParametersType>
    void
    dirichlet (
      AnyPositionTransitionParametersType const & counts,
      Random & random
    );

    /**
     * Return the largest value in the contained distributions.
     */
    ProbabilityType
    maximumValue () const;

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream & is,
      PositionTransitionParameters & pos
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      pos.readPositionTransitionParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // operator>> ( istream &, PositionTransitionParameters & )

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readPositionTransitionParameters (
      std::istream& is
    );

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readParameterCollection (
      std::istream & is
    )
    {
      readPositionTransitionParameters( is );
    } // readParameterCollection( istream & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      PositionTransitionParameters const& pos
    )
    {
      os << "[ ";
      pos.writePositionTransitionParameters( os );
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, PositionTransitionParameters const& )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
    template<class CharT, class Traits>
    void
    writePositionTransitionParameters (
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
      writePositionTransitionParameters( os );
    } // writeParameterCollection( basic_ostream & ) const

  }; // End class PositionTransitionParameters

  /**
   * Each profile position stores certain position-specific parameters.  These
   * are they.  Right now this is just a trivial subclass of
   * MatchEmissionParameters.
   */
  template <typename ResidueType, typename ProbabilityType>
  class PositionSpecificParameters :
    public MatchEmissionParameters<ResidueType, ProbabilityType>
  {
    // Boost serialization
  private:
    typedef MatchEmissionParameters<ResidueType, ProbabilityType> MatchEmissionParameters_ProbabilityType;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      // save/load base class information
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( MatchEmissionParameters_ProbabilityType );
    } // serialize( Archive &, const unsigned int )

  public:
    PositionSpecificParameters ();

    template <typename AnyProbabilityType>
    PositionSpecificParameters (
      PositionSpecificParameters<ResidueType, AnyProbabilityType> const & copy_from
    );

    PositionSpecificParameters (
      PositionSpecificParameters const & copy_from
    );

    // operators
#define GALOSH_OPERATORS_PositionSpecificParameters( probability_type ) \
    GALOSH_OPERATORS_MatchEmissionParameters( probability_type )

    // operators
    GALOSH_OPERATORS_PositionSpecificParameters( ProbabilityType )

    template <typename AnyProbabilityType>
    PositionSpecificParameters &
    operator= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    PositionSpecificParameters &
    operator= ( PositionSpecificParameters const & other_pos );

    template <typename AnyProbabilityType>
    void
    copyFrom ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Add to each contained distribution the values in the given other
     * ProfilePosition.  Note that this violates the rule that the
     * probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PositionSpecificParameters &
    operator+= ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Add to each contained distribution the values in the given other
     * ProfilePosition.  Note that this violates the rule that the
     * probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PositionSpecificParameters &
    operator+= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      MatchEmissionParameters<ResidueType, ProbabilityType>::operator+=( other_pos );

      return *this;
    } // operator+=( MatchEmissionParameters<ResidueType, AnyProbabilityType> const & )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).  Note that
     * the cross entropy is non-symmetric (calling other_pos.crossEntropy(
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
    template <typename AnyProbabilityType>
    double
    crossEntropy (
      PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos
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
    template <typename AnyProbabilityType,
              typename AnyPositionSpecificParametersType>
    double
    crossEntropy (
      PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPositionSpecificParametersType const * const weights
    ) const;

    /**
     * Calculate and return the Euclidean distance between this profile
     * position and another profile position (treating every probability as an
     * orthogonal dimension).
     */
    double
    euclideanDistance (
      PositionSpecificParameters const& other_pos
    ) const;

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile position and another profile position (treating every
     * probability as an orthogonal dimension).
     */
    double
    euclideanDistanceSquared (
      PositionSpecificParameters const& other_pos
    ) const;

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const;

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream & is,
      PositionSpecificParameters & pos
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      pos.readPositionSpecificParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // operator>> ( istream &, PositionSpecificParameters & )

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readPositionSpecificParameters (
      std::istream& is
    );

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readParameterCollection (
      std::istream & is
    )
    {
      readPositionSpecificParameters( is );
    } // readParameterCollection( istream & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      PositionSpecificParameters const& pos
    )
    {
      os << "[ ";
      pos.writePositionSpecificParameters( os );
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, PositionSpecificParameters const& )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
    template<class CharT, class Traits>
    void
    writePositionSpecificParameters (
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
      writePositionSpecificParameters( os );
    } // writeParameterCollection( basic_ostream & ) const

  }; // End class PositionSpecificParameters

  /**
   * These are the pre-align transition parameters and the pre-align
   * insertion distribution parameters, and maybe the del-in transition too.
   */
  template <typename ResidueType, typename ProbabilityType>
  class PreAlignParameters
  {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_NVP( m_PreAlign_Distribution )
        & BOOST_SERIALIZATION_NVP( m_Begin_Distribution );
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      ar & BOOST_SERIALIZATION_NVP( m_PreAlign_Emission_Distribution );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
#ifdef USE_DEL_IN_DEL_OUT
      ar & BOOST_SERIALIZATION_NVP( m_DeletionIn_Distribution );
#endif // USE_DEL_IN_DEL_OUT
    } // serialize( Archive &, const unsigned int )

  public:
// Accessors and declarations for transitions
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, PreAlign )
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, Begin )
#ifdef USE_DEL_IN_DEL_OUT
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, DeletionIn )
#endif // USE_DEL_IN_DEL_OUT

// Accessors and declarations for emissions
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
    GALOSH_OPERATOR_DECLARE_EMISSION( PreAlign, ProbabilityType )

// Accessors for derived classes
#ifdef USE_DEL_IN_DEL_OUT
#define GALOSH_OPERATORS_PreAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, PreAlign )  \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, Begin ) \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, DeletionIn ) \
    GALOSH_OPERATOR_EMISSION( PreAlignParameters, PreAlign, probability_type )
#else
#define GALOSH_OPERATORS_PreAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, PreAlign )  \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, Begin ) \
    GALOSH_OPERATOR_EMISSION( PreAlignParameters, PreAlign, probability_type )
#endif // USE_DEL_IN_DEL_OUT .. else ..

// Accessors for delegating classes
#ifdef USE_DEL_IN_DEL_OUT
#define GALOSH_OPERATORS_DELEGATING_PreAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PreAlign )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, Begin ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, DeletionIn ) \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PreAlign, probability_type )
#else
#define GALOSH_OPERATORS_DELEGATING_PreAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PreAlign )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, Begin ) \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PreAlign, probability_type )
#endif // USE_DEL_IN_DEL_OUT .. else ..

#else // !USE_FLANKING_EMISSION_DISTRIBUTIONS
    GALOSH_OPERATOR_DECLARE_EMISSION_USEINSERTION( PreAlign, ProbabilityType )

// Accessors for derived classes
#ifdef USE_DEL_IN_DEL_OUT
#define GALOSH_OPERATORS_PreAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, PreAlign )  \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, DeletionIn )  \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, Begin )
#else
#define GALOSH_OPERATORS_PreAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, PreAlign )  \
    GALOSH_OPERATOR_TRANSITION( PreAlignParameters, probability_type, Begin )
#endif // USE_DEL_IN_DEL_OUT .. else ..

// Accessors for delegating classes
#ifdef USE_DEL_IN_DEL_OUT
#define GALOSH_OPERATORS_DELEGATING_PreAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PreAlign )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, DeletionIn )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, Begin )
#else
#define GALOSH_OPERATORS_DELEGATING_PreAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PreAlign )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, Begin )
#endif // USE_DEL_IN_DEL_OUT .. else ..

#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS .. else ..

  public:

    PreAlignParameters ();

    PreAlignParameters (
      PreAlignParameters const & copy_from
    );

    void
    reinitialize ();

    template <typename AnyProbabilityType>
    PreAlignParameters<ResidueType, ProbabilityType> &
    operator= ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_node );

    PreAlignParameters<ResidueType, ProbabilityType> &
    operator= ( PreAlignParameters<ResidueType, ProbabilityType> const& other_node );

    template <typename AnyProbabilityType>
    void
    copyFrom ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PreAlignParameters &
    operator/= ( AnyProbabilityType const& denominator );

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PreAlignParameters<ResidueType, ProbabilityType> &
    operator*= ( AnyProbabilityType const& scalar );

    /**
     * Add to each contained distribution the values in the given other
     * PreAlignParameters.  Note that this may violate the rule
     * that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PreAlignParameters<ResidueType, ProbabilityType> &
    operator+= ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_node );

    /**
     * Subtract from each contained distribution the values in the given other
     * PreAlignParameters.  Note that this may violate the
     * rule that the probabilities are greater than 0.
     */
    template <typename AnyProbabilityType>
    PreAlignParameters<ResidueType, ProbabilityType> &
    operator-= ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_node );

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).  Note that
     * the cross entropy is non-symmetric (calling other_pos.crossEntropy(
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
    template <typename AnyProbabilityType>
    double
    crossEntropy (
      PreAlignParameters<ResidueType, AnyProbabilityType> const& other_pos
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
    template <typename AnyProbabilityType,
              typename AnyPreAlignParametersType>
    double
    crossEntropy (
      PreAlignParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPreAlignParametersType const * const weights
    ) const;

    /**
     * Calculate and return the Euclidean distance between this set of
     * transition parameters and another set (treating every probability as an
     * orthogonal dimension).
     */
    double
    euclideanDistance (
      PreAlignParameters const& other_pos
    ) const;

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of transition parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    // Templated this because ProfileTreeInternalNode is not a
    // PreAlignParameters but supports its interface.
    template <typename PreAlignParametersType>
    double
    euclideanDistanceSquared (
      PreAlignParametersType const& other_node
    ) const;

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const;

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ();

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
    void
    even ();

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const;

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
    normalize ( ProbabilityType const & min );

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random );

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
    template <typename AnyPreAlignParametersType>
    void
    dirichlet (
      AnyPreAlignParametersType const & counts,
      Random & random
    );

    /**
     * Return the largest value in the contained distributions.
     */
    ProbabilityType
    maximumValue () const;

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream & is,
      PreAlignParameters & pos
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      pos.readPreAlignParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // operator>> ( istream &, PreAlignParameters & )

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readPreAlignParameters (
      std::istream& is
    );

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readParameterCollection (
      std::istream & is
    )
    {
      readPreAlignParameters( is );
    } // readParameterCollection( istream & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      PreAlignParameters const& pos
    )
    {
      os << "[ ";
      pos.writePreAlignParameters( os );
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, PreAlignParameters const& )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
    template<class CharT, class Traits>
    void
    writePreAlignParameters (
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
      writePreAlignParameters( os );
    } // writeParameterCollection( basic_ostream & ) const

  }; // End class PreAlignParameters

  /**
   * These are the post-align transition parameters and the post-align
   * insertion distribution parameters, and maybe the del-out parameters.
   */
  template <typename ResidueType, typename ProbabilityType>
  class PostAlignParameters
  {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
#ifdef USE_END_DISTRIBUTION
      ar & BOOST_SERIALIZATION_NVP( m_End_Distribution );
#endif // USE_END_DISTRIBUTION
      ar & BOOST_SERIALIZATION_NVP( m_PostAlign_Distribution );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      ar & BOOST_SERIALIZATION_NVP( m_DeletionOut_Distribution );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      ar & BOOST_SERIALIZATION_NVP( m_PostAlign_Emission_Distribution );
#endif //USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // serialize( Archive &, const unsigned int )

  public:

// Accessors and declarations for transitions
#ifdef USE_END_DISTRIBUTION
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, End )
#endif // USE_END_DISTRIBUTION
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, PostAlign )
#ifdef USE_DEL_IN_DEL_OUT
    GALOSH_OPERATOR_DECLARE_TRANSITION( ProbabilityType, DeletionOut )
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
// Accessors and declarations for emissions
    GALOSH_OPERATOR_DECLARE_EMISSION( PostAlign, ProbabilityType )

// Accessors for derived classes

#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, End ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign )  \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, DeletionOut )  \
    GALOSH_OPERATOR_EMISSION( PostAlignParameters, PostAlign, probability_type )
#else
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign )  \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, DeletionOut )  \
    GALOSH_OPERATOR_EMISSION( PostAlignParameters, PostAlign, probability_type )
#endif // USE_END_DISTRIBUTION
#else
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, End ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign )  \
    GALOSH_OPERATOR_EMISSION( PostAlignParameters, PostAlign, probability_type )
#else
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign )  \
    GALOSH_OPERATOR_EMISSION( PostAlignParameters, PostAlign, probability_type )
#endif // USE_END_DISTRIBUTION
#endif // USE_DEL_IN_DEL_OUT .. else ..

// Accessors for delegating classes
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, End ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, DeletionOut )  \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PostAlign, probability_type )
#else
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, DeletionOut )  \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PostAlign, probability_type )
#endif // USE_END_DISTRIBUTION
#else
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, End ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign )  \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PostAlign, probability_type )
#else
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign )  \
    GALOSH_OPERATOR_DELEGATING_EMISSION( obj, PostAlign, probability_type )
#endif // USE_END_DISTRIBUTION
#endif // USE_DEL_IN_DEL_OUT .. else ..

#else // !USE_FLANKING_EMISSION_DISTRIBUTIONS
// Accessors and declarations for emissions
    GALOSH_OPERATOR_DECLARE_EMISSION_USEINSERTION( PostAlign, ProbabilityType )

// Accessors for derived classes
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, End ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign )  \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, DeletionOut )
#else
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, DeletionOut )
#endif // USE_END_DISTRIBUTION
#else
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, End ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign )
#else
#define GALOSH_OPERATORS_PostAlignParameters( probability_type ) \
    GALOSH_OPERATOR_TRANSITION( PostAlignParameters, probability_type, PostAlign )
#endif // USE_END_DISTRIBUTION
#endif // USE_DEL_IN_DEL_OUT .. else ..

// Accessors for delegating classes
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, End ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign )  \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, DeletionOut )
#else
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, DeletionOut )
#endif // USE_END_DISTRIBUTION
#else
#ifdef USE_END_DISTRIBUTION
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, End ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign )
#else
#define GALOSH_OPERATORS_DELEGATING_PostAlignParameters( obj, probability_type ) \
    GALOSH_OPERATOR_DELEGATING_TRANSITION( obj, probability_type, PostAlign )
#endif // USE_END_DISTRIBUTION
#endif // USE_DEL_IN_DEL_OUT .. else ..

#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

  public:

    PostAlignParameters ();

    PostAlignParameters (
      PostAlignParameters const & copy_from
    );

    void
    reinitialize ();

    template <typename AnyProbabilityType>
    PostAlignParameters<ResidueType, ProbabilityType> &
    operator= ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_node );

    PostAlignParameters<ResidueType, ProbabilityType> &
    operator= ( PostAlignParameters<ResidueType, ProbabilityType> const& other_node );

    template <typename AnyProbabilityType>
    void
    copyFrom ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PostAlignParameters &
    operator/= ( AnyProbabilityType const& denominator );

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PostAlignParameters<ResidueType, ProbabilityType> &
    operator*= ( AnyProbabilityType const& scalar );

    /**
     * Add to each contained distribution the values in the given other
     * PostAlignParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    PostAlignParameters<ResidueType, ProbabilityType> &
    operator+= ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_node );

    /**
     * Subtract from each contained distribution the values in the given other
     * PostAlignParameters.  Note that this may violate the
     * rule that the probabilities are greater than 0.
     */
    template <typename AnyProbabilityType>
    PostAlignParameters<ResidueType, ProbabilityType> &
    operator-= ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_node );

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).  Note that
     * the cross entropy is non-symmetric (calling other_pos.crossEntropy(
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
    template <typename AnyProbabilityType>
    double
    crossEntropy (
      PostAlignParameters<ResidueType, AnyProbabilityType> const& other_pos
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
    template <typename AnyProbabilityType,
              typename AnyPostAlignParametersType>
    double
    crossEntropy (
      PostAlignParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPostAlignParametersType const * const weights
    ) const;

    /**
     * Calculate and return the Euclidean distance between this set of
     * transition parameters and another set (treating every probability as an
     * orthogonal dimension).
     */
    double
    euclideanDistance (
      PostAlignParameters const& other_pos
    ) const;

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of transition parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    // Templated this because ProfileTreeInternalNode is not a
    // PostAlignParameters but supports its interface.
    template <typename PostAlignParametersType>
    double
    euclideanDistanceSquared (
      PostAlignParametersType const& other_node
    ) const;

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const;

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ();

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
    void
    even ();

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const;

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
    normalize ( ProbabilityType const & min );

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random );


    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
    template <typename AnyPostAlignParametersType>
    void
    dirichlet (
      AnyPostAlignParametersType const & counts,
      Random & random
    );

    /**
     * Return the largest value in the contained distributions.
     */
    ProbabilityType
    maximumValue () const;

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream & is,
      PostAlignParameters & pos
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      pos.readPostAlignParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // operator>> ( istream &, PostAlignParameters & )

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readPostAlignParameters (
      std::istream& is
    );

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readParameterCollection (
      std::istream & is
    )
    {
      readPostAlignParameters( is );
    } // readParameterCollection( istream & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      PostAlignParameters const& pos
    )
    {
      os << "[ ";
      pos.writePostAlignParameters( os );
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, PostAlignParameters const& )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
    template<class CharT, class Traits>
    void writePostAlignParameters (
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
      writePostAlignParameters( os );
    } // writeParameterCollection( basic_ostream & ) const

  }; // End class PostAlignParameters

  /**
   * A Profile is a data structure for the model parameters.  Conceptually, for
   * every position of the profile there are a bunch of parameters for each
   * kind of Plan 7 profile transition/emission.  We make this a vector of maps
   * from ProfileKeys (eg. Transition) to maps from ProfileKey-specific params
   * (eg. Transition_M_to_D) to the values of those params.
   */
  template <typename ResidueType, typename ProbabilityType>
  class ProfilePosition :
    public PositionSpecificParameters<ResidueType, ProbabilityType>
  {
    // Boost serialization
  private:
    typedef PositionSpecificParameters<ResidueType, ProbabilityType> PositionSpecificParameters_ProbabilityType;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      // save/load base class information
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PositionSpecificParameters_ProbabilityType );
      ar & BOOST_SERIALIZATION_NVP( m_root );
    } // serialize( Archive &, const unsigned int )

  public://protected:
    ProfileTreeRoot<ResidueType, ProbabilityType> * m_root;

  public:
    ProfilePosition ();

    ProfilePosition ( ProfileTreeRoot<ResidueType, ProbabilityType> * const pt_root );

    ProfilePosition ( ProfilePosition<ResidueType, ProbabilityType> const & copy_from );

    void
    reinitialize () { reinitialize( 0 ); }

    void
    reinitialize ( ProfileTreeRoot<ResidueType, ProbabilityType> * const pt_root );

    void
    setProfileTreeRoot ( ProfileTreeRoot<ResidueType, ProbabilityType> * const pt_root );

    ProfileTreeRoot<ResidueType, ProbabilityType> *
    getProfileTreeRoot () const;

    // operators
    GALOSH_OPERATORS_MatchEmissionParameters( ProbabilityType )
    GALOSH_OPERATORS_DELEGATING_PositionTransitionParameters( ( *m_root ), ProbabilityType )
    GALOSH_OPERATORS_DELEGATING_InsertionEmissionParameters( ( *m_root ), ProbabilityType )

    template <typename AnyProbabilityType>
    void
    copyFrom ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    template <typename AnyProbabilityType>
    void
    copyFrom ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos );

    template <typename AnyProbabilityType>
    void
    copyFrom ( ProfilePosition<ResidueType, AnyProbabilityType> const& other_pos );

    template <typename AnyProbabilityType>
    ProfilePosition<ResidueType, ProbabilityType> &
    operator= ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos );

    template <typename AnyProbabilityType>
    ProfilePosition<ResidueType, ProbabilityType> &
    operator= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos );

    template <typename AnyProbabilityType>
    ProfilePosition<ResidueType, ProbabilityType> &
    operator= ( ProfilePosition<ResidueType, AnyProbabilityType> const& other_pos );

    ProfilePosition<ResidueType, ProbabilityType> &
    operator= ( ProfilePosition<ResidueType, ProbabilityType> const& other_pos );

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    ProfilePosition &
    operator/= ( AnyProbabilityType const& denominator );

    /**
     * Add to each contained distribution the values in the given other
     * ProfilePosition.  Note that this violates the rule that the
     * probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    ProfilePosition &
    operator+= ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos );

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ();

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
    void
    even ();

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const;

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
    normalize ( ProbabilityType const & min );

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random );

    /**
     * Return the largest value in the contained distributions.
     */
    ProbabilityType
    maximumValue () const;

  }; // End class ProfilePosition

  template <typename ResidueType, typename ProbabilityType>
  class GlobalParameters :
    //public ProfileTreeRoot<ResidueType, ProbabilityType>
    public PositionTransitionParameters<ResidueType, ProbabilityType>,
    public InsertionEmissionParameters<ResidueType, ProbabilityType>,
    public PreAlignParameters<ResidueType, ProbabilityType>,
    public PostAlignParameters<ResidueType, ProbabilityType>
  {
    // Boost serialization
  private:
    typedef PositionTransitionParameters<ResidueType, ProbabilityType> PositionTransitionParameters_ProbabilityType;
    typedef InsertionEmissionParameters<ResidueType, ProbabilityType> InsertionEmissionParameters_ProbabilityType;
    typedef PreAlignParameters<ResidueType, ProbabilityType> PreAlignParameters_ProbabilityType;
    typedef PostAlignParameters<ResidueType, ProbabilityType> PostAlignParameters_ProbabilityType;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PositionTransitionParameters_ProbabilityType );
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( InsertionEmissionParameters_ProbabilityType );
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PreAlignParameters_ProbabilityType );
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( PostAlignParameters_ProbabilityType );
    } // serialize( Archive &, const unsigned int )

  public:
    GlobalParameters ();
  
    GlobalParameters (
      GlobalParameters const & copy_from
    );

    /**
     * Reinitialize contained distributions.
     */
    void
    reinitialize ();
  
    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const;

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
    normalize ( ProbabilityType const& min );
  
    /**
     * Set all values to 0.
     */
    void
    zero ();

    /**
     * Set all values such that each distrubution is uniformly distributed.
     */
    void
    even ();

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random );

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
    template <typename AnyGlobalParametersType>
    void
    dirichlet (
      AnyGlobalParametersType const & counts,
      Random & random
    );

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream & is,
      GlobalParameters & pos
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      pos.readGlobalParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // operator>> ( istream &, GlobalParameters & )

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readGlobalParameters (
      std::istream& is
    );

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
    void
    readParameterCollection (
      std::istream & is
    )
    {
      readGlobalParameters( is );
    } // readParameterCollection( istream & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      GlobalParameters const& node
    )
    {
      os << "[ ";
      node.writeGlobalParameters( os );
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, GlobalParameters const& )
  
    /**
     * Write a comma-separated list of the parameters to the stream.
     */
    template<class CharT, class Traits>
    void
    writeGlobalParameters (
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
      writeGlobalParameters( os );
    } // writeParameterCollection( basic_ostream & ) const

    /**
     * Return the largest value stored in this ProfileTreeRoot.
     */
    ProbabilityType
    maximumValue () const;

    // operators
#define GALOSH_OPERATORS_GlobalParameters( probability_type ) \
    GALOSH_OPERATORS_PositionTransitionParameters( probability_type ) \
    GALOSH_OPERATORS_InsertionEmissionParameters( probability_type ) \
    GALOSH_OPERATORS_PreAlignParameters( probability_type ) \
    GALOSH_OPERATORS_PostAlignParameters( probability_type )

    GALOSH_OPERATORS_GlobalParameters( ProbabilityType )

    template <typename AnyProbabilityType>
    void
    copyFrom ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node );

    template <typename AnyProbabilityType>
    GlobalParameters<ResidueType, ProbabilityType> &
    operator= ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node );

    GlobalParameters<ResidueType, ProbabilityType> &
    operator= ( GlobalParameters<ResidueType, ProbabilityType> const& other_node );

    /**
     * Multiply each profile value by the given scalar.  Note that this may
     * violate the rule that the probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    GlobalParameters<ResidueType, ProbabilityType> &
    operator*= ( AnyProbabilityType const& scalar );

    /**
     * Divide each non-position-specific profile value by the given
     * denominator.  Note that this may violate the rule that the
     * probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    GlobalParameters<ResidueType, ProbabilityType> &
    operator/= ( AnyProbabilityType const& denominator );

    /**
     * Add to each contained distribution the corresponding values in the given
     * other node.  Note that this may violate the rule that the probabilities
     * sum to 1.
     */
    template <typename AnyProbabilityType>
    GlobalParameters<ResidueType, ProbabilityType> &
    operator+= ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node );

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).  Note that
     * the cross entropy is non-symmetric (calling other_pos.crossEntropy(
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
    template <typename AnyProbabilityType>
    double
    crossEntropy (
      GlobalParameters<ResidueType, AnyProbabilityType> const& other_pos
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
    template <typename AnyProbabilityType,
              typename AnyGlobalParametersType>
    double
    crossEntropy (
      GlobalParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyGlobalParametersType const * const weights
    ) const;

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension).
     */
    template <typename AnyProbabilityType>
    double
    euclideanDistanceSquared (
      GlobalParameters<ResidueType, AnyProbabilityType> const& other_node
    ) const;

    /**
     * How many free parameters are there?  This is the sum of the free
     * parameters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const;

  }; // End class GlobalParameters

  /**
   * A Profile is a data structure for the model parameters.  Conceptually, for
   * every position of the profile there are a bunch of parameters for each
   * kind of Plan 7 profile transition/emission.  We make this a vector of maps
   * from ProfileKeys (eg. Transition) to maps from ProfileKey-specific params
   * (eg. Transition_M_to_D) to the values of those params.
   *
   * This is the basic implementation of the Profile.  It is called TreeRoot
   * because it supports variations as "child profiles" (ProfileTreeInternalNode
   * objects).
   */
  template <typename ResidueType, typename ProbabilityType>
  class ProfileTreeRoot :
    public GlobalParameters<ResidueType, ProbabilityType>,
    public vector<ProfilePosition<ResidueType, ProbabilityType> >
  {
    // Boost serialization
  private:
    typedef GlobalParameters<ResidueType, ProbabilityType> GlobalParameters_ProbabilityType;
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( GlobalParameters_ProbabilityType );

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
         ar & boost::serialization::make_nvp(tmp_str.c_str(), ( *this )[ pos_i ]);
       }
      } // End if length > 0, serialize positions
    } // serialize( Archive &, const unsigned int )

  protected:
    typedef ProfilePosition<ResidueType, ProbabilityType> Position;

  protected:
    GALOSH_OPERATORS_GlobalParameters( ProbabilityType )

  public:
    uint32_t m_profileTreeVertex;

    typedef ResidueType ProfileResidueType;

    friend class ProfilePosition<ResidueType, ProbabilityType>;

    ProfileTreeRoot () :
      GlobalParameters<ResidueType, ProbabilityType>()
      //,m_Loop_Distribution()
    {
      //cout << "ProfileTreeRoot () initializer called: this is " << this << endl;
      //this->reinitialize( 0, 0 );
    } // <init>()

//    ProfileTreeRoot ( uint32_t length ) :
//      GlobalParameters<ResidueType, ProbabilityType>(),
//      vector<Position>()
//    {
//      // TODO: REMOVE
//      //cout << "ProfileTreeRoot ( uint32_t ) initializer called: this is " << this << endl;
//      this->reinitialize( length, 0 );
//    } // <init>( uint32_t )

    ProfileTreeRoot ( 
      uint32_t length
    ) :
      GlobalParameters<ResidueType, ProbabilityType>(),
      vector<Position>()
    {
      // TODO: REMOVE
      //cout << "ProfileTreeRoot ( uint32_t ) initializer called: this is " << this << endl;

      this->reinitialize( length );
    } // <init>( uint32_t )

    ProfileTreeRoot ( 
      ProfileTreeRoot<ResidueType, ProbabilityType> const & copy_from
    ) :
      GlobalParameters<ResidueType, ProbabilityType>(),
      vector<Position>()
    {
      // TODO: REMOVE
      //cout << "ProfileTreeRoot ( ProfileTreeRoot const& ) initializer called: this is " << this << endl;

      this->copyFrom( copy_from );
    } // <init>( ProfileTreeRoot const& )

    template <class AnyProbabilityType>
    ProfileTreeRoot ( 
      ProfileTreeRoot<ResidueType, AnyProbabilityType> const & copy_from
    ) :
      GlobalParameters<ResidueType, ProbabilityType>(),
      vector<Position>()
    {
      // TODO: REMOVE
      //cout << "ProfileTreeRoot ( ProfileTreeRoot<ResidueType, AnyProbabilityType> const& ) initializer called: this is " << this << endl;

      this->copyFrom( copy_from );
    } // <init>( ProfileTreeRoot<ResidueType, AnyProbabilityType> const& )

    // TODO: REMOVE.  HACK to get it to have the same initializer as the
    // internalnode...
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeRoot ( 
      AnyInternalNodeOrRoot const * copy_length_from,
      vector<int> const & unused_variations
    )
    {
      reinitialize( *copy_length_from );
    } // <init>( AnyInternalNodeOrRoot const *, vector<int> const & )

    // TODO: REMOVE.  HACK to get it to have the same initializer as the
    // internalnode...
    /*
    //template <typename AnyProbabilityType>
    ProfileTreeRoot ( 
      ProfileTreeRoot<ResidueType, ProbabilityType> const * copy_from,
      vector<int> const & unused_variations
    ) :
      GlobalParameters<ResidueType, ProbabilityType>(),
      vector<Position>()
    {
      // TODO: REMOVE
      //cout << "ProfileTreeRoot ( AnyRootOrRoot const *, vector<int> const & ) initializer called: this is " << this << endl;

      this->copyFrom( &copy_from );
    } // <init>( ProfileTreeRoot<ResidueType, AnyProbabilityType> const *, vector<int> const & )
    // TODO: REMOVE.  HACK to get it to have the same initializer as the
    // internalnode...
    //template <typename AnyProbabilityType>
    ProfileTreeRoot ( 
      ProfileTreeInternalNode<ResidueType, ProbabilityType> * copy_from,
      vector<int> const & unused_variations
    ) :
      GlobalParameters<ResidueType, ProbabilityType>(),
      vector<Position>()
    {
      // TODO: REMOVE
      //cout << "ProfileTreeRoot ( AnyInternalNodeOrRoot const *, vector<int> const & ) initializer called: this is " << this << endl;

      this->copyFrom( &copy_from );
    } // <init>( ProfileTreeInternalNode<ResidueType, AnyProbabilityType> const *, vector<int> const & )
    */

    // TODO: REMOVE.  HACK to get it to have the same initializer as the
    // internalnode...
    template <typename AnyInternalNodeOrRoot>
    void
    reinitialize ( 
      AnyInternalNodeOrRoot const * copy_length_from,
      vector<int> const & unused_variations
    )
    {
      reinitialize( *copy_length_from );
    }

    template <typename AnyInternalNodeOrRoot>
    void
    reinitialize ( 
      AnyInternalNodeOrRoot const& copy_length_from
    )
    {
      ProfileTreeRoot<ResidueType, ProbabilityType>::reinitialize( copy_length_from.length() );

      // Also copy the profile tree vertex.
      setProfileTreeVertex( copy_length_from.getProfileTreeVertex() );
    } // reinitialize( AnyInternalNodeOrRoot const& )

    void
    reinitialize (
      uint32_t length = 0
    )
    {
      GlobalParameters<ResidueType, ProbabilityType>::reinitialize();

      if( vector<Position>::size() != length ) {
        vector<Position>::resize( length );
      }

      // Initialize and even() positions.
      if( length != 0 ) {
        uint32_t last_pos = length - 1;
        uint32_t pos_i;
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          this->vector<Position>::operator[]( pos_i ).reinitialize( this );
        }
      }
    } // reinitialize( uint32_t )

    /**
     * Set the values of the the positions of this profile to assign
     * probability 1.0 to the given sequence (and to be of the same length).
     * Does not effect global/transition parameters.
     */
    template <typename SequenceResidueType>
    void
    fromSequence (
      Sequence<SequenceResidueType> const & sequence
    )
    {
      reinitialize( sequence.length() );

      // Initialize and even() positions.
      if( length() != 0 ) {
        uint32_t last_pos = length() - 1;
        uint32_t pos_i;

        // First set everything to zero...
        zeroPositions();

        // Now set the relevant probs to 1.0.
        for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
          ( *this )[ pos_i ][ Emission::Match ][ sequence[ pos_i ] ] = 1.0;
        } // End foreach pos..
      } // End if length != 0
    } // fromSequence( Sequence<SequenceResidueType> const & )

    /**
     * Returns a pointer to the root of the tree, which is *this.
     */
    ProfileTreeRoot<ResidueType, ProbabilityType> * const
    getProfileTreeRoot ()
    {
      return this;
    } // getProfileTreeRoot()

    /**
     * Sometimes, for reasons I have yet to figure out, profile positions point
     * to the wrong root.  Until I figure out why/where, I use this as a
     * workaround.  It just goes through and sets each position's root to be
     * this node's root (which is this).
     */
    void
    ensurePositionsKnowTheirRoot ()
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = vector<Position>::size() - 1;
      uint32_t pos_i; // ( here pos_i refers to position into the internal datastore ).
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
#ifndef NDEBUG
        if( this->vector<Position>::operator[]( pos_i ).getProfileTreeRoot() != this ) {
          std::cerr << "WARNING: ensurePositionsKnowTheirRoot() has found that position " << pos_i << " has the wrong root ( " << this->vector<Position>::operator[]( pos_i ).getProfileTreeRoot() << " instead of " << this << " )." << endl;
        }
#endif // !NDEBUG
        this->vector<Position>::operator[]( pos_i ).setProfileTreeRoot( this );
      }
    } // ensurePositionsKnowTheirRoot()

    uint32_t
    getProfileTreeVertex () const
    {
      //return 0; // Root always has 0 vertex in the ProfileTree.
      return m_profileTreeVertex;
    } // getProfileTreeVertex() const

    // TODO: Should use ProfileTree.hpp's vertex_t.
    void
    setProfileTreeVertex ( uint32_t const & new_vertex )
    {
      // TODO: REMOVE
      //cout << "setProfileTreeVertex( " << new_vertex << " )" << endl;
      m_profileTreeVertex = new_vertex;
    } // setProfileTreeVertex( uint32_t const & )

#ifdef USE_DEL_IN_DEL_OUT
    /**
     * Create in the given Match distribution a version of this Profile's Match
     * distribution, scaled for the DeletionOut paths out of the Match state at
     * the given position.
     */
    template <typename ScaledMatchDistributionProbabilityType>
    void
    createScaledMatchDistributionForPosition (
      uint32_t const & pos_i,
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> & match_distribution
    ) const
    {
      // There are no transitions out of the last position.  (or rather, the
      // only transition is to the PostAlign state).
      assert( pos_i < ( this->length() - 1 ) );

      match_distribution =
        ( *this )[ Transition::fromMatch ];

      if( ( pos_i + 2 ) < this->length() ) {
        // then there's del-out extensions to consider (from this pos)
      
        ScaledMatchDistributionProbabilityType del_out_extensions_prob =
          pow( ( *this )[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ], static_cast<ScaledMatchDistributionProbabilityType>( this->length() - ( pos_i + 2 ) ) );
      
        match_distribution[ TransitionFromMatch::toDeletionOut ] *=
          del_out_extensions_prob;
      } // End if there are del-out extensions to consider for this pos
      // And the del-out has to end.
      match_distribution[ TransitionFromMatch::toDeletionOut ] *=
        ( *this )[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ];

      //// TODO: REMOVE.  TESTING (part of debugging USE_DEL_IN_DEL_OUT).  This is an equivalent alternative to the above, and I had wanted to see what the effect of the numeric drift was, so I wanted to do it this way, which is exactly how it is calculated .. partly by backward_calculateRow(..), then by updateGlobalEntenteForSequence(..)
      ScaledMatchDistributionProbabilityType foo =
        ( *this )[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toEnd ];
      for( uint32_t subsequent_del_out_ext_pos_i = ( pos_i + 2 ); ( subsequent_del_out_ext_pos_i < this->length() ); subsequent_del_out_ext_pos_i++ ) {
        foo *=
          ( *this )[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ];
      }
      match_distribution[ TransitionFromMatch::toDeletionOut ] *= foo;

        //// TODO: REMOVE
        //ScaledMatchDistributionProbabilityType other_way = pow( ( *this )[ Transition::fromDeletionOut ][ TransitionFromDeletionOut::toDeletionOut ], static_cast<ScaledMatchDistributionProbabilityType>( this->length() - ( pos_i + 2 ) ) );
        //cout << "del_out_extensions_prob: " << del_out_extensions_prob << endl;
        //cout << "other_way: " << del_out_extensions_prob << endl;
        //assert( del_out_extensions_prob == other_way );





#ifdef USE_END_DISTRIBUTION
      match_distribution[ TransitionFromMatch::toDeletionOut ] *=
         ( *this )[
           Transition::fromEnd
         ][
           TransitionFromEnd::toPostAlign
         ];
#endif // USE_END_DISTRIBUTION
      // Now rescale the other probs
      ScaledMatchDistributionProbabilityType old_scale, new_scale, scale_ratio;
      old_scale = 1.0;
      old_scale -=
        ( *this )[ Transition::fromMatch ][ TransitionFromMatch::toDeletionOut ];
      new_scale = 1.0;
      new_scale -= match_distribution[ TransitionFromMatch::toDeletionOut ];
      scale_ratio = ( new_scale / old_scale );

      match_distribution[ TransitionFromMatch::toMatch ] *= scale_ratio;
      match_distribution[ TransitionFromMatch::toInsertion ] *= scale_ratio;
      match_distribution[ TransitionFromMatch::toDeletion ] *= scale_ratio;

      return;
    } // createScaledMatchDistributionForPosition ( uint32_t const &, MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, ScaledMatchDistributionProbabilityType> & ) const
#endif // USE_DEL_IN_DEL_OUT


    template <typename AnyInternalNodeOrRoot>
    void
    copyFrom ( AnyInternalNodeOrRoot const& other_node )
    {
      // TODO: REMOVE
      //cout << "Root.copyFrom(..)" << endl;

      // Set length, initialize.
      this->reinitialize( other_node.length() );

      // Copy the values of the other profile into this one.
      copyExceptPositions( other_node );

      // And each position.
      copyPositions( other_node );

      // Also copy the profile tree vertex.
      setProfileTreeVertex( other_node.getProfileTreeVertex() );
    } // copyFrom( AnyInternalNodeOrRoot const& )

    ProfileTreeRoot &
    operator= ( ProfileTreeRoot const & other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( ProfileTreeRoot const& )

    template <typename AnyInternalNodeOrRoot>
    ProfileTreeRoot &
    operator= ( AnyInternalNodeOrRoot const& other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( AnyInternalNodeOrRoot const& )

    /**
     * Copy all values from the other profile's multinomial distributions into
     * this profile's, except for those that are position-specific.
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeRoot<ResidueType, ProbabilityType> &
    copyExceptPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      GlobalParameters<ResidueType, ProbabilityType>::copyFrom( other_node );

      return *this;
    } // copyExceptPositions( AnyInternalNodeOrRoot const& )

    /**
     * Copy the position-specific values from the other profile's multinomial
     * distributions into this profile's.  This copies as many positions as
     * possible (the min of the two lengths); you probably only want to use
     * this for profiles of the same length only).
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeRoot<ResidueType, ProbabilityType> &
    copyPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      uint32_t last_pos = min( this->length(), other_node.length() );
      if( last_pos == 0 ) {
        return *this;
      } else {
        last_pos -= 1;
      }
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ) = other_node[ pos_i ];
        this->vector<Position>::operator[]( pos_i ).setProfileTreeRoot( this );
      }
      // TODO: REMOVE
      //cout << "Root.copyPositions(..): this is " << this << endl;

      return *this;
    } // copyPositions( AnyInternalNodeOrRoot const& )

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
    void
    zero ()
    {
      zeroExceptPositions();
      zeroPositions();
    } // zero()

    /**
     * Set all values to 0, except position-specific values.  Note that this
     * violates the rule that the values sum to 1.
     */
    void
    zeroExceptPositions ()
    {
      GlobalParameters<ResidueType, ProbabilityType>::zero();
    } // zeroExceptPositions()

    /**
     * Set all position-specific values to 0.  Note that this violates the rule
     * that the values sum to 1.
     */
    void
    zeroPositions ()
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ).zero();
      }
    } // zeroPositions()

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const
    {
      ProbabilityType t =
        totalExceptPositions();
      t +=
        totalPositions();
      return t;
    } // total() const

    /**
     * Calculate the total of all contained values (except position-specific
     * values).
     */
    ProbabilityType
    totalExceptPositions () const
    {
      return
        GlobalParameters<ResidueType, ProbabilityType>::total();
    } // totalExceptPositions() const

    /**
     * Calculate the total of all contained position-specific values.
     */
    ProbabilityType
    totalPositions () const
    {
      ProbabilityType t( 0 );
      if( length() == 0 ) {
        return t;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        t +=
          this->vector<Position>::operator[]( pos_i ).total();
      }
      return t;
    } // totalPositions() const

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
    void
    even ()
    {
      evenExceptPositions();
      evenPositions();
    } // even()

    /**
     * Set all values (except position-specific values) such that each
     * distrubution is evenly distributed.
     */
    void
    evenExceptPositions ()
    {
      GlobalParameters<ResidueType, ProbabilityType>::even();
    } // evenExceptPositions()

    /**
     * Set allposition-specific values such that each distrubution is evenly
     * distributed.
     */
    void
    evenPositions ()
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ).even();
      }
    } // evenPositions()

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
    template <typename other_type>
    void
    normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
    void
    normalize ( ProbabilityType const & min )
    {
      normalizeExceptPositions( min );
      normalizePositions( min );
    } // normalize( ProbabilityType const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Do all distributions
     * except for those stored on a position-specific basis.
     */
    template <typename other_type>
    void
    normalizeExceptPositions ( other_type const & min )
    {
      //normalizeExceptPositions( min );
      normalizeExceptPositions( static_cast<ProbabilityType>( min ) );
    } // normalizeExceptPositions( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Do all distributions
     * except for those stored on a position-specific basis.
     *
     * NOTE HACK: "End" transition distribution is normalized without enforcing
     * minimum.
     */
    void
    normalizeExceptPositions ( ProbabilityType const & min )
    {
      GlobalParameters<ResidueType, ProbabilityType>::normalize( min );
    } // normalizeExceptPositions( ProbabilityType )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Do only those
     * distributions that are stored on a position-specific basis.
     */
    template <typename other_type>
    void
    normalizePositions ( other_type const & min )
    {
      //normalizePositions( min );
      normalizePositions( static_cast<ProbabilityType>( min ) );
    } // normalizePositions( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Do only those
     * distributions that are stored on a position-specific basis.
     */
    void
    normalizePositions ( ProbabilityType const & min )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ).normalize( min );
      }
    } // normalizePositions( ProbabilityType const & )

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
    void
    uniform ( Random & random )
    {
      uniformExceptPositions( random );
      uniformPositions( random );
    } // uniform( Random & )

    /**
     * Set all values (except position-specific values) such that each
     * distribution is randomly distributed.
     */
    void
    uniformExceptPositions ( Random & random )
    {
      GlobalParameters<ResidueType, ProbabilityType>::uniform( random );
    } // uniformExceptPositions( Random & )

    /**
     * Set all position-specific values such that each distribution is randomly
     * distributed.
     */
    void
    uniformPositions ( Random & random )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ).uniform( random );
      }
    } // uniformPositions( Random & )

    /**
     * Set all values such that each distrubution is distributed proportional
     * to the given counts.
     */
    template <typename CountsType>
    void
    dirichlet (
      CountsType const & counts,
      Random & random
    )
    {
      dirichletExceptPositions( counts, random );
      dirichletPositions( counts, random );
    } // dirichlet( CountsType const &, Random & )

    /**
     * Set all values (except position-specific values) such that each
     * distrubution is distributed proportional to the given counts.
     */
    template <typename CountsType>
    void
    dirichletExceptPositions (
      CountsType const & counts,
      Random & random
    )
    {
      GlobalParameters<ResidueType, ProbabilityType>::dirichlet( counts, random );
    } // dirichletExceptPositions( CountsType const &, Random & )

    /**
     * Set all position-specific values such that each distrubution is
     * distributed proportional to the given counts.
     */
    template <typename CountsType>
    void
    dirichletPositions (
      CountsType const & counts,
      Random & random
    )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ).dirichlet( counts, random );
      }
    } // dirichletPositions( CountsType const &, Random & )

    /**
     * Set all values such that each distrubution is distributed according to
     * the given DirichletMixture distributions (which should be a
     * DynamicProgramming::DirichletMixtureMatchEmissionPrior and a
     * DynamicProgramming::DirichletMixtureGlobalPrior).
     */
    template <typename DirichletMixtureMatchEmissionPriorType,
              typename DirichletMixtureGlobalPriorType>
    void
    dirichletMixture (
      //DirichletMixtureMatchEmissionPrior<DirichletParameterType> const & me_prior,
      DirichletMixtureMatchEmissionPriorType const & me_prior,
      //DirichletMixtureGlobalPrior<DirichletParameterType> const & global_prior,
      DirichletMixtureGlobalPriorType const & global_prior,
      Random & random
    )
    {
      dirichletMixtureExceptPositions( global_prior, random );
      dirichletMixturePositions( me_prior, random );
    } // dirichletMixture( DirichletMixtureMatchEmissionPrior const &, DirichletMixtureGlobalPrior const &, Random & )

    /**
     * Set all values (except position-specific values) such that each
     * distrubution is distributed according to the given
     * DirichletMixtureGlobalPrior.
     */
    template <typename DirichletMixtureGlobalPriorType>
    void
    dirichletMixtureExceptPositions (
      DirichletMixtureGlobalPriorType const & global_prior,
      Random & random
    )
    {
      if( global_prior.size() == 1 ) {
        this->GlobalParameters<ResidueType, ProbabilityType>::dirichlet( global_prior[ 0 ], random );
      } else {
        // First draw from the mixture distribution.
        double u = random.nextUniform();
        uint32_t num_components = global_prior.size();
        ProbabilityType total_so_far( 0.0 );
        for( uint32_t component_i = 0; component_i < num_components; component_i++ ) {
          total_so_far += global_prior.m_mixingProbs[ component_i ];
          if( u < total_so_far ) {
            this->GlobalParameters<ResidueType, ProbabilityType>::dirichlet( global_prior[ component_i ], random );            
            break;
          }
        }
        // This point should never be reached!
        cout << "ERROR in ProfileTreeRoot::dirichletMixtureExceptPositions(..): m_mixingProbs don't add to 1!" << endl;
      }
    } // dirichletMixtureExceptPositions( DirichletMixtureGlobalPrior const &, Random & )

    /**
     * Set all values position-specific such that each distrubution is
     * distributed according to the given DirichletMixtureMatchEmissionPrior.
     */
    template <typename DirichletMixtureMatchEmissionPriorType>
    void
    dirichletMixturePositions (
      DirichletMixtureMatchEmissionPriorType const & me_prior,
      Random & random
    )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ).dirichletMixture(
          me_prior,
          random
        );
      }
    } // dirichletMixturePositions( DirichletMixtureMatchEmissionPrior const &, Random & )

    /**
     * Stream reader.
     */
    friend std::istream &
    operator>> (
      std::istream & is,
      ProfileTreeRoot<ResidueType, ProbabilityType> & prof
    )
    {
      ProfileTreeRoot<ResidueType, ProbabilityType> globals;
      globals.readExceptPositions( is );
      // First count the positions
      std::streampos streampos_before_positions = is.tellg();
      uint32_t length = 0;
      while( !is.eof() && !is.fail() ) {
        length += 1;
        is.ignore( 100000, '\n' );
      }
      is.clear();
      if( length > 0 ) {
        prof.reinitialize( length - 1 );
      }
      prof.copyExceptPositions( globals );
      // Reset the stream reader position
      is.seekg( streampos_before_positions );

      prof.readPositions( is );

      return is;
    } // friend operator>> ( istream, ProfileTreeRoot & )

    void
    readExceptPositions (
      std::istream & is
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      this->readGlobalParameters( is );
      //is >> " ]";
      is.ignore( 2 );
      //is >> endl;
      // TODO: Is endl possibly more than one char?
      is.ignore( 1 );
    } // readExceptPositions( std::istream & )

    // Set the length first!
    void
    readPositions (
      std::istream & is
    )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        //is >> vector<Position>::operator[]( pos_i );
        is.ignore( 2 ); // "[ "
        ( *this )[ pos_i ].readPositionSpecificParameters( is );
        is.ignore( 2 ); // " ]"
        //is >> endl;
        // TODO: Is endl possibly more than one char?
        is.ignore( 1 );
      }
    } // readPositions( std::istream & os )

    /**
     * String reader.  Clobbers existing profile, replacing it with what it
     * gleans from the given string.
     */
    void
    fromString ( string const& str )
    {
      this->clear();

      std::istringstream strm(( str ));

      operator>>( strm, *this );
    } // fromString( string const& )

    /**
     * File reader.  Clobbers existing profile, replacing it with what it
     * gleans from the file with the given filename.
     */
    bool
    fromFile ( string const & filename )
    {
      return fromFile( filename.c_str() );
    } // fromFile( string )

    /**
     * File reader.  Clobbers existing profile, replacing it with what it
     * gleans from the file with the given filename.
     */
    bool
    fromFile ( const char * filename )
    {
      std::ifstream fs ( filename );

      if( !fs.is_open() ) {
        // TODO: ?
        cerr << "The profile file '" << filename << "' could not be opened." << endl;
        return false;
      } else {
        operator>>( fs, *this );
        fs.close();
      }
      return true;
    } // fromFile( const char * )

    /**
     * File reader.  Clobbers existing fasta, replacing it with what it
     * gleans from the given file.
     */
    void
    fromFile ( std::ifstream fs )
    {
      this->clear();

      if( !fs.is_open() ) {
        cerr << "The profile file is not open." << endl;
      } else {
        operator>>( fs, *this );
      }

    } // fromFile( ifstream )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      ProfileTreeRoot<ResidueType, ProbabilityType> const& prof )
    {
      prof.writeExceptPositions( os );
      prof.writePositions( os );

      return os;
    } // friend operator<< ( basic_ostream, ProfileTreeRoot const&)

    // for cout, std::cerr, common cases:
    void
    writeExceptPositions (
      std::basic_ostream<char, std::char_traits<char> >& os
    ) const
    {
      // WARNING: Exact duplicate of templated function!!!!
      os << "[ ";
      this->writeGlobalParameters( os );
      os << " ]";
      os << endl;
    } // writeExceptPositions( std::basic_ostream<char, std::char_traits<char> >& )

    template<class CharT, class Traits>
    void
    writeExceptPositions (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      // WARNING: Exact duplicate of non-templated function!!!!
      os << "[ ";
      this->writeGlobalParameters( os );
      os << " ]";
      os << endl;
    } // writeExceptPositions( basic_ostream & )

    // for cout, std::cerr, common cases:
    void
    writePositions (
      std::basic_ostream<char, std::char_traits<char> >& os
    ) const
    {
      // WARNING: Exact duplicate of templated function!!!!
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        os << vector<Position>::operator[]( pos_i ) << endl;
      }
    } // writePositions( std::basic_ostream<char, std::char_traits<char> >& os )

    template<class CharT, class Traits>
    void
    writePositions (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      // WARNING: Exact duplicate of non-templated function!!!!
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        os << vector<Position>::operator[]( pos_i ) << endl;
      }
    } // writePositions ( basic_ostream )

    uint32_t
    length () const
    {
      return vector<Position>::size();
    } // length()

    /**
     * Return the largest value stored in this ProfileTreeRoot.
     */
    ProbabilityType
    maximumValue () const
    {
      ProbabilityType largest_value( maximumValueExceptPositions() );
      ProbabilityType largest_value2( maximumValuePositions() );

      if( largest_value > largest_value2 ) {
        return largest_value;
      } else {
        return largest_value2;
      }
    } // maximumValue() const

    /**
     * Return the largest value stored in this ProfileTreeRoot (except
     * position-specific values).
     */
    ProbabilityType
    maximumValueExceptPositions () const
    {
      return GlobalParameters<ResidueType, ProbabilityType>::maximumValue();
    } // maximumValueExceptPositions() const

    /**
     * Return the largest value stored in this ProfileTreeRoot (only
     * position-specific values), or 0 if there are no positions.
     */
    ProbabilityType
    maximumValuePositions () const
    {
      // And each position.
      if( length() == 0 ) {
        return 0;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      ProbabilityType largest_value(
        vector<Position>::operator[]( 0 ).maximumValue()
      );
      ProbabilityType largest_value_tmp;
      for( pos_i = 1; pos_i <= last_pos; pos_i++ ) {
        largest_value_tmp =
          vector<Position>::operator[]( pos_i ).maximumValue();
        if( largest_value_tmp > largest_value ) {
          largest_value = largest_value_tmp;
        }
      }
      return largest_value;
    } // maximumValuePositions () const

    /**
     * Multiply each profile value by the given scalar.  Note that this may
     * violate the rule that the probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    ProfileTreeRoot<ResidueType, ProbabilityType> &
    operator*= ( AnyProbabilityType const& scalar )
    {
      multiplyByExceptPositions( scalar );
      multiplyByPositions( scalar );

      return *this;
    } // operator*= ( AnyProbabilityType const& )

    /**
     * Multiply each non-position-specific profile value by the given
     * scalar.  Note that this will violate the rule that the
     * probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    void
    multiplyByExceptPositions ( AnyProbabilityType const& scalar )
    {
      GlobalParameters<ResidueType, ProbabilityType>::operator*=( scalar );
    } // multiplyByExceptPositions( AnyProbabilityType const& )

    /**
     * Multiply each position-specific profile value by the given scalar.
     * Note that this will violate the rule that the probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    void
    multiplyByPositions ( AnyProbabilityType const& scalar )
    {
      // And each position.
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        vector<Position>::operator[]( pos_i ) *= scalar;
      }
    } // multiplyByPositions ( AnyProbabilityType const& )

    /**
     * Divide each profile value by the given denominator.  Note that this will
     * violate the rule that the probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    ProfileTreeRoot<ResidueType, ProbabilityType> &
    operator/= ( AnyProbabilityType const& denominator )
    {
      divideByExceptPositions( denominator );
      divideByPositions( denominator );

      return *this;
    } // operator/= ( AnyProbabilityType const& )

    /**
     * Divide each non-position-specific profile value by the given
     * denominator.  Note that this may violate the rule that the
     * probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    void
    divideByExceptPositions ( AnyProbabilityType const& denominator )
    {
      GlobalParameters<ResidueType, ProbabilityType>::operator/=( denominator );
    } // divideByExceptPositions( AnyProbabilityType const& )

    /**
     * Divide each position-specific profile value by the given denominator.
     * Note that this will violate the rule that the probabilities add to 1.
     */
    template <typename AnyProbabilityType>
    void
    divideByPositions ( AnyProbabilityType const& denominator )
    {
      // And each position.
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        vector<Position>::operator[]( pos_i ) /= denominator;
      }
    } // divideByPositions ( AnyProbabilityType const& )
    
    /**
     * Add to each contained distribution the corresponding values in the given
     * other Profile (which must have the same number of positions as this
     * profile).  Note that this may violate the rule that the probabilities
     * sum to 1.
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeRoot<ResidueType, ProbabilityType> &
    operator+= ( AnyInternalNodeOrRoot const& other_node )
    {
      incrementByExceptPositions( other_node );
      incrementByPositions( other_node );

      return *this;
    } // operator+=( AnyInternalNodeOrRoot const& )

    /**
     * Add to each contained distribution the corresponding values in the given
     * other Profile (which must have the same number of positions as this
     * profile).  Note that this may violate the rule that the probabilities
     * sum to 1.
     */
    template <typename AnyProbabilityType>
    ProfileTreeRoot<ResidueType, ProbabilityType> &
    operator+= ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      incrementByExceptPositions( other_node );

      return *this;
    } // operator+=( GlobalParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Add to each contained non-position-specific distribution the
     * corresponding values in the given other Profile.  Note that this
     * violates the rule that the probabilities sum to 1.
     */
    template <typename AnyProbabilityType>
    void
    incrementByExceptPositions ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      GlobalParameters<ResidueType, ProbabilityType>::operator+=( other_node );
    } // incrementByExceptPositions( GlobalParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Add to each contained position-specific distribution the corresponding
     * values in the given other Profile (which must have the same number of
     * positions).  Note that this violates the rule that the probabilities sum
     * to 1.
     */
    template <typename AnyInternalNodeOrRoot>
    void
    incrementByPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      // ASSERTION: this->length() == other_prof.length()
      // TODO: Check assertion.

      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ) += other_node[ pos_i ];
      }
    } // incrementByPositions( AnyInternalNodeOrRoot const& )

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
    template <typename AnyInternalNodeOrRoot>
    double
    crossEntropy (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      return crossEntropy( other_node, ( ProfileTreeRoot<ResidueType, double> const * )0 );
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
    template <typename AnyInternalNodeOrRoot,
              typename AnyInternalNodeOrRoot_Weights>
    double
    crossEntropy (
      AnyInternalNodeOrRoot const& other_node,
      AnyInternalNodeOrRoot_Weights const * const weights
    ) const
    {
      double cross_entropy = 0.0;
      cross_entropy +=
        crossEntropyExceptPositions( other_node, weights );
      cross_entropy +=
        crossEntropyPositions( other_node, weights );
      return cross_entropy;
    } // crossEntropy( AnyInternalNodeOrRoot const &, AnyInternalNodeOrRoot_Weights const * ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_node)), except don't include position-specific parameters.
     * The weights (if non-null) may be of any type convertible to a double,
     * and the cross entropy will be computed as E(-log(weights*other_node)).
     * Note that the cross entropy is non-symmetric (calling
     * other_node.crossEntropy( *this, weights ) will return a different
     * value).
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
    template <typename AnyInternalNodeOrRoot,
              typename AnyInternalNodeOrRoot_Weights>
    double
    crossEntropyExceptPositions (
      AnyInternalNodeOrRoot const& other_node,
      AnyInternalNodeOrRoot_Weights const * const weights
    ) const
    {
      return
        GlobalParameters<ResidueType, ProbabilityType>::crossEntropy( other_node, weights );
    } // crossEntropyExceptPositions( AnyInternalNodeOrRoot const &, AnyInternalNodeOrRoot_Weights const * ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_node)), except include only the position-specific
     * parameters.  The weights (if non-null) may be of any type convertible to
     * a double, and the cross entropy will be computed as
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
    template <typename AnyInternalNodeOrRoot,
              typename AnyInternalNodeOrRoot_Weights>
    double
    crossEntropyPositions (
      AnyInternalNodeOrRoot const& other_node,
      AnyInternalNodeOrRoot_Weights const * const weights
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
      //cout << "crossEntropyPositions(..):" << endl;
      //double ce;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        cross_entropy +=
          vector<Position>::operator[]( pos_i ).crossEntropy(
            other_node[ pos_i ],
            ( weights ? &( *weights )[ pos_i ] : ( ProfilePosition<ResidueType, double> const * )0 )
        );
        //ce = 
        //  vector<Position>::operator[]( pos_i ).crossEntropy(
        //    other_node[ pos_i ],
        //    ( weights ? &( *weights )[ pos_i ] : ( ProfilePosition<ResidueType, double> const * )0 )
        //);
        //cout << "[ " << pos_i << " ] " << ( ce - vector<Position>::operator[]( pos_i ).crossEntropy(
        //  ( *this )[ pos_i ],
        //    ( weights ? &( *weights )[ pos_i ] : ( ProfilePosition<ResidueType, double> const * )0 )
        //) )
        //<< endl;
        //cross_entropy += ce;
      }
      return cross_entropy;
    } // crossEntropyPositions( AnyInternalNodeOrRoot const &, AnyInternalNodeOrRoot_Weights const * ) const

    /**
     * Calculate and return the Euclidean distance between this profile and
     * another profile (treating every probability as an orthogonal dimension).
     * The other profile must have the same number of positions as this
     * profile.
     */
    template <typename AnyInternalNodeOrRoot>
    double
    euclideanDistance (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      return sqrt( euclideanDistanceSquared( other_node ) );
    } // euclideanDistance( AnyInternalNodeOrRoot const& ) const

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension).  The other profile must have the same number of positions as
     * this profile.
     */
    template <typename AnyInternalNodeOrRoot>
    double euclideanDistanceSquared (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      double squared_euclidean_distance = 0.0;
      squared_euclidean_distance +=
        euclideanDistanceSquaredExceptPositions( other_node );
      squared_euclidean_distance +=
        euclideanDistanceSquaredPositions( other_node );
      return squared_euclidean_distance;
    } // euclideanDistanceSquared( AnyInternalNodeOrRoot const& ) const

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension), but only include the non-position-specific parameters.
     */
    template <typename AnyInternalNodeOrRoot>
    double euclideanDistanceSquaredExceptPositions (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      return GlobalParameters<ResidueType, ProbabilityType>::euclideanDistanceSquared( other_node );
    } // euclideanDistanceSquaredExceptPositions( AnyInternalNodeOrRoot const& ) const

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension), but only include the position-specific parameters.  The
     * other profile must have the same number of positions as this profile.
     */
    template <typename AnyInternalNodeOrRoot>
    double
    euclideanDistanceSquaredPositions (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      assert( this->length() == other_node.length() );

      double squared_euclidean_distance = 0.0;
      // And each position.
      if( length() == 0 ) {
        return squared_euclidean_distance;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        squared_euclidean_distance +=
          vector<Position>::operator[]( pos_i ).euclideanDistanceSquared(
            other_node[ pos_i ]
        );
      }
      return squared_euclidean_distance;
    } // euclideanDistanceSquaredPositions( AnyInternalNodeOrRoot const& ) const

    /**
     * An alias for euclideanDistanceSquaredPositions.
     */
    template <typename AnyInternalNodeOrRoot>
    double
    euclideanDistanceSquaredUniquePositions (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      return euclideanDistanceSquaredPositions( other_node );
    } // euclideanDistanceSquaredUniquePositions( ProfileTreeInternalNode const& ) const

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const
    {
      uint32_t free_params = 0;
      free_params +=
        freeParameterCountExceptPositions();
      free_params +=
        freeParameterCountPositions();

      return free_params;
    } // freeParameterCount() const

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCountExceptPositions () const
    {
      return GlobalParameters<ResidueType, ProbabilityType>::freeParameterCount();
    } // freeParameterCountExceptPositions() const

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCountPositions () const
    {
      // And each position.
      if( length() == 0 ) {
        return 0;
      }

      uint32_t free_params = 0;

      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        free_params +=
          vector<Position>::operator[]( pos_i ).freeParameterCount();
      }

      return free_params;
    } // freeParameterCountPositions() const

    /**
     * An alias for freeParameterCountPositions().
     */
    uint32_t
    freeParameterCountUniquePositions () const
    {
      return freeParameterCountPositions();
    }

    /**
     * Since I am overloading [] to accept global parameters, I need to
     * explicitly tell the compiler to still use the vector [ int ] definition.
     */
    Position &
    operator[] ( uint32_t pos_i )
    {
      return vector<Position>::operator[]( pos_i );
    } // operator[] ( uint32_t )

    /**
     * Since I am overloading [] to accept global parameters, I need to
     * explicitly tell the compiler to still use the vector [ int ] definition.
     */
    Position const&
    operator[] ( uint32_t pos_i ) const
    {
      return vector<Position>::operator[]( pos_i );
    } // const& operator[] ( uint32_t ) const

  }; // End class ProfileTreeRoot

  /**
   * A Profile is a data structure for the model parameters.  Conceptually, for
   * every position of the profile there are a bunch of parameters for each
   * kind of Plan 7 profile transition/emission.  We make this a vector of maps
   * from ProfileKeys (eg. Transition) to maps from ProfileKey-specific params
   * (eg. Transition_M_to_D) to the values of those params.
   *
   * This is a Profile that is a variation on another Profile (its parent).  It
   * can also have children, which are other ProfileTreeInternalNodes.  Its
   * parent might be a ProfileTreeInternalNode or it might be the
   * ProfileTreeRoot.
   */
// TODO: REMOVE?  For now we make trees of ProfileTreeRoots.
// TODO: Add a copy constructor, and in it make sure to also copy the m_profileTreeVertex.
  template <typename ResidueType, typename ProbabilityType>
  class ProfileTreeInternalNode :
    public vector<ProfilePosition<ResidueType, ProbabilityType> >
  {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_NVP( m_parent_when_parent_is_internal_node );
      ar & BOOST_SERIALIZATION_NVP( m_root );
      uint32_t size = this->size();
      ar & BOOST_SERIALIZATION_NVP( size );
      if( size != this->size() ) {
        this->resize( size );
      }
      ar & BOOST_SERIALIZATION_NVP( m_length )
         & BOOST_SERIALIZATION_NVP( m_indexMap )
         & BOOST_SERIALIZATION_NVP( m_parentPositionVariations );

      // Serialize node-specific positions
      if( size > 0 ) {
       uint32_t last_pos = size - 1;
       uint32_t pos_i;
       string tmp_str;
       for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
         tmp_str = "node_specific_pos_" + boost::lexical_cast<std::string>(pos_i);
         ar & boost::serialization::make_nvp( tmp_str.c_str(), this->vector<Position>::operator[]( pos_i ) );
       }
      } // End if size > 0, serialize node-specific positions
    } // serialize( Archive &, const unsigned int )

  protected:
    typedef ProfilePosition<ResidueType, ProbabilityType> Position;

  protected:
    ProfileTreeRoot<ResidueType, ProbabilityType> * m_root;
    ProfileTreeInternalNode<ResidueType, ProbabilityType> * m_parent_when_parent_is_internal_node;
  

    /**
     * The index map maps the indices of this ProfileTreeInternalNode to either indices
     * of the parent profile or of the underlying datastore.  If the index is 0
     * or positive, it refers to the corresponding index into the underlying
     * datastore.  If the number is negative, it refers to an index into the
     * parent: -1 is 0, -2, is 1, etc (-x => ( x-1 ) ).
     *
     * The length of m_indexMap is this.length() (NOT this.size()).
     */
    vector<int> m_indexMap;
    uint32_t m_length;

    uint32_t m_profileTreeVertex;

    // operators
    GALOSH_OPERATORS_DELEGATING_PositionTransitionParameters( ( *m_root ), ProbabilityType )
    GALOSH_OPERATORS_DELEGATING_InsertionEmissionParameters( ( *m_root ), ProbabilityType )
    GALOSH_OPERATORS_DELEGATING_PreAlignParameters( ( *m_root ), ProbabilityType )
    GALOSH_OPERATORS_DELEGATING_PostAlignParameters( ( *m_root ), ProbabilityType )

  public:

  vector<int> m_parentPositionVariations;

    ProfileTreeInternalNode () :
      m_parent_when_parent_is_internal_node( NULL ),
      m_root( NULL ),
      m_length( 0 ),
      m_parentPositionVariations(),
      m_profileTreeVertex( 0 )
    {
      // Do nothing else.
    } // <init>()

    ProfileTreeInternalNode ( ProfileTreeInternalNode<ResidueType, ProbabilityType> const * parent ) :
      m_parent_when_parent_is_internal_node( parent ),
      m_root( ((ProfileTreeRoot<ResidueType,ProbabilityType> * const)parent)->getProfileTreeRoot() ),
      m_length( 0 ),
      m_parentPositionVariations(),
      m_profileTreeVertex( 0 )
    {
      // Do nothing else.
    } // <init>( ProfileTreeInternalNode const * const parent )

    ProfileTreeInternalNode ( ProfileTreeRoot<ResidueType, ProbabilityType> * const parent_is_root ) :
      m_parent_when_parent_is_internal_node( NULL ),
      m_root( parent_is_root ),
      m_length( 0 ),
      m_parentPositionVariations(),
      m_profileTreeVertex( 0 )
    {
      // Do nothing else.
    } // <init>( ProfileTreeRoot * parent )

    /**
     * The parent_position_variations argument is a vector of integers of
     * length one greater than the length of the parent profile.  The number at
     * each position i of this vector indicates how this child modifies the
     * parent: if the number is <= 0, then the parent position ( i + 1 ) is
     * deleted in the child, and if it is <= -1, then that position is replaced
     * by 1 or more positions in the child.  If the number is >= 1, then the
     * parent position is preserved in the child, and if it is some number k >=
     * 2, then there are (k-1) additional positions inserted after the parent
     * position ( i + 1 ).  The absolute value of the integer at the 0th index
     * indicates how many positions the child adds before the parent's first
     * position.
     */
    ProfileTreeInternalNode (
      ProfileTreeInternalNode * parent,
      vector<int> const& parent_position_variations
    ) :
      m_parent_when_parent_is_internal_node( parent ),
      m_root( parent->getProfileTreeRoot() ),
      m_length( 0 ),
      m_parentPositionVariations(), // see reinitializeIndexMap(..)
      m_profileTreeVertex( 0 )
    {
      // TODO: REMOVE
      //cout << "<init>( ProfileTreeInternalNode *, vector<int> )" << endl;
      //cout << "Calling reinitializeIndexMap(..).  parent is root.  parentLength() returns " << parentLength() << endl;
      reinitializeIndexMap( parent_position_variations );
    } // <init>( ProfileTreeInternalNode *, vector<int> )

    /**
     * The parent_position_variations argument is a vector of integers of
     * length one greater than the length of the parent profile.  The number at
     * each position i of this vector indicates how this child modifies the
     * parent: if the number is <= 0, then the parent position ( i + 1 ) is
     * deleted in the child, and if it is <= -1, then that position is replaced
     * by 1 or more positions in the child.  If the number is >= 1, then the
     * parent position is preserved in the child, and if it is some number k >=
     * 2, then there are (k-1) additional positions inserted after the parent
     * position ( i + 1 ).  The absolute value of the integer at the 0th index
     * indicates how many positions the child adds before the parent's first
     * position.
     */
    ProfileTreeInternalNode (
      ProfileTreeRoot<ResidueType, ProbabilityType> * const parent_is_root,
      vector<int> const& parent_position_variations
    ) :
      m_parent_when_parent_is_internal_node( NULL ),
      m_root( parent_is_root ),
      m_length( 0 ),
      m_parentPositionVariations(), // see reinitializeIndexMap(..)
      m_profileTreeVertex( 0 )
    {
      // TODO: REMOVE
      //cout << "Calling reinitializeIndexMap(..).  parent is root.  parentLength() returns " << parentLength() << endl;
      reinitializeIndexMap( parent_position_variations );
    } // <init>( ProfileTreeRoot *, vector<int> )

    void
    reinitialize ( ProfileTreeInternalNode * parent )
    {
      m_parent_when_parent_is_internal_node = parent;
      m_root = parent->getProfileTreeRoot();
      m_length = 0;
      reinitializeIndexMap();
    } // reinitialize( ProfileTreeInternalNode * parent )

    void
    reinitialize ( ProfileTreeRoot<ResidueType, ProbabilityType> * const parent_is_root )
    {
      m_parent_when_parent_is_internal_node = NULL;
      m_root = parent_is_root;
      m_length = 0;
      reinitializeIndexMap();
    } // reinitialize( ProfileTreeRoot * )

    /**
     * The parent_position_variations argument is a vector of integers of
     * length one greater than the length of the parent profile.  The number at
     * each position i of this vector indicates how this child modifies the
     * parent: if the number is <= 0, then the parent position ( i + 1 ) is
     * deleted in the child, and if it is <= -1, then that position is replaced
     * by 1 or more positions in the child.  If the number is >= 1, then the
     * parent position is preserved in the child, and if it is some number k >=
     * 2, then there are (k-1) additional positions inserted after the parent
     * position ( i + 1 ).  The absolute value of the integer at the 0th index
     * indicates how many positions the child adds before the parent's first
     * position.
     */
    void
    reinitialize (
      ProfileTreeInternalNode * parent,
      vector<int> const& parent_position_variations
    )
    {
      m_parent_when_parent_is_internal_node = parent;
      m_root = parent->getProfileTreeRoot();
      reinitializeIndexMap( parent_position_variations );
    } // reinitialize( ProfileTreeInternalNode * parent, vector<int> )

    /**
     * The parent_position_variations argument is a vector of integers of
     * length one greater than the length of the parent profile.  The number at
     * each position i of this vector indicates how this child modifies the
     * parent: if the number is <= 0, then the parent position ( i + 1 ) is
     * deleted in the child, and if it is <= -1, then that position is replaced
     * by 1 or more positions in the child.  If the number is >= 1, then the
     * parent position is preserved in the child, and if it is some number k >=
     * 2, then there are (k-1) additional positions inserted after the parent
     * position ( i + 1 ).  The absolute value of the integer at the 0th index
     * indicates how many positions the child adds before the parent's first
     * position.
     */
    void
    reinitialize (
      ProfileTreeRoot<ResidueType, ProbabilityType> * const parent_is_root,
      vector<int> const& parent_position_variations
    )
    {
      m_parent_when_parent_is_internal_node = NULL;
      m_root = parent_is_root;
      // TODO: REMOVE
      //cout << "reinitialize( ProfileTreeRoot *, vector<int> )" << endl;
      //cout << "parentLength() returns " << parentLength() << endl;
      reinitializeIndexMap( parent_position_variations );
    } // reinitialize( ProfileTreeRoot *, vector<int> )

    void
    reinitialize ( ProfileTreeInternalNode<ResidueType, ProbabilityType> const& copy_parent_and_parent_position_variations_from )
    {
      if( copy_parent_and_parent_position_variations_from.m_parent_when_parent_is_internal_node == NULL ) {
        this->reinitialize(
          copy_parent_and_parent_position_variations_from.m_root,
          copy_parent_and_parent_position_variations_from.m_parentPositionVariations
        );
      } else {
        this->reinitialize(
          copy_parent_and_parent_position_variations_from.m_parent_when_parent_is_internal_node,
          copy_parent_and_parent_position_variations_from.m_parentPositionVariations
        );
      }
    } // reinitialize( ProfileTreeInternalNode<ResidueType, ProbabilityType> const& )

//    /**
//     * The parent_position_variations argument is a vector of integers of
//     * length one greater than the length of the parent profile.  The number at
//     * each position i of this vector indicates how this child modifies the
//     * parent: if the number is <= 0, then the parent position ( i + 1 ) is
//     * deleted in the child, and if it is <= -1, then that position is replaced
//     * by 1 or more positions in the child.  If the number is >= 1, then the
//     * parent position is preserved in the child, and if it is some number k >=
//     * 2, then there are (k-1) additional positions inserted after the parent
//     * position ( i + 1 ).  The absolute value of the integer at the 0th index
//     * indicates how many positions the child adds before the parent's first
//     * position.
//     */
//    void reinitialize (
//      ParentType * parent,
//      vector<int> const& parent_position_variations
//    )
//    {
//      m_parent_when_parent_is_internal_node = parent;
//      reinitializeIndexMap( parent_position_variations );
//    } // reinitialize( ParentType * parent, vector<int> )

    /**
     * This makes this ProfileChildNode do nothing: it just delegates
     * everything to its parent.
     */
    void reinitializeIndexMap ()
    {
      if( vector<Position>::size() != 0 ) {
        vector<Position>::resize( 0 );
      }
      if( m_indexMap.size() != 0 ) {
        m_indexMap.resize( 0 );
      }
      if( m_parentPositionVariations.size() != 0 ) {
        m_parentPositionVariations.resize( 0 );
      }
      m_length = 0;
    } // reinitializeIndexMap()

    /**
     * The parent_position_variations argument is a vector of integers of
     * length one greater than the length of the parent profile.  The number at
     * each position i of this vector indicates how this child modifies the
     * parent: if the number is <= 0, then the parent position ( i + 1 ) is
     * deleted in the child, and if it is <= -1, then that position is replaced
     * by 1 or more positions in the child.  If the number is >= 1, then the
     * parent position is preserved in the child, and if it is some number k >=
     * 2, then there are (k-1) additional positions inserted after the parent
     * position ( i + 1 ).  The absolute value of the integer at the 0th index
     * indicates how many positions the child adds before the parent's first
     * position.
     */
    void reinitializeIndexMap (
      vector<int> const& parent_position_variations
    )
    {
      // TODO: REMOVE
      //cout << "reinitializeIndexMap( ";
      //for( uint32_t i = 0; i < parent_position_variations.size(); i++ ) {
      //  if( i > 0 ) {
      //    cout << ", ";
      //  }
      //  cout << parent_position_variations[ i ];
      //}
      //cout << " )" << endl;

      if( parent_position_variations.size() != ( 1 + parentLength() ) ) {
        // Uh-oh.
        // TODO: Throw an exception.
        cout << "ERROR: Assertion failed.  ProfileTreeInternalNode.reinitializeIndexMap(..): the argument vector should be of length 1 + this.parentLength(), but it is not (its length is " << parent_position_variations.size() << ", but we were expecting " << ( 1 + parentLength() ) << ")." << endl;
        return;
      }
      // First count how many new positions we'll need.
      uint32_t new_position_count = 0;
      m_length = parentLength();
      uint32_t ppv_i;
      for( ppv_i = 0; ppv_i < parent_position_variations.size(); ppv_i++ ) {
        if( ppv_i == 0 ) {
          new_position_count += abs( parent_position_variations[ ppv_i ] );
        } else {
          if( parent_position_variations[ ppv_i ] <= 0 ) {
            m_length--;
            if( parent_position_variations[ ppv_i ] <= -1 ) {
              new_position_count += ( -parent_position_variations[ ppv_i ] );
            }
          } else if( parent_position_variations[ ppv_i ] >= 2 ) {
            new_position_count += ( parent_position_variations[ ppv_i ] - 1 );
          }
        }
      }
      // Our length is:
      m_length += new_position_count;

      // We can now allocate our underlying array of ProfilePositions.
      if( vector<Position>::size() != new_position_count ) {
        vector<Position>::resize( new_position_count );
      }

      // TODO: REMOVE
      //cout << "m_length is " << m_length << endl;
      //cout << "new_position_count is " << new_position_count << endl;;

      // And initialize them.
      uint32_t new_pos_i = 0; // this is the next index into the underlying array
      for( new_pos_i = 0; new_pos_i < new_position_count; new_pos_i++ ) {
        this->vector<Position>::operator[]( new_pos_i ).reinitialize( getProfileTreeRoot() );
      }

      // And we can allocate our indexMap.
      if( m_indexMap.size() != m_length ) {
        m_indexMap.resize( m_length );
      }
      // We also need to set up the index map.
      uint32_t pos_i = 0; // this is the next index as the outside world sees it
      new_pos_i = 0; // this is the next index into the underlying array
      uint32_t pos_i_loop_until;
      // For ppv_i == 0:
      for( pos_i = 0, new_pos_i = 0; pos_i < abs( parent_position_variations[ 0 ] ); pos_i++, new_pos_i++ ) {
        m_indexMap[ pos_i ] = new_pos_i; // pos_i == new_pos_i...
      }
      for( ppv_i = 1; ppv_i < parent_position_variations.size(); ppv_i++ ) {
        // Remember: parent position is ppv_i - 1.
        // First, if we are supposed to preserve the parent position, do so.
        if( parent_position_variations[ ppv_i ] > 0 ) {
          m_indexMap[ pos_i ] = -( ( ppv_i - 1 ) + 1 ); // -( parent_pos_i + 1 )
          pos_i += 1;
        }
        // Now handle new positions.
        if( ( parent_position_variations[ ppv_i ] > 1 ) ||
            ( parent_position_variations[ ppv_i ] < 0 ) ) {
          pos_i_loop_until =
            ( ( parent_position_variations[ ppv_i ] > 1 ) ?
              ( pos_i + ( parent_position_variations[ ppv_i ] - 1 ) - 1 ) :
              ( pos_i + ( -parent_position_variations[ ppv_i ] ) - 1 ) );
          for( ; pos_i <= pos_i_loop_until; pos_i++, new_pos_i++ ) {
            m_indexMap[ pos_i ] = new_pos_i;
          }
          pos_i = pos_i_loop_until + 1;
        }
      } // End foreach ppv_i, update the indexMap.

      // TODO: REMOVE
      //cout << "m_indexMap is [ ";
      //for( uint32_t i = 0; i < m_indexMap.size(); i++ ) {
      //  if( i > 0 ) {
      //    cout << ", ";
      //  }
      //  cout << m_indexMap[ i ];
      //}
      //cout << " ]" << endl;

      // And save the parentPositionVariations for later.
      m_parentPositionVariations = parent_position_variations;
    } // reinitializeIndexMap( vector<int> const& )

    uint32_t
    getProfileTreeVertex () const
    {
      // TODO: REMOVE
      //cout << "getProfileTreeVertex(): returning " << m_profileTreeVertex << "." << endl;
      return m_profileTreeVertex;
    } // getProfileTreeVertex() const

    // TODO: Should use ProfileTree.hpp's vertex_t.
    void
    setProfileTreeVertex ( uint32_t const & new_vertex )
    {
      // TODO: REMOVE
      //cout << "setProfileTreeVertex( " << new_vertex << " )" << endl;
      m_profileTreeVertex = new_vertex;
    } // setProfileTreeVertex( uint32_t const & )

    /**
     * Start using the parent position at the given index.  Uses the parent
     * position with the same index as that given.
     *
     * NOTE for now this does not update m_parentPositionVariations
     */
    // TODO: Something fancier to handle the possibility that they are
    // different lengths, and to allow using a different parent position...
    void useParentPosition ( uint32_t const & pos_i )
    {
      m_indexMap[ pos_i ] = -pos_i - 1;
      // TODO: update m_parentPositionVariations
    } // useParentPosition( uint32_t const & )

    /**
     * Returns a pointer to the root of the tree (by delegating up).
     */
    ProfileTreeRoot<ResidueType, ProbabilityType> * const
    getProfileTreeRoot ()
    {
      return m_root;
    } // getProfileTreeRoot() const

    /**
     * Sometimes, for reasons I have yet to figure out, profile positions point
     * to the wrong root.  Until I figure out why/where, I use this as a
     * workaround.  It just goes through and sets each position's root to be
     * this node's root.
     */
    void
    ensurePositionsKnowTheirRoot ()
    {
      uint32_t last_pos = vector<Position>::size() - 1;
      uint32_t pos_i; // ( here pos_i refers to position into the internal datastore ).
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ).setProfileTreeRoot( getProfileTreeRoot() );
      }
    } // ensurePositionsKnowTheirRoot()

    /**
     * Affects the whole tree (delegates to root).
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    operator= ( AnyInternalNodeOrRoot const& other_node )
    {
      m_root->operator=( other_node );

      return *this;
    } // operator=( AnyInternalNodeOrRoot const& )

    /**
     * Overridden to affect just this node.
     *
     * This sets the parent of this ProfileTreeInternalNode to be that of the
     * given one, and sets all other values to be identical to those of the
     * other node as well.
     */
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    operator= ( ProfileTreeInternalNode<ResidueType, ProbabilityType> const& other_node )
    {
      copyPositions( other_node );

      return *this;
    } // operator=( ProfileTreeInternalNode<ResidueType, ProbabilityType> const& )

    /**
     * Copy all values from the other profile's multinomial distributions into
     * this profile's, except for those that are position-specific.
     *
     * Affects the whole tree (delegates to root).
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    copyExceptPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      m_root->copyExceptPositions( other_node );

      return *this;
    } // copyExceptPositions( AnyInternalNodeOrRoot const& )

    /**
     * Copy the position-specific values from the other profile's multinomial
     * distributions into this profile's.
     *
     * Affects the whole tree (delegates to root).
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    copyPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      m_root->copyPositions( other_node );

      return *this;
    } // copyPositions( AnyInternalNodeOrRoot const& )

    /**
     * Overridden to affect just this node.
     *
     * This sets the parent of this ProfileTreeInternalNode to be that of the
     * given one, and sets all positional variations to be identical to those
     * of the other node as well.
     */
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    copyPositions ( ProfileTreeInternalNode<ResidueType, ProbabilityType> const& other_node )
    {
      m_parent_when_parent_is_internal_node =
        other_node.m_parent_when_parent_is_internal_node;
      m_root = other_node.m_root;

      if( other_node.size() != vector<Position>::size() ) {
        vector<Position>::resize( other_node.size() );
      }

      // TODO: REMOVE
      //cout << "WARNING: RESETTING INDEX MAP!!!" << endl;

      // NOTE: it should work to use the vector<int> copy constructor here.
      m_indexMap = other_node.m_indexMap;
      //if( other_node.m_indexMap.size() != m_indexMap.size() ) {
      //  m_indexMap.resize( other_node.m_indexMap.size() );
      //}
      if( other_node.length() != m_length ) {
        m_length = other_node.length();
      }

      // Now copy the additional positions.
      uint32_t last_pos = vector<Position>::size() - 1;
      uint32_t pos_i; // ( here pos_i refers to position into the internal datastore ).
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->vector<Position>::operator[]( pos_i ) = other_node[ pos_i ];
        this->vector<Position>::operator[]( pos_i ).setProfileTreeRoot( getProfileTreeRoot() );
      }
      // And copy the index map values.
      // ( here pos_i refers to the position as the outside world sees it ).
      // NOTE: This has been replaced by using the copy constructor (above).
      //last_pos = m_length - 1;
      //for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
      //  m_indexMap[ pos_i ] = other_node.m_indexMap[ pos_i ];
      //}

      return *this;
    } // copyPositions( ProfileTreeInternalNode<ResidueType, ProbabilityType> const& )

    /**
     * Copy the position-specific values from the other profile's multinomial
     * distributions into this profile's.
     *
     * May affect the whole tree.  Copies every node, including those that are
     * unique to this node and those that are not.
     *
     * The assumption is that the other node has the same length() as this one.
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    copyAllPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      // Note we assume the lengths are the same.  JIC I take the min.
      uint32_t last_pos = min( length(), other_node.length() ) - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        ( *this )[ pos_i ] = other_node[ pos_i ];
        // Redundant
        //( *this )[ pos_i ].setProfileTreeRoot( getProfileTreeRoot() );
      }

      // TODO: REMOVE
      //cout << "copyAllPositions(): " << endl;
      //cout << "\tlength() is " << length() << endl;
      //cout << "\tsize() is " << vector<ProfilePosition<ResidueType, ProbabilityType> >::size() << endl;
      //// TODO: REMOVE
      //cout << "\tm_indexMap is [ ";
      //for( uint32_t i = 0; i < m_indexMap.size(); i++ ) {
      //  if( i > 0 ) {
      //    cout << ", ";
      //  }
      //  cout << m_indexMap[ i ];
      //}
      //cout << " ]" << endl;
      //cout << "\t*getProfileTreeRoot() returns:" << endl;
      //cout << *getProfileTreeRoot() << endl;
      //cout << "\tgetProfileTreeRoot()->length() returns:" << endl;
      //cout << getProfileTreeRoot()->length() << endl;

      return *this;
    } // copyAllPositions( AnyInternalNodeOrRoot const& )

    // NOTE:
    // The zero(), even(), etc. functions (incl +=, etc) affect the "whole
    // profile" -- that is, the positions of this profile as seen from the
    // outside world.  Since that might not be what is meant, in future we may
    // need to add another set or two of these functions (for zeroing the
    // entire tree, or for zeroing only the additional positions, or the tree
    // from here up, down, etc.).

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will zero() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.  It will also zero the ProfileTreeRoot's global distributions.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions zeroed, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.  It will also cause the entire tree to have
     * zeroed global distributions.
     */
    void zero ()
    {
      zeroExceptPositions();
      zeroPositions();
    } // zero()

    /**
     * Set all values to 0, except position-specific values.  Note that this
     * violates the rule that the values sum to 1.
     *
     * Affects the whole tree (delegates to root).
     */
    void zeroExceptPositions ()
    {
      m_root->zeroExceptPositions();
    } // zeroExceptPositions()

    /**
     * Set all position-specific values to 0.  Note that this violates the rule
     * that the values sum to 1.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will zero() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions zeroed, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    void zeroPositions ()
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->operator[]( pos_i ).zero();
      }
    } // zeroPositions()

    /**
     * Set all values such that each distrubution is evenly distributed.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will even() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.  It will also even the ProfileTreeRoot's global distributions.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions evened, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.  It will also cause the entire tree to have
     * evened global distributions.
     */
    void
    even ()
    {
      evenExceptPositions();
      evenPositions();
    } // even()

    /**
     * Set all values (except position-specific values) such that each
     * distrubution is evenly distributed.
     *
     * Affects the whole tree (delegates to root).
     */
    void
    evenExceptPositions ()
    {
      m_root->evenExceptPositions();
    } // evenExceptPositions()

    /**
     * Set allposition-specific values such that each distrubution is evenly
     * distributed.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will even() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions evened, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    void
    evenPositions ()
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->operator[]( pos_i ).even();
      }
    } // evenPositions()

    /**
     * Calculate the total of all contained values.
     */
    ProbabilityType
    total () const
    {
      ProbabilityType t(
        totalExceptPositions()
      );
      t +=
        totalPositions();
      return t;
    } // total() const

    /**
     * Calculate the total of all contained values (except position-specific
     * values).
     */
    ProbabilityType
    totalExceptPositions () const
    {
      return
        m_root->totalExceptPositions();
    } // totalExceptPositions() const

    /**
     * Calculate the total of all contained position-specific values.
     */
    ProbabilityType
    totalPositions () const
    {
      ProbabilityType t( 0 );
      if( length() == 0 ) {
        return t;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        t +=
          this->operator[]( pos_i ).total();
      }
      return t;
    } // totalPositions() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
    template <typename other_type>
    void normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will normalize() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.  It will also normalize the ProfileTreeRoot's global distributions.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions normalized, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.  It will also cause the entire tree to have
     * normalized global distributions.
     */
    void normalize ( ProbabilityType const & min )
    {
      normalizeExceptPositions( min );
      normalizePositions( min );
    } // normalize( ProbabilityType const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Do all distributions
     * except for those stored on a position-specific basis.
     */
    template <typename other_type>
    void normalizeExceptPositions ( other_type const & min )
    {
      //normalizeExceptPositions( min );
      normalizeExceptPositions( static_cast<ProbabilityType>( min ) );
    } // normalizeExceptPositions( other_type const & )

    /**
     * Set all values (except position-specific values) such that each
     * distrubution is normalized.
     *
     * Affects the whole tree (delegates to root).
     */
    void normalizeExceptPositions ( ProbabilityType min )
    {
      m_root->normalizeExceptPositions( min );
    } // normalizeExceptPositions( ProbabilityType )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Do only those
     * distributions that are stored on a position-specific basis.
     */
    template <typename other_type>
    void normalizePositions ( other_type const & min )
    {
      //normalizePositions( min );
      normalizePositions( static_cast<ProbabilityType>( min ) );
    } // normalizePositions( other_type const & )

    /**
     * Set all position-specific values such that each distrubution is
     * normalized.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will normalize() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions normalized, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    void normalizePositions ( ProbabilityType min )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->operator[]( pos_i ).normalize( min );
      }
    } // normalizePositions( ProbabilityType )

    /**
     * Set all values such that each distrubution is randomly distributed.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will uniform() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.  It will also uniform the ProfileTreeRoot's global distributions.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions uniformed, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.  It will also cause the entire tree to have
     * uniformed global distributions.
     */
    void uniform ( Random & random )
    {
      uniformExceptPositions( random );
      uniformPositions( random );
    } // uniform( Random & )

    /**
     * Set all values (except position-specific values) such that each
     * distrubution is randomly distributed.
     *
     * Affects the whole tree (delegates to root).
     */
    void uniformExceptPositions ( Random & random )
    {
      m_root->uniformExceptPositions( random );
    } // uniformExceptPositions( Random & )

    /**
     * Set all position-specific values such that each distrubution is randomly
     * distributed.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will uniform() the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions uniformed, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    void uniformPositions ( Random & random )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->operator[]( pos_i ).uniform( random );
      }
    } // uniformPositions( Random & )

    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      ProfileTreeInternalNode<ResidueType, ProbabilityType> const& prof )
    {
      prof.writeExceptPositions( os );
      prof.writePositions( os );

      return os;
    } // friend operator<< ( basic_ostream, ProfileTreeInternalNode const&)

    /**
     * delegates to root.
     */
    void writeExceptPositions (
      std::basic_ostream<char, std::char_traits<char> >& os
    ) const
    {
      m_root->writeExceptPositions( os );
    } // writeExceptPositions( std::basic_ostream<char, std::char_traits<char> >& os )

    /**
     * delegates to root.
     */
    template<class CharT, class Traits>
    void writeExceptPositions (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      m_root->writeExceptPositions( os );
    } // writeExceptPositions( basic_ostream & )

    void writePositions (
      std::basic_ostream<char, std::char_traits<char> >& os
    ) const
    {
      // WARNING: Exact duplicate of templated function!
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        os << ( *this )[ pos_i ] << endl;
      }
    } // writePositions( std::basic_ostream<char, std::char_traits<char> >& )

    template<class CharT, class Traits>
    void writePositions (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      // WARNING: Exact duplicate of non-templated function!
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        os << this->operator[]( pos_i ) << endl;
      }
    } // writePositions ( basic_ostream )

    uint32_t length () const
    {
      return m_length;
    }

    uint32_t parentLength () const
    {
      if( m_parent_when_parent_is_internal_node == NULL ) {
        // TODO: REMOVE
        //cout << "parentLength() parent is root" << endl;
        return m_root->length();
      } else {
        // TODO: REMOVE
        //cout << "parentLength() parent is internal node" << endl;
        return m_parent_when_parent_is_internal_node->length();
      }
    } // parentLength()

    /**
     * Divide each profile value by the given denominator.  Note that this will
     * violate the rule that the probabilities add to 1.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will /= the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions altered, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    template <typename AnyProbabilityType>
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    operator/= ( AnyProbabilityType const& denominator )
    {
      divideByExceptPositions( denominator );
      divideByPositions( denominator );

      return *this;
    } // operator/= ( AnyProbabilityType const& )

    /**
     * Divide each non-position-specific profile value by the given
     * denominator.  Note that this will violate the rule that the
     * probabilities add to 1.
     *
     * Affects the whole tree (delegates to root).
     */
    template <typename AnyProbabilityType>
    void
    divideByExceptPositions ( AnyProbabilityType const& denominator )
    {
      m_root->divideByExceptPositions( denominator );
    } // divideByExceptPositions( AnyProbabilityType const& )

    /**
     * Divide each position-specific profile value by the given denominator.
     * Note that this will violate the rule that the probabilities add to 1.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will /= the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions altered, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    template <typename AnyProbabilityType>
    void
    divideByPositions ( AnyProbabilityType const& denominator )
    {
      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->operator[]( pos_i ) /= denominator;
      }
    } // divideByPositions ( AnyProbabilityType const& )
    
    /**
     * Add to each contained distribution the corresponding values in the given
     * other Profile (which must have the same number of positions as this
     * profile).  Note that this violates the rule that the probabilities sum
     * to 1.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will += the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions altered, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    template <typename AnyInternalNodeOrRoot>
    ProfileTreeInternalNode<ResidueType, ProbabilityType> &
    operator+= ( AnyInternalNodeOrRoot const& other_node )
    {
      incrementByExceptPositions( other_node );
      incrementByPositions( other_node );

      return *this;
    } // operator+=( AnyInternalNodeOrRoot const& )

    /**
     * Add to each contained non-position-specific distribution the
     * corresponding values in the given other Profile.  Note that this
     * violates the rule that the probabilities sum to 1.
     *
     * Affects the whole tree (delegates to root).
     */
    template <typename AnyInternalNodeOrRoot>
    void
    incrementByExceptPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      m_root->incrementByExceptPositions( other_node );
    } // incrementByExceptPositions( AnyInternalNodeOrRoot const& )

    /**
     * Add to each contained position-specific distribution the corresponding
     * values in the given other Profile (which must have the same number of
     * positions).  Note that this violates the rule that the probabilities sum
     * to 1.
     *
     * This affects the entirety of this profile as seen from the outside
     * world.  That is, it will += the additional positions and all
     * positions of the parent profile that are not deleted or replaced by this
     * child.
     *
     * WARNING: This may leave the parent profile (and other ancestors) with
     * some positions altered, some not.  Since other descendents of the
     * ancestor profiles might use those positions, this will have a complex
     * effect on the whole tree.
     */
    template <typename AnyInternalNodeOrRoot>
    void
    incrementByPositions ( AnyInternalNodeOrRoot const& other_node )
    {
      // ASSERTION: this->length() == other_prof.length()
      // TODO: Check assertion.

      if( length() == 0 ) {
        return;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        this->operator[]( pos_i ) += other_node[ pos_i ];
      }
    } // incrementByPositions( AnyInternalNodeOrRoot const& )

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
    template <typename AnyInternalNodeOrRoot>
    double
    crossEntropy (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      return crossEntropy( other_node, ( ProfileTreeInternalNode<ResidueType, double> const * )0 );
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
    template <typename AnyInternalNodeOrRoot,
              typename AnyInternalNodeOrRoot_Weights>
    double
    crossEntropy (
      AnyInternalNodeOrRoot const& other_node,
      AnyInternalNodeOrRoot_Weights const * const weights
    ) const
    {
      double cross_entropy = 0.0;
      cross_entropy +=
        crossEntropyExceptPositions( other_node, weights );
      cross_entropy +=
        crossEntropyPositions( other_node, weights );
      return cross_entropy;
    } // crossEntropy( AnyInternalNodeOrRoot const &, AnyInternalNodeOrRoot_Weights const * ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_node)), except don't include position-specific parameters.
     * The weights (if non-null) may be of any type convertible to a double,
     * and the cross entropy will be computed as E(-log(weights*other_node)).
     * Note that the cross entropy is non-symmetric (calling
     * other_node.crossEntropy( *this, weights ) will return a different
     * value).
     *
     * This can be used to calculate the (possibly weighted) (self-)entropy by
     * calling this->crossEntropy( *this, weights ).  It can also be used to
     * calculate the KL divergence by taking the difference of the
     * cross-entropy and the self-entropy, or the symmeterized KL divergence by
     * summing the KL divergences computed both ways.  Note, though, that
     * weights should all be the same within any Multinomial distribution if
     * this is to be used to calculate weighted entropies or KL divergences.
     * Otherwise the usual properties of these metrics will be violated.
     *
     * delegates to root.
     */
    template <typename AnyInternalNodeOrRoot,
              typename AnyInternalNodeOrRoot_Weights>
    double
    crossEntropyExceptPositions (
      AnyInternalNodeOrRoot const& other_node,
      AnyInternalNodeOrRoot_Weights const * const weights
    ) const
    {
      return
        m_root->crossEntropyExceptPositions( other_node, weights );
    } // crossEntropyExceptPositions( AnyInternalNodeOrRoot const &, AnyInternalNodeOrRoot_Weights const * ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_node)), except include only the position-specific
     * parameters.  The weights (if non-null) may be of any type convertible to
     * a double, and the cross entropy will be computed as
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
    template <typename AnyInternalNodeOrRoot,
              typename AnyInternalNodeOrRoot_Weights>
    double
    crossEntropyPositions (
      AnyInternalNodeOrRoot const& other_node,
      AnyInternalNodeOrRoot_Weights const * const weights
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
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        cross_entropy +=
          vector<Position>::operator[]( pos_i ).crossEntropy(
            other_node[ pos_i ],
            ( weights ? &( *weights )[ pos_i ] : ( ProfilePosition<ResidueType, double> const * )0 )
        );
      }
      return cross_entropy;
    } // crossEntropyPositions( AnyInternalNodeOrRoot const &, AnyInternalNodeOrRoot_Weights const * ) const

    /**
     * Calculate and return the Euclidean distance between this profile and
     * another profile (treating every probability as an orthogonal dimension).
     * The other profile must have the same number of positions as this
     * profile.
     */
    template <typename AnyInternalNodeOrRoot>
    double
    euclideanDistance (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      return sqrt( euclideanDistanceSquared( other_node ) );
    } // euclideanDistance( AnyInternalNodeOrRoot const& ) const

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension).  The other profile must have the same number of positions as
     * this profile.
     */
    template <typename AnyInternalNodeOrRoot>
    double euclideanDistanceSquared (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      double squared_euclidean_distance = 0.0;
      squared_euclidean_distance +=
        euclideanDistanceSquaredExceptPositions( other_node );
      squared_euclidean_distance +=
        euclideanDistanceSquaredPositions( other_node );
      return squared_euclidean_distance;
    } // euclideanDistanceSquared( AnyInternalNodeOrRoot const& ) const
    
    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension), but only include the non-position-specific parameters.
     *
     * delegates to root.
     */
    template <typename AnyInternalNodeOrRoot>
    double euclideanDistanceSquaredExceptPositions (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      return m_root->euclideanDistanceSquaredExceptPositions( other_node );
    } // euclideanDistanceSquaredExceptPositions( AnyInternalNodeOrRoot const& ) const

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension), but only include the position-specific parameters.  The
     * other profile must have the same number of positions as this profile.
     */
    template <typename AnyInternalNodeOrRoot>
    double euclideanDistanceSquaredPositions (
      AnyInternalNodeOrRoot const& other_node
    ) const
    {
      // ASSERTION: this->length() == other_node.length()
      // TODO: Check assertion.

      double squared_euclidean_distance = 0.0;
      // And each position.
      if( length() == 0 ) {
        return squared_euclidean_distance;
      }
      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        squared_euclidean_distance +=
          this->operator[]( pos_i ).euclideanDistanceSquared(
            other_node[ pos_i ]
        );
      }
      return squared_euclidean_distance;
    } // euclideanDistanceSquaredPositions( AnyInternalNodeOrRoot const& ) const

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension), but only include the position-specific parameters.  This
     * method assumes that the given profile is a child of the same parent as
     * this one, with the same pattern of replaced positions.  Just those
     * positions are used for the calculation.
     */
    double euclideanDistanceSquaredUniquePositions (
      ProfileTreeInternalNode const& other_node
    ) const
    {
      // ASSERTION: parents are the same
      // ASSERTION: this->size() == other_node.size()
      // ASSERTION: m_indexMaps are the same
      // TODO: Check assertions.

      double squared_euclidean_distance = 0.0;
      // And each position.
      if( this->size() == 0 ) {
        return squared_euclidean_distance;
      }
      uint32_t last_pos = this->size() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        squared_euclidean_distance +=
          vector<Position>::operator[]( pos_i ).euclideanDistanceSquared(
            other_node[ pos_i ]
        );
      }
      return squared_euclidean_distance;
    } // euclideanDistanceSquaredUniquePositions( ProfileTreeInternalNode const& ) const

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCount () const
    {
      uint32_t free_params = 0;
      free_params +=
        freeParameterCountExceptPositions();
      free_params +=
        freeParameterCountPositions();

      return free_params;
    } // freeParameterCount() const

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
    uint32_t
    freeParameterCountExceptPositions () const
    {
      return m_root->freeParameterCountExceptPositions();
    } // freeParameterCountExceptPositions() const

    /**
     * How many free parameters are there?  This is the sum of the free
     * parameters in the contained distributions.
     */
    uint32_t
    freeParameterCountPositions () const
    {
      // And each position.
      if( length() == 0 ) {
        return 0;
      }

      uint32_t free_params = 0;

      uint32_t last_pos = length() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        free_params +=
          this->operator[]( pos_i ).freeParameterCount();
      }

      return free_params;
    } // freeParameterCountPositions() const

    /**
     * How many free parameters are there?  This is the sum of the free
     * parameters in the contained distributions, including only those
     * positions that differ between this node and its parent.
     */
    uint32_t
    freeParameterCountUniquePositions () const
    {
      // And each position.
      if( this->size() == 0 ) {
        return 0;
      }

      uint32_t free_params = 0;

      uint32_t last_pos = this->size() - 1;
      uint32_t pos_i;
      for( pos_i = 0; pos_i <= last_pos; pos_i++ ) {
        free_params +=
          vector<Position>::operator[]( pos_i ).freeParameterCount();
      }

      return free_params;
    } // freeParameterCountUniquePositions() const

    /**
     * Since I am overloading [] to accept global parameters, I need to
     * explicitly tell the compiler to still use the vector [ int ] definition.
     *
     * This sometimes delegates to parent.
     */
    Position &
    operator[] ( uint32_t pos_i )
    {
      // ASSERT: pos_i < length() and >0.
      if( m_indexMap[ pos_i ] < 0 ) {
        if( m_parent_when_parent_is_internal_node == NULL ) {
          return m_root->operator[]( -( m_indexMap[ pos_i ] + 1 ) );
        } else {
          return m_parent_when_parent_is_internal_node->operator[]( -( m_indexMap[ pos_i ] + 1 ) );
        }
      } else {
        return this->vector<Position>::operator[]( m_indexMap[ pos_i ] );
      }
    } // operator[] ( uint32_t )

    /**
     * Since I am overloading [] to accept global parameters, I need to
     * explicitly tell the compiler to still use the vector [ int ] definition.
     *
     * This sometimes delegates to parent.
     */
    Position const&
    operator[] ( uint32_t pos_i ) const
    {
      if( m_indexMap[ pos_i ] < 0 ) {
        if( m_parent_when_parent_is_internal_node == NULL ) {
          return m_root->operator[]( -( m_indexMap[ pos_i ] + 1 ) );
        } else {
          return m_parent_when_parent_is_internal_node->operator[]( -( m_indexMap[ pos_i ] + 1 ) );
        }
      } else {
        return this->vector<Position>::operator[]( m_indexMap[ pos_i ] );
      }
    } // const& operator[] const ( uint32_t )

  }; // End class ProfileTreeInternalNode

  //======//// traits classes ////========//

  template <typename UnusedParameterCollectionType>
  struct has_m_scalar : boost::false_type
  {}; // generic has_m_scalar

  template <typename ScoreType, typename ParameterCollectionType>
  inline
  ScoreType const &
  getScalarImplConst ( ParameterCollectionType const & pc, boost::false_type const & )
  {
    static const ScoreType one = 1.0;
    return one;
  } // getScalarImplConst( ParameterCollectionType const &, boost::false_type const & )

  // TODO: REMOVE/FIX.  This should never happen.
  template <typename ScoreType, typename ParameterCollectionType>
  inline
  ScoreType &
  getScalarImpl ( ParameterCollectionType & pc, boost::false_type const & )
  {
    // TODO: ?
    //boost::static_assert<boost::false_type>();
    assert( false );
    static ScoreType one = 1.0;
    return one;
  } // getScalarImpl( ParameterCollectionType const &, boost::false_type const & )

  template <typename ScoreType, typename ParameterCollectionType>
  inline
  ScoreType &
  getScalarImpl ( ParameterCollectionType & pc, boost::true_type const & )
  {
    return pc.m_scalar;
  } // getScalarImpl( ParameterCollectionType &, boost::true_type const & )

  template <typename ScoreType, typename ParameterCollectionType>
  inline
  ScoreType const &
  getScalarImplConst ( ParameterCollectionType const & pc, boost::true_type const & )
  {
    return pc.m_scalar;
  } // getScalarImplConst( ParameterCollectionType const &, boost::true_type const & )

  template <typename ScoreType, typename ParameterCollectionType>
  inline
  ScoreType &
  getScalar ( ParameterCollectionType & pc )
  {
    return getScalarImpl<ScoreType, ParameterCollectionType>( pc, has_m_scalar<ParameterCollectionType>() );
  } // getScalar( ParameterCollectionType & )

  template <typename ScoreType, typename ParameterCollectionType>
  inline
  ScoreType const &
  getScalar ( ParameterCollectionType const & pc )
  {
    return getScalarImplConst<ScoreType, ParameterCollectionType>( pc, has_m_scalar<ParameterCollectionType>() );
  } // getScalar( ParameterCollectionType const & )

  template <typename UnusedParameterCollectionProbabilityType>
  struct parameter_collection_traits
  {
  }; // generic parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<MatchEmissionParameters<ParameterCollectionResidueType, ParameterCollectionProbabilityType> >
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // MatchEmissionParameters parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<InsertionEmissionParameters<ParameterCollectionResidueType, ParameterCollectionProbabilityType> > 
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // InsertionEmissionParameters parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<PositionTransitionParameters<ParameterCollectionResidueType, ParameterCollectionProbabilityType> > 
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // PositionTransitionParameters parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<PositionSpecificParameters<ParameterCollectionResidueType, ParameterCollectionProbabilityType> > 
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // PositionSpecificParameters parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<PreAlignParameters<ParameterCollectionResidueType, ParameterCollectionProbabilityType> > 
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // PreAlignParameters parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<PostAlignParameters<ParameterCollectionResidueType, ParameterCollectionProbabilityType> > 
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // PostAlignParameters parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<GlobalParameters<ParameterCollectionResidueType, ParameterCollectionProbabilityType> > 
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // GlobalParameters parameter_collection_traits

  template <typename ParameterCollectionResidueType, typename ParameterCollectionProbabilityType>
  struct parameter_collection_traits<ProfilePosition<ParameterCollectionResidueType, ParameterCollectionProbabilityType> > 
  {
    typedef ParameterCollectionResidueType ResidueType;
    typedef ParameterCollectionProbabilityType ProbabilityType;
  }; // ProfilePosition parameter_collection_traits

  template <typename UnusedProfileType>
  struct profile_traits
  {
  }; // generic profile_traits

  template <typename ProfileResidueType, typename ProfileProbabilityType>
  struct profile_traits<ProfileTreeRoot<ProfileResidueType, ProfileProbabilityType> > 
  {
    typedef ProfileResidueType ResidueType;
    typedef ProfileProbabilityType ProbabilityType;
  }; // ProfileTreeRoot profile_traits

  template <typename ProfileResidueType, typename ProfileProbabilityType>
  struct profile_traits<ProfileTreeInternalNode<ProfileResidueType, ProfileProbabilityType> >
  {
    typedef ProfileResidueType ResidueType;
    typedef ProfileProbabilityType ProbabilityType;
  }; // ProfileTreeInternalNode profile_traits

  //======//// Generic Extensions ////========//

  template <typename ScoreType,
            typename ParameterCollectionType>
  class ScalableParameterCollection :
    public ParameterCollectionType
  {
  public:
    typedef typename parameter_collection_traits<ParameterCollectionType >::ProbabilityType ProbabilityType;

    /**
     * Sometimes parameter values get tiny (eg when multiplied repeatedly).
     * Since those values can be very small, there is a real danger of
     * numerical underflow.  To avoid this, all parameter collections that are
     * being added together can be scaled by the same value (the operators
     * won't account for different scales).  This is that scalar; its type
     * should be able to handle very small values, eg. the ScoreType that is
     * used by the DynamicProgramming class.  Divide by this to return to the
     * true parameter values.
     */
    ScoreType m_scalar;
  
    ScalableParameterCollection () :
      ParameterCollectionType(),
      m_scalar( 1.0 )
    {
      // Do nothing else.
    } // <init>()
  
    template <typename AnyParameterCollectionType>
    ScalableParameterCollection ( ScalableParameterCollection<ScoreType, AnyParameterCollectionType> const & copy_from ) :
      ParameterCollectionType( copy_from ),
      m_scalar( copy_from.m_scalar )
    {
      // Do nothing else.
    } // <init>( ScalableParameterCollection<ScoreType, AnyParameterCollectionType> const& )

    template <typename other_type>
    ScalableParameterCollection ( other_type const & other_obj ) :
      ParameterCollectionType(),
      m_scalar( 1.0 )
    {
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::<init>(other_type const &)" << endl;

      *this = other_obj;

      // This should be redundant.
      if( has_m_scalar<other_type>() ) {
        m_scalar = getScalar<ScoreType, other_type>( other_obj );
      }
    } // <init>( other_type const & )
  
    void
    reinitialize ()
    {
      ParameterCollectionType::reinitialize();
  
      m_scalar = 1.0;
    } // reinitialize()

    ParameterCollectionType
    createNonscalableCopy () const
    {
      return *this;
    } // createNonscalableCopy() const
  
    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Also sets m_scalar
     * to 1.
     */
    template <typename other_type>
    void
    normalize ( other_type const& min )
    {
      ParameterCollectionType::normalize( static_cast<ProbabilityType>( min ) );

      m_scalar = 1.0;
    } // normalize( ParameterValueType const& )
  
    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.  Also sets m_scalar
     * to 1.
     */
    void
    normalize ( ProbabilityType const& min )
    {
      ParameterCollectionType::normalize( min );

      m_scalar = 1.0;
    } // normalize( ProbabilityType const& )
  
    /**
     * Set all values to 0, and the scalar to 1.
     */
    void
    zero ()
    {
      ParameterCollectionType::zero();

      m_scalar = 1.0;
    } // zero()

    /**
     * Adjust each distribution's values such that the largest value is 1.0,
     * utilizing the m_scalar value.
     */
    void
    rescale ()
    {
      rescale(
        static_cast<ProbabilityType>( 1.0 )
      );
    } // rescale()
  
    /**
     * Adjust each distribution's values such that the largest value is the
     * given maximum, utilizing the m_scalar value.
     */
    virtual
    void
    rescale (
      ProbabilityType const & max
    )
    {

      // TODO: PUT BACK.  TESTING.
      //ProbabilityType largest_value( ParameterCollectionType::maximumValue() );
      ScoreType largest_value( ParameterCollectionType::maximumValue() );

      // TODO: REMOVE
      //cout << "in rescale(): largest_value is " << largest_value << endl;

      // TODO: REMOVE
      // This can be caused by underflow in the MatrixValueType in
      // dp calcs such as calculateAlignmentProfilePosition() in
      // DynamicProgramming.hpp.  A workaround is to switch to a type that can
      // represent a larger range, like the bfloat.
      if( isinf( largest_value ) ) {
        // TODO: REMOVE
        std::cerr << "WARNING: Underflow detected (maximumValue() returns infinity).  You should try setting the MatrixValueType to bfloat, which is the type with the largest range available." << endl;
        std::cerr << "\t *this is " << *this << endl;
        assert( !isinf( largest_value ) );
        return;
      }
      // TODO: REMOVE
      // This can be caused by underflow in the MatrixValueType in
      // dp calcs such as calculateAlignmentProfilePosition() in
      // DynamicProgramming.hpp.  A workaround is to switch to a type that can
      // represent a larger range, like the bfloat.
      if( largest_value <= 0.0 ) {
        // TODO: REMOVE
        std::cerr << "WARNING: Underflow detected (maximumValue() returns 0).  You should try setting the MatrixValueType to bfloat, which is the type with the largest range available." << endl;
        std::cerr << "\t *this is " << *this << endl;
      ////// TODO: PUT BACK!!!!!!!
        //assert( largest_value > 0.0 );
        return;
      }

      //ProbabilityType scale_factor( max );
      ScoreType scale_factor( max );
      scale_factor /= largest_value;

      // TODO: REMOVE
      // This can be caused by underflow in the MatrixValueType in
      // dp calcs such as calculateAlignmentProfilePosition() in
      // DynamicProgramming.hpp.  A workaround is to switch to a type that can
      // represent a larger range, like the bfloat.
      if( isinf( scale_factor ) ) {
        std::cerr << "WARNING: Underflow detected: ( " << max << " / " << largest_value << " ) returns infinity.  You should try setting the MatrixValueType to bfloat, which is the type with the largest range available." << endl;
        std::cerr << "\t *this is " << *this << endl;
        assert( !isinf( scale_factor ) );
        return;
      } // End if isinf( scale_factor )

      // TODO: REMOVE
      //cout << "in rescale(): scale_factor is " << scale_factor << endl;

      // TODO: REMOVE
      //cout << "in rescale(): m_scalar (before) is " << m_scalar << endl;

      m_scalar *= scale_factor;
      // TODO: REMOVE
      //cout << "in rescale(): m_scalar (after) is " << m_scalar << endl;

      // TODO: REMOVE
      assert( !isinf( m_scalar ) );
      assert( m_scalar != 0.0 );

      *this *= scale_factor;
    } // rescale( ProbabilityType const & )
  
    /**
     * Divide each value by m_scalar and set m_scalar to 1.0.
     */
    void
    unscale ()
    {
      if( m_scalar == 1.0 ) {
        return;
      }
      // TODO: REMOVE
      //cout << "in unscale(); about to call *this /= m_scalar.  BEFORE: " << *this << endl;
      *this /= m_scalar;
      //cout << "in unscale(); just called *this /= m_scalar.  AFTER: " << *this << endl;
      m_scalar = 1.0;
    } // unscale()

    /**
     * Convenience method; returns an unscaled (but scalable) copy.
     */
    ScalableParameterCollection<ScoreType, ParameterCollectionType>
    createUnscaledCopy () const
    {
      ScalableParameterCollection<ScoreType, ParameterCollectionType> unscaled_copy( *this );
      assert( unscaled_copy.m_scalar == m_scalar );
      //cout << "about to unscale unscaled_copy.  BEFORE: " << unscaled_copy << endl;
      unscaled_copy.unscale();
      // TODO: REMOVE
      //cout << "just unscaled unscaled_copy.  AFTER: " << unscaled_copy << endl;
      return unscaled_copy;
    } // createUnscaledCopy() const
  
    /**
     * Stream reader.
     */
    friend std::istream &
    operator>> (
      std::istream & is,
      ScalableParameterCollection<ScoreType, ParameterCollectionType> & spc
    )
    {
      //is >> "[ ";
      is.ignore( 2 );
      spc.readParameterCollection( is );
      if( is.peek() == ',' ) {
        //is >> ", scalar:";
        is.ignore( 9 );
        is >> spc.m_scalar;
      }
      //is >> " ]";
      is.ignore( 2 );
      return is;
    } // friend operator>> ( istream &, ScalableParameterCollection & )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      ScalableParameterCollection<ScoreType, ParameterCollectionType> const& spc
    )
    {
      os << "[ ";
      spc.writeParameterCollection( os );
      if( spc.m_scalar != 1.0 ) {
        os << ", scalar:" << spc.m_scalar;
      }
      os << " ]";
      return os;
    } // friend operator<< ( basic_ostream &, ScalableParameterCollection const& )
  
    /**
     * Copies values using ParameterCollectionType's operator=, and sets
     * m_scalar to getScalar( other_obj ).
     */
    template <typename other_type>
    ScalableParameterCollection<ScoreType, ParameterCollectionType> &
    operator= ( other_type const & other_obj )
    {
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator=(..)" << endl;

      ParameterCollectionType::operator=( other_obj );

      if( has_m_scalar<other_type>() ) {
        m_scalar = getScalar<ScoreType, other_type>( other_obj );

        //cout << "ScalableParameterCollection::operator=(..): m_scalar is now " << m_scalar << endl;
      }

      return *this;
    } // operator= ( other_type const & )

    /**
     * Multiplies values using ParameterCollectionType's operator*=, assumes
     * other value is already scaled by m_scalar if it is not scalable, or that
     * it has the same m_scalar if it is scalable.
     */
    template <typename other_type>
    ScalableParameterCollection<ScoreType, ParameterCollectionType> &
    operator*= ( other_type & other_obj )
    {
      if( has_m_scalar<other_type>() ) {
        ScoreType other_obj_m_scalar =
          getScalar<ScoreType, other_type>( other_obj );
        assert( other_obj_m_scalar == m_scalar );
      }

      ParameterCollectionType::operator*=( other_obj );

      return *this;
    } // operator*= ( other_type const & )

    /**
     * Divides values using ParameterCollectionType's operator/=, assumes
     * other value is already scaled by m_scalar if it is not scalable, or that
     * it has the same m_scalar if it is scalable.
     */
    template <typename other_type>
    ScalableParameterCollection<ScoreType, ParameterCollectionType> &
    operator/= ( other_type & other_obj )
    {
      if( has_m_scalar<other_type>() ) {
        ScoreType other_obj_m_scalar =
          getScalar<ScoreType, other_type>( other_obj );
        assert( other_obj_m_scalar == m_scalar );
      }

      // TODO: REMOVE
      //cout << *this << ": ScalableParameterCollection<..>::operator/=( " << other_obj << " )" << endl;
      ParameterCollectionType::operator/=( other_obj );
      // TODO: REMOVE
      //cout << "after: " << *this << endl;

      return *this;
    } // operator/= ( other_type const & )

    /**
     * Increments values using ParameterCollectionType's operator+=, assumes
     * other value is already scaled by m_scalar if it is not scalable.
     * If it is scalable, accounts correctly for the scale difference (if any).
     */
    template <typename other_type>
    ScalableParameterCollection<ScoreType, ParameterCollectionType> &
    operator+= ( other_type const & other_obj )
    {
      if( !has_m_scalar<other_type>() ) {
        this->ParameterCollectionType::operator+=( other_obj );
        return *this;
      }

      ScoreType other_obj_m_scalar =
        getScalar<ScoreType, other_type>( other_obj );

      // TODO: REMOVE 
      //cout << "ScalableParameterCollection::operator+=(ScalableParameterCollection)" << endl;

      if( other_obj_m_scalar == 0 ) {
        // Nothing to do here.
        return *this;
      }

      // Note this copies just the relevant portion of the other spc.
      ScalableParameterCollection<ScoreType, ParameterCollectionType> other_obj_copy( other_obj );

      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator+=(..): other_obj_copy is " << other_obj_copy << endl;
      //cout << "ScalableParameterCollection::operator+=(..): unscaled other_obj_copy is " << other_obj_copy.createUnscaledCopy() << endl;
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator+=(..): ( *this ).m_scalar is " << m_scalar << "; other_obj_copy.m_scalar is " << other_obj_copy.m_scalar << endl;

      if( m_scalar == 0.0 ) { // Special signal to copy directly.
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator+=(..): First update; copying directly." << endl;

        // If this is the first update, then just copy it directly (including its scalar).
        *this = other_obj_copy;
      } else if( m_scalar == other_obj_copy.m_scalar ) {
        // Scales are the same.
        this->ParameterCollectionType::operator+=( other_obj_copy );
      } else {
        // Scales are different; we need to account for that difference..

        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator+=(..): ( *this ) (before update) is " << ( *this ) << "; m_scalar is " << ( *this ).m_scalar << endl;
        // First account for the different scales
        ScoreType scale_ratio_st = this->m_scalar;
        scale_ratio_st /= other_obj_copy.m_scalar;
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator+=(..): scale_ratio is " << toDouble( scale_ratio_st ) << " = ( " << ( *this ).m_scalar << " / " << other_obj_copy.m_scalar << " )" << endl;
        other_obj_copy *= scale_ratio_st;
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator+=(..): other_obj_copy after multiplying by scale_ratio is " << other_obj_copy << endl;
        this->ParameterCollectionType::operator+=( other_obj_copy );
  
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator+=(..): ( *this ) after adding other_obj_copy is " << ( *this ) << endl;
      } // End if this is the first update .. else ..
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator+=(..): ( *this ) (after update) is " << ( *this ) << "; m_scalar is " << ( *this ).m_scalar << endl;

      this->rescale();
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator+=(..): ( *this ) after rescaling is " << ( *this ) << endl;

      return *this;
    } // operator+= ( ScalableParameterCollection<ScoreType,AnyParameterCollectionType> const & )

    /**
     * Decrements values using ParameterCollectionType's operator-=, assumes
     * other value is already scaled by m_scalar if it is not scalable.
     * If it is scalable, accounts correctly for the scale difference (if any).
     */
    template <typename other_type>
    ScalableParameterCollection<ScoreType, ParameterCollectionType> &
    operator-= ( other_type const & other_obj )
    {
      if( !has_m_scalar<other_type>() ) {
        // TODO: REMOVE 
        //cout << "ScalableParameterCollection::operator-=(non-scalable other type): other_type is " << typeid(other_obj).name() << endl;

        this->ParameterCollectionType::operator-=( other_obj );
        return *this;
      }

      ScoreType other_obj_m_scalar =
        getScalar<ScoreType, other_type>( other_obj );

      // TODO: REMOVE 
      //cout << "ScalableParameterCollection::operator-=(ScalableParameterCollection)" << endl;

      if( other_obj_m_scalar == 0 ) {
        // Nothing to do here.
        return *this;
      }

      // Note this copies just the relevant portion of the other spc.
      ScalableParameterCollection<ScoreType, ParameterCollectionType> other_obj_copy( other_obj );

      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator-=(..): other_obj_copy is " << other_obj_copy << endl;
      //cout << "ScalableParameterCollection::operator-=(..): unscaled other_obj_copy is " << other_obj_copy.createUnscaledCopy() << endl;
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator-=(..): ( *this ).m_scalar is " << m_scalar << "; other_obj_copy.m_scalar is " << other_obj_copy.m_scalar << endl;

      if( m_scalar == 0.0 ) { // Special signal to copy directly.
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator-=(..): First update; copying directly." << endl;

        // If this is the first update, then just copy it directly (including its scalar).
        *this = other_obj_copy;
      } else if( m_scalar == other_obj_copy.m_scalar ) {
        // Scales are the same.
        this->ParameterCollectionType::operator-=( other_obj_copy );
      } else {
        // Scales are different; we need to account for that difference..

        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator-=(..): ( *this ) (before update) is " << ( *this ) << "; m_scalar is " << ( *this ).m_scalar << endl;
        // First account for the different scales
        ScoreType scale_ratio_st = this->m_scalar;
        scale_ratio_st /= other_obj_copy.m_scalar;
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator-=(..): scale_ratio is " << toDouble( scale_ratio_st ) << " = ( " << ( *this ).m_scalar << " / " << other_obj_copy.m_scalar << " )" << endl;
        other_obj_copy *= scale_ratio_st;
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator-=(..): other_obj_copy after multiplying by scale_ratio is " << other_obj_copy << endl;
        this->ParameterCollectionType::operator-=( other_obj_copy );
  
        // TODO: REMOVE
        //cout << "ScalableParameterCollection::operator-=(..): ( *this ) after adding other_obj_copy is " << ( *this ) << endl;
      } // End if this is the first update .. else ..
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator-=(..): ( *this ) (after update) is " << ( *this ) << "; m_scalar is " << ( *this ).m_scalar << endl;

      this->rescale();
      // TODO: REMOVE
      //cout << "ScalableParameterCollection::operator-=(..): ( *this ) after rescaling is " << ( *this ) << endl;

      return *this;
    } // operator-= ( ScalableParameterCollection<ScoreType,AnyParameterCollectionType> const & )

  }; // End class ScalableParameterCollection

  template <typename ScoreType,
            typename ParameterCollectionType>
  struct has_m_scalar<ScalableParameterCollection<ScoreType, ParameterCollectionType> > : boost::true_type
  {}; // ScalableParameterCollection has_m_scalar

  template <typename ScoreType,
            typename ParameterCollectionType>
  struct parameter_collection_traits<ScalableParameterCollection<ScoreType, ParameterCollectionType> >
  {
    typedef typename parameter_collection_traits<ParameterCollectionType>::ResidueType ResidueType;
    typedef typename parameter_collection_traits<ParameterCollectionType>::ProbabilityType ProbabilityType;
  }; // ScalableParameterCollection parameter_collection_traits


  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::MatchEmissionParameters ////
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  MatchEmissionParameters<ResidueType, ProbabilityType>::
      MatchEmissionParameters () :
      m_Match_Emission_Distribution()
    {
      // Do nothing else.
    } // <init>()

  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_INIT
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    MatchEmissionParameters (
      MatchEmissionParameters<ResidueType, AnyProbabilityType> const & copy_from
    ) :
      m_Match_Emission_Distribution( copy_from.m_Match_Emission_Distribution )
    {
      // Do nothing else.
      //cout << "MEP::<init>( " << copy_from << " )" << endl;
    } // <init>( MatchEmissionParameters<ResidueType, AnyProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    MatchEmissionParameters (
      MatchEmissionParameters<ResidueType, ProbabilityType> const & copy_from
    ) :
      m_Match_Emission_Distribution( copy_from.m_Match_Emission_Distribution )
    {
      // Do nothing else.
      //cout << "MEP::<init>( " << copy_from << " )" << endl;
    } // <init>( MatchEmissionParameters<ResidueType, ProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
  reinitialize ()
    {
      m_Match_Emission_Distribution.reinitialize();
    } // reinitialize ()

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
  copyFrom ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      // TODO: REMOVE
      //cout << "MEP::copyFrom( " << other_pos << " )" << endl;
      m_Match_Emission_Distribution = other_pos.m_Match_Emission_Distribution;
    } // copyFrom( MatchEmissionParameters<ResidueType, AnyProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
    MatchEmissionParameters<ResidueType, ProbabilityType> &
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    operator= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      copyFrom( other_pos );

      return *this;
    } // operator=( MatchEmissionParameters<ResidueType, AnyProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_COPY
    MatchEmissionParameters<ResidueType, ProbabilityType> &
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    operator= ( MatchEmissionParameters<ResidueType, ProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      copyFrom( other_pos );

      return *this;
    } // operator=( MatchEmissionParameters<ResidueType, ProbabilityType>& )

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  MatchEmissionParameters<ResidueType, ProbabilityType> &
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    operator/= ( AnyProbabilityType const& denominator )
    {
      // TODO: REMOVE!
      //cout << "MEP::operator/=( " << denominator << " )" << endl;
      //cout << "MEP::operator/=( " << denominator << " ): about to do " << m_Match_Emission_Distribution << " /= " << denominator << endl;
      m_Match_Emission_Distribution /= denominator;
      //cout << "\t after, it is " << m_Match_Emission_Distribution << endl

      return *this;
    } // operator/=( AnyProbabilityType& )

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    MatchEmissionParameters<ResidueType, ProbabilityType> &
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    operator*= ( AnyProbabilityType const& scalar )
    {
      m_Match_Emission_Distribution *= scalar;

      return *this;
    } // operator*=( AnyProbabilityType& )

    /**
     * Add to each contained distribution the values in the given other
     * MatchEmissionParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    MatchEmissionParameters<ResidueType, ProbabilityType> &
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    operator+= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      m_Match_Emission_Distribution += other_pos.m_Match_Emission_Distribution;

      return *this;
    } // operator+=( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Subtract from each contained distribution the values in the given other
     * MatchEmissionParameters.  Note that this may violate the rule that the
     * probabilities are greater than 0.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    MatchEmissionParameters<ResidueType, ProbabilityType> &
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    operator-= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      m_Match_Emission_Distribution -= other_pos.m_Match_Emission_Distribution;

      return *this;
    } // operator-=( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos
    ) const
    {
      return crossEntropy( other_pos, ( MatchEmissionParameters<ResidueType, double> const * )0 );
    } // crossEntropy( MatchEmissionParameters const & ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_pos)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType,
            typename AnyMatchEmissionParametersType>
  GALOSH_INLINE_PROFILE_CROSSENTROPY
    double
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyMatchEmissionParametersType const * const weights
    ) const
  {
      double negative_cross_entropy = 0.0;
      MultinomialDistribution<ResidueType, double> tmp_dist;
      // log-transform it.
      other_pos.m_Match_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos.m_Match_Emission_Distribution << " is " << tmp_dist << "." << endl;

      if( weights != 0 ) {
        tmp_dist *= weights->m_Match_Emission_Distribution;
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos.m_Match_Emission_Distribution) is " << tmp_dist << endl;
      }
      negative_cross_entropy +=
        m_Match_Emission_Distribution.calculateExpectedValue(
          tmp_dist
        );

//        other_pos.m_Insertion_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );
//        if( weights != 0 ) {
//          tmp_dist *= weights->m_Insertion_Emission_Distribution;
//        }
//        negative_cross_entropy +=
//          m_Insertion_Emission_Distribution.calculateExpectedValue(
//            tmp_dist
//          );

      return ( 0.0 - negative_cross_entropy );
  } // crossEntropy( MatchEmissionParameters const&, AnyMatchEmissionParametersType const * ) const

    /**
     * Calculate and return the Euclidean distance between this set of emission
     * parameters and another set (treating every probability as an orthogonal
     * dimension).
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    euclideanDistance (
      MatchEmissionParameters const& other_pos
    ) const {
      return sqrt( euclideanDistanceSquared( other_pos ) );
    } // euclideanDistance( MatchEmissionParameters const& other_pos )

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of emission parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ
    double
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    euclideanDistanceSquared (
      MatchEmissionParameters const& other_pos
    ) const {
      double squared_euclidean_distance = 0.0;
        squared_euclidean_distance +=
          m_Match_Emission_Distribution.euclideanDistanceSquared(
            other_pos.m_Match_Emission_Distribution
          );
//        squared_euclidean_distance +=
//          m_Insertion_Emission_Distribution.euclideanDistanceSquared(
//            other_pos.m_Insertion_Emission_Distribution
//          );

      return squared_euclidean_distance;
    } // euclideanDistanceSquared( MatchEmissionParameters const& )

    /**
     * How many free parameters are there?  This is the sum of the free
     * parameters in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_FREEPARAMCOUNT
    uint32_t
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    freeParameterCount () const
    {
      static uint32_t free_params = 0;
      if( free_params > 0 ) {
        return free_params;
      }
        free_params +=
          m_Match_Emission_Distribution.freeParameterCount();

      return free_params;
    } // freeParameterCount() const

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
  void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    zero ()
    {
        m_Match_Emission_Distribution.zero();
    } // zero()

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    even ()
    {
      m_Match_Emission_Distribution.even();
    } // even()

    /**
     * Calculate the total of all contained values.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  ProbabilityType
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    total () const
    {
      return
        m_Match_Emission_Distribution.total();
    } // total() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )


    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    normalize ( ProbabilityType const & min )
    {
      m_Match_Emission_Distribution.normalize( min );
    } // normalize( ProbabilityType const & )


    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    uniform ( Random & random )
    {
      m_Match_Emission_Distribution.uniform( random );
    } // uniform( Random & )

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyMatchEmissionParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
  dirichlet (
    AnyMatchEmissionParametersType const & counts,
    Random & random
  ) {
    m_Match_Emission_Distribution.fromDirichlet(
      counts.m_Match_Emission_Distribution,
      random
    );
  } // dirichlet( AnyMatchEmissionParametersType const &, Random & )

    /**
     * Set all values such that each distrubution is distributed according to
     * the given DirichletMixture distribution.
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename DirichletMixtureMatchEmissionPriorType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    dirichletMixture (
      DirichletMixtureMatchEmissionPriorType const & me_prior,
      Random & random
    )
    {
      double u;
      uint32_t num_components;
      ProbabilityType total_so_far;
      uint32_t component_i;
      if( me_prior.size() == 1 ) {
        this->dirichlet( me_prior[ 0 ], random );
      } else {
        // First draw from the mixture distribution.
        u = random.nextUniform();
        num_components = me_prior.size();
        total_so_far = 0.0;
        for( component_i = 0; component_i < num_components; component_i++ ) {
          total_so_far += me_prior.m_mixingProbs[ component_i ];
          if( u < total_so_far ) {
            this->dirichlet( me_prior[ component_i ], random );            
            break;
          }
        }
        // This point should never be reached!
        cout << "ERROR in MatchEmissionParameters::dirichletMixture(..): m_mixingProbs don't add to 1!" << endl;
      }
    } // dirichletMixture( DirichletMixtureMatchEmissionPrior const &, Random & )

    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ProbabilityType
  MatchEmissionParameters<ResidueType, ProbabilityType>::
    maximumValue () const
    {
      return m_Match_Emission_Distribution.maximumValue();
    } // maximumValue() const

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_OSTREAM
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
  readMatchEmissionParameters (
      std::istream & is
    )
    {
      //is >> "M:";
      is.ignore( 2 );
      is >> m_Match_Emission_Distribution;
    } // readMatchEmissionParameters( istream & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
    void
  MatchEmissionParameters<ResidueType, ProbabilityType>::
  writeMatchEmissionParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      os << "M:";
      os << m_Match_Emission_Distribution;
    } // writeMatchEmissionParameters( basic_ostream & )

  ////// Class galosh::InsertionEmissionParameters ////
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_INIT
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    InsertionEmissionParameters () :
      m_Insertion_Emission_Distribution()
    {
      // Do nothing else.
    } // <init>()

  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_INIT
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    InsertionEmissionParameters (
      InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> const & copy_from
    ) :
      m_Insertion_Emission_Distribution( copy_from.m_Insertion_Emission_Distribution )
    {
      // Do nothing else.
      //cout << "IEP::<init>( " << copy_from << " )" << endl;
    } // <init>( InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> const & )

  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_REINITIALIZE
    void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
  reinitialize ()
    {
      m_Insertion_Emission_Distribution.reinitialize();
    } // reinitialize ()

  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
  operator= ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( InsertionEmissionParameters<ResidueType, AnyProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_TRIVIAL
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    operator= ( InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>& )

  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    copyFrom ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      m_Insertion_Emission_Distribution = other_pos.m_Insertion_Emission_Distribution;

      return;
    } // copyFrom( InsertionEmissionParameters<ResidueType, AnyProbabilityType>& )


    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    operator/= ( AnyProbabilityType const& denominator )
    {
      m_Insertion_Emission_Distribution /= denominator;

      return *this;
    } // operator/=( AnyProbabilityType& )

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    operator*= ( AnyProbabilityType const& scalar )
    {
      m_Insertion_Emission_Distribution *= scalar;

      return *this;
    } // operator*=( AnyProbabilityType& )

    /**
     * Add to each contained distribution the values in the given other
     * InsertionEmissionParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    operator+= ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos )
    {
      m_Insertion_Emission_Distribution += other_pos.m_Insertion_Emission_Distribution;

      return *this;
    } // operator+=( InsertionEmissionParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Subtract from each contained distribution the values in the given other
     * InsertionEmissionParameters.  Note that this may violate the rule that the
     * probabilities are greater than 0.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType> &
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    operator-= ( InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos )
    {
      m_Insertion_Emission_Distribution -= other_pos.m_Insertion_Emission_Distribution;

      return *this;
    } // operator-=( InsertionEmissionParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    crossEntropy (
      InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos
    ) const
    {
      return crossEntropy( other_pos, ( InsertionEmissionParameters<ResidueType, double> const * )0 );
    } // crossEntropy( InsertionEmissionParameters const & ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_pos)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  template <typename AnyProbabilityType,
            typename AnyInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_CROSSENTROPY
    double
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    crossEntropy (
      InsertionEmissionParameters<ResidueType, AnyProbabilityType, IsActualInsertionEmissionParametersType> const& other_pos,
      AnyInsertionEmissionParametersType const * const weights
    ) const
  {
      double negative_cross_entropy = 0.0;
      MultinomialDistribution<ResidueType, double> tmp_dist;
      // log-transform it.
      other_pos.m_Insertion_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );

      if( weights != 0 ) {
        tmp_dist *= weights->m_Insertion_Emission_Distribution;
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos.m_Insertion_Emission_Distribution) is " << tmp_dist << endl;
      }
      negative_cross_entropy +=
        m_Insertion_Emission_Distribution.calculateExpectedValue(
          tmp_dist
        );

//      other_pos.m_Insertion_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );
//      if( weights != 0 ) {
//        tmp_dist *= weights->m_Insertion_Emission_Distribution;
//      }
//      negative_cross_entropy +=
//        m_Insertion_Emission_Distribution.calculateExpectedValue(
//          tmp_dist
//        );

      // TODO: REMOVE
      //cout << "InsertionEmissionParameters::crossEntropy(..): negative_cross_entropy is " << negative_cross_entropy << endl;

      return ( 0.0 - negative_cross_entropy );
  } // crossEntropy( InsertionEmissionParameters const&, AnyInsertionEmissionParametersType const * ) const

    /**
     * Calculate and return the Euclidean distance between this set of emission
     * parameters and another set (treating every probability as an orthogonal
     * dimension).
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_TRIVIAL
    double
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    euclideanDistance (
      InsertionEmissionParameters const& other_pos
    ) const {
      return sqrt( euclideanDistanceSquared( other_pos ) );
    } // euclideanDistance( InsertionEmissionParameters const& other_pos )

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of emission parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ
    double
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    euclideanDistanceSquared (
      InsertionEmissionParameters const& other_pos
    ) const
    {
      return
        m_Insertion_Emission_Distribution.euclideanDistanceSquared(
          other_pos.m_Insertion_Emission_Distribution
        );
    } // euclideanDistanceSquared( InsertionEmissionParameters const& )


    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_FREEPARAMCOUNT
    uint32_t
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    freeParameterCount () const
    {
      static uint32_t free_params = 0;
      if( free_params > 0 ) {
        return free_params;
      }
      free_params =
        m_Insertion_Emission_Distribution.freeParameterCount();

      return free_params;
    } // freeParameterCount() const

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_REINITIALIZE
    void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    zero ()
    {
      m_Insertion_Emission_Distribution.zero();
    } // zero()

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_REINITIALIZE
    void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    even ()
    {
      m_Insertion_Emission_Distribution.even();
    } // even()

    /**
     * Calculate the total of all contained values.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  ProbabilityType
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    total () const
    {
      return
        m_Insertion_Emission_Distribution.total();
    } // total() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
  void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
  normalize ( other_type const & min )
    {
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    normalize ( ProbabilityType const & min )
    {
      m_Insertion_Emission_Distribution.normalize( min );
    } // normalize( ProbabilityType const & )


    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    uniform ( Random & random )
    {
      m_Insertion_Emission_Distribution.uniform( random );
    } // uniform( Random & )

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  template <typename AnyInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
  dirichlet (
    AnyInsertionEmissionParametersType const & counts,
    Random & random
  ) {
    m_Insertion_Emission_Distribution.fromDirichlet(
      counts.m_Insertion_Emission_Distribution,
      random
    );
  } // dirichlet( AnyInsertionEmissionParametersType const &, Random & )

    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ProbabilityType
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    maximumValue () const
    {
      return m_Insertion_Emission_Distribution.maximumValue();
    } // maximumValue() const



    /**
     * Read a comma-separated list of the parameters from the stream.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
  GALOSH_INLINE_OSTREAM
    void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
  readInsertionEmissionParameters (
      std::istream & is
    )
    {
      //is >> "I:";
      is.ignore( 2 );
      is >> m_Insertion_Emission_Distribution;
    } // readInsertionEmissionParameters( istream & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ProbabilityType, typename IsActualInsertionEmissionParametersType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  InsertionEmissionParameters<ResidueType, ProbabilityType, IsActualInsertionEmissionParametersType>::
    writeInsertionEmissionParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      os << "I:";
      os << m_Insertion_Emission_Distribution;
    } // writeInsertionEmissionParameters( basic_ostream & )

  ////// Class galosh::PositionTransitionParameters ////
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    PositionTransitionParameters () :
      m_Match_Distribution(),
      m_Insertion_Distribution(),
      m_Deletion_Distribution()
    {
      // Do nothing else.
    } // <init>()

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    PositionTransitionParameters (
      PositionTransitionParameters<ResidueType, ProbabilityType> const & copy_from
    ) :
      m_Match_Distribution( copy_from.m_Match_Distribution ),
      m_Insertion_Distribution( copy_from.m_Insertion_Distribution ),
      m_Deletion_Distribution( copy_from.m_Deletion_Distribution )
    {
      // Do nothing else.
      //cout << "PTP::<init>( " << copy_from << " )" << endl;
    } // <init>( PositionTransitionParameters<ResidueType, ProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    reinitialize ()
    {
      m_Match_Distribution.reinitialize();
      m_Insertion_Distribution.reinitialize();
      m_Deletion_Distribution.reinitialize();
    } // reinitialize ()

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  PositionTransitionParameters<ResidueType, ProbabilityType> &
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    operator= ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( PositionTransitionParameters<ResidueType, AnyProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_COPY
  PositionTransitionParameters<ResidueType, ProbabilityType> &
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    operator= ( PositionTransitionParameters<ResidueType, ProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( PositionTransitionParameters<ResidueType, ProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    copyFrom ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.

      m_Match_Distribution = other_pos.m_Match_Distribution;
      m_Insertion_Distribution = other_pos.m_Insertion_Distribution;
      m_Deletion_Distribution = other_pos.m_Deletion_Distribution;

      return;
    } // copyFrom( PositionTransitionParameters<ResidueType, AnyProbabilityType>& )

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  PositionTransitionParameters<ResidueType, ProbabilityType> &
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    operator/= ( AnyProbabilityType const& denominator )
    {
      m_Match_Distribution /= denominator;
      m_Insertion_Distribution /= denominator;
      m_Deletion_Distribution /= denominator;

      return *this;
    } // operator/=( AnyProbabilityType& )

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PositionTransitionParameters<ResidueType, ProbabilityType> &
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    operator*= ( AnyProbabilityType const& scalar )
    {
      m_Match_Distribution *= scalar;
      m_Insertion_Distribution *= scalar;
      m_Deletion_Distribution *= scalar;

      return *this;
    } // operator*=( AnyProbabilityType& )

    /**
     * Add to each contained distribution the values in the given other
     * PositionTransitionParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PositionTransitionParameters<ResidueType, ProbabilityType> &
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    operator+= ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      m_Match_Distribution += other_pos.m_Match_Distribution;
      m_Insertion_Distribution += other_pos.m_Insertion_Distribution;
      m_Deletion_Distribution += other_pos.m_Deletion_Distribution;

      return *this;
    } // operator+=( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Subtract from each contained distribution the values in the given other
     * PositionTransitionParameters.  Note that this may violate the
     * rule that the probabilities are greater than 0.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_INIT
    PositionTransitionParameters<ResidueType, ProbabilityType> &
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    operator-= ( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      m_Match_Distribution -= other_pos.m_Match_Distribution;
      m_Insertion_Distribution -= other_pos.m_Insertion_Distribution;
      m_Deletion_Distribution -= other_pos.m_Deletion_Distribution;

      return *this;
    } // operator-=( PositionTransitionParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos
    ) const
    {
      return crossEntropy( other_pos, ( PositionTransitionParameters<ResidueType, double> const * )0 );
    } // crossEntropy( PositionTransitionParameters const & ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_pos)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType,
            typename AnyPositionTransitionParametersType>
  GALOSH_INLINE_PROFILE_CROSSENTROPY
    double
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PositionTransitionParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPositionTransitionParametersType const * const weights
    ) const
  {
      double negative_cross_entropy = 0.0;

      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::MatchStateLabel, galosh::Plan7>::Type, double> tmp_match_dist;
      // log-transform it.
      other_pos[ profile_Match_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_match_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_Match_distribution_tag() ] << " is " << tmp_match_dist << "." << endl;
      if( weights != 0 ) {
        tmp_match_dist *= weights->operator[]( profile_Match_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_Match_distribution_tag() ]) is " << tmp_match_dist << endl;
      }
      negative_cross_entropy +=
        m_Match_Distribution.calculateExpectedValue(
          tmp_match_dist
        );

      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::InsertionStateLabel, galosh::Plan7>::Type, double> tmp_insertion_dist;
      // log-transform it.
      other_pos[ profile_Insertion_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_insertion_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_Insertion_distribution_tag() ] << " is " << tmp_insertion_dist << "." << endl;
      if( weights != 0 ) {
        tmp_insertion_dist *= weights->operator[]( profile_Insertion_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_Insertion_distribution_tag() ]) is " << tmp_insertion_dist << endl;
      }
      negative_cross_entropy +=
        m_Insertion_Distribution.calculateExpectedValue(
          tmp_insertion_dist
        );

      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::DeletionStateLabel, galosh::Plan7>::Type, double> tmp_deletion_dist;
      // log-transform it.
      other_pos[ profile_Deletion_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_deletion_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_Deletion_distribution_tag() ] << " is " << tmp_deletion_dist << "." << endl;
      if( weights != 0 ) {
        tmp_deletion_dist *= weights->operator[]( profile_Deletion_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_Deletion_distribution_tag() ]) is " << tmp_deletion_dist << endl;
      }
      negative_cross_entropy +=
        m_Deletion_Distribution.calculateExpectedValue(
          tmp_deletion_dist
        );

      // TODO: REMOVE
      //cout << "PositionTransitionParameters::crossEntropy(..): negative_cross_entropy is " << negative_cross_entropy << endl;

      return ( 0.0 - negative_cross_entropy );
  } // crossEntropy( PositionTransitionParameters const&, AnyPositionTransitionParametersType const * ) const

    /**
     * Calculate and return the Euclidean distance between this set of
     * transition parameters and another set (treating every probability as an
     * orthogonal dimension).
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    euclideanDistance (
      PositionTransitionParameters const& other_pos
    ) const {
      return sqrt( euclideanDistanceSquared( other_pos ) );
    } // euclideanDistance( PositionTransitionParameters const& other_pos )

    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of transition parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    // Templated this because ProfileTreeInternalNode is not a
    // PositionTransitionParameters but supports its interface.
  template <typename ResidueType, typename ProbabilityType>
    template <typename TransitionParametersType>
  GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ
    double
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    euclideanDistanceSquared (
      TransitionParametersType const& other_pos
    ) const {
      double squared_euclidean_distance = 0.0;
      squared_euclidean_distance +=
        m_Match_Distribution.euclideanDistanceSquared(
          other_pos[ profile_Match_distribution_tag() ]
        );
      squared_euclidean_distance +=
        m_Insertion_Distribution.euclideanDistanceSquared(
          other_pos[ profile_Insertion_distribution_tag() ]
        );
      squared_euclidean_distance +=
        m_Deletion_Distribution.euclideanDistanceSquared(
          other_pos[ profile_Deletion_distribution_tag() ]
        );

      return squared_euclidean_distance;
    } // euclideanDistanceSquared( TransitionParametersType const& )


    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_FREEPARAMCOUNT
    uint32_t
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    freeParameterCount () const
    {
      static uint32_t free_params = 0;
      if( free_params > 0 ) {
        return free_params;
      }
      free_params +=
        m_Match_Distribution.freeParameterCount();
      free_params +=
        m_Insertion_Distribution.freeParameterCount();
      free_params +=
        m_Deletion_Distribution.freeParameterCount();

      return free_params;
    } // freeParameterCount() const

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    zero ()
    {
      m_Match_Distribution.zero();
      m_Insertion_Distribution.zero();
      m_Deletion_Distribution.zero();
    } // zero()


    /**
     * Set all values such that each distrubution is evenly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    even ()
    {
      m_Match_Distribution.even();
      m_Insertion_Distribution.even();
      m_Deletion_Distribution.even();
    } // even()

    /**
     * Calculate the total of all contained values.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  ProbabilityType
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    total () const
    {
      ProbabilityType t =
        m_Match_Distribution.total();
      t +=
        m_Insertion_Distribution.total();
      t +=
        m_Deletion_Distribution.total();
      return t;
    } // total() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    normalize ( ProbabilityType const & min )
    {
      m_Match_Distribution.normalize( min );
      m_Insertion_Distribution.normalize( min );
      m_Deletion_Distribution.normalize( min );
    } // normalize( ProbabilityType const & )

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    uniform ( Random & random )
    {
      m_Match_Distribution.uniform( random );
      m_Insertion_Distribution.uniform( random );
      m_Deletion_Distribution.uniform( random );
    } // uniform( Random & )

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyPositionTransitionParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
  dirichlet (
    AnyPositionTransitionParametersType const & counts,
    Random & random
  ) {
    m_Match_Distribution.fromDirichlet(
      counts.m_Match_Distribution,
      random
    );
    m_Insertion_Distribution.fromDirichlet(
      counts.m_Insertion_Distribution,
      random
    );
    m_Deletion_Distribution.fromDirichlet(
      counts.m_Deletion_Distribution,
      random
    );
  } // dirichlet( AnyPositionTransitionParametersType const &, Random & )

    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ProbabilityType
  PositionTransitionParameters<ResidueType, ProbabilityType>::
    maximumValue () const
    {
      ProbabilityType largest_value;
      ProbabilityType largest_value_tmp;
      largest_value =
        m_Match_Distribution.maximumValue();
      largest_value_tmp =
        m_Insertion_Distribution.maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
      largest_value_tmp =
        m_Deletion_Distribution.maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
      return largest_value;
    } // maximumValue() const

    /**
     * Read a comma-separated list of the parameters from the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_OSTREAM
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
  readPositionTransitionParameters (
      std::istream & is
    )
    {
#ifdef NDEBUG
      static const bool do_extra_debugging = false;
#else
      static const bool do_extra_debugging = true;
#endif // NDEBUG .. else ..

      //is >> "M->";
      if( do_extra_debugging ) {
        assert( is.get() == 'M' );
        assert( is.get() == '-' );
        assert( is.get() == '>' );
      } else {
        is.ignore( 3 );
      }
      // TODO: REMOVE
      //cout << "READING insertion distribution." << endl;
      is >> m_Match_Distribution;
      // TODO: REMOVE
      //cout << "GOT " << m_Match_Distribution << endl;
      //is >> ", ";
      if( do_extra_debugging ) {
        // TODO: REMOVE
        if( is.peek() != ',' ) {
          cout << "Expecting ',' but getting '" << (char)is.peek() << "'" << endl;
        }
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      //is >> "I->";
      if( do_extra_debugging ) {
        assert( is.get() == 'I' );
        assert( is.get() == '-' );
        assert( is.get() == '>' );
      } else {
        is.ignore( 3 );
      }
      // TODO: REMOVE
      //cout << "READING insertion distribution." << endl;
      is >> m_Insertion_Distribution;
      // TODO: REMOVE
      //cout << "GOT " << m_Insertion_Distribution << endl;
      //is >> ", ";
      if( do_extra_debugging ) {
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      //is >> "D->";
      if( do_extra_debugging ) {
        assert( is.get() == 'D' );
        assert( is.get() == '-' );
        assert( is.get() == '>' );
      } else {
        is.ignore( 3 );
      }
      // TODO: REMOVE
      //cout << "READING deletion distribution." << endl;
      is >> m_Deletion_Distribution;
      // TODO: REMOVE
      //cout << "GOT " << m_Deletion_Distribution << endl;
    } // readPositionTransitionParameters( istream & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
    void
  PositionTransitionParameters<ResidueType, ProbabilityType>::
  writePositionTransitionParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      os << "M->";
      os << m_Match_Distribution;
      os << ", ";
      os << "I->";
      os << m_Insertion_Distribution;
      os << ", ";
      os << "D->";
      os << m_Deletion_Distribution;
    } // writePositionTransitionParameters( basic_ostream & )

  ////// Class galosh::PositionSpecificParameters ////
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    PositionSpecificParameters () :
      MatchEmissionParameters<ResidueType, ProbabilityType>()
    {
      // Do nothing else.
    } // <init>()

  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_INIT
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    PositionSpecificParameters (
      PositionSpecificParameters<ResidueType, AnyProbabilityType> const & copy_from
    ) :
      MatchEmissionParameters<ResidueType, ProbabilityType>( copy_from )
    {
      // Do nothing else.
    } // <init>( PositionSpecificParameters<ResidueType, ProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    PositionSpecificParameters (
      PositionSpecificParameters<ResidueType, ProbabilityType> const & copy_from
    ) :
      MatchEmissionParameters<ResidueType, ProbabilityType>( copy_from )
    {
      // Do nothing else.
    } // <init>( PositionSpecificParameters<ResidueType, ProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    copyFrom ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      MatchEmissionParameters<ResidueType, ProbabilityType>::copyFrom( other_pos );
    } // copyFrom( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
  PositionSpecificParameters<ResidueType, ProbabilityType> &
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    operator= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
  PositionSpecificParameters<ResidueType, ProbabilityType> &
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    operator= ( PositionSpecificParameters<ResidueType, ProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( PositionSpecificParamters<ResidueType, ProbabilityType> const& )

    /**
     * Add to each contained distribution the values in the given other
     * ProfilePosition.  Note that this violates the rule that the
     * probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  PositionSpecificParameters<ResidueType, ProbabilityType> &
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    operator+= ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      MatchEmissionParameters<ResidueType, ProbabilityType>::operator+=( other_pos );

      return *this;
    } // operator+=( PositionSpecificParameters const& )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos
    ) const
    {
      return crossEntropy( other_pos, ( PositionSpecificParameters<ResidueType, double> const * )0 );
    } // crossEntropy( PositionSpecificParameters const & ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_pos)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType,
            typename AnyPositionSpecificParametersType>
  GALOSH_INLINE_PROFILE_CROSSENTROPY
    double
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPositionSpecificParametersType const * const weights
    ) const
  {
    return
      MatchEmissionParameters<ResidueType, ProbabilityType>::crossEntropy( other_pos, weights );
  } // crossEntropy( PositionSpecificParameters const&, AnyPositionSpecificParametersType const * ) const

    /**
     * Calculate and return the Euclidean distance between this profile
     * position and another profile position (treating every probability as an
     * orthogonal dimension).
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    euclideanDistance (
      PositionSpecificParameters const& other_pos
    ) const {
      return sqrt( euclideanDistanceSquared( other_pos ) );
    } // euclideanDistance( PositionSpecificParameters const& other_pos )

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile position and another profile position (treating every
     * probability as an orthogonal dimension).
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ
    double
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    euclideanDistanceSquared (
      PositionSpecificParameters const& other_pos
    ) const {
      double squared_euclidean_distance = 0.0;

      squared_euclidean_distance =
        MatchEmissionParameters<ResidueType, ProbabilityType>::euclideanDistanceSquared( other_pos );

      return squared_euclidean_distance;
    } // euclideanDistanceSquared( PositionSpecificParameters const & )

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_FREEPARAMCOUNT
    uint32_t
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    freeParameterCount () const
    {
      static uint32_t free_params = 0;
      if( free_params > 0 ) {
        return free_params;
      }
      free_params +=
        MatchEmissionParameters<ResidueType, ProbabilityType>::freeParameterCount();

      return free_params;
    } // freeParameterCount() const


    /**
     * Read a comma-separated list of the parameters from the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_OSTREAM
    void
  PositionSpecificParameters<ResidueType, ProbabilityType>::
  readPositionSpecificParameters (
      std::istream & is
    )
    {
      this->readMatchEmissionParameters( is );
    } // readPositionSpecificParameters( istream & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
    template<class CharT, class Traits>
  GALOSH_INLINE_INIT
  void
  PositionSpecificParameters<ResidueType, ProbabilityType>::
    writePositionSpecificParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      this->writeMatchEmissionParameters( os );
    } // writePositionSpecificParameters( basic_ostream & )

  ////// Class galosh::PreAlignParameters ////
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PreAlignParameters<ResidueType, ProbabilityType>::
    PreAlignParameters () :
      m_PreAlign_Distribution(),
      m_Begin_Distribution()
#ifdef USE_DEL_IN_DEL_OUT
      ,
      m_DeletionIn_Distribution()
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      ,
      m_PreAlign_Emission_Distribution()
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    {
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution.zero();
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS
      // Do nothing else.
    } // <init>()

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PreAlignParameters<ResidueType, ProbabilityType>::
    PreAlignParameters (
      PreAlignParameters<ResidueType, ProbabilityType> const & copy_from
    ) :
      m_PreAlign_Distribution( copy_from.m_PreAlign_Distribution ),
      m_Begin_Distribution( copy_from.m_Begin_Distribution )
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      ,
      m_DeletionIn_Distribution()
#else
      ,
      m_DeletionIn_Distribution( copy_from.m_DeletionIn_Distribution )
#endif // USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      ,
      m_PreAlign_Emission_Distribution( copy_from.m_PreAlign_Emission_Distribution )
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    {
      // Do nothing else.
      //cout << "PreAP::<init>( " << copy_from << " )" << endl;
    } // <init>( PreAlignParameters<ResidueType, ProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
    reinitialize ()
    {
      m_PreAlign_Distribution.reinitialize();
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution.zero();
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution.reinitialize();
#ifdef USE_DEL_IN_DEL_OUT
      m_DeletionIn_Distribution.reinitialize();
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution.reinitialize();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // reinitialize ()

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
    PreAlignParameters<ResidueType, ProbabilityType> &
  PreAlignParameters<ResidueType, ProbabilityType>::
    operator= ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( PreAlignParameters<ResidueType, AnyProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_COPY
    PreAlignParameters<ResidueType, ProbabilityType> &
  PreAlignParameters<ResidueType, ProbabilityType>::
    operator= ( PreAlignParameters<ResidueType, ProbabilityType> const& other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( PreAlignParameters<ResidueType, ProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  PreAlignParameters<ResidueType, ProbabilityType>::
  copyFrom ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      // Copy the values of the other position into this one.
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution.zero();
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#else
      m_PreAlign_Distribution =
        other_node[ Transition::fromPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution =
        other_node[ Transition::fromBegin ];
#if USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionIn_Distribution = 1;
#else
      m_DeletionIn_Distribution =
        other_node[ Transition::fromDeletionIn ];
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution =
        other_node.m_PreAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return;
    } // copyFrom( PreAlignParameters<ResidueType, AnyProbabilityType>& )


    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  PreAlignParameters<ResidueType, ProbabilityType> &
  PreAlignParameters<ResidueType, ProbabilityType>::
    operator/= ( AnyProbabilityType const& denominator )
    {
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution /= denominator;
#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution /= denominator;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      m_DeletionIn_Distribution /= denominator;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution /= denominator;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator/=( AnyProbabilityType& )

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PreAlignParameters<ResidueType, ProbabilityType> &
  PreAlignParameters<ResidueType, ProbabilityType>::
    operator*= ( AnyProbabilityType const& scalar )
    {
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution *= scalar;
#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution *= scalar;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      m_DeletionIn_Distribution *= scalar;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution *= scalar;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator*=( AnyProbabilityType& )

    /**
     * Add to each contained distribution the values in the given other
     * PreAlignParameters.  Note that this may violate the rule
     * that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PreAlignParameters<ResidueType, ProbabilityType> &
  PreAlignParameters<ResidueType, ProbabilityType>::
    operator+= ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution += other_node[ Transition::fromPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution += other_node[ Transition::fromBegin ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      m_DeletionIn_Distribution += other_node[ Transition::fromDeletionIn ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution += other_node.m_PreAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator+=( PreAlignParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Subtract from each contained distribution the values in the given other
     * PreAlignParameters.  Note that this may violate the
     * rule that the probabilities are greater than 0.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PreAlignParameters<ResidueType, ProbabilityType> &
  PreAlignParameters<ResidueType, ProbabilityType>::
    operator-= ( PreAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution -= other_node[ Transition::fromPreAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution -= other_node[ Transition::fromBegin ];
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      m_DeletionIn_Distribution -= other_node[ Transition::fromDeletionIn ];
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution -= other_node.m_PreAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator-=( PreAlignParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PreAlignParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PreAlignParameters<ResidueType, AnyProbabilityType> const& other_pos
    ) const
    {
      return crossEntropy( other_pos, ( PreAlignParameters<ResidueType, double> const * )0 );
    } // crossEntropy( PreAlignParameters const & ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_pos)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType,
            typename AnyPreAlignParametersType>
  GALOSH_INLINE_PROFILE_CROSSENTROPY
    double
  PreAlignParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PreAlignParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPreAlignParametersType const * const weights
    ) const
  {
      double negative_cross_entropy = 0.0;

#ifndef DISALLOW_FLANKING_TRANSITIONS
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::PreAlignStateLabel, galosh::Plan7>::Type, double> tmp_prealign_dist;
      // log-transform it.
      other_pos[ profile_PreAlign_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_prealign_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_PreAlign_distribution_tag() ] << " is " << tmp_prealign_dist << "." << endl;
      if( weights != 0 ) {
        tmp_prealign_dist *= weights->operator[]( profile_PreAlign_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_PreAlign_distribution_tag() ]) is " << tmp_prealign_dist << endl;
      }
      negative_cross_entropy +=
        m_PreAlign_Distribution.calculateExpectedValue(
          tmp_prealign_dist
        );
#endif // !DISALLOW_FLANKING_TRANSITIONS

      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::BeginStateLabel, galosh::Plan7>::Type, double> tmp_begin_dist;
      // log-transform it.
      other_pos[ profile_Begin_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_begin_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_Begin_distribution_tag() ] << " is " << tmp_begin_dist << "." << endl;
      if( weights != 0 ) {
        tmp_begin_dist *= weights->operator[]( profile_Begin_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_Begin_distribution_tag() ]) is " << tmp_begin_dist << endl;
      }
      negative_cross_entropy +=
        m_Begin_Distribution.calculateExpectedValue(
          tmp_begin_dist
        );

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::DeletionInStateLabel, galosh::Plan7>::Type, double> tmp_deletionIn_dist;
      // log-transform it.
      other_pos[ profile_DeletionIn_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_deletionIn_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_DeletionIn_distribution_tag() ] << " is " << tmp_deletionIn_dist << "." << endl;
      if( weights != 0 ) {
        tmp_deletionIn_dist *= weights->operator[]( profile_DeletionIn_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_DeletionIn_distribution_tag() ]) is " << tmp_deletionIn_dist << endl;
      }
      negative_cross_entropy +=
        m_DeletionIn_Distribution.calculateExpectedValue(
          tmp_deletionIn_dist
        );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      MultinomialDistribution<ResidueType, double> tmp_dist;
      // log-transform it.
      other_pos.m_PreAlign_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos.m_PreAlign_Emission_Distribution << " is " << tmp_dist << "." << endl;

      if( weights != 0 ) {
        tmp_dist *= weights->m_PreAlign_Emission_Distribution;
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos.m_PreAlign_Emission_Distribution) is " << tmp_dist << endl;
      }
      negative_cross_entropy +=
        m_PreAlign_Emission_Distribution.calculateExpectedValue(
          tmp_dist
        );
//        other_pos.m_PreAlign_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );
//        if( weights != 0 ) {
//          tmp_dist *= weights->m_PreAlign_Emission_Distribution;
//        }
//        negative_cross_entropy +=
//          m_PreAlign_Emission_Distribution.calculateExpectedValue(
//            tmp_dist
//          );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      // TODO: REMOVE
      //cout << "PreAlignParameters::crossEntropy(..): negative_cross_entropy is " << negative_cross_entropy << endl;

      return ( 0.0 - negative_cross_entropy );
  } // crossEntropy( PreAlignParameters const&, AnyPreAlignParametersType const * ) const

    /**
     * Calculate and return the Euclidean distance between this set of
     * transition parameters and another set (treating every probability as an
     * orthogonal dimension).
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PreAlignParameters<ResidueType, ProbabilityType>::
    euclideanDistance (
      PreAlignParameters const& other_pos
    ) const {
      return sqrt( euclideanDistanceSquared( other_pos ) );
    } // euclideanDistance( PreAlignParameters const& other_pos )


    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of transition parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    // Templated this because ProfileTreeInternalNode is not a
    // PreAlignParameters but supports its interface.
  template <typename ResidueType, typename ProbabilityType>
    template <typename PreAlignParametersType>
  GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ
  double
  PreAlignParameters<ResidueType, ProbabilityType>::
  euclideanDistanceSquared (
      PreAlignParametersType const& other_node
    ) const {
      double squared_euclidean_distance = 0.0;
#ifndef DISALLOW_FLANKING_TRANSITIONS
      squared_euclidean_distance +=
        m_PreAlign_Distribution.euclideanDistanceSquared(
          other_node[ Transition::fromPreAlign ]
        );
#endif // !DISALLOW_FLANKING_TRANSITIONS
      squared_euclidean_distance +=
        m_Begin_Distribution.euclideanDistanceSquared(
          other_node[ Transition::fromBegin ]
        );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      squared_euclidean_distance +=
        m_DeletionIn_Distribution.euclideanDistanceSquared(
          other_node[ Transition::fromDeletionIn ]
        );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      squared_euclidean_distance +=
        m_PreAlign_Emission_Distribution.euclideanDistanceSquared(
          other_node.m_PreAlign_Emission_Distribution
        );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return squared_euclidean_distance;
    } // euclideanDistanceSquared( PreAlignParametersType const& )

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_FREEPARAMCOUNT
    uint32_t
  PreAlignParameters<ResidueType, ProbabilityType>::
    freeParameterCount () const
    {
      static uint32_t free_params = 0;
      if( free_params > 0 ) {
        return free_params;
      }
#ifndef DISALLOW_FLANKING_TRANSITIONS
      free_params +=
        m_PreAlign_Distribution.freeParameterCount();
#endif // !DISALLOW_FLANKING_TRANSITIONS
      free_params +=
        m_Begin_Distribution.freeParameterCount();
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      free_params +=
        m_DeletionIn_Distribution.freeParameterCount();
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      free_params +=
        m_PreAlign_Emission_Distribution.freeParameterCount();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return free_params;
    } // freeParameterCount() const


    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
    zero ()
    {
      m_PreAlign_Distribution.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution.zero();
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionIn_Distribution = 1; // This is a hack which makes m_DeletionOut_Distribution not actually a distribution.. but doing this effectively implements USE_SWENTRY_SWENTRY by hijacking the existing USE_DEL_IN_DEL_OUT code.
#else
      m_DeletionIn_Distribution.zero();
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution.zero();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // zero()

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
    even ()
    {
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution.zero();
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#else
      m_PreAlign_Distribution.even();
#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution.even();
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionIn_Distribution = 1;
#else
      m_DeletionIn_Distribution.even();
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution.even();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // even()

    /**
     * Calculate the total of all contained values.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  ProbabilityType
  PreAlignParameters<ResidueType, ProbabilityType>::
    total () const
    {
      ProbabilityType t( 0 );
#ifndef DISALLOW_FLANKING_TRANSITIONS
      t +=
        m_PreAlign_Distribution.total();
#endif // !DISALLOW_FLANKING_TRANSITIONS
      t +=
        m_Begin_Distribution.total();
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      t +=
        m_DeletionIn_Distribution.total();
#endif // USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      t +=
        m_PreAlign_Emission_Distribution.total();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      return t;
    } // total() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
    normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
    normalize ( ProbabilityType const & min )
    {
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution.zero();
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#else
      m_PreAlign_Distribution.normalize( min );
#endif // !DISALLOW_FLANKING_TRANSITIONS
//#ifndef DISALLOW_FLANKING_TRANSITIONS
//      m_PreAlign_Distribution.normalize( min );
//#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution.normalize( min );
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionIn_Distribution = 1;
#else
      m_DeletionIn_Distribution.normalize( min );
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution.normalize( min );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // normalize( ProbabilityType const & )

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
    uniform ( Random & random )
    {
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution.zero();
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#else
      m_PreAlign_Distribution.uniform( random );
#endif // !DISALLOW_FLANKING_TRANSITIONS
      m_Begin_Distribution.uniform( random );
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionIn_Distribution = 1;
#else
      m_DeletionIn_Distribution.uniform( random );
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PreAlign_Emission_Distribution.uniform( random );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // uniform( Random & )

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyPreAlignParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
  dirichlet (
    AnyPreAlignParametersType const & counts,
    Random & random
  ) {
#ifdef DISALLOW_FLANKING_TRANSITIONS
    m_PreAlign_Distribution.zero();
    m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#else
    m_PreAlign_Distribution.fromDirichlet(
      counts.m_PreAlign_Distribution,
      random
    );
#endif // !DISALLOW_FLANKING_TRANSITIONS
    m_Begin_Distribution.fromDirichlet(
      counts.m_Begin_Distribution,
      random
    );
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    m_DeletionIn_Distribution.fromDirichlet(
      counts.m_DeletionIn_Distribution,
      random
    );
#endif // USE_DEL_IN_DEL_OUT && !USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
    m_PreAlign_Emission_Distribution.fromDirichlet(
      counts.m_PreAlign_Emission_Distribution,
      random
    );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
  } // dirichlet( AnyPreAlignParametersType const &, Random & )

    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ProbabilityType
  PreAlignParameters<ResidueType, ProbabilityType>::
    maximumValue () const
    {
      ProbabilityType largest_value;
      ProbabilityType largest_value_tmp;

      // We're allowing the maximumValue to include the Pre- and/or Post- align transition probs, since we freak out if the maximum value is 0, and that can happen in the GlobalEntente when really the only possible transition out of the first row is C->T.  For consistency I'm allowing it at the post- end, too.
#ifdef DISALLOW_FLANKING_TRANSITIONS
      largest_value = m_Begin_Distribution.maximumValue();
#else
      largest_value =
        m_PreAlign_Distribution.maximumValue();
      largest_value_tmp = m_Begin_Distribution.maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
#endif // !DISALLOW_FLANKING_TRANSITIONS
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      largest_value_tmp = m_DeletionIn_Distribution.maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifndef DISALLOW_FLANKING_TRANSITIONS
      // See note above about allowing this even when DISALLOW_FLANKING_TRANSITIONS is defined.
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      largest_value_tmp = 
        m_PreAlign_Emission_Distribution.maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
#endif // DISALLOW_FLANKING_TRANSITIONS
      return largest_value;
    } // maximumValue() const


    /**
     * Read a comma-separated list of the parameters from the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_OSTREAM
    void
  PreAlignParameters<ResidueType, ProbabilityType>::
  readPreAlignParameters (
      std::istream & is
    )
    {
#ifdef NDEBUG
      static const bool do_extra_debugging = false;
#else
      static const bool do_extra_debugging = true;
#endif // NDEBUG .. else ..
      //is >> "N->";
      if( do_extra_debugging ) {
        assert( is.get() == 'N' );
        assert( is.get() == '-' );
        assert( is.get() == '>' );
      } else {
        is.ignore( 3 );
      }
      // TODO: REMOVE
      //cout << "READING preAlign distribution." << endl;
      is >> m_PreAlign_Distribution;
      // TODO: REMOVE
      //cout << "GOT " << m_PreAlign_Distribution << endl;
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PreAlign_Distribution.zero();
      m_PreAlign_Distribution[ TransitionFromPreAlign::toBegin ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      //is >> ", ";
      if( do_extra_debugging ) {
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      //is >> "N:";
      if( do_extra_debugging ) {
        assert( is.get() == 'N' );
        assert( is.get() == ':' );
      } else {
        is.ignore( 2 );
      }
      // TODO: REMOVE
      //cout << "READING preAlign Emission distribution." << endl;
      is >> m_PreAlign_Emission_Distribution;
      // TODO: REMOVE
      //cout << "GOT " << m_PreAlign_Emission_Distribution << endl;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      //is >> ", ";
      if( do_extra_debugging ) {
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      //is >> "B->";
      if( do_extra_debugging ) {
        assert( is.get() == 'B' );
        assert( is.get() == '-' );
        assert( is.get() == '>' );
      } else {
        is.ignore( 3 );
      }
      // TODO: REMOVE
      //cout << "READING Begin distribution." << endl;
      is >> m_Begin_Distribution;
      // TODO: REMOVE
      //cout << "GOT " << m_Begin_Distribution << endl;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      //is >> ", ";
      if( do_extra_debugging ) {
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      //is >> "Z->";
      if( do_extra_debugging ) {
        assert( is.get() == 'Z' );
        assert( is.get() == '-' );
        assert( is.get() == '>' );
      } else {
        is.ignore( 3 );
      }
      // TODO: REMOVE
      //cout << "READING DeletionIn distribution." << endl;
      is >> m_DeletionIn_Distribution;
      // TODO: REMOVE
      //cout << "GOT " << m_DeletionIn_Distribution << endl;
#endif // USE_DEL_IN_DEL_OUT && !USE_DEL_IN_DEL_OUT
    } // readPreAlignParameters( istream & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  PreAlignParameters<ResidueType, ProbabilityType>::
    writePreAlignParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      os << "N->";
      os << m_PreAlign_Distribution;
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      os << ", ";
      os << "N:";
      os << m_PreAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      os << ", ";
      os << "B->";
      os << m_Begin_Distribution;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      os << ", ";
      os << "Z->";
      os << m_DeletionIn_Distribution;
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT
    } // writePreAlignParameters( basic_ostream & )

  ////// Class galosh::PostAlignParameters ////
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PostAlignParameters<ResidueType, ProbabilityType>::
    PostAlignParameters () :
#ifdef USE_END_DISTRIBUTION
      m_End_Distribution(),
#endif // USE_END_DISTRIBUTION
      m_PostAlign_Distribution()
#ifdef USE_DEL_IN_DEL_OUT
      ,
      m_DeletionOut_Distribution()
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      ,
      m_PostAlign_Emission_Distribution()
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
#endif // USE_END_DISTRIBUTION

#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution.zero();
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS
      // Do nothing else.
    } // <init>()

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  PostAlignParameters<ResidueType, ProbabilityType>::
    PostAlignParameters (
      PostAlignParameters<ResidueType, ProbabilityType> const & copy_from
    ) :
#ifdef USE_END_DISTRIBUTION
      m_End_Distribution( copy_from.m_End_Distribution ),
#endif // USE_END_DISTRIBUTION
      m_PostAlign_Distribution( copy_from.m_PostAlign_Distribution )
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      ,
      m_DeletionOut_Distribution()
#else
      ,
      m_DeletionOut_Distribution(  copy_from.m_DeletionOut_Distribution )
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      ,
      m_PostAlign_Emission_Distribution( copy_from.m_PostAlign_Emission_Distribution )
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
#endif // USE_END_DISTRIBUTION
      // Do nothing else.
      //cout << "PostAP::<init>( " << copy_from << " )" << endl;
    } // <init>( PostAlignParameters<ResidueType, ProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
    reinitialize ()
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution.reinitialize();
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
#endif // USE_END_DISTRIBUTION
      m_PostAlign_Distribution.reinitialize();
#ifdef USE_DEL_IN_DEL_OUT
      m_DeletionOut_Distribution.reinitialize();
#endif // USE_DEL_IN_DEL_OUT
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution.zero();
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution.reinitialize();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // reinitialize ()

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
    PostAlignParameters<ResidueType, ProbabilityType> &
  PostAlignParameters<ResidueType, ProbabilityType>::
    operator= ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( PostAlignParameters<ResidueType, AnyProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_COPY
    PostAlignParameters<ResidueType, ProbabilityType> &
  PostAlignParameters<ResidueType, ProbabilityType>::
    operator= ( PostAlignParameters<ResidueType, ProbabilityType> const& other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( PostAlignParameters<ResidueType, ProbabilityType>& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  PostAlignParameters<ResidueType, ProbabilityType>::
    copyFrom ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      // Copy the values of the other position into this one.
      // TODO: REMOVE
      //cout << "PostAP::copyFrom( " << other_node << " )" << endl;
      // Copy the values of the other position into this one.

      // Note that for now we don't use the End distribution, but I'll copy it
      // anyway.
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
      //m_End_Distribution =
      //  other_node[ Transition::fromEnd ];
#endif // USE_END_DISTRIBUTION
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution.zero();
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#else
      m_PostAlign_Distribution =
        other_node[ Transition::fromPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS

#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution = 1; // This is a hack which makes m_DeletionOut_Distribution not actually a distribution.. but doing this effectively implements USE_SWENTRY_SWEXIT by hijacking the existing USE_DEL_IN_DEL_OUT code.
#else
      m_DeletionOut_Distribution =
        other_node[ Transition::fromDeletionOut ];
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution =
        other_node.m_PostAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return;
    } // copyFrom( PostAlignParameters<ResidueType, AnyProbabilityType>& )

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  PostAlignParameters<ResidueType, ProbabilityType> &
  PostAlignParameters<ResidueType, ProbabilityType>::
    operator/= ( AnyProbabilityType const& denominator )
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //m_End_Distribution /= denominator;
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution /= denominator;
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution /= denominator;
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution /= denominator;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator/=( AnyProbabilityType& )

    /**
     * Multiply each contained distrubution value by scalar.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PostAlignParameters<ResidueType, ProbabilityType> &
  PostAlignParameters<ResidueType, ProbabilityType>::
    operator*= ( AnyProbabilityType const& scalar )
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //m_End_Distribution *= scalar;
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution *= scalar;
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution *= scalar;
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution *= scalar;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator*=( AnyProbabilityType& )


    /**
     * Add to each contained distribution the values in the given other
     * PostAlignParameters.  Note that this violates the rule
     * that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PostAlignParameters<ResidueType, ProbabilityType> &
  PostAlignParameters<ResidueType, ProbabilityType>::
    operator+= ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //m_End_Distribution += other_node[ Transition::fromEnd ];
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution += other_node[ Transition::fromPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution += other_node[ Transition::fromDeletionOut ];
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution += other_node.m_PostAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator+=( PostAlignParameters<ResidueType, AnyProbabilityType> const& )


    /**
     * Subtract from each contained distribution the values in the given other
     * PostAlignParameters.  Note that this may violate the
     * rule that the probabilities are greater than 0.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    PostAlignParameters<ResidueType, ProbabilityType> &
  PostAlignParameters<ResidueType, ProbabilityType>::
    operator-= ( PostAlignParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //m_End_Distribution -= other_node[ Transition::fromEnd ];
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution -= other_node[ Transition::fromPostAlign ];
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution -= other_node[ Transition::fromDeletionOut ];
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution -= other_node.m_PostAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return *this;
    } // operator-=( PostAlignParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PostAlignParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PostAlignParameters<ResidueType, AnyProbabilityType> const& other_pos
    ) const
    {
      return crossEntropy( other_pos, ( PostAlignParameters<ResidueType, double> const * )0 );
    } // crossEntropy( PostAlignParameters const & ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_pos)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType,
            typename AnyPostAlignParametersType>
  GALOSH_INLINE_PROFILE_CROSSENTROPY
    double
  PostAlignParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      PostAlignParameters<ResidueType, AnyProbabilityType> const& other_pos,
      AnyPostAlignParametersType const * const weights
    ) const
  {
      double negative_cross_entropy = 0.0;

#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::EndStateLabel, galosh::Plan7>::Type, double> tmp_end_dist;
      //// log-transform it.
      //other_pos[ profile_End_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_end_dist );
      //// TODO: REMOVE
      ////cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_End_distribution_tag() ] << " is " << tmp_end_dist << "." << endl;
      //if( weights != 0 ) {
      //  tmp_end_dist *= weights->operator[]( profile_End_distribution_tag() );
      //  // TODO: REMOVE
      //  //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_End_distribution_tag() ]) is " << tmp_end_dist << endl;
      //}
      //negative_cross_entropy +=
      //  m_End_Distribution.calculateExpectedValue(
      //    tmp_end_dist
      //  );
#endif // USE_END_DISTRIBUTION

#ifndef DISALLOW_FLANKING_TRANSITIONS
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::PostAlignStateLabel, galosh::Plan7>::Type, double> tmp_postalign_dist;
      // log-transform it.
      other_pos[ profile_PostAlign_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_postalign_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_PostAlign_distribution_tag() ] << " is " << tmp_postalign_dist << "." << endl;
      if( weights != 0 ) {
        tmp_postalign_dist *= weights->operator[]( profile_PostAlign_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_PostAlign_distribution_tag() ]) is " << tmp_postalign_dist << endl;
      }
      negative_cross_entropy +=
        m_PostAlign_Distribution.calculateExpectedValue(
          tmp_postalign_dist
        );
#endif // !DISALLOW_FLANKING_TRANSITIONS
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      MultinomialDistribution<galosh::StateLabelTransitionTargets<galosh::DeletionOutStateLabel, galosh::Plan7>::Type, double> tmp_deletionOut_dist;
      // log-transform it.
      other_pos[ profile_DeletionOut_distribution_tag() ].toBoltzmannGibbs( 1.0, tmp_deletionOut_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos[ profile_DeletionOut_distribution_tag() ] << " is " << tmp_deletionOut_dist << "." << endl;
      if( weights != 0 ) {
        tmp_deletionOut_dist *= weights->operator[]( profile_DeletionOut_distribution_tag() );
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos[ profile_DeletionOut_distribution_tag() ]) is " << tmp_deletionOut_dist << endl;
      }
      negative_cross_entropy +=
        m_DeletionOut_Distribution.calculateExpectedValue(
          tmp_deletionOut_dist
        );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      MultinomialDistribution<ResidueType, double> tmp_dist;
      // log-transform it.
      other_pos.m_PostAlign_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );
      // TODO: REMOVE
      //cout << "in crossEntropy(..): the log transformation of " << other_pos.m_PostAlign_Emission_Distribution << " is " << tmp_dist << "." << endl;

      if( weights != 0 ) {
        tmp_dist *= weights->m_PostAlign_Emission_Distribution;
        // TODO: REMOVE
        //cout << "in crossEntropy(..): after applying the weights, the thing of which we'll take the expected value (the weighted log-transformed other_pos.m_PostAlign_Emission_Distribution) is " << tmp_dist << endl;
      }
      negative_cross_entropy +=
        m_PostAlign_Emission_Distribution.calculateExpectedValue(
          tmp_dist
        );

//        other_pos.m_PostAlign_Emission_Distribution.toBoltzmannGibbs( 1.0, tmp_dist );
//        if( weights != 0 ) {
//          tmp_dist *= weights->m_PostAlign_Emission_Distribution;
//        }
//        negative_cross_entropy +=
//          m_PostAlign_Emission_Distribution.calculateExpectedValue(
//            tmp_dist
//          );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      // TODO: REMOVE
      //cout << "PostAlignParameters::crossEntropy(..): negative_cross_entropy is " << negative_cross_entropy << endl;

      return ( 0.0 - negative_cross_entropy );
  } // crossEntropy( PostAlignParameters const&, AnyPostAlignParametersType const * ) const

    /**
     * Calculate and return the Euclidean distance between this set of
     * transition parameters and another set (treating every probability as an
     * orthogonal dimension).
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  PostAlignParameters<ResidueType, ProbabilityType>::
    euclideanDistance (
      PostAlignParameters const& other_pos
    ) const {
      return sqrt( euclideanDistanceSquared( other_pos ) );
    } // euclideanDistance( PostAlignParameters const& other_pos )


    /**
     * Calculate and return the square of the Euclidean distance between this
     * set of transition parameters and another set (treating every probability
     * as an orthogonal dimension).
     */
    // Templated this because ProfileTreeInternalNode is not a
    // PostAlignParameters but supports its interface.
  template <typename ResidueType, typename ProbabilityType>
    template <typename PostAlignParametersType>
  GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ
    double
  PostAlignParameters<ResidueType, ProbabilityType>::
    euclideanDistanceSquared (
      PostAlignParametersType const& other_node
    ) const {
      double squared_euclidean_distance = 0.0;
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //squared_euclidean_distance +=
      //  m_End_Distribution.euclideanDistanceSquared(
      //    other_node[ Transition::fromEnd ]
      //  );
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
      squared_euclidean_distance +=
        m_PostAlign_Distribution.euclideanDistanceSquared(
          other_node[ Transition::fromPostAlign ]
        );
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      squared_euclidean_distance +=
        m_DeletionOut_Distribution.euclideanDistanceSquared(
          other_node[ Transition::fromDeletionOut ]
        );
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      squared_euclidean_distance +=
        m_PostAlign_Emission_Distribution.euclideanDistanceSquared(
          other_node.m_PostAlign_Emission_Distribution
        );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return squared_euclidean_distance;
    } // euclideanDistanceSquared( PostAlignParametersType const& )

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_FREEPARAMCOUNT
    uint32_t
  PostAlignParameters<ResidueType, ProbabilityType>::
    freeParameterCount () const
    {
      static uint32_t free_params = 0;
      if( free_params > 0 ) {
        return free_params;
      }
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //free_params +=
      //  m_End_Distribution.freeParameterCount();
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
      free_params +=
        m_PostAlign_Distribution.freeParameterCount();
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      free_params +=
        m_DeletionOut_Distribution.freeParameterCount();
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      free_params +=
        m_PostAlign_Emission_Distribution.freeParameterCount();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS

      return free_params;
    } // freeParameterCount() const


    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
    zero ()
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
      //m_End_Distribution.zero();
#endif // USE_END_DISTRIBUTION
      m_PostAlign_Distribution.zero();
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution = 1;
#else
      m_DeletionOut_Distribution.zero();
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution.zero();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // zero()

    /**
     * Set all values such that each distrubution is evenly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
    even ()
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
      //m_End_Distribution.even();
#endif // USE_END_DISTRIBUTION
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution.zero();
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#else
      m_PostAlign_Distribution.even();
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution = 1;
#else
      m_DeletionOut_Distribution.even();
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution.even();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // even()

    /**
     * Calculate the total of all contained values.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  ProbabilityType
  PostAlignParameters<ResidueType, ProbabilityType>::
  total () const
    {
      ProbabilityType t( 0 );
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //t +=
      //  m_End_Distribution.total();
#endif // USE_END_DISTRIBUTION
#ifndef DISALLOW_FLANKING_TRANSITIONS
      t +=
        m_PostAlign_Distribution.total();
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifndef USE_SWENTRY_SWEXIT
      t +=
        m_DeletionOut_Distribution.total();
#endif // !USE_SWENTRY_SWEXIT
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      t +=
        m_PostAlign_Emission_Distribution.total();
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
      return t;
    } // total() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
    normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
    normalize ( ProbabilityType const & min )
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
      //m_End_Distribution.normalize( min );
      //m_End_Distribution.normalize(); // Don't enforce min for End
#endif // USE_END_DISTRIBUTION
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution.zero();
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#else
      m_PostAlign_Distribution.normalize( min );
#endif // !DISALLOW_FLANKING_TRANSITIONS
//#ifndef DISALLOW_FLANKING_TRANSITIONS
//      m_PostAlign_Distribution.normalize( min );
//#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution = 1;
#else
      m_DeletionOut_Distribution.normalize( min );
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution.normalize( min );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // normalize( ProbabilityType const & )

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
    uniform ( Random & random )
    {
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
      //m_End_Distribution.uniform( random );
#endif // USE_END_DISTRIBUTION
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution.zero();
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#else
      m_PostAlign_Distribution.uniform( random );
#endif // !DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_DEL_IN_DEL_OUT
#ifdef USE_SWENTRY_SWEXIT
      m_DeletionOut_Distribution = 1;
#else
      m_DeletionOut_Distribution.uniform( random );
#endif // USE_SWENTRY_SWEXIT .. else ..
#endif // USE_DEL_IN_DEL_OUT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      m_PostAlign_Emission_Distribution.uniform( random );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // uniform( Random & )

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyPostAlignParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
  dirichlet (
    AnyPostAlignParametersType const & counts,
    Random & random
  ) {
#ifdef USE_END_DISTRIBUTION
    // For now we we never allow loops, so the End distribution is constant.
    m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
    m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
    //m_End_Distribution.fromDirichlet(
    //  counts.m_End_Distribution,
    //  random
    //);
#endif // USE_END_DISTRIBUTION
#ifdef DISALLOW_FLANKING_TRANSITIONS
    m_PostAlign_Distribution.zero();
    m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#else
    m_PostAlign_Distribution.fromDirichlet(
      counts.m_PostAlign_Distribution,
      random
    );
#endif // !DISALLOW_FLANKING_TRANSITIONS
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
    m_DeletionOut_Distribution.fromDirichlet(
      counts.m_DeletionOut_Distribution,
      random
    );
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
    m_PostAlign_Emission_Distribution.fromDirichlet(
      counts.m_PostAlign_Emission_Distribution,
      random
    );
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
  } // dirichlet( AnyPostAlignParametersType const &, Random & )

    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ProbabilityType
  PostAlignParameters<ResidueType, ProbabilityType>::
    maximumValue () const
    {
      ProbabilityType largest_value;
      ProbabilityType largest_value_tmp;
#ifdef USE_END_DISTRIBUTION
      // For now we we never allow loops, so the End distribution is constant.
      //largest_value =
      //  m_End_Distribution.maximumValue();
#endif // USE_END_DISTRIBUTION
      // We're allowing the maximumValue to include the Pre- and/or Post- align transition probs, since we freak out if the maximum value is 0, and that can happen in the GlobalEntente when really the only possible transition out of the first row is C->T.  For consistency I'm allowing it at the post- end, too.
#ifdef DISALLOW_FLANKING_TRANSITIONS
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS // & DISALLOW_FLANKING_TRANSITIONS
      largest_value_tmp =
        m_PostAlign_Emission_Distribution.maximumValue();
      //if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      //}
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS ( & DISALLOW_FLANKING_TRANSITIONS )
#else // !DISALLOW_FLANKING_TRANSITIONS
      largest_value_tmp = m_PostAlign_Distribution.maximumValue();
      //if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      //}

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS // & !DISALLOW_FLANKING_TRANSITIONS
      largest_value_tmp = 
        m_PostAlign_Emission_Distribution.maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS ( & !DISALLOW_FLANKING_TRANSITIONS )
#endif // DISALLOW_FLANKING_TRANSITIONS .. else ..

#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      largest_value_tmp = m_DeletionOut_Distribution.maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
#endif // USE_DEL_IN_DEL_OUT && !USE_SWENTRY_SWEXIT

      return largest_value;
    } // maximumValue() const



    /**
     * Read a comma-separated list of the parameters from the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_OSTREAM
    void
  PostAlignParameters<ResidueType, ProbabilityType>::
  readPostAlignParameters (
      std::istream & is
    )
    {
#ifdef USE_END_DISTRIBUTION
      //is >> "E->";
      is.ignore( 3 );
      is >> m_End_Distribution;
      // For now we we never allow loops, so the End distribution is constant.
      m_End_Distribution[ TransitionFromEnd::toPostAlign ] = 1.0;
      m_End_Distribution[ TransitionFromEnd::toLoop ] = 0.0;
      //is >> ", ";
      is.ignore( 2 );
#endif // USE_END_DISTRIBUTION
      //is >> "C->";
      is.ignore( 3 );
      is >> m_PostAlign_Distribution;
#ifdef DISALLOW_FLANKING_TRANSITIONS
      m_PostAlign_Distribution.zero();
      m_PostAlign_Distribution[ TransitionFromPostAlign::toTerminal ] = 1.0;
#endif // DISALLOW_FLANKING_TRANSITIONS
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      //is >> ", ";
      is.ignore( 2 );
      //is >> "W->";
      is.ignore( 3 );
      is >> m_DeletionOut_Distribution;
#endif // USE_DEL_IN_DEL_OUT & !USE_SWENTRY_SWEXIT

#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      //is >> ", ";
      is.ignore( 2 );
      //is >> "C:";
      is.ignore( 2 );
      is >> m_PostAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // readPostAlignParameters( istream & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  PostAlignParameters<ResidueType, ProbabilityType>::
    writePostAlignParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
#ifdef USE_END_DISTRIBUTION
      os << "E->";
      os << m_End_Distribution;
      os << ", ";
#endif // USE_END_DISTRIBUTION
      os << "C->";
      os << m_PostAlign_Distribution;
#if defined( USE_DEL_IN_DEL_OUT ) && !defined( USE_SWENTRY_SWEXIT )
      os << ", ";
      os << "W->";
      os << m_DeletionOut_Distribution;
#endif // USE_DEL_IN_DEL_OUT & !USE_SWENTRY_SWEXIT
#ifdef USE_FLANKING_EMISSION_DISTRIBUTIONS
      os << ", ";
      os << "C:";
      os << m_PostAlign_Emission_Distribution;
#endif // USE_FLANKING_EMISSION_DISTRIBUTIONS
    } // writePostAlignParameters( basic_ostream & )

  ////// Class galosh::ProfilePosition ////
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  ProfilePosition<ResidueType, ProbabilityType>::
    ProfilePosition () :
      PositionSpecificParameters<ResidueType, ProbabilityType>(),
      m_root( 0 )
    {
      // Do nothing else.
    } // <init>()

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  ProfilePosition<ResidueType, ProbabilityType>::
    ProfilePosition ( ProfileTreeRoot<ResidueType, ProbabilityType> * const pt_root ) :
      PositionSpecificParameters<ResidueType, ProbabilityType>(),
      m_root( pt_root )
    {
      // Do nothing else.
    } // <init>( ProfileTreeRoot<ResidueType, ProbabilityType> * profile )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  ProfilePosition<ResidueType, ProbabilityType>::
    ProfilePosition ( ProfilePosition<ResidueType, ProbabilityType> const & copy_from ) :
      PositionSpecificParameters<ResidueType, ProbabilityType>( copy_from ),
      m_root( copy_from.m_root )
    {
      // Do nothing else.
    } // <init>( ProfilePosition<ResidueType, ProbabilityType> const & )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  ProfilePosition<ResidueType, ProbabilityType>::
    copyFrom ( ProfilePosition<ResidueType, AnyProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      PositionSpecificParameters<ResidueType, ProbabilityType>::copyFrom( other_pos );

      // TODO: Put back?  I think this is why we keep having to call ensurePositionsKnowTheirRoot(..)...
      //m_root = other_pos.m_root;
    } // copyFrom( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  ProfilePosition<ResidueType, ProbabilityType>::
    copyFrom ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      MatchEmissionParameters<ResidueType, ProbabilityType>::copyFrom( other_pos );
    } // copyFrom( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
  void
  ProfilePosition<ResidueType, ProbabilityType>::
    copyFrom ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      // Copy the values of the other position into this one.
      PositionSpecificParameters<ResidueType, ProbabilityType>::copyFrom( other_pos );
    } // copyFrom( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
  ProfilePosition<ResidueType, ProbabilityType> &
  ProfilePosition<ResidueType, ProbabilityType>::
    operator= ( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( MatchEmissionParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
  ProfilePosition<ResidueType, ProbabilityType> &
  ProfilePosition<ResidueType, ProbabilityType>::
    operator= ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
  ProfilePosition<ResidueType, ProbabilityType> &
  ProfilePosition<ResidueType, ProbabilityType>::
    operator= ( ProfilePosition<ResidueType, AnyProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( ProfilePosition<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
  ProfilePosition<ResidueType, ProbabilityType> &
  ProfilePosition<ResidueType, ProbabilityType>::
    operator= ( ProfilePosition<ResidueType, ProbabilityType> const& other_pos )
    {
      copyFrom( other_pos );

      return *this;
    } // operator=( ProfilePosition<ResidueType, ProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  ProfilePosition<ResidueType, ProbabilityType>::
    reinitialize ( ProfileTreeRoot<ResidueType, ProbabilityType> * const pt_root )
    {
      PositionSpecificParameters<ResidueType, ProbabilityType>::reinitialize();
      setProfileTreeRoot( pt_root );
    } // reinitialize ( ProfileTreeRoot<ResidueType, ProbabilityType> * const )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
    void
  ProfilePosition<ResidueType, ProbabilityType>::
    setProfileTreeRoot ( ProfileTreeRoot<ResidueType, ProbabilityType> * const pt_root )
    {
      m_root = pt_root;
    } // setProfileTreeRoot()

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_ACCESSOR
  ProfileTreeRoot<ResidueType, ProbabilityType> * 
  ProfilePosition<ResidueType, ProbabilityType>::
    getProfileTreeRoot () const
    {
      return m_root;
    } // getProfileTreeRoot() const

    /**
     * Divide each contained distrubution value by denominator.  Note that
     * this violates the rule that the probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  ProfilePosition<ResidueType, ProbabilityType> &
  ProfilePosition<ResidueType, ProbabilityType>::
    operator/= ( AnyProbabilityType const& denominator )
    {
      PositionSpecificParameters<ResidueType, ProbabilityType>::operator/=( denominator );

      return *this;
    } // operator/=( AnyProbabilityType const& )

    /**
     * Add to each contained distribution the values in the given other
     * ProfilePosition.  Note that this violates the rule that the
     * probabilities sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
  ProfilePosition<ResidueType, ProbabilityType> &
  ProfilePosition<ResidueType, ProbabilityType>::
    operator+= ( PositionSpecificParameters<ResidueType, AnyProbabilityType> const& other_pos )
    {
      PositionSpecificParameters<ResidueType, ProbabilityType>::operator+=( other_pos );

      return *this;
    } // operator+=( PositionSpecificParameters const& )

    /**
     * Set all values to 0.  Note that this violates the rule that the values
     * sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  ProfilePosition<ResidueType, ProbabilityType>::
    zero ()
    {
      MatchEmissionParameters<ResidueType, ProbabilityType>::zero();
    } // zero()


    /**
     * Set all values such that each distrubution is evenly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  ProfilePosition<ResidueType, ProbabilityType>::
    even ()
    {
      MatchEmissionParameters<ResidueType, ProbabilityType>::even();
    } // even()

    /**
     * Calculate the total of all contained values.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
  ProbabilityType
  ProfilePosition<ResidueType, ProbabilityType>::
    total () const
    {
      return
        MatchEmissionParameters<ResidueType, ProbabilityType>::total();
    } // total() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
    void
  ProfilePosition<ResidueType, ProbabilityType>::
    normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  ProfilePosition<ResidueType, ProbabilityType>::
    normalize ( ProbabilityType const & min )
    {
      MatchEmissionParameters<ResidueType, ProbabilityType>::normalize( min );
    } // normalize( ProbabilityType const & )

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  ProfilePosition<ResidueType, ProbabilityType>::
    uniform ( Random & random )
    {
      MatchEmissionParameters<ResidueType, ProbabilityType>::uniform( random );
    } // uniform( Random & )


    /**
     * Return the largest value in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ProbabilityType
  ProfilePosition<ResidueType, ProbabilityType>::
    maximumValue () const
    {
      //bool have_set_largest_value = false;
      ProbabilityType largest_value;
      //ProbabilityType largest_value_tmp;
      //largest_value =
      return
        MatchEmissionParameters<ResidueType, ProbabilityType>::maximumValue();
      //if( !have_set_largest_value ) {
      //  // TODO: ?!
      //  largest_value = 0;
      //}
      //return largest_value;
    } // maximumValue() const


  ////// Class galosh::GlobalParameters ////
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  GlobalParameters<ResidueType, ProbabilityType>::
    GlobalParameters () :
      PositionTransitionParameters<ResidueType, ProbabilityType>(),
      InsertionEmissionParameters<ResidueType, ProbabilityType>(),
      PreAlignParameters<ResidueType, ProbabilityType>(),
      PostAlignParameters<ResidueType, ProbabilityType>()
    {
      // Do nothing else
    } // <init>()

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_INIT
  GlobalParameters<ResidueType, ProbabilityType>::
    GlobalParameters (
      GlobalParameters<ResidueType, ProbabilityType> const & copy_from
    ) :
      PositionTransitionParameters<ResidueType, ProbabilityType>( copy_from ),
      InsertionEmissionParameters<ResidueType, ProbabilityType>( copy_from ),
      PreAlignParameters<ResidueType, ProbabilityType>( copy_from ),
      PostAlignParameters<ResidueType, ProbabilityType>( copy_from )
    {
      // Do nothing else.
      //cout << "GP::<init>( " << copy_from << " )" << endl;
    } // <init>( GlobalParameters<ResidueType, ProbabilityType> const & )

    /**
     * Reinitialize contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    reinitialize ()
    {
      PositionTransitionParameters<ResidueType, ProbabilityType>::reinitialize();
      InsertionEmissionParameters<ResidueType, ProbabilityType>::reinitialize();
      PreAlignParameters<ResidueType, ProbabilityType>::reinitialize();
      PostAlignParameters<ResidueType, ProbabilityType>::reinitialize();
    } // reinitialize()

  /**
   * Calculate the total of all contained values.
   */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_TRIVIAL
    ProbabilityType
  GlobalParameters<ResidueType, ProbabilityType>::
    total () const
    {
      ProbabilityType t( 0 );
      t +=
        PositionTransitionParameters<ResidueType, ProbabilityType>::total();
      t +=
        InsertionEmissionParameters<ResidueType, ProbabilityType>::total();
      t +=
        PreAlignParameters<ResidueType, ProbabilityType>::total();
      t +=
        PostAlignParameters<ResidueType, ProbabilityType>::total();
    } // total() const

    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename other_type>
  GALOSH_INLINE_TRIVIAL
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    normalize ( other_type const & min )
    {
      //normalize( min );
      normalize( static_cast<ProbabilityType>( min ) );
    } // normalize( other_type const & )


    /**
     * Adjust each distribution's values such that they sum to one, ensuring
     * that no value is less than the specified minimum.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    normalize ( ProbabilityType const& min )
    {
      PositionTransitionParameters<ResidueType, ProbabilityType>::normalize( min );
      InsertionEmissionParameters<ResidueType, ProbabilityType>::normalize( min );
      PreAlignParameters<ResidueType, ProbabilityType>::normalize( min );
      PostAlignParameters<ResidueType, ProbabilityType>::normalize( min );
    } // normalize( ProbabilityType const& )

    /**
     * Set all values to 0.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    zero ()
    {
      PositionTransitionParameters<ResidueType, ProbabilityType>::zero();
      InsertionEmissionParameters<ResidueType, ProbabilityType>::zero();
      PreAlignParameters<ResidueType, ProbabilityType>::zero();
      PostAlignParameters<ResidueType, ProbabilityType>::zero();
    } // zero()

    /**
     * Set all values such that each distrubution is uniformly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    even ()
    {
      PositionTransitionParameters<ResidueType, ProbabilityType>::even();
      InsertionEmissionParameters<ResidueType, ProbabilityType>::even();
      PreAlignParameters<ResidueType, ProbabilityType>::even();
      PostAlignParameters<ResidueType, ProbabilityType>::even();
    } // even()

    /**
     * Set all values such that each distrubution is randomly distributed.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_REINITIALIZE
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    uniform ( Random & random )
    {
      PreAlignParameters<ResidueType, ProbabilityType>::uniform( random );
      PostAlignParameters<ResidueType, ProbabilityType>::uniform( random );
      InsertionEmissionParameters<ResidueType, ProbabilityType>::uniform( random );
      PositionTransitionParameters<ResidueType, ProbabilityType>::uniform( random );
    } // uniform( Random & )

    /**
     * Change the probabilities of the contained Multinomials to values drawn
     * from dirichlet distributions with the given counts.  AnyCountType can be
     * a Real type (anything coercible to a double using toDouble( count )).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyGlobalParametersType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    dirichlet (
      AnyGlobalParametersType const & counts,
      Random & random
    )
    {
      PreAlignParameters<ResidueType, ProbabilityType>::dirichlet( counts, random );
      PostAlignParameters<ResidueType, ProbabilityType>::dirichlet( counts, random );
      InsertionEmissionParameters<ResidueType, ProbabilityType>::dirichlet( counts, random );
      PositionTransitionParameters<ResidueType, ProbabilityType>::dirichlet( counts, random );
    } // dirichlet( AnyGlobalParametersType const &, Random & )


    /**
     * Read a comma-separated list of the parameters from the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_OSTREAM
    void
  GlobalParameters<ResidueType, ProbabilityType>::
  readGlobalParameters (
      std::istream & is
    )
    {
#ifdef NDEBUG
      static const bool do_extra_debugging = false;
#else
      static const bool do_extra_debugging = true;
#endif // NDEBUG .. else ..
      this->readPositionTransitionParameters( is );
      //is >> ", ";
      if( do_extra_debugging ) {
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      // TODO: REMOVE
      //cout << "READING insertion emission parameters." << endl;
      this->readInsertionEmissionParameters( is );
      // TODO: REMOVE
      //cout << "GOT "; this->writeInsertionEmissionParameters( cout ); cout << endl;
      //is >> ", ";
      if( do_extra_debugging ) {
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      // TODO: REMOVE
      //cout << "READING preAlign parameters." << endl;
      this->readPreAlignParameters( is );
      // TODO: REMOVE
      //cout << "GOT "; this->writePreAlignParameters( cout ); cout << endl;
      //is >> ", ";
      if( do_extra_debugging ) {
        assert( is.get() == ',' );
        assert( is.get() == ' ' );
      } else {
        is.ignore( 2 );
      }
      // TODO: REMOVE
      //cout << "READING postAlign parameters." << endl;
      this->readPostAlignParameters( is );
      // TODO: REMOVE
      //cout << "GOT "; this->writePostAlignParameters( cout ); cout << endl;
    } // readGlobalParameters( istream & )

    /**
     * Write a comma-separated list of the parameters to the stream.
     */
  template <typename ResidueType, typename ProbabilityType>
    template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    writeGlobalParameters (
      std::basic_ostream<CharT,Traits>& os
    ) const
    {
      this->writePositionTransitionParameters( os );
      os << ", ";
      this->writeInsertionEmissionParameters( os );
      os << ", ";
      this->writePreAlignParameters( os );
      os << ", ";
      this->writePostAlignParameters( os );
    } // writeGlobalParameters( basic_ostream & )

    /**
     * Return the largest value stored in this ProfileTreeRoot.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR
    ProbabilityType
  GlobalParameters<ResidueType, ProbabilityType>::
    maximumValue () const
    {
      ProbabilityType largest_value;
      ProbabilityType largest_value_tmp;
      largest_value =
        PreAlignParameters<ResidueType, ProbabilityType>::maximumValue();
      largest_value_tmp =
        PostAlignParameters<ResidueType, ProbabilityType>::maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
      largest_value_tmp =
        InsertionEmissionParameters<ResidueType, ProbabilityType>::maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
      largest_value_tmp =
        PositionTransitionParameters<ResidueType, ProbabilityType>::maximumValue();
      if( largest_value_tmp > largest_value ) {
        largest_value = largest_value_tmp;
      }
      return largest_value;
    } // maximumValue() const

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
    void
  GlobalParameters<ResidueType, ProbabilityType>::
    copyFrom ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      // TODO: REMOVE
      //cout << "GlobalParameters.copyFrom(..)" << endl;
      InsertionEmissionParameters<ResidueType, ProbabilityType>::operator=( other_node );
      PreAlignParameters<ResidueType, ProbabilityType>::operator=( other_node );
      PostAlignParameters<ResidueType, ProbabilityType>::operator=( other_node );
      PositionTransitionParameters<ResidueType, ProbabilityType>::operator=( other_node );
    } // copyFrom( GlobalParameters const& )

  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_COPY
    GlobalParameters<ResidueType, ProbabilityType> &
  GlobalParameters<ResidueType, ProbabilityType>::
    operator= ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( GlobalParameters<ResidueType, AnyProbabilityType> const& )

  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_COPY
    GlobalParameters<ResidueType, ProbabilityType> &
  GlobalParameters<ResidueType, ProbabilityType>::
    operator= ( GlobalParameters<ResidueType, ProbabilityType> const& other_node )
    {
      copyFrom( other_node );

      return *this;
    } // operator=( GlobalParameters<ResidueType, ProbabilityType> const& )

    /**
     * Multiply each profile value by the given scalar.  Note that this may
     * violate the rule that the probabilities add to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    GlobalParameters<ResidueType, ProbabilityType> &
  GlobalParameters<ResidueType, ProbabilityType>::
    operator*= ( AnyProbabilityType const& scalar )
    {
      PreAlignParameters<ResidueType, ProbabilityType>::operator*=( scalar );
      PostAlignParameters<ResidueType, ProbabilityType>::operator*=( scalar );
      InsertionEmissionParameters<ResidueType, ProbabilityType>::operator*=( scalar );
      PositionTransitionParameters<ResidueType, ProbabilityType>::operator*=( scalar );

      return *this;
    } // operator*= ( AnyProbabilityType const& )

    /**
     * Divide each non-position-specific profile value by the given
     * denominator.  Note that this may violate the rule that the
     * probabilities add to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    GlobalParameters<ResidueType, ProbabilityType> &
  GlobalParameters<ResidueType, ProbabilityType>::
    operator/= ( AnyProbabilityType const& denominator )
    {
      PreAlignParameters<ResidueType, ProbabilityType>::operator/=( denominator );
      PostAlignParameters<ResidueType, ProbabilityType>::operator/=( denominator );
      InsertionEmissionParameters<ResidueType, ProbabilityType>::operator/=( denominator );
      PositionTransitionParameters<ResidueType, ProbabilityType>::operator/=( denominator );

      return *this;
    } // operator/= ( AnyProbabilityType const& )

    /**
     * Add to each contained distribution the corresponding values in the given
     * other node.  Note that this may violate the rule that the probabilities
     * sum to 1.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_ARITHMETIC
    GlobalParameters<ResidueType, ProbabilityType> &
  GlobalParameters<ResidueType, ProbabilityType>::
    operator+= ( GlobalParameters<ResidueType, AnyProbabilityType> const& other_node )
    {
      PreAlignParameters<ResidueType, ProbabilityType>::operator+=( other_node );
      PostAlignParameters<ResidueType, ProbabilityType>::operator+=( other_node );
      InsertionEmissionParameters<ResidueType, ProbabilityType>::operator+=( other_node );
      PositionTransitionParameters<ResidueType, ProbabilityType>::operator+=( other_node );

      return *this;
    } // operator+=( GlobalParameters<ResidueType, AnyProbabilityType> const& )

    /**
     * Calculate and return the cross entropy E(-log(other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType>
  GALOSH_INLINE_TRIVIAL
    double
  GlobalParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      GlobalParameters<ResidueType, AnyProbabilityType> const& other_node
    ) const
    {
      return crossEntropy( other_node, ( GlobalParameters<ResidueType, double> const * )0 );
    } // crossEntropy( GlobalParameters const & ) const

    /**
     * Calculate and return the (possibly weighted) cross entropy
     * E(-log(other_pos)).  The weights (if non-null) may be of any type
     * convertible to a double, and the cross entropy will be computed as
     * E(-log(weights*other_pos)).
     */
  template <typename ResidueType, typename ProbabilityType>
  template <typename AnyProbabilityType,
            typename AnyGlobalParametersType>
  GALOSH_INLINE_PROFILE_CROSSENTROPY
    double
  GlobalParameters<ResidueType, ProbabilityType>::
    crossEntropy (
      GlobalParameters<ResidueType, AnyProbabilityType> const& other_node,
      AnyGlobalParametersType const * const weights
    ) const
  {
      double cross_entropy = 0.0;

      cross_entropy +=
        PreAlignParameters<ResidueType, ProbabilityType>::crossEntropy( other_node, weights );
      cross_entropy +=
        PostAlignParameters<ResidueType, ProbabilityType>::crossEntropy( other_node, weights );
      cross_entropy +=
        InsertionEmissionParameters<ResidueType, ProbabilityType>::crossEntropy( other_node, weights );
      cross_entropy +=
        PositionTransitionParameters<ResidueType, ProbabilityType>::crossEntropy( other_node, weights );

      return cross_entropy;
  } // crossEntropy( GlobalParameters const&, AnyGlobalParametersType const * ) const

    /**
     * Calculate and return the square of the Euclidean distance between this
     * profile and another profile (treating every probability as an orthogonal
     * dimension).  The other profile must have the same number of positions as
     * this profile.
     */
  template <typename ResidueType, typename ProbabilityType>
    template <typename AnyProbabilityType>
  GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ
    double
  GlobalParameters<ResidueType, ProbabilityType>::
    euclideanDistanceSquared (
      GlobalParameters<ResidueType, AnyProbabilityType> const& other_node
    ) const
    {
      double squared_euclidean_distance = 0.0;

      squared_euclidean_distance +=
        PreAlignParameters<ResidueType, ProbabilityType>::euclideanDistanceSquared( other_node );
      squared_euclidean_distance +=
        PostAlignParameters<ResidueType, ProbabilityType>::euclideanDistanceSquared( other_node );
      squared_euclidean_distance +=
        InsertionEmissionParameters<ResidueType, ProbabilityType>::euclideanDistanceSquared( other_node );
      squared_euclidean_distance +=
        PositionTransitionParameters<ResidueType, ProbabilityType>::euclideanDistanceSquared( other_node );

      return squared_euclidean_distance;
    } // euclideanDistanceSquared( GlobalParameters<ResidueType, AnyProbabilityType> const& ) const

    /**
     * How many free parameters are there?  This is the sum of the free
     * paramters in the contained distributions.
     */
  template <typename ResidueType, typename ProbabilityType>
  GALOSH_INLINE_PROFILE_FREEPARAMCOUNT
    uint32_t
  GlobalParameters<ResidueType, ProbabilityType>::
    freeParameterCount () const
    {
      static uint32_t free_params = 0;
      if( free_params > 0 ) {
        return free_params;
      }
      free_params +=
        PreAlignParameters<ResidueType, ProbabilityType>::freeParameterCount();
      free_params +=
        PostAlignParameters<ResidueType, ProbabilityType>::freeParameterCount();
      free_params +=
        InsertionEmissionParameters<ResidueType, ProbabilityType>::freeParameterCount();
      free_params +=
        PositionTransitionParameters<ResidueType, ProbabilityType>::freeParameterCount();

      return free_params;
    } // freeParameterCount() const

} // End namespace galosh

#endif // __GALOSH_PROFILE_HPP__
