/*---------------------------------------------------------------------------##
##  Library:
##      galosh::prolific
##  File:
##      ProlificParameters.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      The galosh::Parameters decendent for programs using the
##      prolific library (inherits from
##      galosh::DynamicProgramming::Parameters).  Adds parameters relating to
##      priors and initial values for Profile HMMs.
#*
#******************************************************************************
#*
#*    This file is part of prolific, a library of useful C++ classes for
#*    working with genomic sequence data and Profile HMMs.  Please see the
#*    document CITING, which should have been included with this file.  You may
#*    use at will, subject to the license (Apache v2.0), but *please cite the
#*    relevant papers* in your documentation and publications associated with
#*    uses of this library.  Thank you!
#*
#*    Copyright (C) 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
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
#*****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROLIFICPARAMETERS_HPP__
#define __GALOSH_PROLIFICPARAMETERS_HPP__

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "DynamicProgramming.hpp"
using galosh::DynamicProgramming;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

namespace galosh {

template <class ResidueType,
          class ProbabilityType,
          class ScoreType,
          class MatrixValueType>
  class ProlificParameters {
    public:

    class Parameters :
    public DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters
    {
      // Boost serialization
    private:
      typedef typename DynamicProgramming<ResidueType, ProbabilityType,ScoreType,MatrixValueType>::Parameters dynamic_programming_parameters_t;
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( dynamic_programming_parameters_t );

        ar & BOOST_SERIALIZATION_NVP( useDeletionsForInsertionsParameters );
        ar & BOOST_SERIALIZATION_NVP( expectedDeletionsCounts );
        ar & BOOST_SERIALIZATION_NVP( expectedInsertionsCounts );
        ar & BOOST_SERIALIZATION_NVP( expectedDeletionLengthAsProfileLengthFractions );
        ar & BOOST_SERIALIZATION_NVP( expectedInsertionLengthAsProfileLengthFractions );
        ar & BOOST_SERIALIZATION_NVP( minExpectedDeletionLength );
        ar & BOOST_SERIALIZATION_NVP( minExpectedInsertionLength );
        ar & BOOST_SERIALIZATION_NVP( preAlignInsertion );
        ar & BOOST_SERIALIZATION_NVP( postAlignInsertion );
        ar & BOOST_SERIALIZATION_NVP( priorStrength );
        ar & BOOST_SERIALIZATION_NVP( priorStrength_internal_transitions );
        ar & BOOST_SERIALIZATION_NVP( priorMtoM );
        ar & BOOST_SERIALIZATION_NVP( priorMtoI );
        ar & BOOST_SERIALIZATION_NVP( priorMtoD );
        ar & BOOST_SERIALIZATION_NVP( priorItoM );
        ar & BOOST_SERIALIZATION_NVP( priorItoI );
        ar & BOOST_SERIALIZATION_NVP( priorDtoM );
        ar & BOOST_SERIALIZATION_NVP( priorDtoD );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_scalar );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_maxNtoN );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_maxBtoD );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_maxMtoI );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_maxMtoD );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_maxItoI );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_maxDtoD );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformGlobals_maxCtoC );
        ar & BOOST_SERIALIZATION_NVP( startWithUniformPositions );
        ar & BOOST_SERIALIZATION_NVP( startWithGlobalsDrawnFromPrior );
        ar & BOOST_SERIALIZATION_NVP( startWithPositionsDrawnFromPrior );
        ar & BOOST_SERIALIZATION_NVP( profileProfileIndelOpenCost );
        ar & BOOST_SERIALIZATION_NVP( profileProfileIndelExtensionCost );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /**
       * Lock the indel parameters of the true profile to be the same for
       * insertions as for deletions?  This makes the expectedInsertionsCounts,
       * expectedInsertionLengthAsProfileLengthFractions, and
       * minExpectedInsertionLength unused, since the corresponding deletion
       * values will be used instead.  It also reduces the number of tests by
       * reducing the number of possible combinations (since deletions and
       * insertions will go in lock step).
       */
      bool useDeletionsForInsertionsParameters;
  #define DEFAULT_useDeletionsForInsertionsParameters true

      /**
       * The deletionOpen value of the true profile will be set to (
       * expectedDeletionsCount / profileLength ).  If
       * useDeletionsForInsertionsParameters is true, the insertionOpen value
       * of the true profile will also be set to ( expectedDeletionsCount /
       * profileLength ).
       *
       * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
       * expected_deletions_count in expectedDeletionCounts.  If it is NULL,
       * { 1.0 } will be used (this is the default).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> * expectedDeletionsCounts;
  #define DEFAULT_expectedDeletionsCounts NULL

      /**
       * If useDeletionsForInsertionsParameters is false, the insertionOpen
       * value of the true profile will be set to ( expectedInsertionsCount /
       * profileLength ).
       *
       * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
       * expected_insertions_count in expectedInsertionCounts.  If it is NULL,
       * { 1.0 } will be used (this is the default).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> * expectedInsertionsCounts;
  #define DEFAULT_expectedInsertionsCounts NULL

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFraction *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).  If
       * useDeletionsForInsertionsParameters is true, the insertionExtension
       * value of the true profile will also be set to be the minimum of ( 1.0
       * / ( expectedDeletionLengthAsProfileLengthFraction * profileLength ) )
       * and ( 1.0 / minExpectedDeletionLength ).
       *
       * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
       * expected_deletion_length_as_profile_length_fraction in
       * expectedDeletionLengthAsProfileLengthFraction.  If it is NULL, { 0.1 }
       * will be used (this is the default).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> * expectedDeletionLengthAsProfileLengthFractions;
  #define DEFAULT_expectedDeletionLengthAsProfileLengthFractions NULL

      /**
       * If useDeletionsForInsertionsParameters is false, the
       * insertionExtension value of the true profile will be the minimum of (
       * 1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength
       * ) ) and ( 1.0 / minExpectedInsertionLength ).
       *
       * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
       * expected_insertion_length_as_profile_length_fraction in
       * expectedInsertionLengthAsProfileLengthFraction.  If it is NULL, { 0.1 }
       * will be used (this is the default).
       *
       * @see useDeletionsForInsertionsParameters
       */
      vector<double> * expectedInsertionLengthAsProfileLengthFractions;
  #define DEFAULT_expectedInsertionLengthAsProfileLengthFractions NULL

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFraction *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).  If
       * useDeletionsForInsertionsParameters is true, the insertionExtension
       * value of the true profile will also be the minimum of ( 1.0 / (
       * expectedDeletionLengthAsProfileLengthFraction * profileLength ) ) and
       * ( 1.0 / minExpectedDeletionLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      double minExpectedDeletionLength;
  #define DEFAULT_minExpectedDeletionLength 1.25

      /**
       * If useDeletionsForInsertionsParameters is false, the
       * insertionExtension value of the true profile will be the minimum of (
       * 1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength
       * ) ) and ( 1.0 / minExpectedInsertionLength ).
       *
       * @see useDeletionsForInsertionsParameters
       */
      double minExpectedInsertionLength;
  #define DEFAULT_minExpectedInsertionLength 1.25

      /**
       * The preAlignInsertion value of the true profile.
       */
      double preAlignInsertion;
  #define DEFAULT_preAlignInsertion .01

      /**
       * The postAlignInsertion value of the true profile.
       */
      double postAlignInsertion;
  #define DEFAULT_postAlignInsertion .01

      /**
       * The effective number of sequences "observed" a priori.  Note that we
       * use a different prior strength for main-model transitions: see
       * priorStrength_internal_transitions.
       */
      float priorStrength;
  #define DEFAULT_priorStrength 1.0f

      /**
       * The effective number of sequences "observed" a priori, for main-model
       * transitions.
       */
      float priorStrength_internal_transitions;
  #define DEFAULT_priorStrength_internal_transitions 10.0f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorMtoM;
  #define DEFAULT_priorMtoM .95f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorMtoI;
  #define DEFAULT_priorMtoI .025f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorMtoD;
  #define DEFAULT_priorMtoD .025f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorItoM;
  #define DEFAULT_priorItoM .05f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorItoI;
  #define DEFAULT_priorItoI .95f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorDtoM;
  #define DEFAULT_priorDtoM .95f

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      float priorDtoD;
  #define DEFAULT_priorDtoD .05f

      /**
       * If startWithGlobalsDrawnFromPrior is not true, and
       * if startWithUniformGlobals is true, then we set the global values of
       * the startingProfile to random values between 0 and
       * min(startWithUniformGlobals_scalar times the true
       * values,startWithUniformGlobals_maxXtoY).  If it is false, we start
       * with the known, true globals.
       *
       * @see startWithUniformGlobals_scalar
       */
      bool startWithUniformGlobals;
  #define DEFAULT_startWithUniformGlobals false

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_scalar;
  #define DEFAULT_startWithUniformGlobals_scalar 2.0

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxNtoN;
  #define DEFAULT_startWithUniformGlobals_maxNtoN .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxBtoD;
  #define DEFAULT_startWithUniformGlobals_maxBtoD .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxMtoI;
  #define DEFAULT_startWithUniformGlobals_maxMtoI .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxMtoD;
  #define DEFAULT_startWithUniformGlobals_maxMtoD .2

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxItoI;
  #define DEFAULT_startWithUniformGlobals_maxItoI .5

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxDtoD;
  #define DEFAULT_startWithUniformGlobals_maxDtoD .5

      /**
       * @see startWithUniformGlobals
       */
      double startWithUniformGlobals_maxCtoC;
  #define DEFAULT_startWithUniformGlobals_maxCtoC .2

      /**
       * If startWithUniformPositions is true, then we set the
       * position-specific values of the startingProfile to random values
       * between 0 and 1.  If it is false, we start with the known, true
       * parameter values.  Note that if startWithPositionsDrawnFromPrior is
       * also true, then the first half of the starting profiles will start
       * with positions drawn from the prior and the second half will start
       * with uniform() positions (possibly excluding the index-0 starting
       * profile, if alsoStartWithEvenPositions is true).
       *
       * @see startWithPositionsDrawnFromPrior
       * @see alsoStartWithEvenPositions
       */
      bool startWithUniformPositions;
  #define DEFAULT_startWithUniformPositions false

      /**
       * If startWithGlobalsDrawnFromPrior is true, the
       * global values of the starting profile will be drawn from the prior.
       *
       * @see startWithUniformGlobals
       */
      bool startWithGlobalsDrawnFromPrior;
  #define DEFAULT_startWithGlobalsDrawnFromPrior false

      /**
       * If startWithPositionsDrawnFromPrior is true, the
       * position-specific values of the starting profile will be drawn from
       * the prior... but see the notes in startWithUniformPositions.
       *
       * @see startWithUniformPositions
       * @see alsoStartWithEvenPositions
       */
      bool startWithPositionsDrawnFromPrior;
  #define DEFAULT_startWithPositionsDrawnFromPrior false

      /**
       * The cost of a gap open when performing SKL profile-profile alignements.
       */
      double profileProfileIndelOpenCost;
  #define DEFAULT_profileProfileIndelOpenCost .25

      /**
       * The cost of a gap extension when performing SKL profile-profile alignements.
       */
      double profileProfileIndelExtensionCost;
  #define DEFAULT_profileProfileIndelExtensionCost .25

      Parameters ();
      virtual ~Parameters () {};
    
      // Copy constructor
      template <class AnyParameters>
      Parameters ( const AnyParameters & copy_from );
    
      // Copy constructor/operator
      template <class AnyParameters>
      Parameters &
      operator= (
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

      virtual
      void
      copyFrom ( const Parameters & copy_from );
    
      virtual
      void
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
      public DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::template ParametersModifierTemplate<ParametersType>
    {
      typedef typename DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t; 

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
        ar & BOOST_SERIALIZATION_NVP( isModified_useDeletionsForInsertionsParameters );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedDeletionsCounts );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedInsertionsCounts );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedDeletionLengthAsProfileLengthFractions );
        ar & BOOST_SERIALIZATION_NVP( isModified_expectedInsertionLengthAsProfileLengthFractions );
        ar & BOOST_SERIALIZATION_NVP( isModified_minExpectedDeletionLength );
        ar & BOOST_SERIALIZATION_NVP( isModified_minExpectedInsertionLength );
        ar & BOOST_SERIALIZATION_NVP( isModified_preAlignInsertion );
        ar & BOOST_SERIALIZATION_NVP( isModified_postAlignInsertion );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorStrength );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorStrength_internal_transitions );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorMtoM );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorMtoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorMtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorItoM );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorItoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorDtoM );
        ar & BOOST_SERIALIZATION_NVP( isModified_priorDtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_scalar );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxNtoN );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxBtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxMtoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxMtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxItoI );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxDtoD );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformGlobals_maxCtoC );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithUniformPositions );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithGlobalsDrawnFromPrior );
        ar & BOOST_SERIALIZATION_NVP( isModified_startWithPositionsDrawnFromPrior );
        ar & BOOST_SERIALIZATION_NVP( isModified_profileProfileIndelOpenCost );
        ar & BOOST_SERIALIZATION_NVP( isModified_profileProfileIndelExtensionCost );
      } // serialize( Archive &, const unsigned int )

    public:
  
      /// isModified flags for Parameters
      /**
       * Lock the indel parameters of the true profile to be the same for
       * insertions as for deletions?  This makes the expectedInsertionsCounts,
       * expectedInsertionLengthAsProfileLengthFractions, and
       * minExpectedInsertionLength unused, since the corresponding deletion
       * values will be used instead.  It also reduces the number of tests by
       * reducing the number of possible combinations (since deletions and
       * insertions will go in lock step).
       */
      bool isModified_useDeletionsForInsertionsParameters;

      /**
       * The deletionOpen value of the true profile will be set to (
       * expectedDeletionsCounts / profileLength ).
       */
      bool isModified_expectedDeletionsCounts;

      /**
       * The insertionOpen value of the true profile will be set to (
       * expectedInsertionsCounts / profileLength ).
       */
      bool isModified_expectedInsertionsCounts;

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).
       */
      bool isModified_expectedDeletionLengthAsProfileLengthFractions;

      /**
       * The insertionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedInsertionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedInsertionLength ).
       */
      bool isModified_expectedInsertionLengthAsProfileLengthFractions;

      /**
       * The deletionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedDeletionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedDeletionLength ).
       */
      bool isModified_minExpectedDeletionLength;

      /**
       * The insertionExtension value of the true profile will be the minimum of
       * ( 1.0 / ( expectedInsertionLengthAsProfileLengthFractions *
       * profileLength ) ) and ( 1.0 / minExpectedInsertionLength ).
       */
      bool isModified_minExpectedInsertionLength;

      /**
       * The preAlignInsertion value of the true profile.
       */
      bool isModified_preAlignInsertion;

      /**
       * The postAlignInsertion value of the true profile.
       */
      bool isModified_postAlignInsertion;

      /**
       * The effective number of sequences "observed" a priori.
       */
      bool isModified_priorStrength;

      /**
       * The effective number of sequences "observed" a priori.
       */
      bool isModified_priorStrength_internal_transitions;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorMtoM;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorMtoI;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorMtoD;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorItoM;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * I->I transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorItoI;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorDtoM;

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * D->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
      bool isModified_priorDtoD;

      bool isModified_startWithUniformGlobals;

      bool isModified_startWithUniformGlobals_scalar;

      bool isModified_startWithUniformGlobals_maxNtoN;

      bool isModified_startWithUniformGlobals_maxBtoD;

      bool isModified_startWithUniformGlobals_maxMtoI;

      bool isModified_startWithUniformGlobals_maxMtoD;

      bool isModified_startWithUniformGlobals_maxItoI;

      bool isModified_startWithUniformGlobals_maxDtoD;

      bool isModified_startWithUniformGlobals_maxCtoC;

      bool isModified_startWithUniformPositions;

      bool isModified_startWithGlobalsDrawnFromPrior;

      bool isModified_startWithPositionsDrawnFromPrior;

      bool isModified_profileProfileIndelOpenCost;

      bool isModified_profileProfileIndelExtensionCost;

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

      template <class AnyParameters>
      void
      applyModifications ( AnyParameters & target_parameters );

    }; // End inner class ParametersModifierTemplate

    typedef ParametersModifierTemplate<typename ProlificParameters::Parameters> ParametersModifier;
  }; // End class ProlificParameters

  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::ProlificParameters::Parameters ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      Parameters ()
      {
        if( DEFAULT_debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::Parameters::<init>()" << endl;
        } // End if DEBUG_All
        resetToDefaults();
      } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      // Copy constructor
      Parameters ( const AnyParameters & copy_from )
      {
        //if( static_cast<galosh::Parameters>( copy_from ).debug >= DEBUG_All ) {
        //  cout << "[debug] ProlificParameters::Parameters::<init>( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_TRIVIAL
  typename ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters &
      // Copy constructor/operator
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      operator= (
        const AnyParameters & copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::Parameters::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_TRIVIAL
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProlificParameters::Parameters::copyFromNonVirtual( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtualDontDelegate( copy_from );
      } // copyFromNonVirtual( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      )
      {
        useDeletionsForInsertionsParameters =                             copy_from.useDeletionsForInsertionsParameters;
        expectedDeletionsCounts =                          copy_from.expectedDeletionsCounts;
        expectedInsertionsCounts =                          copy_from.expectedInsertionsCounts;
        expectedDeletionLengthAsProfileLengthFractions =                          copy_from.expectedDeletionLengthAsProfileLengthFractions;
        expectedInsertionLengthAsProfileLengthFractions =                          copy_from.expectedInsertionLengthAsProfileLengthFractions;
        minExpectedDeletionLength =                          copy_from.minExpectedDeletionLength;
        minExpectedInsertionLength =                          copy_from.minExpectedInsertionLength;
        preAlignInsertion =                          copy_from.preAlignInsertion;
        postAlignInsertion =                          copy_from.postAlignInsertion;
        priorStrength =                          copy_from.priorStrength;
        priorStrength_internal_transitions =                          copy_from.priorStrength_internal_transitions;
        priorMtoM =                          copy_from.priorMtoM;
        priorMtoI =                          copy_from.priorMtoI;
        priorMtoD =                          copy_from.priorMtoD;
        priorItoM =                          copy_from.priorItoM;
        priorItoI =                          copy_from.priorItoI;
        priorDtoM =                          copy_from.priorDtoM;
        priorDtoD =                          copy_from.priorDtoD;
        startWithUniformGlobals =                          copy_from.startWithUniformGlobals;
        startWithUniformGlobals_scalar =                          copy_from.startWithUniformGlobals_scalar;
        startWithUniformGlobals_maxNtoN =                          copy_from.startWithUniformGlobals_maxNtoN;
        startWithUniformGlobals_maxBtoD =                          copy_from.startWithUniformGlobals_maxBtoD;
        startWithUniformGlobals_maxMtoI =                          copy_from.startWithUniformGlobals_maxMtoI;
        startWithUniformGlobals_maxMtoD =                          copy_from.startWithUniformGlobals_maxMtoD;
        startWithUniformGlobals_maxItoI =                          copy_from.startWithUniformGlobals_maxItoI;
        startWithUniformGlobals_maxDtoD =                          copy_from.startWithUniformGlobals_maxDtoD;
        startWithUniformGlobals_maxCtoC =                          copy_from.startWithUniformGlobals_maxCtoC;
        startWithUniformPositions =                          copy_from.startWithUniformPositions;
        startWithGlobalsDrawnFromPrior =                          copy_from.startWithGlobalsDrawnFromPrior;
        startWithPositionsDrawnFromPrior =                          copy_from.startWithPositionsDrawnFromPrior;
        profileProfileIndelOpenCost =   copy_from.profileProfileIndelOpenCost;
        profileProfileIndelExtensionCost =   copy_from.profileProfileIndelExtensionCost;
      } // copyFromNonVirtualDontDelegate( AnyParameters const & )


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  GALOSH_INLINE_TRIVIAL
      void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFrom ( const Parameters & copy_from )
      {
        copyFromNonVirtual( copy_from );
      } // copyFrom( Parameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      resetToDefaults ()
      {
        DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::Parameters::resetToDefaults();
        // TODO: Why isn't the compiler finding "debug" in galosh::Parameters?
        //if( debug >= DEBUG_All ) {
        //  cout << "[debug] ProlificParameters::Parameters::resetToDefaults()" << endl;
        //} // End if DEBUG_All

        useDeletionsForInsertionsParameters =                             DEFAULT_useDeletionsForInsertionsParameters;
        expectedDeletionsCounts =                          DEFAULT_expectedDeletionsCounts;
        expectedInsertionsCounts =                          DEFAULT_expectedInsertionsCounts;
        expectedDeletionLengthAsProfileLengthFractions =                          DEFAULT_expectedDeletionLengthAsProfileLengthFractions;
        expectedInsertionLengthAsProfileLengthFractions =                          DEFAULT_expectedInsertionLengthAsProfileLengthFractions;
        minExpectedDeletionLength =                          DEFAULT_minExpectedDeletionLength;
        minExpectedInsertionLength =                          DEFAULT_minExpectedInsertionLength;
        preAlignInsertion =                          DEFAULT_preAlignInsertion;
        postAlignInsertion =                          DEFAULT_postAlignInsertion;
        priorStrength =                          DEFAULT_priorStrength;
        priorStrength_internal_transitions =                          DEFAULT_priorStrength_internal_transitions;
        priorMtoM =                          DEFAULT_priorMtoM;
        priorMtoI =                          DEFAULT_priorMtoI;
        priorMtoD =                          DEFAULT_priorMtoD;
        priorItoM =                          DEFAULT_priorItoM;
        priorItoI =                          DEFAULT_priorItoI;
        priorDtoM =                          DEFAULT_priorDtoM;
        priorDtoD =                          DEFAULT_priorDtoD;
        startWithUniformGlobals =                          DEFAULT_startWithUniformGlobals;
        startWithUniformGlobals_scalar =                          DEFAULT_startWithUniformGlobals_scalar;
        startWithUniformGlobals_maxNtoN =                          DEFAULT_startWithUniformGlobals_maxNtoN;
        startWithUniformGlobals_maxBtoD =                          DEFAULT_startWithUniformGlobals_maxBtoD;
        startWithUniformGlobals_maxMtoI =                          DEFAULT_startWithUniformGlobals_maxMtoI;
        startWithUniformGlobals_maxMtoD =                          DEFAULT_startWithUniformGlobals_maxMtoD;
        startWithUniformGlobals_maxItoI =                          DEFAULT_startWithUniformGlobals_maxItoI;
        startWithUniformGlobals_maxDtoD =                          DEFAULT_startWithUniformGlobals_maxDtoD;
        startWithUniformGlobals_maxCtoC =                          DEFAULT_startWithUniformGlobals_maxCtoC;
        startWithUniformPositions =                          DEFAULT_startWithUniformPositions;
        startWithGlobalsDrawnFromPrior =                          DEFAULT_startWithGlobalsDrawnFromPrior;
        startWithPositionsDrawnFromPrior =                          DEFAULT_startWithPositionsDrawnFromPrior;
        profileProfileIndelOpenCost =   DEFAULT_profileProfileIndelOpenCost;
        profileProfileIndelExtensionCost =   DEFAULT_profileProfileIndelExtensionCost;
      } // resetToDefaults()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      writeParameters (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::Parameters::writeParameters( os );
        os << endl;

        os << "[ProlificParameters]" << endl;

        os << "useDeletionsForInsertionsParameters = " <<                         useDeletionsForInsertionsParameters << endl;
        if( expectedDeletionsCounts == NULL ) {
          os << "expectedDeletionsCounts = NULL" << endl;
              } else {
          os << "expectedDeletionsCounts = { ";
          for( uint32_t cr_i = 0; cr_i < expectedDeletionsCounts->size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << ( *expectedDeletionsCounts )[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if expectedDeletionsCounts == NULL .. else ..

        if( expectedInsertionsCounts == NULL ) {
          os << "expectedInsertionsCounts = NULL" << endl;
              } else {
          os << "expectedInsertionsCounts = { ";
          for( uint32_t cr_i = 0; cr_i < expectedInsertionsCounts->size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << ( *expectedInsertionsCounts )[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if expectedInsertionsCounts == NULL .. else ..

        if( expectedDeletionLengthAsProfileLengthFractions == NULL ) {
          os << "expectedDeletionLengthAsProfileLengthFractions = NULL" << endl;
              } else {
          os << "expectedDeletionLengthAsProfileLengthFractions = { ";
          for( uint32_t cr_i = 0; cr_i < expectedDeletionLengthAsProfileLengthFractions->size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << ( *expectedDeletionLengthAsProfileLengthFractions )[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if expectedDeletionLengthAsProfileLengthFractions == NULL .. else ..

        if( expectedInsertionLengthAsProfileLengthFractions == NULL ) {
          os << "expectedInsertionLengthAsProfileLengthFractions = NULL" << endl;
              } else {
          os << "expectedInsertionLengthAsProfileLengthFractions = { ";
          for( uint32_t cr_i = 0; cr_i < expectedInsertionLengthAsProfileLengthFractions->size(); cr_i++ ) {
            if( cr_i > 0 ) {
              os << ", ";
            }
            os << ( *expectedInsertionLengthAsProfileLengthFractions )[ cr_i ];
          } // End foreach conservation rate..
          os << "}" << endl;
              } // End if expectedInsertionLengthAsProfileLengthFractions == NULL .. else ..

        os << "minExpectedDeletionLength = " <<                         minExpectedDeletionLength << endl;
        os << "minExpectedInsertionLength = " <<                         minExpectedInsertionLength << endl;
        os << "preAlignInsertion = " <<                         preAlignInsertion << endl;
        os << "postAlignInsertion = " <<                         postAlignInsertion << endl;
        os << "priorStrength = " <<                         priorStrength << endl;
        os << "priorStrength_internal_transitions = " <<                         priorStrength_internal_transitions << endl;
        os << "priorMtoM = " <<                         priorMtoM << endl;
        os << "priorMtoI = " <<                         priorMtoI << endl;
        os << "priorMtoD = " <<                         priorMtoD << endl;
        os << "priorItoM = " <<                         priorItoM << endl;
        os << "priorItoI = " <<                         priorItoI << endl;
        os << "priorDtoM = " <<                         priorDtoM << endl;
        os << "priorDtoD = " <<                         priorDtoD << endl;
        os << "startWithUniformGlobals = " <<                         startWithUniformGlobals << endl;
        os << "startWithUniformGlobals_scalar = " <<                         startWithUniformGlobals_scalar << endl;
        os << "startWithUniformGlobals_maxNtoN = " <<                         startWithUniformGlobals_maxNtoN << endl;
        os << "startWithUniformGlobals_maxBtoD = " <<                         startWithUniformGlobals_maxBtoD << endl;
        os << "startWithUniformGlobals_maxMtoI = " <<                         startWithUniformGlobals_maxMtoI << endl;
        os << "startWithUniformGlobals_maxMtoD = " <<                         startWithUniformGlobals_maxMtoD << endl;
        os << "startWithUniformGlobals_maxItoI = " <<                         startWithUniformGlobals_maxItoI << endl;
        os << "startWithUniformGlobals_maxDtoD = " <<                         startWithUniformGlobals_maxDtoD << endl;
        os << "startWithUniformGlobals_maxCtoC = " <<                         startWithUniformGlobals_maxCtoC << endl;
        os << "startWithUniformPositions = " <<                         startWithUniformPositions << endl;
        os << "startWithGlobalsDrawnFromPrior = " <<                         startWithGlobalsDrawnFromPrior << endl;
        os << "startWithPositionsDrawnFromPrior = " <<                         startWithPositionsDrawnFromPrior << endl;
        os << "profileProfileIndelOpenCost = " <<  profileProfileIndelOpenCost << endl;
        os << "profileProfileIndelExtensionCost = " <<  profileProfileIndelExtensionCost << endl;
      } // writeParameters ( basic_ostream & )

  ////// Class galosh::ProlificParameters::ParametersModifierTemplate ////
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        if( base_parameters_modifier_t::parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::ParametersModifierTemplate::<init>()" << endl;
        } // End if DEBUG_All
        isModified_reset();
      } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::ParametersModifierTemplate::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_TRIVIAL
  typename ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template ParametersModifierTemplate<ParametersType> &
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor/operator
  operator= (
        const AnyParametersModifierTemplate & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::ParametersModifierTemplate::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
        } // End if DEBUG_All

        isModified_copyFromNonVirtual( copy_from );

        base_parameters_modifier_t::parameters.copyFromNonVirtual( copy_from.parameters );
      } // copyFromNonVirtual( AnyParametersModifierTemplate const & )


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        base_parameters_modifier_t::isModified_copyFromNonVirtual( copy_from );

        isModified_useDeletionsForInsertionsParameters =                             copy_from.isModified_useDeletionsForInsertionsParameters;
        isModified_expectedDeletionsCounts =                          copy_from.isModified_expectedDeletionsCounts;
        isModified_expectedInsertionsCounts =                          copy_from.isModified_expectedInsertionsCounts;
        isModified_expectedDeletionLengthAsProfileLengthFractions =                          copy_from.isModified_expectedDeletionLengthAsProfileLengthFractions;
        isModified_expectedInsertionLengthAsProfileLengthFractions =                          copy_from.isModified_expectedInsertionLengthAsProfileLengthFractions;
        isModified_minExpectedDeletionLength =                          copy_from.isModified_minExpectedDeletionLength;
        isModified_minExpectedInsertionLength =                          copy_from.isModified_minExpectedInsertionLength;
        isModified_preAlignInsertion =                          copy_from.isModified_preAlignInsertion;
        isModified_postAlignInsertion =                          copy_from.isModified_postAlignInsertion;
        isModified_priorStrength =                          copy_from.isModified_priorStrength;
        isModified_priorStrength_internal_transitions =                          copy_from.isModified_priorStrength_internal_transitions;
        isModified_priorMtoM =                          copy_from.isModified_priorMtoM;
        isModified_priorMtoI =                          copy_from.isModified_priorMtoI;
        isModified_priorMtoD =                          copy_from.isModified_priorMtoD;
        isModified_priorItoM =                          copy_from.isModified_priorItoM;
        isModified_priorItoI =                          copy_from.isModified_priorItoI;
        isModified_priorDtoM =                          copy_from.isModified_priorDtoM;
        isModified_priorDtoD =                          copy_from.isModified_priorDtoD;
        isModified_startWithUniformGlobals =                          copy_from.isModified_startWithUniformGlobals;
        isModified_startWithUniformGlobals_scalar =                          copy_from.isModified_startWithUniformGlobals_scalar;
        isModified_startWithUniformGlobals_maxNtoN =                          copy_from.isModified_startWithUniformGlobals_maxNtoN;
        isModified_startWithUniformGlobals_maxBtoD =                          copy_from.isModified_startWithUniformGlobals_maxBtoD;
        isModified_startWithUniformGlobals_maxMtoI =                          copy_from.isModified_startWithUniformGlobals_maxMtoI;
        isModified_startWithUniformGlobals_maxMtoD =                          copy_from.isModified_startWithUniformGlobals_maxMtoD;
        isModified_startWithUniformGlobals_maxItoI =                          copy_from.isModified_startWithUniformGlobals_maxItoI;
        isModified_startWithUniformGlobals_maxDtoD =                          copy_from.isModified_startWithUniformGlobals_maxDtoD;
        isModified_startWithUniformGlobals_maxCtoC =                          copy_from.isModified_startWithUniformGlobals_maxCtoC;
        isModified_startWithUniformPositions =                          copy_from.isModified_startWithUniformPositions;
        isModified_startWithGlobalsDrawnFromPrior =                          copy_from.isModified_startWithGlobalsDrawnFromPrior;
        isModified_startWithPositionsDrawnFromPrior =                          copy_from.isModified_startWithPositionsDrawnFromPrior;
        isModified_profileProfileIndelOpenCost =   copy_from.isModified_profileProfileIndelOpenCost;
        isModified_profileProfileIndelExtensionCost =   copy_from.isModified_profileProfileIndelExtensionCost;
      } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_TRIVIAL
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      reset ()
      {
        isModified_reset();
        base_parameters_modifier_t::parameters.resetToDefaults();
      } // reset()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        base_parameters_modifier_t::isModified_reset();

        isModified_useDeletionsForInsertionsParameters = false;
        isModified_expectedDeletionsCounts = false;
        isModified_expectedInsertionsCounts = false;
        isModified_expectedDeletionLengthAsProfileLengthFractions = false;
        isModified_expectedInsertionLengthAsProfileLengthFractions = false;
        isModified_minExpectedDeletionLength = false;
        isModified_minExpectedInsertionLength = false;
        isModified_preAlignInsertion = false;
        isModified_postAlignInsertion = false;
        isModified_priorStrength = false;
        isModified_priorStrength_internal_transitions = false;
        isModified_priorMtoM = false;
        isModified_priorMtoI = false;
        isModified_priorMtoD = false;
        isModified_priorItoM = false;
        isModified_priorItoI = false;
        isModified_priorDtoM = false;
        isModified_priorDtoD = false;
        isModified_startWithUniformGlobals = false;
        isModified_startWithUniformGlobals_scalar = false;
        isModified_startWithUniformGlobals_maxNtoN =  false;
        isModified_startWithUniformGlobals_maxBtoD =  false;
        isModified_startWithUniformGlobals_maxMtoI =  false;
        isModified_startWithUniformGlobals_maxMtoD =  false;
        isModified_startWithUniformGlobals_maxItoI =  false;
        isModified_startWithUniformGlobals_maxDtoD =  false;
        isModified_startWithUniformGlobals_maxCtoC =  false;
        isModified_startWithUniformPositions =        false;
        isModified_startWithGlobalsDrawnFromPrior =   false;
        isModified_startWithPositionsDrawnFromPrior = false;
        isModified_profileProfileIndelOpenCost = false;
        isModified_profileProfileIndelExtensionCost = false;
      } // isModified_reset()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::basic_ostream<CharT,Traits>& os
      ) const
      {
        //base_parameters_modifier_t::operator<<( os, parameters_modifier );
        base_parameters_modifier_t::writeParametersModifier( os );
        os << endl;

        os << "[ProlificParameters]" << endl;
        if( isModified_useDeletionsForInsertionsParameters ) {
          os << "useDeletionsForInsertionsParameters = " <<                         base_parameters_modifier_t::parameters.useDeletionsForInsertionsParameters << endl;
        }
        if( isModified_expectedDeletionsCounts ) {
          if( base_parameters_modifier_t::parameters.expectedDeletionsCounts == NULL ) {
            os << "expectedDeletionsCounts = NULL" << endl;
          } else {
            os << "expectedDeletionsCounts = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedDeletionsCounts->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedDeletionsCounts )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedDeletionsCounts == NULL .. else ..
        }
        if( isModified_expectedInsertionsCounts ) {
          if( base_parameters_modifier_t::parameters.expectedInsertionsCounts == NULL ) {
            os << "expectedInsertionsCounts = NULL" << endl;
          } else {
            os << "expectedInsertionsCounts = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedInsertionsCounts->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedInsertionsCounts )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedInsertionsCounts == NULL .. else ..
        }
        if( isModified_expectedDeletionLengthAsProfileLengthFractions ) {
          if( base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions == NULL ) {
            os << "expectedDeletionLengthAsProfileLengthFractions = NULL" << endl;
          } else {
            os << "expectedDeletionLengthAsProfileLengthFractions = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedDeletionLengthAsProfileLengthFractions == NULL .. else ..
        }
        if( isModified_expectedInsertionLengthAsProfileLengthFractions ) {
          if( base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions == NULL ) {
            os << "expectedInsertionLengthAsProfileLengthFractions = NULL" << endl;
          } else {
            os << "expectedInsertionLengthAsProfileLengthFractions = { ";
            for( uint32_t cr_i = 0; cr_i < base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions->size(); cr_i++ ) {
              if( cr_i > 0 ) {
                os << ", ";
              }
              os << ( *base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions )[ cr_i ];
            } // End foreach conservation rate..
            os << "}" << endl;
          } // End if expectedInsertionLengthAsProfileLengthFractions == NULL .. else ..
        }
        if( isModified_minExpectedDeletionLength ) {
          os << "minExpectedDeletionLength = " <<                         base_parameters_modifier_t::parameters.minExpectedDeletionLength << endl;
        }
        if( isModified_minExpectedInsertionLength ) {
          os << "minExpectedInsertionLength = " <<                         base_parameters_modifier_t::parameters.minExpectedInsertionLength << endl;
        }
        if( isModified_preAlignInsertion ) {
          os << "preAlignInsertion = " <<                         base_parameters_modifier_t::parameters.preAlignInsertion << endl;
        }
        if( isModified_postAlignInsertion ) {
          os << "postAlignInsertion = " <<                         base_parameters_modifier_t::parameters.postAlignInsertion << endl;
        }
        if( isModified_priorStrength ) {
          os << "priorStrength = " <<                         base_parameters_modifier_t::parameters.priorStrength << endl;
        }
        if( isModified_priorStrength_internal_transitions ) {
          os << "priorStrength_internal_transitions = " <<                         base_parameters_modifier_t::parameters.priorStrength_internal_transitions << endl;
        }
        if( isModified_priorMtoM ) {
          os << "priorMtoM = " <<                         base_parameters_modifier_t::parameters.priorMtoM << endl;
        }
        if( isModified_priorMtoI ) {
          os << "priorMtoI = " <<                         base_parameters_modifier_t::parameters.priorMtoI << endl;
        }
        if( isModified_priorMtoD ) {
          os << "priorMtoD = " <<                         base_parameters_modifier_t::parameters.priorMtoD << endl;
        }
        if( isModified_priorItoM ) {
          os << "priorItoM = " <<                         base_parameters_modifier_t::parameters.priorItoM << endl;
        }
        if( isModified_priorItoI ) {
          os << "priorItoI = " <<                         base_parameters_modifier_t::parameters.priorItoI << endl;
        }
        if( isModified_priorDtoM ) {
          os << "priorDtoM = " <<                         base_parameters_modifier_t::parameters.priorDtoM << endl;
        }
        if( isModified_priorDtoD ) {
          os << "priorDtoD = " <<                         base_parameters_modifier_t::parameters.priorDtoD << endl;
        }
        if( isModified_startWithUniformGlobals ) {
          os << "startWithUniformGlobals = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals << endl;
        }
        if( isModified_startWithUniformGlobals_scalar ) {
          os << "startWithUniformGlobals_scalar = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_scalar << endl;
        }
        if( isModified_startWithUniformGlobals_maxNtoN ) {
          os << "startWithUniformGlobals_maxNtoN = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxNtoN << endl;
        }
        if( isModified_startWithUniformGlobals_maxBtoD ) {
          os << "startWithUniformGlobals_maxBtoD = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxBtoD << endl;
        }
        if( isModified_startWithUniformGlobals_maxMtoI ) {
          os << "startWithUniformGlobals_maxMtoI = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoI << endl;
        }
        if( isModified_startWithUniformGlobals_maxMtoD ) {
          os << "startWithUniformGlobals_maxMtoD = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoD << endl;
        }
        if( isModified_startWithUniformGlobals_maxItoI ) {
          os << "startWithUniformGlobals_maxItoI = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxItoI << endl;
        }
        if( isModified_startWithUniformGlobals_maxDtoD ) {
          os << "startWithUniformGlobals_maxDtoD = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxDtoD << endl;
        }
        if( isModified_startWithUniformGlobals_maxCtoC ) {
          os << "startWithUniformGlobals_maxCtoC = " <<                         base_parameters_modifier_t::parameters.startWithUniformGlobals_maxCtoC << endl;
        }
        if( isModified_startWithUniformPositions ) {
          os << "startWithUniformPositions = " <<                         base_parameters_modifier_t::parameters.startWithUniformPositions << endl;
        }
        if( isModified_startWithGlobalsDrawnFromPrior ) {
          os << "startWithGlobalsDrawnFromPrior = " <<                         base_parameters_modifier_t::parameters.startWithGlobalsDrawnFromPrior << endl;
        }
        if( isModified_startWithPositionsDrawnFromPrior ) {
          os << "startWithPositionsDrawnFromPrior = " <<                         base_parameters_modifier_t::parameters.startWithPositionsDrawnFromPrior << endl;
        }
        if( isModified_profileProfileIndelOpenCost ) {
          os << "profileProfileIndelOpenCost = " <<  base_parameters_modifier_t::parameters.profileProfileIndelOpenCost << endl;
        }
        if( isModified_profileProfileIndelExtensionCost ) {
          os << "profileProfileIndelExtensionCost = " <<  base_parameters_modifier_t::parameters.profileProfileIndelExtensionCost << endl;
        }
      } // writeParametersModifier ( basic_ostream & ) const


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParameters>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( AnyParameters & target_parameters )
      {
        base_parameters_modifier_t::applyModifications( target_parameters );

        if( isModified_useDeletionsForInsertionsParameters ) {
          target_parameters.useDeletionsForInsertionsParameters =
            base_parameters_modifier_t::parameters.useDeletionsForInsertionsParameters;
        }
        if( isModified_expectedDeletionsCounts ) {
          target_parameters.expectedDeletionsCounts =
            base_parameters_modifier_t::parameters.expectedDeletionsCounts;
        }
        if( isModified_expectedInsertionsCounts ) {
          target_parameters.expectedInsertionsCounts =
            base_parameters_modifier_t::parameters.expectedInsertionsCounts;
        }
        if( isModified_expectedDeletionLengthAsProfileLengthFractions ) {
          target_parameters.expectedDeletionLengthAsProfileLengthFractions =
            base_parameters_modifier_t::parameters.expectedDeletionLengthAsProfileLengthFractions;
        }
        if( isModified_expectedInsertionLengthAsProfileLengthFractions ) {
          target_parameters.expectedInsertionLengthAsProfileLengthFractions =
            base_parameters_modifier_t::parameters.expectedInsertionLengthAsProfileLengthFractions;
        }
        if( isModified_minExpectedDeletionLength ) {
          target_parameters.minExpectedDeletionLength =
            base_parameters_modifier_t::parameters.minExpectedDeletionLength;
        }
        if( isModified_minExpectedInsertionLength ) {
          target_parameters.minExpectedInsertionLength =
            base_parameters_modifier_t::parameters.minExpectedInsertionLength;
        }
        if( isModified_preAlignInsertion ) {
          target_parameters.preAlignInsertion =
            base_parameters_modifier_t::parameters.preAlignInsertion;
        }
        if( isModified_postAlignInsertion ) {
          target_parameters.postAlignInsertion =
            base_parameters_modifier_t::parameters.postAlignInsertion;
        }
        if( isModified_priorStrength ) {
          target_parameters.priorStrength =
            base_parameters_modifier_t::parameters.priorStrength;
        }
        if( isModified_priorStrength_internal_transitions ) {
          target_parameters.priorStrength_internal_transitions =
            base_parameters_modifier_t::parameters.priorStrength_internal_transitions;
        }
        if( isModified_priorMtoM ) {
          target_parameters.priorMtoM =
            base_parameters_modifier_t::parameters.priorMtoM;
        }
        if( isModified_priorMtoI ) {
          target_parameters.priorMtoI =
            base_parameters_modifier_t::parameters.priorMtoI;
        }
        if( isModified_priorMtoD ) {
          target_parameters.priorMtoD =
            base_parameters_modifier_t::parameters.priorMtoD;
        }
        if( isModified_priorItoM ) {
          target_parameters.priorItoM =
            base_parameters_modifier_t::parameters.priorItoM;
        }
        if( isModified_priorItoI ) {
          target_parameters.priorItoI =
            base_parameters_modifier_t::parameters.priorItoI;
        }
        if( isModified_priorDtoM ) {
          target_parameters.priorDtoM =
            base_parameters_modifier_t::parameters.priorDtoM;
        }
        if( isModified_priorDtoD ) {
          target_parameters.priorDtoD =
            base_parameters_modifier_t::parameters.priorDtoD;
        }
        if( isModified_startWithUniformGlobals ) {
          target_parameters.startWithUniformGlobals =
            base_parameters_modifier_t::parameters.startWithUniformGlobals;
        }
        if( isModified_startWithUniformGlobals_scalar ) {
          target_parameters.startWithUniformGlobals_scalar =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_scalar;
        }
        if( isModified_startWithUniformGlobals_maxNtoN ) {
          target_parameters.startWithUniformGlobals_maxNtoN =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxNtoN;
        }
        if( isModified_startWithUniformGlobals_maxBtoD ) {
          target_parameters.startWithUniformGlobals_maxBtoD =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxBtoD;
        }
        if( isModified_startWithUniformGlobals_maxMtoI ) {
          target_parameters.startWithUniformGlobals_maxMtoI =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoI;
        }
        if( isModified_startWithUniformGlobals_maxMtoD ) {
          target_parameters.startWithUniformGlobals_maxMtoD =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxMtoD;
        }
        if( isModified_startWithUniformGlobals_maxItoI ) {
          target_parameters.startWithUniformGlobals_maxItoI =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxItoI;
        }
        if( isModified_startWithUniformGlobals_maxDtoD ) {
          target_parameters.startWithUniformGlobals_maxDtoD =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxDtoD;
        }
        if( isModified_startWithUniformGlobals_maxCtoC ) {
          target_parameters.startWithUniformGlobals_maxCtoC =
            base_parameters_modifier_t::parameters.startWithUniformGlobals_maxCtoC;
        }
        if( isModified_startWithUniformPositions ) {
          target_parameters.startWithUniformPositions =
            base_parameters_modifier_t::parameters.startWithUniformPositions;
        }
        if( isModified_startWithGlobalsDrawnFromPrior ) {
          target_parameters.startWithGlobalsDrawnFromPrior =
            base_parameters_modifier_t::parameters.startWithGlobalsDrawnFromPrior;
        }
        if( isModified_startWithPositionsDrawnFromPrior ) {
          target_parameters.startWithPositionsDrawnFromPrior =
            base_parameters_modifier_t::parameters.startWithPositionsDrawnFromPrior;
        }
        if( isModified_profileProfileIndelOpenCost ) {
          target_parameters.profileProfileIndelOpenCost =
            base_parameters_modifier_t::parameters.profileProfileIndelOpenCost;
        }
        if( isModified_profileProfileIndelExtensionCost ) {
          target_parameters.profileProfileIndelExtensionCost =
            base_parameters_modifier_t::parameters.profileProfileIndelExtensionCost;
        }
      } // applyModifications( Parameters & )

} // End namespace galosh

#endif // __GALOSH_PROLIFICPARAMETERS_HPP__
