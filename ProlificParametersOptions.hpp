/**
 * \file ProlificParametersOptions.hpp
 * \author Ted Holzman (tholzman@scharp.org)
 * \par Library:
 * Galosh prolific
 * \brief Options for all programs using the ProlificParameters header
 * \copyright &copy; 2008, 2011, 2012, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************/
/**
 * \note Note the lack of "protective" \#define symbols in this file.  We may
 * want to include it several times, redefining GALOSH_DEF_OPT as needed.
 */

///Useful for modifying GALOSH_DEF_OPT behavior.  See vector-valued options below
#ifndef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
#endif

 /**
  * A scalar to use to multiply the matrix values.  TODO: DISABLE?
  */
GALOSH_DEF_OPT(matrixRowScaleFactor,double,1.0,"A scalar by which to multiply matrix rows; best to leave this one alone.");

/**
 * Lock the indel parameters of the true profile to be the same for
 * insertions as for deletions?  This makes the expectedInsertionsCounts,
 * expectedInsertionLengthAsProfileLengthFractions, and
 * minExpectedInsertionLength unused, since the corresponding deletion
 * values will be used instead.  It also reduces the number of tests by
 * reducing the number of possible combinations (since deletions and
 * insertions will go in lock step).
 */
GALOSH_DEF_OPT(useDeletionsForInsertionsParameters,bool,true,"Should insertion prob be forced equal to deletion prob?");

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
GALOSH_DEF_OPT(minExpectedDeletionLength,double,1.25,"The deletionExtension value of the true profile will be the minimum of  ( 1.0 / ( expectedDeletionLengthAsProfileLengthFraction * profileLength ) ) and ( 1.0 / minExpectedDeletionLength )");

 /**
  * If useDeletionsForInsertionsParameters is false, the
  * insertionExtension value of the true profile will be the minimum of (
  * 1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength
  * ) ) and ( 1.0 / minExpectedInsertionLength ).
  *
  * @see useDeletionsForInsertionsParameters
  */
GALOSH_DEF_OPT(minExpectedInsertionLength,double,1.25,"If useDeletionsForInsertionsParameters is false, the insertionExtension value of the true profile will be the minimum of (1.0 / ( expectedInsertionLengthAsProfileLengthFraction * profileLength) ) and ( 1.0 / minExpectedInsertionLength ).");

 /**
  * The effective number of sequences "observed" a priori.  Note that we
  * use a different prior strength for main-model transitions: see
  * priorStrength_internal_transitions.
  */
GALOSH_DEF_OPT(priorStrength,float,1.0f,"The effective number of sequences \"observed\" a priori");

 /**
  * The effective number of sequences "observed" a priori, for main-model
  * transitions.
  */
GALOSH_DEF_OPT(priorStrength_internal_transitions,float,10.0f,"The effective number of sequences \"observed\" a priori, for main-model transitions");

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * B->D transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
#ifdef USE_DEL_IN_DEL_OUT
GALOSH_DEF_OPT(priorBtoD,float,0.01f,"The prior contribution (per 'a priori sequence': see priorStrength) of B->D transitions");
#else
GALOSH_DEF_OPT(priorBtoD,float,0.01f,"The prior contribution (per 'a priori sequence': see priorStrength) of B->D transitions");
#endif

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * B->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
#ifdef USE_DEL_IN_DEL_OUT
GALOSH_DEF_OPT(priorBtoM,float,0.09f,"The prior contribution (per 'a priori sequence': see priorStrength) of B->M transitions");
#else
GALOSH_DEF_OPT(priorBtoM,float,0.99f,"The prior contribution (per 'a priori sequence': see priorStrength) of B->M transitions");
#endif

#ifdef USE_DEL_IN_DEL_OUT
      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * B->Z transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
GALOSH_DEF_OPT(priorBtoZ,float,0.9f,"The prior contribution (per 'a priori sequence': see priorStrength) of B->Z transitions");
#endif

#ifndef DISALLOW_FLANKING_TRANSITIONS
      /**
       * The prior contribution  of
       * C->C transitions. 
       */
GALOSH_DEF_OPT(priorCtoC,float,0.01f,"The prior contribution of C->C transitions");

      /**
       * The prior contribution  of
       * C->T transitions. 
       */
GALOSH_DEF_OPT(priorCtoT,float,0.99f,"The prior contribution of C->T transitions");

      /**
       * The prior contribution  of
       * N->N transitions. 
       */
GALOSH_DEF_OPT(priorNtoN,float,0.01f,"The prior contribution of N->N transitions");

      /**
       * The prior contribution  of
       * N->B transitions. 
       */
GALOSH_DEF_OPT(priorNtoB,float,0.99f,"The prior contribution of N->B transitions");

#endif /* !DISALLOW_FLANKING_TRANSITIONS */

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * M->M transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
#ifdef USE_DEL_IN_DEL_OUT
GALOSH_DEF_OPT(priorMtoM,float,0.095f,"The prior contribution (per 'a priori sequence': see priorStrength) of M->M transitions");
#else
GALOSH_DEF_OPT(priorMtoM,float,0.95f,"The prior contribution (per \"a priori sequence\": see priorStrength) of M->M transitions");
#endif

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * M->I transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
#ifdef USE_DEL_IN_DEL_OUT
GALOSH_DEF_OPT(priorMtoI,float,0.0025f,"The prior contribution (per 'a priori sequence': see priorStrength) of M->I transitions");
#else
GALOSH_DEF_OPT(priorMtoI,float,0.025f,"The prior contribution (per \"a priori sequence\": see priorStrength) of M->I transitions");
#endif

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * M->D transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
#ifdef USE_DEL_IN_DEL_OUT
GALOSH_DEF_OPT(priorMtoD,float,0.0025f,"The prior contribution (per 'a priori sequence': see priorStrength) of M->D transitions");
#else
GALOSH_DEF_OPT(priorMtoD,float,0.025f,"The prior contribution (per \"a priori sequence\": see priorStrength) of M->D transitions");
#endif


#ifdef USE_DEL_IN_DEL_OUT
      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * M->W transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
GALOSH_DEF_OPT(priorMtoW,float,0.9f,"The prior contribution (per 'a priori sequence': see priorStrength) of M->W transitions");

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * W->W transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
GALOSH_DEF_OPT(priorWtoW,float,0.9999f,"The prior contribution (per 'a priori sequence': see priorStrength) of W->W transitions");

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * W->E transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
GALOSH_DEF_OPT(priorWtoE,float,0.0001f,"The prior contribution (per 'a priori sequence': see priorStrength) of W->E transitions");


      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * Z->Z transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
GALOSH_DEF_OPT(priorZtoZ,float,0.9999f,"The prior contribution (per 'a priori sequence': see priorStrength) of Z->Z transitions");

      /**
       * The prior contribution (per "a priori sequence": see priorStrength) of
       * Z->M transitions.  This will be multiplied by the profile length and
       * by the priorStrength when setting up the global prior.
       */
GALOSH_DEF_OPT(priorZtoM,float,0.0001f,"The prior contribution (per 'a priori sequence': see priorStrength) of Z->M transitions");
#endif /* USE_DEL_IN_DEL_OUT */

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * I->M transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorItoM,float,.05f,"The prior contribution (per \"a priori sequence\": see priorStrength) of I->M transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * I->I transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorItoI,float,.95f,"The prior contribution (per \"a priori sequence\": see priorStrength) of I->I transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * D->M transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorDtoM,float,.95f," The prior contribution (per \"a priori sequence\": see priorStrength) of D->M transitions");

 /**
  * The prior contribution (per "a priori sequence": see priorStrength) of
  * D->D transitions.  This will be multiplied by the profile length and
  * by the priorStrength when setting up the global prior.
  */
GALOSH_DEF_OPT(priorDtoD,float,.05f,"The prior contribution (per \"a priori sequence\": see priorStrength) of D->D transitions");

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
GALOSH_DEF_OPT(startWithUniformGlobals,bool,false,"If startWithGlobalsDrawnFromPrior is not true, and if startWithUniformGlobals is true, then we set the global values of the startingProfile to random values between 0 and min(startWithUniformGlobals_scalar times the true values,startWithUniformGlobals_maxXtoY).  If it is false, we start with the known, true globals.");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_scalar,double,2.0,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxNtoN,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxBtoD,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxMtoI,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxMtoD,double,.2,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxItoI,double,.5,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxDtoD,double,.5,"See \"startWithUniformGlobals\"");

 /**
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithUniformGlobals_maxCtoC,double,.2,"See \"startWithUniformGlobals\"");

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
GALOSH_DEF_OPT(startWithUniformPositions,bool,false,"If startWithUniformPositions is true, then we set the position-specific values of the startingProfile to random values between 0 and 1.  If it is false, we start with the known, true parameter values.");

 /**
  * If startWithGlobalsDrawnFromPrior is true, the
  * global values of the starting profile will be drawn from the prior.
  *
  * @see startWithUniformGlobals
  */
GALOSH_DEF_OPT(startWithGlobalsDrawnFromPrior,bool,false,"If startWithGlobalsDrawnFromPrior is true, the global values of the starting profile will be drawn from the prior");

 /**
  * If startWithPositionsDrawnFromPrior is true, the
  * position-specific values of the starting profile will be drawn from
  * the prior... but see the notes in startWithUniformPositions.
  *
  * @see startWithUniformPositions
  * @see alsoStartWithEvenPositions
  */
GALOSH_DEF_OPT(startWithPositionsDrawnFromPrior,bool,false,"If startWithPositionsDrawnFromPrior is true, the position-specific values of the starting profile will be drawn from the prior");

 /**
  * The cost of a gap open when performing SKL profile-profile alignements.
  */
GALOSH_DEF_OPT(profileProfileIndelOpenCost,double,.25,"The cost of a gap open when performing SKL profile-profile alignements");

 /**
  * The cost of a gap extension when performing SKL profile-profile alignments.
  */
GALOSH_DEF_OPT(profileProfileIndelExtensionCost,double,0.25,"The cost of a gap extension when performing SKL profile-profile alignment");

/**
 * BROKEN, DO NOT USE.
 */
GALOSH_DEF_OPT(useRabinerScaling,bool,false,"Temporarily broken; do not use.");

/**
 * BROKEN, DO NOT USE.
 */
GALOSH_DEF_OPT(rabinerScaling_useMaximumValue,bool,false,"Temporarily broken; do not use.");




///
/// Vector/Array section
///
/// This is how we're handling vectors.  It is a work-around because vectors are handled specially
/// by boost::program_options.  It allows the command line to look something like
///
/// --profileLengths 10 20 30
///
/// The TMP_EXTRA_STUFF must be set to include (at least) the ->multitoken() thing.
/// It should also be unset at the bottom of the vector initializations.

#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF ->multitoken()
///Vector definitions go here

/// Make the expected number of deletions be .5 or 1.0 per sequence.
/// Note that this will apply to the insertions, too, unless
/// m_parameters.useDeletionsForInsertionsParameters is set to false.
GALOSH_DEF_OPT(expectedDeletionsCounts,myVector<double>,myVector<double>(1,1.0) BOOST_PP_COMMA() string("1.0"),"Iterate through this set of expected deletions.  Note that deletions-counts=insertion counts unless useDeletionsForInsertionsParameters is false");

/**
  * If useDeletionsForInsertionsParameters is false, the insertionOpen
  * value of the true profile will be set to ( expectedInsertionsCount /
  * profileLength ).
  *
  * UPDATE: This is now a pointer to a vector.  Tests will be run foreach
  * expected_insertions_count in expectedInsertionCounts.  If it is unspecified,
  * { 1.0 } will be used.
  * @see useDeletionsForInsertionsParameters
  */
GALOSH_DEF_OPT(expectedInsertionsCounts,myVector<double>,myVector<double>(1,1.0) BOOST_PP_COMMA() string("1.0"),"Iterate through this series of expected insertions.  useDeletionsForInsertionsParameters must be false or this parameter is ignored.");

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

// Make the expected length of each deletion be ( profile_length / 20 ) or (
// profile_length / 10 )...
// Note that this will apply to the insertions, too, unless
// m_parameters.useDeletionsForInsertionsParameters is set to false.
GALOSH_DEF_OPT(expectedDeletionLengthAsProfileLengthFractions,myVector<double>,myVector<double>(1,0.1) BOOST_PP_COMMA() string("0.1"),"Expected lengths of deletions. Iterate through all lengths in this list.");

// Make the expected length of each insertion be ( profile_length / 20 ) or (
// profile_length / 10 )...
// Note that this is not used unless
// m_parameters.useDeletionsForInsertionsParameters is set to false.
// ..(or 1.25, whichever is larger).
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
GALOSH_DEF_OPT(expectedInsertionLengthAsProfileLengthFractions,myVector<double>,myVector<double>(1,0.1) BOOST_PP_COMMA() string("0.1"),"Expected lengths of deletions. Iterate through all lengths in this list.");

/** do this after the vector definition section */
#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
