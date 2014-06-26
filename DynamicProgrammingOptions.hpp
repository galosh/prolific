/**
 * \file DynamicProgrammingOptions.hpp
 * \author Ted Holzman (tholzman@scharp.org)
 * \par Library:
 * Galosh prolific
 * \brief Options for all programs using the DynamicProgramming header
 * \copyright &copy; 2008, 2011, 2012, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************/

///Useful for modifying GALOSH_DEF_OPT behavior.  See vector-valued options below
#ifndef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
#endif

/**
 * \note Note the lack of "protective" \#define symbols in this file.  We may
 * want to include it several times, redefining GALOSH_DEF_OPT as needed.
 */

      /**
       * Use (our rotated, state-specific) Rabiner-style scaling to keep
       * forward matrix values from underflowing?  Note that if this is false,
       * you should use a log type for the ScoreType (unless the profile and
       * sequences are very short and there are few sequences, in which case
       * you might not encounter underflow).
       */
GALOSH_DEF_OPT(useRabinerScaling,bool,false,"Use (our rotated, state-specific) Rabiner-style scaling to keep forward matrix values from underflowing?  Note that if this is false, you should use a log type for the ScoreType");

      /// NOTE THAT RIGHT NOW, RABINER SCALING SEEMS TO BE BROKEN.  \todo FIX!


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

GALOSH_DEF_OPT(rabinerScaling_useMaximumValue,bool,false,"We can either scale them by the total of the values there (or by the total of just the Match and Deletion values) or by the maximum value, which lets us use more of the range of the MatrixValueType.  We use the total of all values in the position unless rabinerScaling_useMaximumValue is true.");

      /// \todo Make this depend on the MatrixValueType?  If using bfloats or logs, we don't need to do this.  
      /// We can avoid a few calculations if we have this false.

      /**
       * Vanilla rabiner scaling would divide all matrix row values by the
       * total of all values in that row, but when there is a wide discrepency
       * between the largest and smallest values, the smallest ones will
       * underflow (with a non-log MatrixValueType).  So we can multiply
       * everything in the row by this scale factor to use more of the range of
       * the MatrixValueType.
       */
GALOSH_DEF_OPT(matrixRowScaleFactor,double,1.0,"We can multiply everything in the row by this scale factor to use more of the range of the MatrixValueType");

      /// \todo Make this depend on the MatrixValueType?  If using bfloats or logs, we don't need to scale it.
//#define DEFAULT_matrixRowScaleFactor ( pow( numeric_limits<float>::max(), .25f ) - 1.0f )
//#define DEFAULT_matrixRowScaleFactor ( pow( numeric_limits<double>::max(), .0625 ) - 1.0 )


/// This is how we're handling vectors.  It is a work-around because vectors are handled specially
/// by boost::program_options.  It allows the command line to look something like
///
/// --profileLengths 10 20 30
///
/// The TMP_EXTRA_STUFF must be set to include (at least) the ->multitoken() thing.
/// It should also be unset at the bottom of the vector initializations.

///Vector definition section
#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF ->multitoken()
///Vector definitions go here

/** do this after the vector definition section */
#undef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
