/**
 * \file GaloshOptions.hpp
 * \author Ted Holzman (tholzman@scharp.org)
 * \par Library
 * Galosh
 * \brief 
 * Basic options for "all" (kindof a lot of) Galosh related programs
 * \copyright &copy; 2008, 2011, 2012, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************/

// TODO: I think this is _not_ the final solution -
// I think the options should somehow be chained together later, but I'm not sure how that's done. - DOPTE

/**
 * \note Note the lack of "protective" \#define symbols in this file.  We may
 * want to include it several times, redefining GALOSH_DEF_OPT as needed.
 */

///Useful for modifying GALOSH_DEF_OPT behavior.  See vector-valued options below
#ifndef TMP_EXTRA_STUFF
#define TMP_EXTRA_STUFF BOOST_PP_EMPTY()
#endif

GALOSH_DEF_OPT(debug,int,DEBUG_None,"Controls amount of diagnostic information");
GALOSH_DEF_OPT(verbosity,int,VERBOSITY_Low,"Controls verbosity of informative messages");

