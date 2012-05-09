/**
 * \file Prolific.hpp
 * \author
 *      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      Global \#defines for Prolific.
 * \par Overview:
 *    This file is part of prolific, a library of useful C++ classes for
 *    working with genomic sequence data and Profile HMMs.  Please see the
 *    document CITING, which should have been included with this file.  You may
 *    use at will, subject to the license (Apache v2.0), but *please cite the
 *    relevant papers* in your documentation and publications associated with
 *    uses of this library.  Thank you!
 *
 *  \copyright &copy; 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  \par License:
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

#ifndef __GALOSH_PROLIFIC_HPP__
#define __GALOSH_PROLIFIC_HPP__

#if defined __MSC_VER
#include <boost/math/special_functions/log1p.hpp>
using boost::math::log1p;
#endif

////// What shall we inline?
#define GALOSH_INLINE_INIT inline
#define GALOSH_INLINE_REINITIALIZE inline
#define GALOSH_INLINE_COPY inline
#define GALOSH_INLINE_ACCESSOR inline
#define GALOSH_INLINE_OSTREAM

// "trivial" means it is just a renaming of another fn (so all it does is call that other fn), or perhaps returning simple function of that call.
#define GALOSH_INLINE_TRIVIAL inline

// Probabilities
#ifdef GALOSH_USE_ICSILOG
  #define GALOSH_LOG(x) icsilogv2(x)
#else
  #define GALOSH_LOG(x) ::log(x)
#endif
#if defined __MSC_VER
#define GALOSH_LOG1P(x) boost::math::log1p(x)
#else
#define GALOSH_LOG1P(x) ::log1p(x)
#endif

// MultinomialDistribution
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_ARITHMETIC inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPARE inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_NORMALIZE_ENSUREMIN inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_COMPLEX_ACCESSOR inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_EUCLIDEANDISTSQ inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_DRAW inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_FROMDIRICHLET inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_TO_BOLTZMANNGIBBS inline
#define GALOSH_INLINE_MULTINOMIALDISTRIBUTION_FROM_BOLTZMANNGIBBS inline

// Parameters, DynamicProgramming, ...
#define GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS inline

// Profile
#define GALOSH_INLINE_PROFILE_ARITHMETIC inline
#define GALOSH_INLINE_PROFILE_EUCLIDEANDISTSQ inline
#define GALOSH_INLINE_PROFILE_CROSSENTROPY inline
#define GALOSH_INLINE_PROFILE_FREEPARAMCOUNT inline
#define GALOSH_INLINE_PROFILE_COMPLEX_ACCESSOR inline

// DynamicProgramming
#define GALOSH_INLINE_ALGORITHM inline
#define GALOSH_INLINE_ALGORITHM_INNERLOOP inline

#endif // __GALOSH_PROLIFIC_HPP__
