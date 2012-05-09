/**
 * \file Ambiguous.hpp
 * \author  D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      Metafunctions for describing how one residue type (eg seqan::Dna5 or
 *      seqan::Iupac) is ambiguous over another (eg seqan::Dna).
 * \par Oveview:
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

#ifndef __GALOSH_AMBIGUOUS_HPP__
#define __GALOSH_AMBIGUOUS_HPP__

// TODO: REMOVE
#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h> // for ordValue

#include "Prolific.hpp"
#include "AminoAcid20.hpp"

namespace galosh {

/* ///////////////////////////////////////////////////////////////////////////// */
//IsAmbiguous
/* ////////////////////////////////////////////////////////////////////////////// */

/**
.Metafunction.IsAmbiguous:
..summary:Tests whether one type (the "larger type") is ambiguous over another (the "smaller type") -- that is, whether there exists a multivalued mapping from the elements of the larger (enumerated) type to elements of the smaller (enumerated) type.
..signature:IsAmbiguous<LargerType, SmallerType>::Type
..param.LargerType:The enumerated type being tested for ambiguity over the other type.
..param.SmallerType:The other enumerated type.
..returns.param.Type:@Tag.Logical Values.True@, if $LargerType$ is ambiguous over $SmallerType$, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..remarks:Both types should have a @Metafunction.ValueSize@, and the @Metafunction.ValueSize@ of the former should be greater than or equal to that of the latter.  The @Metafunction.AmbiguousValueSize@ of at least one value in 0..ValueSize<LargerType>::VALUE should be greater than 1, and none should be less than 1.
..see:Metafunction.AmbiguousValueSize
..see:Metafunction.ValueSize
*/

template <typename LargerType, typename SmallerType>
struct _IsAmbiguous {
  typedef seqan::False Type;
};

template <typename LargerType, typename SmallerType>
struct IsAmbiguous:
    public _IsAmbiguous<LargerType, SmallerType> {};
template <typename LargerType, typename SmallerType>
struct IsAmbiguous<LargerType const, SmallerType const>:
    public IsAmbiguous<LargerType, SmallerType> {};

//////////////////////////////////////////////////////////////////////////////
// Dna5 is ambiguous over Dna
template <>
struct IsAmbiguous<seqan::Dna5, seqan::Dna>
{
  typedef seqan::True Type;
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ambiguous_Dna5_2_Dna
{
  static char const VALUE[ 5 ][ 4 ];
};
template <typename T>
char const _Translate_Table_Ambiguous_Dna5_2_Dna<T>::VALUE[ 5 ][ 4 ] =
{
  { 0, 0, 0, 0 },
  { 1, 1, 1, 1 },
  { 2, 2, 2, 2 },
  { 3, 3, 3, 3 },
  { 0, 1, 2, 3 }
};
template <typename T = void>
struct _Size_Table_Ambiguous_Dna5_2_Dna
{
  static char const VALUE[ 5 ];
};
template <typename T>
char const _Size_Table_Ambiguous_Dna5_2_Dna<T>::VALUE[ 5 ] =
{
  1, 1, 1, 1, 4
};

//////////////////////////////////////////////////////////////////////////////

static
inline void ambiguousAssign ( seqan::Dna & c_target, seqan::Dna5 const & source, unsigned int const & which )
{
  SEQAN_CHECKPOINT
  c_target.value =
    galosh::_Translate_Table_Ambiguous_Dna5_2_Dna<>::VALUE[ source.value ][ which ];
}

static
inline size_t ambiguousCount ( seqan::Dna5 const & source, seqan::Dna const & )
{
  return
    galosh::_Size_Table_Ambiguous_Dna5_2_Dna<>::VALUE[ source.value ];
}

/*/////////////////////////////////////////////////////////////////////////////
// Iupac is ambiguous over Dna */
template <>
struct IsAmbiguous<seqan::Iupac, seqan::Dna>
{
  typedef seqan::True Type;
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ambiguous_Iupac_2_Dna
{
  static char const VALUE[ 16 ][ 4 ];
};
template <typename T>
char const _Translate_Table_Ambiguous_Iupac_2_Dna<T>::VALUE[ 16 ][ 4 ] =
{
  // 'U'=0, 'T', 'A', 'W', 'C', 'Y', 'M', 'H', 'G', 'K', 'R', 'D', 'S', 'B', 'V', 'N'=15.
  { 3, 3, 3, 3 }, // U
  { 3, 3, 3, 3 }, // T
  { 0, 0, 0, 0 }, // A
  { 0, 3, 0, 0 }, // W: A or T/U
  { 1, 1, 1, 1 }, // C
  { 1, 3, 1, 1 }, // Y: Pyrimidine (C or T/U)
  { 0, 1, 0, 0 }, // M: A or C
  { 0, 1, 3, 0 }, // H: A or C or T/U
  { 2, 2, 2, 2 }, // G
  { 2, 3, 2, 2 }, // K: G or T/U
  { 0, 2, 0, 0 }, // R: Purine (A or G)
  { 0, 2, 3, 0 }, // D: A or G or T/U
  { 1, 2, 1, 1 }, // S: C or G
  { 1, 2, 3, 1 }, // B: C or G or T/U
  { 0, 1, 2, 0 }, // V: A or C or G
  { 0, 1, 2, 3 }  // N: A or C or G or T/U
};
template <typename T = void>
struct _Size_Table_Ambiguous_Iupac_2_Dna
{
  static char const VALUE[ 16 ];
};
template <typename T>
char const _Size_Table_Ambiguous_Iupac_2_Dna<T>::VALUE[ 16 ] =
{
// U, T, A, W, C, Y, M, H, G, K, R, D, S, B, V, N.
   1, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4
};

//////////////////////////////////////////////////////////////////////////////

static
inline void ambiguousAssign ( seqan::Dna & c_target, seqan::Iupac const & source, unsigned int const & which )
{
  SEQAN_CHECKPOINT
  c_target.value =
    galosh::_Translate_Table_Ambiguous_Iupac_2_Dna<>::VALUE[ source.value ][ which ];
}

static
inline size_t ambiguousCount ( seqan::Iupac const & source, seqan::Dna const & )
{
  return
    galosh::_Size_Table_Ambiguous_Iupac_2_Dna<>::VALUE[ source.value ];
}

// mark
//////////////////////////////////////////////////////////////////////////////
// AminoAcid is ambiguous over AminoAcid20
template <>
struct IsAmbiguous<seqan::AminoAcid, seqan::AminoAcid20>
{
  typedef seqan::True Type;
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ambiguous_AminoAcid_2_AminoAcid20
{
  static char const VALUE[ 24 ][ 20 ];
};
template <typename T>
char const _Translate_Table_Ambiguous_AminoAcid_2_AminoAcid20<T>::VALUE[ 24 ][ 20 ] =
{
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // A
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, // R
  { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 }, // N
  { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 }, // D
  { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 }, // C
  { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 }, // Q
  { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }, // E
  { 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 }, // G
  { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 }, // H
  { 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 }, // I
  { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 }, // L
  { 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 }, // K
  { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 }, // M
  { 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13 }, // F
  { 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 }, // P
  { 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 }, // S
  { 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 }, // T
  { 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17 }, // W
  { 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 }, // Y
  { 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19 }, // V
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }, // (it's B, but we'll map it to all)
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }, // (it's Z, but we'll map it to all)
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }, // (it's X, but we'll map it to all)
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 } // (it's *, but we'll map it to all)
};
template <typename T = void>
struct _Size_Table_Ambiguous_AminoAcid_2_AminoAcid20
{
  static char const VALUE[ 24 ];
};
template <typename T>
char const _Size_Table_Ambiguous_AminoAcid_2_AminoAcid20<T>::VALUE[ 24 ] =
{
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 20, 20, 20, 20
};

//////////////////////////////////////////////////////////////////////////////

static
inline void ambiguousAssign ( seqan::AminoAcid20 & c_target, seqan::AminoAcid const & source, unsigned int const & which )
{
  SEQAN_CHECKPOINT
  c_target.value =
    galosh::_Translate_Table_Ambiguous_AminoAcid_2_AminoAcid20<>::VALUE[ source.value ][ which ];
}

static
inline size_t ambiguousCount ( seqan::AminoAcid const & source, seqan::AminoAcid20 const & )
{
  return
    galosh::_Size_Table_Ambiguous_AminoAcid_2_AminoAcid20<>::VALUE[ source.value ];
}

// endmark

} // End namespace galosh

#endif // __GALOSH_AMBIGUOUS_HPP__
