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
 * \copyright &copy; 2015 by Paul T. Edlefsen, Fred Hutchinson Cancer
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
#include <seqan/stream.h> // now required to avoid "use of undeclared identifier 'directionIterator'" error
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
  //              //TGCA=[hex]            // IUPAC
  { 3, 3, 3, 3 }, //0000=0 = or U         // U
  { 0, 0, 0, 0 }, //0001=1                // A
  { 1, 1, 1, 1 }, //0010=2                // C
  { 0, 1, 0, 0 }, //0011=3 AC             // M: A or C
  { 2, 2, 2, 2 }, //0100=4                // G
  { 0, 2, 0, 0 }, //0101=5 AG (purine)    // R: Purine (A or G)
  { 1, 2, 0, 0 }, //0110=6 CG             // S: C or G
  { 0, 1, 2, 0 }, //0111=7 non-T          // V: A or C or G
  { 3, 3, 3, 3 }, //1000=8                // T
  { 0, 3, 0, 0 }, //1001=9 TA             // W: A or T/U
  { 1, 3, 0, 0 }, //1010=A TC (pyrimidine)// Y: Pyrimidine (C or T/U)
  { 0, 1, 3, 0 }, //1011=B not-G          // H: A or C or T/U
  { 2, 3, 0, 0 }, //1100=C TG             // K: G or T/U
  { 0, 2, 3, 0 }, //1101=D not-C          // D: A or G or T/U
  { 1, 2, 3, 1 }, //1110=E non-A          // B: C or G or T/U
  { 0, 1, 2, 3 }  //1111=F any            // N: A or C or G or T/U
};
template <typename T = void>
struct _Size_Table_Ambiguous_Iupac_2_Dna
{
  static char const VALUE[ 16 ];
};
template <typename T>
char const _Size_Table_Ambiguous_Iupac_2_Dna<T>::VALUE[ 16 ] =
{
// U, A, C, M, G, R, S, V, T, W, Y, H, K, D, B, N.
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
  static char const VALUE[ 27 ][ 20 ];
};
template <typename T>
char const _Translate_Table_Ambiguous_AminoAcid_2_AminoAcid20<T>::VALUE[ 27 ][ 20 ] =
{
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // A
  { 2, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // (it's B = D or N)
  { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }, // C
  { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 }, // D
  { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 }, // E
  { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 }, // F
  { 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 }, // G
  { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }, // H
  { 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7 }, // I
  { 7, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // (it's J = L or I)
  { 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8 }, // K
  { 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9 }, // L
  { 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10 }, // M
  { 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11 }, // N
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }, // (it's O, but we'll map it to all)
  { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12 }, // P
  { 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13 }, // Q
  { 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14 }, // R
  { 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15 }, // S
  { 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16 }, // T
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }, // (it's U, but we'll map it to all)
  { 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17 }, // V
  { 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 }, // W
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 }, // (it's X, but we'll map it to all)
  { 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19 }, // Y
  { 3, 13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, // (it's Z = E or Q)
  { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 } // (it's *, but we'll map it to all)
};
template <typename T = void>
struct _Size_Table_Ambiguous_AminoAcid_2_AminoAcid20
{
  static char const VALUE[ 27 ];
};
template <typename T>
char const _Size_Table_Ambiguous_AminoAcid_2_AminoAcid20<T>::VALUE[ 27 ] =
{
  1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 20, 1, 1, 1, 1, 1, 20, 1, 1, 20, 1, 2, 20
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
