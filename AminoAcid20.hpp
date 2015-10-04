/**
 * \file AminoAcid20.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      A 20-char alternative to seqan::AminoAcid (which includes an extra 4 chars).
 * \par Overview:
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

#ifndef __GALOSH_AMINOACID20_HPP__
#define __GALOSH_AMINOACID20_HPP__

#include <seqan/basic.h>
#include "Prolific.hpp"

namespace seqan {
// Copied and modified (Oct 3, 2015) from latest git update of SeqAn (GNU Lesser General Public License v3)

// --------------------------------------------------------------------------
// Amino Acid
// --------------------------------------------------------------------------

template <typename T = void>
struct TranslateTableAminoAcid20ToChar_
{
    static char const VALUE[20];
};
template <typename T>
char const TranslateTableAminoAcid20ToChar_<T>::VALUE[20] =
{
    'A', // Ala Alanine
    'C', // Cys Cystine
    'D', // Asp Aspartic Acid
    'E', // Glu Glutamic Acid
    'F', // Phe Phenylalanine
    'G', // Gly Glycine
    'H', // His Histidine
    'I', // Ile Isoleucine
    'K', // Lys Lysine
    'L', // Leu Leucine
    'M', // Met Methionine
    'N', // Asn Asparagine
    'P', // Pro Proline
    'Q', // Gln Glutamine
    'R', // Arg Arginine
    'S', // Ser Serine
    'T', // Thr Threonine
    'V', // Val Valine
    'W', // Trp Tryptophan
    'Y' // Tyr Tyrosine
};

template <typename T = void>
struct TranslateTableCharToAminoAcid20_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableCharToAminoAcid20_<T>::VALUE[256] =
{
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //0
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //1
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //2
//                                                     *
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //3
    0,   0,   0,   1,   2,   3,   4,   5,   6,   7,   0,  8,  9,  10,  11,  0, //4
//    ,   A,   ,    C,   D,   E,   F,   G,   H,   I,    ,   K,   L,   M,   N,   ,

    12,  13,  14,  15,  16,  0,  17,  18,  0,  19,  0,  0,  0,  0,  0,  0, //5
//   P,   Q,   R,   S,   T,    ,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    ,

    0,   0,   0,   1,   2,   3,   4,   5,   6,   7,   0,  8,  9,  10,  11,  0, //6
//    ,   a,    ,   c,   d,   e,   f,   g,   h,   i,    ,   k,   l,   m,   n,   ,

    12,  13,  14,  15,  16,  0,  17,  18,  0,  19,  0,  0,  0,  0,  0,  0, //7
//   p,   q,   r,   s,   t,   ,   v,   w,   ,  y,   ,    ,    ,    ,    ,    ,

    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //8
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //9
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //10
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //11
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //12
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //13
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //14
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  //15
};

template <typename T = void>
struct TranslateTableByteToAminoAcid20_
{
    static char const VALUE[256];
};

template <typename T>
char const TranslateTableByteToAminoAcid20_<T>::VALUE[256] =
{
    0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,  13,  14,  15, //0
    16,  17,  18,  19,  20,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //1
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //2
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //3
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //4
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //5
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //6
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //7
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //8
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //9
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //10
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //11
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //12
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //13
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, //14
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0  //15
};
  
// ----------------------------------------------------------------------------
// Specialization AminoAcid20
// ----------------------------------------------------------------------------

/*!
 * @class AminoAcid20
 * @extends SimpleType
 * @headerfile <seqan/basic.h>
 * @brief IUPAC code for amino acids.
 * @signature typedef SingleType<unsigned char, AminoAcid20_> AminoAcid20;
 *
 * The ValueSize of <tt>AminoAcid20</tt> is 20.
 *
 * The amino acid symbols are as follows, i.e. they are sorted alphabetically
 * up until the last two symbols:
 *
 * 'A' = 0, 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'
 *
 * These are the 20 AAs, differing from the AminoAcid type in that it does _not_ include the wildcards, rare AAs, or terminator.
 * 'B' is a wildcard for (Aspartic Acid, Asparagine),
 * 'J' for (Leucine, Isoleucine), 'Z' for (Glutamic Acid, Glutamine) and
 * 'X' for "any amino acid".
 *
 * 'O' refers to the rare Pyrrolysine, 'U' refers to the rare Selenocysteine and '*' to the terminator tRNA.
 *
 * Objects of type <tt>AminoAcid20</tt> can be converted to <tt>char</tt> and vice versa.  Unknown values are converted to
 * <tt>'A'</tt>.
 *
 * @see FiniteOrderedAlphabetConcept#ValueSize
 * @see PeptideIterator
 * @see Peptide
 */

struct AminoAcid20_ {};
typedef SimpleType<unsigned char, AminoAcid20_> AminoAcid20;

template <> struct ValueSize<AminoAcid20>
{
    typedef __uint8 Type;
    static const Type VALUE = 20;
};

template <> struct BitsPerValue<AminoAcid20>
{
    typedef __uint8 Type;
    static const Type VALUE = 5;
};

inline AminoAcid20
unknownValueImpl(AminoAcid20 *)
{
    static const AminoAcid20 _result = AminoAcid20('A');
    return _result;
}

inline void assign(char & c_target, AminoAcid20 const & source)
{
    c_target = TranslateTableAminoAcid20ToChar_<>::VALUE[source.value];
}

template <>
struct CompareTypeImpl<AminoAcid20, __uint8>
{
    typedef AminoAcid20 Type;
};

inline void assign(AminoAcid20 & target, __uint8 c_source)
{
    target.value = TranslateTableByteToAminoAcid20_<>::VALUE[c_source];
}

template <>
struct CompareTypeImpl<AminoAcid20, char>
{
    typedef AminoAcid20 Type;
};

inline void assign(AminoAcid20 & target, char c_source)
{
    target.value = TranslateTableCharToAminoAcid20_<>::VALUE[(unsigned char) c_source];
}


//////////////////////////////////////////////////////////////////////////////
} // End namespace seqan

#endif // __GALOSH_AMINOACID20_HPP__
