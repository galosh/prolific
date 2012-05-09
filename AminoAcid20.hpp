/**
 * \file AminoAcid20.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      A 20-char alternative to seqan::AminoAcid (which includes a 21st char for
 *      the stop codon).
 * \par Overview:
 *    This file is part of prolific, a library of useful C++ classes for
 *    working with genomic sequence data and Profile HMMs.  Please see the
 *    document CITING, which should have been included with this file.  You may
 *    use at will, subject to the license (Apache v2.0), but *please cite the
 *    relevant papers* in your documentation and publications associated with
 *    uses of this library.  Thank you!
 *
 * \copyright &copy; 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
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

// Copied and modified from SeqAn (GNU Lesser General Public License v3) : $Id: basic_alphabet_simple.h 4412 2009-06-08 16:21:28Z rausch@PCPOOL.MI.FU-BERLIN.DE $

//____________________________________________________________________________

/**
.Spec.AminoAcid20:
..cat:Alphabets
..summary:Iupac code for amino acids.
..general:Class.SimpleType
..signature:AminoAcid20
..remarks:
...text:The @Metafunction.ValueSize@ of $AminoAcid20$ is 20. 
...text:The amino acids are enumerated from 0 to 19 in this order: 
...text:'A'=0, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=19.
...text:Objects of type $AminoAcid20$ can be converted to $char$ and vice versa. 
Unkown values are converted to $'A'$.
...text:$AminoAcid20$ is typedef for $SimpleType<char,_AminoAcid20>$, while $_AminoAcid20$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct _AminoAcid20 {};
typedef SimpleType<unsigned char, _AminoAcid20> AminoAcid20;

template <> struct ValueSize< AminoAcid20 > { enum { VALUE = 20 }; };
template <> struct BitsPerValue< AminoAcid20 > { enum { VALUE = 5 }; };

inline void assign(Ascii & c_target, AminoAcid20 const & source)
{
SEQAN_CHECKPOINT
	c_target = _Translate_Table_AA_2_Ascii<>::VALUE[source.value];
}

//////////////////////////////////////////////////////////////////////////////
//Amino Acid (5 bits)

template <>
struct CompareType<AminoAcid20, Byte> { typedef AminoAcid20 Type; };
inline void assign(AminoAcid20 & target, Byte c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Byte_2_AA<>::VALUE[c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<AminoAcid20, Ascii> { typedef AminoAcid20 Type; };
inline void assign(AminoAcid20 & target, Ascii c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source];
}
//____________________________________________________________________________

template <>
struct CompareType<AminoAcid20, Unicode> { typedef AminoAcid20 Type; };
inline void assign(AminoAcid20 & target, Unicode c_source)
{
SEQAN_CHECKPOINT
	target.value = _Translate_Table_Ascii_2_AA<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////
} // End namespace seqan

#endif // __GALOSH_AMINOACID20_HPP__
