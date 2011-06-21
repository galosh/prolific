/*---------------------------------------------------------------------------##
##  Library:
##      galosh::prolific
##  File:
##      AminoAcid20.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      A 20-char alternative to seqan::AminoAcid (which includes a 21st char for
##      the stop codon).
##
#******************************************************************************
#*
#*    This file is part of prolific, a library of useful C++ classes for
#*    working with genomic sequence data and Profile HMMs.  Please see the
#*    document CITING, which should have been included with this file.  You may
#*    use at will, subject to the license (LGPL v3), but *please cite the
#*    relevant papers* in your documentation and publications associated with
#*    uses of this library.  Thank you!
#*
#*    Copyright (C) 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
#*    prolific is free software: you can redistribute it and/or modify it under
#*    the terms of the GNU Lesser Public License as published by the Free
#*    Software Foundation, either version 3 of the License, or (at your option)
#*    any later version.
#*
#*    prolific is distributed in the hope that it will be useful, but WITHOUT
#*    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#*    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for
#*    more details.
#*
#*    You should have received a copy of the GNU Lesser Public License along
#*    with prolific.  If not, see <http://www.gnu.org/licenses/>.
#*****************************************************************************/

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
