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
 * \copyright &copy; 2012 by Paul T. Edlefsen, Fred Hutchinson Cancer
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
  // Copied and modified (Nov 21, 2012) from latest svn update of SeqAn (GNU Lesser General Public License v3): seqan-trunk/core/include/seqan/basic/alphabet_residue.h 

// ----------------------------------------------------------------------------
// Specialization AminoAcid20
// ----------------------------------------------------------------------------

/**
.Spec.AminoAcid20:
..cat:Alphabets
..summary:Iupac code for amino acids with exactly 20 categories.
..general:Class.SimpleType
..signature:AminoAcid
..remarks:
...text:The @Metafunction.ValueSize@ of $AminoAcid$ is 20. 
...text:The amino acids are enumerated from 0 to 19 in this order: 
...text:'A'=0, 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'=19.
...text:Objects of type $AminoAcid$ can be converted to $char$ and vice versa. 
Unknown values are converted to $'A'$.
...text:$AminoAcid$ is typedef for $SimpleType<char,AminoAcid20_>$, while $AminoAcid20_$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
..include:AminoAcid20.h
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


//____________________________________________________________________________
inline void assign(char & c_target, AminoAcid20 const & source)
{
    c_target = TranslateTableAAToAscii_<>::VALUE[source.value];
}

//____________________________________________________________________________
// ---------------------------------------------------------------------------
// Amino Acid 20
// ---------------------------------------------------------------------------

template <>
struct CompareType<AminoAcid20, __uint8>
{
    typedef AminoAcid20 Type;
};

inline void assign(AminoAcid20 & target, __uint8 c_source)
{
    target.value = TranslateTableByteToAA_<>::VALUE[c_source];
}

template <>
struct CompareType<AminoAcid20, char>
{
    typedef AminoAcid20 Type;
};

inline void assign(AminoAcid20 & target, char c_source)
{
    target.value = TranslateTableAsciiToAA_<>::VALUE[(unsigned char) c_source];
}

template <>
struct CompareType<AminoAcid20, Unicode>
{
    typedef AminoAcid20 Type;
};

inline void assign(AminoAcid20 & target, Unicode c_source)
{
    target.value = TranslateTableAsciiToAA_<>::VALUE[(unsigned char) c_source];
}

//////////////////////////////////////////////////////////////////////////////
} // End namespace seqan

#endif // __GALOSH_AMINOACID20_HPP__
