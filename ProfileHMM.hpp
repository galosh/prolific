/**
 * \file ProfileHMM.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      Class definition for the Galosh Profile HMM types, based on the
 *      "Plan 7" model used by HMMer and SAM.
 * \details
 *      More about Plan7 can be found at
 *      http://www.csb.yale.edu/userguides/seq/hmmer/docs/node11.html
 *      (or any other HMMER docs mirror, node 11)
 * \par Overview:
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

#ifndef __GALOSH_PROFILE_HMM_HPP__
#define __GALOSH_PROFILE_HMM_HPP__

#include "Prolific.hpp"

#include <seqan/basic.h>

namespace galosh {

/* ///////////////////////////////////////////////////////////////////////////// */

/**
.Spec.StateLabel:
..cat:Galosh.Profile
..summary:State labels of the Profile.  Note that the states of the profile are (label, position) pairs.
..general:Class.SimpleType
..signature:StateLabel
..remarks:
...text:The @Metafunction.ValueSize@ of $StateLabel$ is 12.
The state labels are enumerated this way: $'S' = 0, 'N' = 1, 'B' = 2, 'M' = 3, 'I' = 4, 'D' = 5, 'E' = 6, 'J' = 7, 'C' = 8, 'T' = 9, 'Z' = 10, 'W' = 11$.
...text:Objects of type $StateLabel$ can be converted to various other types and vice versa. 
An object that has a value not in ${'S', 'N', 'B', 'M', 'I', 'D', 'E', 'J', 'C', 'T', 'Z', 'W'}$ is converted to $'T'$.
...text:$StateLabel$ is typedef for $SimpleType<unsigned char,_StateLabel>$, while $_StateLables$ is a helper
specialization tag class.
..see:Metafunction.ValueSize
*/
struct _StateLabel {};
typedef seqan::SimpleType<unsigned char,_StateLabel> StateLabel;

/* ///////////////////////////////////////////////////////////////////////////// */

template <typename T = void>
struct _Translate_Table_StateLabel_2_Ascii
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_Ascii<T>::VALUE[ 12 ] =
  {'S', 'N', 'B', 'M', 'I', 'D', 'E', 'J', 'C', 'T', 'Z', 'W'};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_StateLabel
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_StateLabel<T>::VALUE[ 256 ] = 
{
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //0
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //1
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //2
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //3

	9,   9,   2,   8,   5,   6,   9,   9,   9,   4,   7,   9,   9,   3,   1,   9, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	9,   9,   9,   0,   9,   9,   9,   11,  9,   9,   10,  9,   9,   9,   9,   9, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	9,   9,   2,   8,   5,   6,   9,   9,   9,   4,   7,   9,   9,   3,   1,   9, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	9,   9,   9,   0,   9,   9,   9,   11,  9,   9,   10,  9,   9,   9,   9,   9, //7
//      p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //8
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //9
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //10
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //11
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //12
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //13
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //14
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_StateLabel
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_StateLabel<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  9,   9,   9,   9, //0
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //1
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //2
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //3
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //4
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //5
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //6
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //7
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //8
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //9
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //10
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //11
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //12
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //13
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9, //14
	9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9,   9  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::StateLabel > { enum { VALUE = 12 }; };
template <> struct BitsPerValue< galosh::StateLabel > { enum { VALUE = 4 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::StateLabel const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_StateLabel_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//StateLabel (10 labels)

template <>
struct CompareType<galosh::StateLabel, Byte> {
  typedef galosh::StateLabel Type;
};
inline void assign ( galosh::StateLabel & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_StateLabel<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::StateLabel, Ascii> {
  typedef galosh::StateLabel Type;
};
inline void assign ( galosh::StateLabel & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_StateLabel<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::StateLabel, Unicode> {
  typedef galosh::StateLabel Type;
};
inline void assign ( galosh::StateLabel & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_StateLabel<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//____________________________________________________________________________
/**
.Spec.Plan7:
..cat:Galosh.Profile
..summary:specialization tag indicating the "Plan7" architecture.
..signature:Plan7
..remarks:
..see:Spec.Plan9
*/
struct Plan7 {};

//____________________________________________________________________________
/**
.Spec.Plan9:
..cat:Galosh.Profile
..summary:specialization tag indicating the "Plan9" architecture, which is just like Plan7 except that D->I and I->D transitions are allowed.
..signature:Plan9
..remarks:
..see:Spec.Plan7
*/
struct Plan9 {};

//////////////////////////////////////////////////////////////////////////////
// METAFUNCTIONS, ETC
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//StateLabelTransitionTargets
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.StateLabelTransitionTargets:
..summary:Type representing the transition targets of a State Label.
..signature:StateLabelTransitionTargets<TStateLabel, TSpec>::Type
..param.TStateLabel:A StateLabel tag.
..param.TSpec:A specialization tag, either Plan7 or Plan9
..returns.param.Type:A SimpleType representing the transition targets of $T$.
..see:Spec.StateLabel
*/
template <typename TStateLabel, typename TSpec>
struct StateLabelTransitionTargets {};

//////////////////////////////////////////////////////////////////////////////
//StateLabelId
//////////////////////////////////////////////////////////////////////////////
/**
.Metafunction.StateLabelId:
..summary:byte representation of a state label.
..signature:StateLabelId<T>::VALUE
..param.T:A StateLabel tag.
..returns.param.VALUE:State Label Id of $T$.
...default:$9$
..see:Spec.StateLabel
*/
template <typename TValue>
struct StateLabelId
{
  enum { VALUE = 9 };
};
template <typename TValue>
struct StateLabelId<TValue const>:
    public StateLabelId<TValue> {};

//////////////////////////////////////////////////////////////////////////////
//IsEmitting
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.IsEmitting:
..summary:Tests state label to be an emitting state.
..signature:IsEmitting<T>::Type
..param.T:State Label that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is an emitting state, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..remarks:An emitting state is a state from which residues are emitted, such as Match and Insert states.
..see:Spec.StateLabel
*/

template <typename T>
struct _IsEmitting {
  typedef seqan::False Type;
};

template <typename T>
struct IsEmitting:
    public _IsEmitting<T> {};
template <typename T>
struct IsEmitting<T const>:
    public IsEmitting<T> {};

//////////////////////////////////////////////////////////////////////////////
//IsAssociatedWithPosition
//////////////////////////////////////////////////////////////////////////////

/**
.Metafunction.IsAssociatedWithPosition:
..summary:Tests state label to determine if it is associated with an ancestral sequence position.
..signature:IsAssociatedWithPosition<T>::Type
..param.T:State Label that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is an emitting state, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..remarks:State labels that are associated with ancestral sequence positions increment the position of the state (which is a (label, position) pair) on transition-in.
..see:Spec.StateLabel
*/

template <typename T>
struct _IsAssociatedWithPosition {
  typedef seqan::False Type;
};

template <typename T>
struct IsAssociatedWithPosition:
    public _IsAssociatedWithPosition<T> {};
template <typename T>
struct IsAssociatedWithPosition<T const>:
    public IsAssociatedWithPosition<T> {};

//////////////////////////////////////////////////////////////////////////////
} // End namespace galosh

//////////////////////////////////////////////////////////////////////////////
// State Labels
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
// Start State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.StartStateLabel:
..cat:Galosh.Profile
..summary:Start state label.
..signature:StartStateLabel
..remarks:
...text:$StartStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $StartStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $StartStateLabel$ is 0. 
...text:The @Metafunction.IsEmitting@ of $StartStateLabel$ is false. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $StartStateLabel$ is false. 
..see:Metafunction.StateLabelId
*/
struct StartStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.StartStateLabel
template <> struct StateLabelId< galosh::StartStateLabel > { enum { VALUE = 0 }; };

///.Metafunction.IsEmitting.param.T.type:Class.StartStateLabel
template <>
struct IsEmitting< galosh::StartStateLabel >
{
  typedef seqan::False Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.StartStateLabel
template <>
struct IsAssociatedWithPosition< galosh::StartStateLabel >
{
  typedef seqan::False Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.StartStateLabel
template <>
struct IsSimple< galosh::StartStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.StartStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Start state.
..general:Class.SimpleType
..signature:StartStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $StartStateTransitionTargets$ is 1. 
The start state transition targets are enumerated this way: $'N' = 0$.
...text:Objects of type $StartStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'N'}$ is converted to $'N'$.
...text:$StartStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_StartStateTransitionTargets>$, while $_StartStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _StartStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_StartStateTransitionTargets> StartStateTransitionTargets;

template <typename TSpec>
struct StateLabelTransitionTargets<StartStateLabel, TSpec> {
  typedef StartStateTransitionTargets Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_StartStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 1 ];
};
template <typename T>
char const _Translate_Table_StartStateTransitionTargets_2_Ascii<T>::VALUE[ 1 ] =
  { 'N' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_StartStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_StartStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StartStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 1 ];
};
template <typename T>
char const _Translate_Table_StartStateTransitionTargets_2_StateLabel<T>::VALUE[ 1 ] =
  { 1 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_StartStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_StartStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_StartStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_StartStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::StartStateTransitionTargets > { enum { VALUE = 1 }; };
template <> struct BitsPerValue< galosh::StartStateTransitionTargets > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::StartStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_StartStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//StartStateTransitionTargets

template <>
struct CompareType<galosh::StartStateTransitionTargets, Byte> {
  typedef galosh::StartStateTransitionTargets Type;
};
inline void assign ( galosh::StartStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_StartStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::StartStateTransitionTargets, Ascii> {
  typedef galosh::StartStateTransitionTargets Type;
};
inline void assign ( galosh::StartStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_StartStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::StartStateTransitionTargets, Unicode> {
  typedef galosh::StartStateTransitionTargets Type;
};
inline void assign ( galosh::StartStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_StartStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// PreAlign State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.PreAlignStateLabel:
..cat:Galosh.Profile
..summary:PreAlign state label.
..signature:PreAlignStateLabel
..remarks:
...text:$PreAlignStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $PreAlignStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $PreAlignStateLabel$ is 1. 
...text:The @Metafunction.IsEmitting@ of $PreAlignStateLabel$ is true. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $PreAlignStateLabel$ is false. 
..see:Metafunction.StateLabelId
*/
struct PreAlignStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.PreAlignStateLabel
template <> struct StateLabelId< galosh::PreAlignStateLabel > { enum { VALUE = 1 }; };

///.Metafunction.IsEmitting.param.T.type:Class.PreAlignStateLabel
template <>
struct IsEmitting< galosh::PreAlignStateLabel >
{
  typedef seqan::True Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.PreAlignStateLabel
template <>
struct IsAssociatedWithPosition< galosh::PreAlignStateLabel >
{
  typedef seqan::False Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.PreAlignStateLabel
template <>
struct IsSimple< galosh::PreAlignStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.PreAlignStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the PreAlign state.
..general:Class.SimpleType
..signature:PreAlignStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $PreAlignStateTransitionTargets$ is 2. 
The PreAlign state transition targets are enumerated this way: $'N' = 0, 'B' = 1$.
...text:Objects of type $PreAlignStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'N', 'B'}$ is converted to $'N'$.
...text:$PreAlignStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_PreAlignStateTransitionTargets>$, while $_PreAlignStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _PreAlignStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_PreAlignStateTransitionTargets> PreAlignStateTransitionTargets;

template <typename TSpec>
struct StateLabelTransitionTargets<PreAlignStateLabel, TSpec> {
  typedef PreAlignStateTransitionTargets Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_PreAlignStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_PreAlignStateTransitionTargets_2_Ascii<T>::VALUE[ 2 ] =
  { 'N', 'B' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_PreAlignStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_PreAlignStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_PreAlignStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_PreAlignStateTransitionTargets_2_StateLabel<T>::VALUE[ 2 ] =
  { 1, 2 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_PreAlignStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_PreAlignStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_PreAlignStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_PreAlignStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::PreAlignStateTransitionTargets > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::PreAlignStateTransitionTargets > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::PreAlignStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_PreAlignStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//PreAlignStateTransitionTargets

template <>
struct CompareType<galosh::PreAlignStateTransitionTargets, Byte> {
  typedef galosh::PreAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PreAlignStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_PreAlignStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::PreAlignStateTransitionTargets, Ascii> {
  typedef galosh::PreAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PreAlignStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_PreAlignStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::PreAlignStateTransitionTargets, Unicode> {
  typedef galosh::PreAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PreAlignStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_PreAlignStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// Begin State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.BeginStateLabel:
..cat:Galosh.Profile
..summary:Begin state label.
..signature:BeginStateLabel
..remarks:
...text:$BeginStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $BeginStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $BeginStateLabel$ is 2. 
...text:The @Metafunction.IsEmitting@ of $BeginStateLabel$ is false. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $BeginStateLabel$ is false. 
..see:Metafunction.StateLabelId
*/
struct BeginStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.BeginStateLabel
template <> struct StateLabelId< galosh::BeginStateLabel > { enum { VALUE = 2 }; };

///.Metafunction.IsEmitting.param.T.type:Class.BeginStateLabel
template <>
struct IsEmitting< galosh::BeginStateLabel >
{
  typedef seqan::False Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.BeginStateLabel
template <>
struct IsAssociatedWithPosition< galosh::BeginStateLabel >
{
  typedef seqan::False Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.BeginStateLabel
template <>
struct IsSimple< galosh::BeginStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BeginStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Begin state.
..general:Class.SimpleType
..signature:BeginStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $BeginStateTransitionTargets$ is 3.
The Begin state transition targets are enumerated this way: $'M' = 0, 'D' = 1, 'Z' = 2$.
...text:Objects of type $BeginStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'D', 'Z'}$ is converted to $'M'$.
...text:$BeginStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_BeginStateTransitionTargets>$, while $_BeginStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _BeginStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_BeginStateTransitionTargets> BeginStateTransitionTargets;

#ifdef USE_DEL_IN_DEL_OUT
template <typename TSpec>
struct StateLabelTransitionTargets<BeginStateLabel, TSpec> {
  typedef BeginStateTransitionTargets Type;
};
// #else replaced by BeginStateTransitionTargetsOld (below)
#endif // USE_DEL_IN_DEL_OUT

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_BeginStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_BeginStateTransitionTargets_2_Ascii<T>::VALUE[ 3 ] =
  { 'M', 'D', 'Z' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_BeginStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_BeginStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_BeginStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_BeginStateTransitionTargets_2_StateLabel<T>::VALUE[ 3 ] =
  { 3, 5, 10 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_BeginStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_BeginStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   2,   0,   0,   0,   0,   0, //7
//      p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_BeginStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_BeginStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::BeginStateTransitionTargets > { enum { VALUE = 3 }; };
template <> struct BitsPerValue< galosh::BeginStateTransitionTargets > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::BeginStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_BeginStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//BeginStateTransitionTargets

template <>
struct CompareType<galosh::BeginStateTransitionTargets, Byte> {
  typedef galosh::BeginStateTransitionTargets Type;
};
inline void assign ( galosh::BeginStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_BeginStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::BeginStateTransitionTargets, Ascii> {
  typedef galosh::BeginStateTransitionTargets Type;
};
inline void assign ( galosh::BeginStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_BeginStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::BeginStateTransitionTargets, Unicode> {
  typedef galosh::BeginStateTransitionTargets Type;
};
inline void assign ( galosh::BeginStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_BeginStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN


// OLD, for compatability with Profuse code (no transition to 'Z')
namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.BeginStateTransitionTargetsOld:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Begin state.
..general:Class.SimpleType
..signature:BeginStateTransitionTargetsOld
..remarks:
...text:The @Metafunction.ValueSize@ of $BeginStateTransitionTargetsOld$ is 2.
The Begin state transition targets are enumerated this way: $'M' = 0, 'D' = 1$.
...text:Objects of type $BeginStateTransitionTargetsOld$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'D'}$ is converted to $'M'$.
...text:$BeginStateTransitionTargetsOld$ is typedef for $SimpleType<unsigned char,_BeginStateTransitionTargetsOld>$, while $_BeginStateTransitionTargetsOld$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _BeginStateTransitionTargetsOld {};
typedef seqan::SimpleType<unsigned char,_BeginStateTransitionTargetsOld> BeginStateTransitionTargetsOld;

#ifndef USE_DEL_IN_DEL_OUT
template <typename TSpec>
struct StateLabelTransitionTargets<BeginStateLabel, TSpec> {
  typedef BeginStateTransitionTargetsOld Type;
};
#endif // !USE_DEL_IN_DEL_OUT

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_BeginStateTransitionTargetsOld_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_BeginStateTransitionTargetsOld_2_Ascii<T>::VALUE[ 2 ] =
  { 'M', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_BeginStateTransitionTargetsOld
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_BeginStateTransitionTargetsOld<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_BeginStateTransitionTargetsOld_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_BeginStateTransitionTargetsOld_2_StateLabel<T>::VALUE[ 2 ] =
  { 3, 5 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_BeginStateTransitionTargetsOld
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_BeginStateTransitionTargetsOld<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_BeginStateTransitionTargetsOld
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_BeginStateTransitionTargetsOld<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::BeginStateTransitionTargetsOld > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::BeginStateTransitionTargetsOld > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::BeginStateTransitionTargetsOld const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_BeginStateTransitionTargetsOld_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//BeginStateTransitionTargetsOld

template <>
struct CompareType<galosh::BeginStateTransitionTargetsOld, Byte> {
  typedef galosh::BeginStateTransitionTargetsOld Type;
};
inline void assign ( galosh::BeginStateTransitionTargetsOld & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_BeginStateTransitionTargetsOld<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::BeginStateTransitionTargetsOld, Ascii> {
  typedef galosh::BeginStateTransitionTargetsOld Type;
};
inline void assign ( galosh::BeginStateTransitionTargetsOld & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_BeginStateTransitionTargetsOld<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::BeginStateTransitionTargetsOld, Unicode> {
  typedef galosh::BeginStateTransitionTargetsOld Type;
};
inline void assign ( galosh::BeginStateTransitionTargetsOld & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_BeginStateTransitionTargetsOld<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// Match State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.MatchStateLabel:
..cat:Galosh.Profile
..summary:Match state label.
..signature:MatchStateLabel
..remarks:
...text:$MatchStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $MatchStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $MatchStateLabel$ is 3. 
...text:The @Metafunction.IsEmitting@ of $MatchStateLabel$ is true. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $MatchStateLabel$ is true.
..see:Metafunction.StateLabelId
*/
struct MatchStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.MatchStateLabel
template <> struct StateLabelId< galosh::MatchStateLabel > { enum { VALUE = 3 }; };

///.Metafunction.IsEmitting.param.T.type:Class.MatchStateLabel
template <>
struct IsEmitting< galosh::MatchStateLabel >
{
  typedef seqan::True Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.MatchStateLabel
template <>
struct IsAssociatedWithPosition< galosh::MatchStateLabel >
{
  typedef seqan::True Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.MatchStateLabel
template <>
struct IsSimple< galosh::MatchStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.MatchStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Match state.
..general:Class.SimpleType
..signature:MatchStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $MatchStateTransitionTargets$ is 4.
The Match state transition targets are enumerated this way: $'M' = 0, 'I' = 1, 'D' = 2, 'W' = 3$.
...text:Objects of type $MatchStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'I', 'D', 'W'}$ is converted to $'M'$.
...text:$MatchStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_MatchStateTransitionTargets>$, while $_MatchStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _MatchStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_MatchStateTransitionTargets> MatchStateTransitionTargets;

#ifdef USE_DEL_IN_DEL_OUT
template <typename TSpec>
struct StateLabelTransitionTargets<MatchStateLabel, TSpec> {
  typedef MatchStateTransitionTargets Type;
};
// #else replaced by MatchStateTransitionTargetsOld (below)
#endif // USE_DEL_IN_DEL_OUT

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_MatchStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const _Translate_Table_MatchStateTransitionTargets_2_Ascii<T>::VALUE[ 4 ] =
  { 'M', 'I', 'D', 'W' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_MatchStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_MatchStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_MatchStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const _Translate_Table_MatchStateTransitionTargets_2_StateLabel<T>::VALUE[ 4 ] =
  { 3, 4, 5, 11 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_MatchStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_MatchStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   3,   0,   0,   0,   0,   0,   0,   0,   0, //7
    //  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_MatchStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_MatchStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::MatchStateTransitionTargets > { enum { VALUE = 4 }; };
template <> struct BitsPerValue< galosh::MatchStateTransitionTargets > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::MatchStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_MatchStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//MatchStateTransitionTargets

template <>
struct CompareType<galosh::MatchStateTransitionTargets, Byte> {
  typedef galosh::MatchStateTransitionTargets Type;
};
inline void assign ( galosh::MatchStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_MatchStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::MatchStateTransitionTargets, Ascii> {
  typedef galosh::MatchStateTransitionTargets Type;
};
inline void assign ( galosh::MatchStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_MatchStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::MatchStateTransitionTargets, Unicode> {
  typedef galosh::MatchStateTransitionTargets Type;
};
inline void assign ( galosh::MatchStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_MatchStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

// OLD, for compatability with Profuse code (no transition to 'W')
namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.MatchStateTransitionTargetsOld:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Match state.
..general:Class.SimpleType
..signature:MatchStateTransitionTargetsOld
..remarks:
...text:The @Metafunction.ValueSize@ of $MatchStateTransitionTargetsOld$ is 3.
The Match state transition targets are enumerated this way: $'M' = 0, 'I' = 1, 'D' = 2$.
...text:Objects of type $MatchStateTransitionTargetsOld$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'I', 'D'}$ is converted to $'M'$.
...text:$MatchStateTransitionTargetsOld$ is typedef for $SimpleType<unsigned char,_MatchStateTransitionTargetsOld>$, while $_MatchStateTransitionTargetsOld$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _MatchStateTransitionTargetsOld {};
typedef seqan::SimpleType<unsigned char,_MatchStateTransitionTargetsOld> MatchStateTransitionTargetsOld;

#ifndef USE_DEL_IN_DEL_OUT
template <typename TSpec>
struct StateLabelTransitionTargets<MatchStateLabel, TSpec> {
  typedef MatchStateTransitionTargetsOld Type;
};
#endif // !USE_DEL_IN_DEL_OUT

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_MatchStateTransitionTargetsOld_2_Ascii
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_MatchStateTransitionTargetsOld_2_Ascii<T>::VALUE[ 3 ] =
  { 'M', 'I', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_MatchStateTransitionTargetsOld
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_MatchStateTransitionTargetsOld<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_MatchStateTransitionTargetsOld_2_StateLabel
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_MatchStateTransitionTargetsOld_2_StateLabel<T>::VALUE[ 3 ] =
  { 3, 4, 5 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_MatchStateTransitionTargetsOld
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_MatchStateTransitionTargetsOld<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_MatchStateTransitionTargetsOld
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_MatchStateTransitionTargetsOld<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::MatchStateTransitionTargetsOld > { enum { VALUE = 3 }; };
template <> struct BitsPerValue< galosh::MatchStateTransitionTargetsOld > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::MatchStateTransitionTargetsOld const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_MatchStateTransitionTargetsOld_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//MatchStateTransitionTargetsOld

template <>
struct CompareType<galosh::MatchStateTransitionTargetsOld, Byte> {
  typedef galosh::MatchStateTransitionTargetsOld Type;
};
inline void assign ( galosh::MatchStateTransitionTargetsOld & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_MatchStateTransitionTargetsOld<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::MatchStateTransitionTargetsOld, Ascii> {
  typedef galosh::MatchStateTransitionTargetsOld Type;
};
inline void assign ( galosh::MatchStateTransitionTargetsOld & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_MatchStateTransitionTargetsOld<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::MatchStateTransitionTargetsOld, Unicode> {
  typedef galosh::MatchStateTransitionTargetsOld Type;
};
inline void assign ( galosh::MatchStateTransitionTargetsOld & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_MatchStateTransitionTargetsOld<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN
// End OLD, for compatability with Profuse code (no transition to 'E')

//////////////////////////////////////////////////////////////////////////////
// Insertion State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.InsertionStateLabel:
..cat:Galosh.Profile
..summary:Insertion state label.
..signature:InsertionStateLabel
..remarks:
...text:$InsertionStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $InsertionStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $InsertionStateLabel$ is 4. 
...text:The @Metafunction.IsEmitting@ of $InsertionStateLabel$ is true. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $InsertionStateLabel$ is false.
..see:Metafunction.StateLabelId
*/
struct InsertionStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.InsertionStateLabel
template <> struct StateLabelId< galosh::InsertionStateLabel > { enum { VALUE = 4 }; };

///.Metafunction.IsEmitting.param.T.type:Class.InsertionStateLabel
template <>
struct IsEmitting< galosh::InsertionStateLabel >
{
  typedef seqan::True Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.InsertionStateLabel
template <>
struct IsAssociatedWithPosition< galosh::InsertionStateLabel >
{
  typedef seqan::False Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.InsertionStateLabel
template <>
struct IsSimple< galosh::InsertionStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.InsertionStateTransitionTargetsPlan9:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Insertion state, for the "Plan 9" architecture.
..general:Class.SimpleType
..signature:InsertionStateTransitionTargetsPlan9
..remarks:
...text:The @Metafunction.ValueSize@ of $InsertionStateTransitionTargetsPlan9$ is 3.
The Insertion state transition targets are enumerated this way: $'M' = 0, 'I' = 1, 'D' = 2$.
...text:Objects of type $InsertionStateTransitionTargetsPlan9$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'I', 'D'}$ is converted to $'M'$.
...text:$InsertionStateTransitionTargetsPlan9$ is typedef for $SimpleType<unsigned char,_InsertionStateTransitionTargetsPlan9>$, while $_InsertionStateTransitionTargetsPlan9$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _InsertionStateTransitionTargetsPlan9 {};
typedef seqan::SimpleType<unsigned char,_InsertionStateTransitionTargetsPlan9> InsertionStateTransitionTargetsPlan9;

template <>
struct StateLabelTransitionTargets<InsertionStateLabel, Plan9> {
  typedef InsertionStateTransitionTargetsPlan9 Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_InsertionStateTransitionTargetsPlan9_2_Ascii
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_InsertionStateTransitionTargetsPlan9_2_Ascii<T>::VALUE[ 3 ] =
  { 'M', 'I', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_InsertionStateTransitionTargetsPlan9
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_InsertionStateTransitionTargetsPlan9<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_InsertionStateTransitionTargetsPlan9_2_StateLabel
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_InsertionStateTransitionTargetsPlan9_2_StateLabel<T>::VALUE[ 3 ] =
  { 3, 4, 5 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan9
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan9<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_InsertionStateTransitionTargetsPlan9
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_InsertionStateTransitionTargetsPlan9<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::InsertionStateTransitionTargetsPlan9 > { enum { VALUE = 3 }; };
template <> struct BitsPerValue< galosh::InsertionStateTransitionTargetsPlan9 > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::InsertionStateTransitionTargetsPlan9 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_InsertionStateTransitionTargetsPlan9_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//InsertionStateTransitionTargetsPlan9

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan9, Byte> {
  typedef galosh::InsertionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan9 & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_InsertionStateTransitionTargetsPlan9<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan9, Ascii> {
  typedef galosh::InsertionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan9 & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan9<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan9, Unicode> {
  typedef galosh::InsertionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan9 & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan9<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.InsertionStateTransitionTargetsPlan7:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Insertion state, for the "Plan 7" architecture.
..general:Class.SimpleType
..signature:InsertionStateTransitionTargetsPlan7
..remarks:
...text:The @Metafunction.ValueSize@ of $InsertionStateTransitionTargetsPlan7$ is 2.
The Insertion state transition targets are enumerated this way: $'M' = 0, 'I' = 1$.
...text:Objects of type $InsertionStateTransitionTargetsPlan7$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'I'}$ is converted to $'M'$.
...text:$InsertionStateTransitionTargetsPlan7$ is typedef for $SimpleType<unsigned char,_InsertionStateTransitionTargetsPlan7>$, while $_InsertionStateTransitionTargetsPlan7$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _InsertionStateTransitionTargetsPlan7 {};
typedef seqan::SimpleType<unsigned char,_InsertionStateTransitionTargetsPlan7> InsertionStateTransitionTargetsPlan7;

template <>
struct StateLabelTransitionTargets<InsertionStateLabel, Plan7> {
  typedef InsertionStateTransitionTargetsPlan7 Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_InsertionStateTransitionTargetsPlan7_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_InsertionStateTransitionTargetsPlan7_2_Ascii<T>::VALUE[ 2 ] =
  { 'M', 'I' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_InsertionStateTransitionTargetsPlan7
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_InsertionStateTransitionTargetsPlan7<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_InsertionStateTransitionTargetsPlan7_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_InsertionStateTransitionTargetsPlan7_2_StateLabel<T>::VALUE[ 2 ] =
  { 3, 4 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan7
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan7<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_InsertionStateTransitionTargetsPlan7
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_InsertionStateTransitionTargetsPlan7<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::InsertionStateTransitionTargetsPlan7 > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::InsertionStateTransitionTargetsPlan7 > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::InsertionStateTransitionTargetsPlan7 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_InsertionStateTransitionTargetsPlan7_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//InsertionStateTransitionTargetsPlan7

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan7, Byte> {
  typedef galosh::InsertionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan7 & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_InsertionStateTransitionTargetsPlan7<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan7, Ascii> {
  typedef galosh::InsertionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan7 & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan7<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan7, Unicode> {
  typedef galosh::InsertionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan7 & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_InsertionStateTransitionTargetsPlan7<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// Deletion State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.DeletionStateLabel:
..cat:Galosh.Profile
..summary:Deletion state label.
..signature:DeletionStateLabel
..remarks:
...text:$DeletionStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $DeletionStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $DeletionStateLabel$ is 5. 
...text:The @Metafunction.IsEmitting@ of $DeletionStateLabel$ is false. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $DeletionStateLabel$ is true.
..see:Metafunction.StateLabelId
*/
struct DeletionStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.DeletionStateLabel
template <> struct StateLabelId< galosh::DeletionStateLabel > { enum { VALUE = 5 }; };

///.Metafunction.IsEmitting.param.T.type:Class.DeletionStateLabel
template <>
struct IsEmitting< galosh::DeletionStateLabel >
{
  typedef seqan::False Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.DeletionStateLabel
template <>
struct IsAssociatedWithPosition< galosh::DeletionStateLabel >
{
  typedef seqan::True Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.DeletionStateLabel
template <>
struct IsSimple< galosh::DeletionStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.DeletionStateTransitionTargetsPlan9:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Deletion state, for the "Plan 9" architecture.
..general:Class.SimpleType
..signature:DeletionStateTransitionTargetsPlan9
..remarks:
...text:The @Metafunction.ValueSize@ of $DeletionStateTransitionTargetsPlan9$ is 4.
The Deletion state transition targets are enumerated this way: $'M' = 0, 'I' = 1, 'D' = 2, 'E' = 3$.
...text:Objects of type $DeletionStateTransitionTargetsPlan9$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'I', 'D', 'E'}$ is converted to $'M'$.
...text:$DeletionStateTransitionTargetsPlan9$ is typedef for $SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan9>$, while $_DeletionStateTransitionTargetsPlan9$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _DeletionStateTransitionTargetsPlan9 {};
typedef seqan::SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan9> DeletionStateTransitionTargetsPlan9;

// Temporarily replaced by DeletionStateTransitionTargetsPlan9Old (below)
//template <>
//struct StateLabelTransitionTargets<DeletionStateLabel, Plan9> {
//  typedef DeletionStateTransitionTargetsPlan9 Type;
//};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan9_2_Ascii
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan9_2_Ascii<T>::VALUE[ 4 ] =
  { 'M', 'I', 'D', 'E' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan9
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan9<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan9_2_StateLabel
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan9_2_StateLabel<T>::VALUE[ 4 ] =
  { 3, 4, 5, 6 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   2,   3,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   2,   3,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan9
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan9<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan9 > { enum { VALUE = 4 }; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan9 > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::DeletionStateTransitionTargetsPlan9 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_DeletionStateTransitionTargetsPlan9_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan9

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9, Byte> {
  typedef galosh::DeletionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9 & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan9<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9, Ascii> {
  typedef galosh::DeletionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9 & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9, Unicode> {
  typedef galosh::DeletionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9 & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.DeletionStateTransitionTargetsPlan7:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Deletion state, for the "Plan 7" architecture.
..general:Class.SimpleType
..signature:DeletionStateTransitionTargetsPlan7
..remarks:
...text:The @Metafunction.ValueSize@ of $DeletionStateTransitionTargetsPlan7$ is 3.
The Deletion state transition targets are enumerated this way: $'M' = 0, 'D' = 1, 'E' = 2$.
...text:Objects of type $DeletionStateTransitionTargetsPlan7$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'D', 'E'}$ is converted to $'M'$.
...text:$DeletionStateTransitionTargetsPlan7$ is typedef for $SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan7>$, while $_DeletionStateTransitionTargetsPlan7$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _DeletionStateTransitionTargetsPlan7 {};
typedef seqan::SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan7> DeletionStateTransitionTargetsPlan7;

// Temporarily replaced by DeletionStateTransitionTargetsPlan7Old (below)
//template <>
//struct StateLabelTransitionTargets<DeletionStateLabel, Plan7> {
//  typedef DeletionStateTransitionTargetsPlan7 Type;
//};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan7_2_Ascii
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan7_2_Ascii<T>::VALUE[ 3 ] =
  { 'M', 'D', 'E' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan7
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan7<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan7_2_StateLabel
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan7_2_StateLabel<T>::VALUE[ 3 ] =
  { 3, 5, 6 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan7
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan7<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan7 > { enum { VALUE = 3 }; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan7 > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::DeletionStateTransitionTargetsPlan7 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_DeletionStateTransitionTargetsPlan7_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan7

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7, Byte> {
  typedef galosh::DeletionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7 & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan7<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7, Ascii> {
  typedef galosh::DeletionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7 & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7, Unicode> {
  typedef galosh::DeletionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7 & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

// OLD, for compatability with Profuse code (no transition to 'W')
namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.DeletionStateTransitionTargetsPlan9Old:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Deletion state, for the "Plan 9" architecture.
..general:Class.SimpleType
..signature:DeletionStateTransitionTargetsPlan9Old
..remarks:
...text:The @Metafunction.ValueSize@ of $DeletionStateTransitionTargetsPlan9Old$ is 3.
The Deletion state transition targets are enumerated this way: $'M' = 0, 'I' = 1, 'D' = 2$.
...text:Objects of type $DeletionStateTransitionTargetsPlan9Old$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'I', 'D'}$ is converted to $'M'$.
...text:$DeletionStateTransitionTargetsPlan9Old$ is typedef for $SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan9Old>$, while $_DeletionStateTransitionTargetsPlan9Old$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _DeletionStateTransitionTargetsPlan9Old {};
typedef seqan::SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan9Old> DeletionStateTransitionTargetsPlan9Old;

template <>
struct StateLabelTransitionTargets<DeletionStateLabel, Plan9> {
  typedef DeletionStateTransitionTargetsPlan9Old Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan9Old_2_Ascii
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan9Old_2_Ascii<T>::VALUE[ 3 ] =
  { 'M', 'I', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan9Old
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan9Old<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan9Old_2_StateLabel
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan9Old_2_StateLabel<T>::VALUE[ 3 ] =
  { 3, 4, 5 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9Old
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9Old<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   2,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan9Old
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan9Old<T>::VALUE[ 256 ] = 
{
	0,   1,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan9Old > { enum { VALUE = 3 }; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan9Old > { enum { VALUE = 2 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::DeletionStateTransitionTargetsPlan9Old const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_DeletionStateTransitionTargetsPlan9Old_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan9Old

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9Old, Byte> {
  typedef galosh::DeletionStateTransitionTargetsPlan9Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9Old & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan9Old<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9Old, Ascii> {
  typedef galosh::DeletionStateTransitionTargetsPlan9Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9Old & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9Old<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9Old, Unicode> {
  typedef galosh::DeletionStateTransitionTargetsPlan9Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9Old & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan9Old<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.DeletionStateTransitionTargetsPlan7Old:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Deletion state, for the "Plan 7" architecture.
..general:Class.SimpleType
..signature:DeletionStateTransitionTargetsPlan7Old
..remarks:
...text:The @Metafunction.ValueSize@ of $DeletionStateTransitionTargetsPlan7Old$ is 2.
The Deletion state transition targets are enumerated this way: $'M' = 0, 'D' = 1$.
...text:Objects of type $DeletionStateTransitionTargetsPlan7Old$ can be converted to various other types and vice versa. 
An object that has a value not in ${'M', 'D'}$ is converted to $'M'$.
...text:$DeletionStateTransitionTargetsPlan7Old$ is typedef for $SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan7Old>$, while $_DeletionStateTransitionTargetsPlan7Old$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _DeletionStateTransitionTargetsPlan7Old {};
typedef seqan::SimpleType<unsigned char,_DeletionStateTransitionTargetsPlan7Old> DeletionStateTransitionTargetsPlan7Old;

template <>
struct StateLabelTransitionTargets<DeletionStateLabel, Plan7> {
  typedef DeletionStateTransitionTargetsPlan7Old Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan7Old_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan7Old_2_Ascii<T>::VALUE[ 2 ] =
  { 'M', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan7Old
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_DeletionStateTransitionTargetsPlan7Old<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_DeletionStateTransitionTargetsPlan7Old_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_DeletionStateTransitionTargetsPlan7Old_2_StateLabel<T>::VALUE[ 2 ] =
  { 3, 5 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7Old
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7Old<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan7Old
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan7Old<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan7Old > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan7Old > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::DeletionStateTransitionTargetsPlan7Old const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_DeletionStateTransitionTargetsPlan7Old_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan7Old

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7Old, Byte> {
  typedef galosh::DeletionStateTransitionTargetsPlan7Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7Old & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_DeletionStateTransitionTargetsPlan7Old<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7Old, Ascii> {
  typedef galosh::DeletionStateTransitionTargetsPlan7Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7Old & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7Old<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7Old, Unicode> {
  typedef galosh::DeletionStateTransitionTargetsPlan7Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7Old & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionStateTransitionTargetsPlan7Old<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN
// End OLD, for compatability with Profuse code (no transition to 'W')

//////////////////////////////////////////////////////////////////////////////
// End State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.EndStateLabel:
..cat:Galosh.Profile
..summary:End state label.
..signature:EndStateLabel
..remarks:
...text:$EndStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $EndStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $EndStateLabel$ is 6. 
...text:The @Metafunction.IsEmitting@ of $EndStateLabel$ is false. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $EndStateLabel$ is false. 
..see:Metafunction.StateLabelId
*/
struct EndStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.EndStateLabel
template <> struct StateLabelId< galosh::EndStateLabel > { enum { VALUE = 6 }; };

///.Metafunction.IsEmitting.param.T.type:Class.EndStateLabel
template <>
struct IsEmitting< galosh::EndStateLabel >
{
  typedef seqan::False Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.EndStateLabel
template <>
struct IsAssociatedWithPosition< galosh::EndStateLabel >
{
  typedef seqan::False Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.EndStateLabel
template <>
struct IsSimple< galosh::EndStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.EndStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the End state.
..general:Class.SimpleType
..signature:EndStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $EndStateTransitionTargets$ is 2.
The End state transition targets are enumerated this way: $'C' = 0, 'J' = 1$.
...text:Objects of type $EndStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'C', 'J'}$ is converted to $'C'$.
...text:$EndStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_EndStateTransitionTargets>$, while $_EndStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _EndStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_EndStateTransitionTargets> EndStateTransitionTargets;

template <typename TSpec>
struct StateLabelTransitionTargets<EndStateLabel, TSpec> {
  typedef EndStateTransitionTargets Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_EndStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_EndStateTransitionTargets_2_Ascii<T>::VALUE[ 2 ] =
  { 'C', 'J' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_EndStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_EndStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_EndStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_EndStateTransitionTargets_2_StateLabel<T>::VALUE[ 2 ] =
  { 7, 8 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_EndStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_EndStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_EndStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_EndStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::EndStateTransitionTargets > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::EndStateTransitionTargets > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::EndStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_EndStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//EndStateTransitionTargets

template <>
struct CompareType<galosh::EndStateTransitionTargets, Byte> {
  typedef galosh::EndStateTransitionTargets Type;
};
inline void assign ( galosh::EndStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_EndStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::EndStateTransitionTargets, Ascii> {
  typedef galosh::EndStateTransitionTargets Type;
};
inline void assign ( galosh::EndStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_EndStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::EndStateTransitionTargets, Unicode> {
  typedef galosh::EndStateTransitionTargets Type;
};
inline void assign ( galosh::EndStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_EndStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// Loop State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.LoopStateLabel:
..cat:Galosh.Profile
..summary:Loop state label.
..signature:LoopStateLabel
..remarks:
...text:$LoopStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $LoopStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $LoopStateLabel$ is 7. 
...text:The @Metafunction.IsEmitting@ of $LoopStateLabel$ is true. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $LoopStateLabel$ is false. 
..see:Metafunction.StateLabelId
*/
struct LoopStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.LoopStateLabel
template <> struct StateLabelId< galosh::LoopStateLabel > { enum { VALUE = 7 }; };

///.Metafunction.IsEmitting.param.T.type:Class.LoopStateLabel
template <>
struct IsEmitting< galosh::LoopStateLabel >
{
  typedef seqan::True Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.LoopStateLabel
template <>
struct IsAssociatedWithPosition< galosh::LoopStateLabel >
{
  typedef seqan::False Type;
};

} // Loop namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.LoopStateLabel
template <>
struct IsSimple< galosh::LoopStateLabel >
{
  typedef True Type;
};
} // Loop namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.LoopStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the Loop state.
..general:Class.SimpleType
..signature:LoopStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $LoopStateTransitionTargets$ is 2.
The Loop state transition targets are enumerated this way: $'J' = 0, 'B' = 1$.
...text:Objects of type $LoopStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'J', 'B'}$ is converted to $'B'$.
...text:$LoopStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_LoopStateTransitionTargets>$, while $_LoopStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _LoopStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_LoopStateTransitionTargets> LoopStateTransitionTargets;

template <typename TSpec>
struct StateLabelTransitionTargets<LoopStateLabel, TSpec> {
  typedef LoopStateTransitionTargets Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_LoopStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_LoopStateTransitionTargets_2_Ascii<T>::VALUE[ 2 ] =
  { 'J', 'B' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_LoopStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_LoopStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_LoopStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_LoopStateTransitionTargets_2_StateLabel<T>::VALUE[ 2 ] =
  { 7, 2 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_LoopStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_LoopStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_LoopStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_LoopStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // Loop namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::LoopStateTransitionTargets > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::LoopStateTransitionTargets > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::LoopStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_LoopStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//LoopStateTransitionTargets

template <>
struct CompareType<galosh::LoopStateTransitionTargets, Byte> {
  typedef galosh::LoopStateTransitionTargets Type;
};
inline void assign ( galosh::LoopStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_LoopStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::LoopStateTransitionTargets, Ascii> {
  typedef galosh::LoopStateTransitionTargets Type;
};
inline void assign ( galosh::LoopStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_LoopStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::LoopStateTransitionTargets, Unicode> {
  typedef galosh::LoopStateTransitionTargets Type;
};
inline void assign ( galosh::LoopStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_LoopStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // Loop namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// PostAlign State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.PostAlignStateLabel:
..cat:Galosh.Profile
..summary:PostAlign state label.
..signature:PostAlignStateLabel
..remarks:
...text:$PostAlignStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $PostAlignStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $PostAlignStateLabel$ is 8. 
...text:The @Metafunction.IsEmitting@ of $PostAlignStateLabel$ is true. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $PostAlignStateLabel$ is false. 
..see:Metafunction.StateLabelId
*/
struct PostAlignStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.PostAlignStateLabel
template <> struct StateLabelId< galosh::PostAlignStateLabel > { enum { VALUE = 8 }; };

///.Metafunction.IsEmitting.param.T.type:Class.PostAlignStateLabel
template <>
struct IsEmitting< galosh::PostAlignStateLabel >
{
  typedef seqan::True Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.PostAlignStateLabel
template <>
struct IsAssociatedWithPosition< galosh::PostAlignStateLabel >
{
  typedef seqan::False Type;
};

} // PostAlign namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.PostAlignStateLabel
template <>
struct IsSimple< galosh::PostAlignStateLabel >
{
  typedef True Type;
};
} // PostAlign namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.PostAlignStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the PostAlign state.
..general:Class.SimpleType
..signature:PostAlignStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $PostAlignStateTransitionTargets$ is 2.
The PostAlign state transition targets are enumerated this way: $'C' = 0, 'T' = 1$.
...text:Objects of type $PostAlignStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'C', 'T'}$ is converted to $'C'$.
...text:$PostAlignStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_PostAlignStateTransitionTargets>$, while $_PostAlignStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _PostAlignStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_PostAlignStateTransitionTargets> PostAlignStateTransitionTargets;

template <typename TSpec>
struct StateLabelTransitionTargets<PostAlignStateLabel, TSpec> {
  typedef PostAlignStateTransitionTargets Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_PostAlignStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_PostAlignStateTransitionTargets_2_Ascii<T>::VALUE[ 2 ] =
  { 'C', 'T' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_PostAlignStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_PostAlignStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_PostAlignStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_PostAlignStateTransitionTargets_2_StateLabel<T>::VALUE[ 2 ] =
  { 8, 9 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_PostAlignStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_PostAlignStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//      p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_PostAlignStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_PostAlignStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // PostAlign namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::PostAlignStateTransitionTargets > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::PostAlignStateTransitionTargets > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::PostAlignStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_PostAlignStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//PostAlignStateTransitionTargets

template <>
struct CompareType<galosh::PostAlignStateTransitionTargets, Byte> {
  typedef galosh::PostAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PostAlignStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_PostAlignStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::PostAlignStateTransitionTargets, Ascii> {
  typedef galosh::PostAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PostAlignStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_PostAlignStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::PostAlignStateTransitionTargets, Unicode> {
  typedef galosh::PostAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PostAlignStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_PostAlignStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // PostAlign namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// Terminal State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.TerminalStateLabel:
..cat:Galosh.Profile
..summary:Terminal state label.
..signature:TerminalStateLabel
..remarks:
...text:$TerminalStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $TerminalStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $TerminalStateLabel$ is 9. 
...text:The @Metafunction.IsEmitting@ of $TerminalStateLabel$ is false. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $TerminalStateLabel$ is false. 
..see:Metafunction.StateLabelId
*/
struct TerminalStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.TerminalStateLabel
template <> struct StateLabelId< galosh::TerminalStateLabel > { enum { VALUE = 9 }; };

///.Metafunction.IsEmitting.param.T.type:Class.TerminalStateLabel
template <>
struct IsEmitting< galosh::TerminalStateLabel >
{
  typedef seqan::False Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.TerminalStateLabel
template <>
struct IsAssociatedWithPosition< galosh::TerminalStateLabel >
{
  typedef seqan::False Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.TerminalStateLabel
template <>
struct IsSimple< galosh::TerminalStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// DeletionIn State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.DeletionInStateLabel:
..cat:Galosh.Profile
..summary:DeletionIn state label.
..signature:DeletionInStateLabel
..remarks:
...text:$DeletionInStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $DeletionInStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $DeletionInStateLabel$ is 10. 
...text:The @Metafunction.IsEmitting@ of $DeletionInStateLabel$ is false. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $DeletionInStateLabel$ is true. 
..see:Metafunction.StateLabelId
*/
struct DeletionInStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.DeletionInStateLabel
template <> struct StateLabelId< galosh::DeletionInStateLabel > { enum { VALUE = 10 }; };

///.Metafunction.IsEmitting.param.T.type:Class.DeletionInStateLabel
template <>
struct IsEmitting< galosh::DeletionInStateLabel >
{
  typedef seqan::False Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.DeletionInStateLabel
template <>
struct IsAssociatedWithPosition< galosh::DeletionInStateLabel >
{
  typedef seqan::True Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.DeletionInStateLabel
template <>
struct IsSimple< galosh::DeletionInStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.DeletionInStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the DeletionIn state.
..general:Class.SimpleType
..signature:DeletionInStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $DeletionInStateTransitionTargets$ is 2.
The DeletionIn state transition targets are enumerated this way: $'Z' = 0, 'M' = 1$.
...text:Objects of type $DeletionInStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'Z', 'M'}$ is converted to $'Z'$.
...text:$DeletionInStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_DeletionInStateTransitionTargets>$, while $_DeletionInStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _DeletionInStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_DeletionInStateTransitionTargets> DeletionInStateTransitionTargets;

template <typename TSpec>
struct StateLabelTransitionTargets<DeletionInStateLabel, TSpec> {
  typedef DeletionInStateTransitionTargets Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_DeletionInStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_DeletionInStateTransitionTargets_2_Ascii<T>::VALUE[ 2 ] =
  { 'Z', 'M' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_DeletionInStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_DeletionInStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_DeletionInStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_DeletionInStateTransitionTargets_2_StateLabel<T>::VALUE[ 2 ] =
  { 10, 3 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_DeletionInStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_DeletionInStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//  p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_DeletionInStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_DeletionInStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::DeletionInStateTransitionTargets > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::DeletionInStateTransitionTargets > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::DeletionInStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_DeletionInStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionInStateTransitionTargets

template <>
struct CompareType<galosh::DeletionInStateTransitionTargets, Byte> {
  typedef galosh::DeletionInStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionInStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_DeletionInStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionInStateTransitionTargets, Ascii> {
  typedef galosh::DeletionInStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionInStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionInStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionInStateTransitionTargets, Unicode> {
  typedef galosh::DeletionInStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionInStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionInStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

//////////////////////////////////////////////////////////////////////////////
// DeletionOut State Label
//////////////////////////////////////////////////////////////////////////////

namespace galosh {
/**
.Spec.DeletionOutStateLabel:
..cat:Galosh.Profile
..summary:DeletionOut state label.
..signature:DeletionOutStateLabel
..remarks:
...text:$DeletionOutStateLabel$ is a tag.
...text:The @Metafunction.IsSimple@ of $DeletionOutStateLabel$ is true. 
...text:The @Metafunction.StateLabelId@ of $DeletionOutStateLabel$ is 11. 
...text:The @Metafunction.IsEmitting@ of $DeletionOutStateLabel$ is false. 
...text:The @Metafunction.IsAssociatedWithPosition@ of $DeletionOutStateLabel$ is true. 
..see:Metafunction.StateLabelId
*/
struct DeletionOutStateLabel {};

///.Metafunction.StateLabelId.param.T.type:Class.DeletionOutStateLabel
template <> struct StateLabelId< galosh::DeletionOutStateLabel > { enum { VALUE = 11 }; };

///.Metafunction.IsEmitting.param.T.type:Class.DeletionOutStateLabel
template <>
struct IsEmitting< galosh::DeletionOutStateLabel >
{
  typedef seqan::False Type;
};

///.Metafunction.IsAssociatedWithPosition.param.T.type:Class.DeletionOutStateLabel
template <>
struct IsAssociatedWithPosition< galosh::DeletionOutStateLabel >
{
  typedef seqan::True Type;
};

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{
///.Metafunction.IsSimple.param.T.type:Class.DeletionOutStateLabel
template <>
struct IsSimple< galosh::DeletionOutStateLabel >
{
  typedef True Type;
};
} // End namespace SEQAN_NAMESPACE_MAIN

namespace galosh {

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.DeletionOutStateTransitionTargets:
..cat:Galosh.Profile
..summary:Valid targets of transitions from the DeletionOut state.
..general:Class.SimpleType
..signature:DeletionOutStateTransitionTargets
..remarks:
...text:The @Metafunction.ValueSize@ of $DeletionOutStateTransitionTargets$ is 2.
The DeletionOut state transition targets are enumerated this way: $'W' = 0, 'E' = 1$.
...text:Objects of type $DeletionOutStateTransitionTargets$ can be converted to various other types and vice versa. 
An object that has a value not in ${'W', 'E'}$ is converted to $'W'$.
...text:$DeletionOutStateTransitionTargets$ is typedef for $SimpleType<unsigned char,_DeletionOutStateTransitionTargets>$, while $_DeletionOutStateTransitionTargets$ is a helper
specialization tag class.
..see:Spec.StateLabel
..see:Metafunction.ValueSize
*/
struct _DeletionOutStateTransitionTargets {};
typedef seqan::SimpleType<unsigned char,_DeletionOutStateTransitionTargets> DeletionOutStateTransitionTargets;

template <typename TSpec>
struct StateLabelTransitionTargets<DeletionOutStateLabel, TSpec> {
  typedef DeletionOutStateTransitionTargets Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T = void>
struct _Translate_Table_DeletionOutStateTransitionTargets_2_Ascii
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_DeletionOutStateTransitionTargets_2_Ascii<T>::VALUE[ 2 ] =
  { 'W', 'E' };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_StateLabel_2_DeletionOutStateTransitionTargets
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const _Translate_Table_StateLabel_2_DeletionOutStateTransitionTargets<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_DeletionOutStateTransitionTargets_2_StateLabel
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const _Translate_Table_DeletionOutStateTransitionTargets_2_StateLabel<T>::VALUE[ 2 ] =
  { 11, 6 };

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Ascii_2_DeletionOutStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Ascii_2_DeletionOutStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3

	0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
//	 ,   A,   B,   C,   D,   E,   D,   G,   H,   I,   J,   K,   L,   M,   N,   O,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
//	P,   Q,   R,   S,   T,   U,   V,   W,   X,   Y,   Z,    ,    ,    ,    ,    

	0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
//       ,   a,   b,   c,   d,   e,   f,   g,   h,   i,   j,   k,   l,   m,   n,   o,

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
//      p,   q,   r,   s,   t,   u,   v,   w,   x,   y,   z,    ,    ,    ,    ,   

	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//____________________________________________________________________________

template <typename T = void>
struct _Translate_Table_Byte_2_DeletionOutStateTransitionTargets
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const _Translate_Table_Byte_2_DeletionOutStateTransitionTargets<T>::VALUE[ 256 ] = 
{
	0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //0
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //1
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //2
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //3
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //4
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //5
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //6
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //7
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //8
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //9
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //10
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //11
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //12
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //13
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0, //14
	0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0  //15
};

//////////////////////////////////////////////////////////////////////////////

} // End namespace galosh

namespace SEQAN_NAMESPACE_MAIN
{

template <> struct ValueSize< galosh::DeletionOutStateTransitionTargets > { enum { VALUE = 2 }; };
template <> struct BitsPerValue< galosh::DeletionOutStateTransitionTargets > { enum { VALUE = 1 }; };

//////////////////////////////////////////////////////////////////////////////
//ASCII

inline void assign ( Ascii & c_target, galosh::DeletionOutStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::_Translate_Table_DeletionOutStateTransitionTargets_2_Ascii<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionOutStateTransitionTargets

template <>
struct CompareType<galosh::DeletionOutStateTransitionTargets, Byte> {
  typedef galosh::DeletionOutStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionOutStateTransitionTargets & target, Byte c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Byte_2_DeletionOutStateTransitionTargets<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionOutStateTransitionTargets, Ascii> {
  typedef galosh::DeletionOutStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionOutStateTransitionTargets & target, Ascii c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionOutStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionOutStateTransitionTargets, Unicode> {
  typedef galosh::DeletionOutStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionOutStateTransitionTargets & target, Unicode c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::_Translate_Table_Ascii_2_DeletionOutStateTransitionTargets<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

#endif // __GALOSH_PROFILE_HMM_HPP__
