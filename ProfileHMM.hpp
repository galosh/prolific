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

typedef unsigned char __uint8;  ///TAH 9/13

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
struct TranslateTableStateLabelToChar_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToChar_<T>::VALUE[ 12 ] =
  {'S', 'N', 'B', 'M', 'I', 'D', 'E', 'J', 'C', 'T', 'Z', 'W'};

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToStateLabel_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToStateLabel_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToStateLabel_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToStateLabel_<T>::VALUE[ 256 ] = 
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

  template <> struct ValueSize< galosh::StateLabel > { typedef __uint8 Type; static const Type VALUE = 12; };
  template <> struct BitsPerValue< galosh::StateLabel > { typedef __uint8 Type; static const Type VALUE = 12; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::StateLabel const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableStateLabelToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//StateLabel (10 labels)

template <>
struct CompareType<galosh::StateLabel, __uint8> {
  typedef galosh::StateLabel Type;
};
inline void assign ( galosh::StateLabel & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToStateLabel_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::StateLabel, char> {
  typedef galosh::StateLabel Type;
};
inline void assign ( galosh::StateLabel & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToStateLabel_<>::VALUE[
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
  typedef __uint8 Type; static const Type VALUE = 9;
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
template <> struct StateLabelId< galosh::StartStateLabel > { typedef __uint8 Type; static const Type VALUE = 0; };

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
struct TranslateTableStartStateTransitionTargetsToChar_
{
  static char const VALUE[ 1 ];
};
template <typename T>
char const TranslateTableStartStateTransitionTargetsToChar_<T>::VALUE[ 1 ] =
  { 'N' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToStartStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToStartStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStartStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 1 ];
};
template <typename T>
char const TranslateTableStartStateTransitionTargetsToStateLabel_<T>::VALUE[ 1 ] =
  { 1 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToStartStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToStartStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToStartStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToStartStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

  template <> struct ValueSize< galosh::StartStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };
  template <> struct BitsPerValue< galosh::StartStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::StartStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableStartStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//StartStateTransitionTargets

template <>
struct CompareType<galosh::StartStateTransitionTargets, __uint8> {
  typedef galosh::StartStateTransitionTargets Type;
};
inline void assign ( galosh::StartStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToStartStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::StartStateTransitionTargets, char> {
  typedef galosh::StartStateTransitionTargets Type;
};
inline void assign ( galosh::StartStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToStartStateTransitionTargets_<>::VALUE[
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
template <> struct StateLabelId< galosh::PreAlignStateLabel > { typedef __uint8 Type; static const Type VALUE = 1; };

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
struct TranslateTablePreAlignStateTransitionTargetsToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTablePreAlignStateTransitionTargetsToChar_<T>::VALUE[ 2 ] =
  { 'N', 'B' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToPreAlignStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToPreAlignStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTablePreAlignStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTablePreAlignStateTransitionTargetsToStateLabel_<T>::VALUE[ 2 ] =
  { 1, 2 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToPreAlignStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToPreAlignStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToPreAlignStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToPreAlignStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::PreAlignStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::PreAlignStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::PreAlignStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTablePreAlignStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//PreAlignStateTransitionTargets

template <>
struct CompareType<galosh::PreAlignStateTransitionTargets, __uint8> {
  typedef galosh::PreAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PreAlignStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToPreAlignStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::PreAlignStateTransitionTargets, char> {
  typedef galosh::PreAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PreAlignStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToPreAlignStateTransitionTargets_<>::VALUE[
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
template <> struct StateLabelId< galosh::BeginStateLabel > { typedef __uint8 Type; static const Type VALUE = 2; };

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
struct TranslateTableBeginStateTransitionTargetsToChar_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableBeginStateTransitionTargetsToChar_<T>::VALUE[ 3 ] =
  { 'M', 'D', 'Z' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToBeginStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToBeginStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableBeginStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableBeginStateTransitionTargetsToStateLabel_<T>::VALUE[ 3 ] =
  { 3, 5, 10 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToBeginStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToBeginStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToBeginStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToBeginStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

  template <> struct ValueSize< galosh::BeginStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 3; };
  template <> struct BitsPerValue< galosh::BeginStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::BeginStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableBeginStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//BeginStateTransitionTargets

template <>
struct CompareType<galosh::BeginStateTransitionTargets, __uint8> {
  typedef galosh::BeginStateTransitionTargets Type;
};
inline void assign ( galosh::BeginStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToBeginStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::BeginStateTransitionTargets, char> {
  typedef galosh::BeginStateTransitionTargets Type;
};
inline void assign ( galosh::BeginStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToBeginStateTransitionTargets_<>::VALUE[
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
struct TranslateTableBeginStateTransitionTargetsOldToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableBeginStateTransitionTargetsOldToChar_<T>::VALUE[ 2 ] =
  { 'M', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToBeginStateTransitionTargetsOld_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToBeginStateTransitionTargetsOld_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableBeginStateTransitionTargetsOldToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableBeginStateTransitionTargetsOldToStateLabel_<T>::VALUE[ 2 ] =
  { 3, 5 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToBeginStateTransitionTargetsOld_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToBeginStateTransitionTargetsOld_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToBeginStateTransitionTargetsOld_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToBeginStateTransitionTargetsOld_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::BeginStateTransitionTargetsOld > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::BeginStateTransitionTargetsOld > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::BeginStateTransitionTargetsOld const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableBeginStateTransitionTargetsOldToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//BeginStateTransitionTargetsOld

template <>
struct CompareType<galosh::BeginStateTransitionTargetsOld, __uint8> {
  typedef galosh::BeginStateTransitionTargetsOld Type;
};
inline void assign ( galosh::BeginStateTransitionTargetsOld & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToBeginStateTransitionTargetsOld_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::BeginStateTransitionTargetsOld, char> {
  typedef galosh::BeginStateTransitionTargetsOld Type;
};
inline void assign ( galosh::BeginStateTransitionTargetsOld & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToBeginStateTransitionTargetsOld_<>::VALUE[
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
template <> struct StateLabelId< galosh::MatchStateLabel > { typedef __uint8 Type; static const Type VALUE = 3; };

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
struct TranslateTableMatchStateTransitionTargetsToChar_
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const TranslateTableMatchStateTransitionTargetsToChar_<T>::VALUE[ 4 ] =
  { 'M', 'I', 'D', 'W' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToMatchStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToMatchStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableMatchStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const TranslateTableMatchStateTransitionTargetsToStateLabel_<T>::VALUE[ 4 ] =
  { 3, 4, 5, 11 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToMatchStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToMatchStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToMatchStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToMatchStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::MatchStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 4; };
template <> struct BitsPerValue< galosh::MatchStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::MatchStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableMatchStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//MatchStateTransitionTargets

template <>
struct CompareType<galosh::MatchStateTransitionTargets, __uint8> {
  typedef galosh::MatchStateTransitionTargets Type;
};
inline void assign ( galosh::MatchStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToMatchStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::MatchStateTransitionTargets, char> {
  typedef galosh::MatchStateTransitionTargets Type;
};
inline void assign ( galosh::MatchStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToMatchStateTransitionTargets_<>::VALUE[
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
struct TranslateTableMatchStateTransitionTargetsOldToChar_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableMatchStateTransitionTargetsOldToChar_<T>::VALUE[ 3 ] =
  { 'M', 'I', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToMatchStateTransitionTargetsOld_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToMatchStateTransitionTargetsOld_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableMatchStateTransitionTargetsOldToStateLabel_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableMatchStateTransitionTargetsOldToStateLabel_<T>::VALUE[ 3 ] =
  { 3, 4, 5 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToMatchStateTransitionTargetsOld_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToMatchStateTransitionTargetsOld_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToMatchStateTransitionTargetsOld_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToMatchStateTransitionTargetsOld_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::MatchStateTransitionTargetsOld > { typedef __uint8 Type; static const Type VALUE = 3; };
template <> struct BitsPerValue< galosh::MatchStateTransitionTargetsOld > { typedef __uint8 Type; static const Type VALUE = 2; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::MatchStateTransitionTargetsOld const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableMatchStateTransitionTargetsOldToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//MatchStateTransitionTargetsOld

template <>
struct CompareType<galosh::MatchStateTransitionTargetsOld, __uint8> {
  typedef galosh::MatchStateTransitionTargetsOld Type;
};
inline void assign ( galosh::MatchStateTransitionTargetsOld & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToMatchStateTransitionTargetsOld_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::MatchStateTransitionTargetsOld, char> {
  typedef galosh::MatchStateTransitionTargetsOld Type;
};
inline void assign ( galosh::MatchStateTransitionTargetsOld & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToMatchStateTransitionTargetsOld_<>::VALUE[
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
template <> struct StateLabelId< galosh::InsertionStateLabel > { typedef __uint8 Type; static const Type VALUE = 4; };

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
struct TranslateTableInsertionStateTransitionTargetsPlan9ToChar_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableInsertionStateTransitionTargetsPlan9ToChar_<T>::VALUE[ 3 ] =
  { 'M', 'I', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToInsertionStateTransitionTargetsPlan9_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToInsertionStateTransitionTargetsPlan9_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableInsertionStateTransitionTargetsPlan9ToStateLabel_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableInsertionStateTransitionTargetsPlan9ToStateLabel_<T>::VALUE[ 3 ] =
  { 3, 4, 5 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToInsertionStateTransitionTargetsPlan9_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToInsertionStateTransitionTargetsPlan9_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToInsertionStateTransitionTargetsPlan9_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToInsertionStateTransitionTargetsPlan9_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::InsertionStateTransitionTargetsPlan9 > { typedef __uint8 Type; static const Type VALUE = 3; };
template <> struct BitsPerValue< galosh::InsertionStateTransitionTargetsPlan9 > { typedef __uint8 Type; static const Type VALUE = 2; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::InsertionStateTransitionTargetsPlan9 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableInsertionStateTransitionTargetsPlan9ToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//InsertionStateTransitionTargetsPlan9

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan9, __uint8> {
  typedef galosh::InsertionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan9 & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToInsertionStateTransitionTargetsPlan9_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan9, char> {
  typedef galosh::InsertionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan9 & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToInsertionStateTransitionTargetsPlan9_<>::VALUE[
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
struct TranslateTableInsertionStateTransitionTargetsPlan7ToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableInsertionStateTransitionTargetsPlan7ToChar_<T>::VALUE[ 2 ] =
  { 'M', 'I' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToInsertionStateTransitionTargetsPlan7_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToInsertionStateTransitionTargetsPlan7_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableInsertionStateTransitionTargetsPlan7ToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableInsertionStateTransitionTargetsPlan7ToStateLabel_<T>::VALUE[ 2 ] =
  { 3, 4 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToInsertionStateTransitionTargetsPlan7_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToInsertionStateTransitionTargetsPlan7_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToInsertionStateTransitionTargetsPlan7_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToInsertionStateTransitionTargetsPlan7_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::InsertionStateTransitionTargetsPlan7 > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::InsertionStateTransitionTargetsPlan7 > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::InsertionStateTransitionTargetsPlan7 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableInsertionStateTransitionTargetsPlan7ToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//InsertionStateTransitionTargetsPlan7

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan7, __uint8> {
  typedef galosh::InsertionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan7 & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToInsertionStateTransitionTargetsPlan7_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::InsertionStateTransitionTargetsPlan7, char> {
  typedef galosh::InsertionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::InsertionStateTransitionTargetsPlan7 & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToInsertionStateTransitionTargetsPlan7_<>::VALUE[
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
template <> struct StateLabelId< galosh::DeletionStateLabel > { typedef __uint8 Type; static const Type VALUE = 5; };

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
struct TranslateTableDeletionStateTransitionTargetsPlan9ToChar_
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan9ToChar_<T>::VALUE[ 4 ] =
  { 'M', 'I', 'D', 'E' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToDeletionStateTransitionTargetsPlan9_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToDeletionStateTransitionTargetsPlan9_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDeletionStateTransitionTargetsPlan9ToStateLabel_
{
  static char const VALUE[ 4 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan9ToStateLabel_<T>::VALUE[ 4 ] =
  { 3, 4, 5, 6 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToDeletionStateTransitionTargetsPlan9_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToDeletionStateTransitionTargetsPlan9_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToDeletionStateTransitionTargetsPlan9_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToDeletionStateTransitionTargetsPlan9_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan9 > { typedef __uint8 Type; static const Type VALUE = 4; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan9 > { typedef __uint8 Type; static const Type VALUE = 2; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::DeletionStateTransitionTargetsPlan9 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableDeletionStateTransitionTargetsPlan9ToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan9

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9, __uint8> {
  typedef galosh::DeletionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9 & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToDeletionStateTransitionTargetsPlan9_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9, char> {
  typedef galosh::DeletionStateTransitionTargetsPlan9 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9 & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToDeletionStateTransitionTargetsPlan9_<>::VALUE[
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
struct TranslateTableDeletionStateTransitionTargetsPlan7ToChar_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan7ToChar_<T>::VALUE[ 3 ] =
  { 'M', 'D', 'E' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToDeletionStateTransitionTargetsPlan7_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToDeletionStateTransitionTargetsPlan7_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDeletionStateTransitionTargetsPlan7ToStateLabel_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan7ToStateLabel_<T>::VALUE[ 3 ] =
  { 3, 5, 6 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToDeletionStateTransitionTargetsPlan7_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToDeletionStateTransitionTargetsPlan7_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToDeletionStateTransitionTargetsPlan7_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToDeletionStateTransitionTargetsPlan7_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan7 > { typedef __uint8 Type; static const Type VALUE = 3; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan7 > { typedef __uint8 Type; static const Type VALUE = 2; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::DeletionStateTransitionTargetsPlan7 const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableDeletionStateTransitionTargetsPlan7ToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan7

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7, __uint8> {
  typedef galosh::DeletionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7 & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToDeletionStateTransitionTargetsPlan7_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7, char> {
  typedef galosh::DeletionStateTransitionTargetsPlan7 Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7 & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToDeletionStateTransitionTargetsPlan7_<>::VALUE[
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
struct TranslateTableDeletionStateTransitionTargetsPlan9OldToChar_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan9OldToChar_<T>::VALUE[ 3 ] =
  { 'M', 'I', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToDeletionStateTransitionTargetsPlan9Old_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToDeletionStateTransitionTargetsPlan9Old_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDeletionStateTransitionTargetsPlan9OldToStateLabel_
{
  static char const VALUE[ 3 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan9OldToStateLabel_<T>::VALUE[ 3 ] =
  { 3, 4, 5 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToDeletionStateTransitionTargetsPlan9Old_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToDeletionStateTransitionTargetsPlan9Old_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToDeletionStateTransitionTargetsPlan9Old_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToDeletionStateTransitionTargetsPlan9Old_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan9Old > { typedef __uint8 Type; static const Type VALUE = 3; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan9Old > { typedef __uint8 Type; static const Type VALUE = 2; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::DeletionStateTransitionTargetsPlan9Old const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableDeletionStateTransitionTargetsPlan9OldToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan9Old

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9Old, __uint8> {
  typedef galosh::DeletionStateTransitionTargetsPlan9Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9Old & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToDeletionStateTransitionTargetsPlan9Old_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan9Old, char> {
  typedef galosh::DeletionStateTransitionTargetsPlan9Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan9Old & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToDeletionStateTransitionTargetsPlan9Old_<>::VALUE[
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
struct TranslateTableDeletionStateTransitionTargetsPlan7OldToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan7OldToChar_<T>::VALUE[ 2 ] =
  { 'M', 'D' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToDeletionStateTransitionTargetsPlan7Old_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToDeletionStateTransitionTargetsPlan7Old_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDeletionStateTransitionTargetsPlan7OldToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableDeletionStateTransitionTargetsPlan7OldToStateLabel_<T>::VALUE[ 2 ] =
  { 3, 5 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToDeletionStateTransitionTargetsPlan7Old_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToDeletionStateTransitionTargetsPlan7Old_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToDeletionStateTransitionTargetsPlan7Old_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToDeletionStateTransitionTargetsPlan7Old_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::DeletionStateTransitionTargetsPlan7Old > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::DeletionStateTransitionTargetsPlan7Old > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::DeletionStateTransitionTargetsPlan7Old const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableDeletionStateTransitionTargetsPlan7OldToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionStateTransitionTargetsPlan7Old

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7Old, __uint8> {
  typedef galosh::DeletionStateTransitionTargetsPlan7Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7Old & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToDeletionStateTransitionTargetsPlan7Old_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionStateTransitionTargetsPlan7Old, char> {
  typedef galosh::DeletionStateTransitionTargetsPlan7Old Type;
};
inline void assign ( galosh::DeletionStateTransitionTargetsPlan7Old & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToDeletionStateTransitionTargetsPlan7Old_<>::VALUE[
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
template <> struct StateLabelId< galosh::EndStateLabel > { typedef __uint8 Type; static const Type VALUE = 6; };

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
struct TranslateTableEndStateTransitionTargetsToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableEndStateTransitionTargetsToChar_<T>::VALUE[ 2 ] =
  { 'C', 'J' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToEndStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToEndStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableEndStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableEndStateTransitionTargetsToStateLabel_<T>::VALUE[ 2 ] =
  { 7, 8 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToEndStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToEndStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToEndStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToEndStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::EndStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::EndStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::EndStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableEndStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//EndStateTransitionTargets

template <>
struct CompareType<galosh::EndStateTransitionTargets, __uint8> {
  typedef galosh::EndStateTransitionTargets Type;
};
inline void assign ( galosh::EndStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToEndStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::EndStateTransitionTargets, char> {
  typedef galosh::EndStateTransitionTargets Type;
};
inline void assign ( galosh::EndStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToEndStateTransitionTargets_<>::VALUE[
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
template <> struct StateLabelId< galosh::LoopStateLabel > { typedef __uint8 Type; static const Type VALUE = 7; };

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
struct TranslateTableLoopStateTransitionTargetsToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableLoopStateTransitionTargetsToChar_<T>::VALUE[ 2 ] =
  { 'J', 'B' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToLoopStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToLoopStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableLoopStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableLoopStateTransitionTargetsToStateLabel_<T>::VALUE[ 2 ] =
  { 7, 2 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToLoopStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToLoopStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToLoopStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToLoopStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::LoopStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::LoopStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::LoopStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableLoopStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//LoopStateTransitionTargets

template <>
struct CompareType<galosh::LoopStateTransitionTargets, __uint8> {
  typedef galosh::LoopStateTransitionTargets Type;
};
inline void assign ( galosh::LoopStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToLoopStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::LoopStateTransitionTargets, char> {
  typedef galosh::LoopStateTransitionTargets Type;
};
inline void assign ( galosh::LoopStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToLoopStateTransitionTargets_<>::VALUE[
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
template <> struct StateLabelId< galosh::PostAlignStateLabel > { typedef __uint8 Type; static const Type VALUE = 8; };

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
struct TranslateTablePostAlignStateTransitionTargetsToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTablePostAlignStateTransitionTargetsToChar_<T>::VALUE[ 2 ] =
  { 'C', 'T' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToPostAlignStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToPostAlignStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTablePostAlignStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTablePostAlignStateTransitionTargetsToStateLabel_<T>::VALUE[ 2 ] =
  { 8, 9 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToPostAlignStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToPostAlignStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToPostAlignStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToPostAlignStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::PostAlignStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::PostAlignStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::PostAlignStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTablePostAlignStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//PostAlignStateTransitionTargets

template <>
struct CompareType<galosh::PostAlignStateTransitionTargets, __uint8> {
  typedef galosh::PostAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PostAlignStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToPostAlignStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::PostAlignStateTransitionTargets, char> {
  typedef galosh::PostAlignStateTransitionTargets Type;
};
inline void assign ( galosh::PostAlignStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToPostAlignStateTransitionTargets_<>::VALUE[
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
template <> struct StateLabelId< galosh::TerminalStateLabel > { typedef __uint8 Type; static const Type VALUE = 9; };

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
template <> struct StateLabelId< galosh::DeletionInStateLabel > { typedef __uint8 Type; static const Type VALUE = 10; };

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
struct TranslateTableDeletionInStateTransitionTargetsToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableDeletionInStateTransitionTargetsToChar_<T>::VALUE[ 2 ] =
  { 'Z', 'M' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToDeletionInStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToDeletionInStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDeletionInStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableDeletionInStateTransitionTargetsToStateLabel_<T>::VALUE[ 2 ] =
  { 10, 3 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToDeletionInStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToDeletionInStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToDeletionInStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToDeletionInStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::DeletionInStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::DeletionInStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::DeletionInStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableDeletionInStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionInStateTransitionTargets

template <>
struct CompareType<galosh::DeletionInStateTransitionTargets, __uint8> {
  typedef galosh::DeletionInStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionInStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToDeletionInStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionInStateTransitionTargets, char> {
  typedef galosh::DeletionInStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionInStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToDeletionInStateTransitionTargets_<>::VALUE[
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
template <> struct StateLabelId< galosh::DeletionOutStateLabel > { typedef __uint8 Type; static const Type VALUE = 11; };

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
struct TranslateTableDeletionOutStateTransitionTargetsToChar_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableDeletionOutStateTransitionTargetsToChar_<T>::VALUE[ 2 ] =
  { 'W', 'E' };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableStateLabelToDeletionOutStateTransitionTargets_
{
  static char const VALUE[ 12 ];
};
template <typename T>
char const TranslateTableStateLabelToDeletionOutStateTransitionTargets_<T>::VALUE[ 12 ] =
  { 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableDeletionOutStateTransitionTargetsToStateLabel_
{
  static char const VALUE[ 2 ];
};
template <typename T>
char const TranslateTableDeletionOutStateTransitionTargetsToStateLabel_<T>::VALUE[ 2 ] =
  { 11, 6 };

//____________________________________________________________________________

template <typename T = void>
struct TranslateTableCharToDeletionOutStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableCharToDeletionOutStateTransitionTargets_<T>::VALUE[ 256 ] = 
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
struct TranslateTableByteToDeletionOutStateTransitionTargets_
{
  static char const VALUE[ 256 ];
};
template <typename T>
char const TranslateTableByteToDeletionOutStateTransitionTargets_<T>::VALUE[ 256 ] = 
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

template <> struct ValueSize< galosh::DeletionOutStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 2; };
template <> struct BitsPerValue< galosh::DeletionOutStateTransitionTargets > { typedef __uint8 Type; static const Type VALUE = 1; };

//////////////////////////////////////////////////////////////////////////////
//Char

inline void assign ( char & c_target, galosh::DeletionOutStateTransitionTargets const & source )
{
SEQAN_CHECKPOINT
  c_target =
    galosh::TranslateTableDeletionOutStateTransitionTargetsToChar_<>::VALUE[ source.value ];
}

//////////////////////////////////////////////////////////////////////////////
//DeletionOutStateTransitionTargets

template <>
struct CompareType<galosh::DeletionOutStateTransitionTargets, __uint8> {
  typedef galosh::DeletionOutStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionOutStateTransitionTargets & target, __uint8 c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableByteToDeletionOutStateTransitionTargets_<>::VALUE[ c_source ];
}
//____________________________________________________________________________

template <>
struct CompareType<galosh::DeletionOutStateTransitionTargets, char> {
  typedef galosh::DeletionOutStateTransitionTargets Type;
};
inline void assign ( galosh::DeletionOutStateTransitionTargets & target, char c_source )
{
SEQAN_CHECKPOINT
  target.value =
    galosh::TranslateTableCharToDeletionOutStateTransitionTargets_<>::VALUE[
      ( unsigned char ) c_source
    ];
}
//____________________________________________________________________________

} // End namespace SEQAN_NAMESPACE_MAIN

#endif // __GALOSH_PROFILE_HMM_HPP__
