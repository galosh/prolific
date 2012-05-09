/**
 * \file Sequence.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      Class definition for the Galosh Sequence class, presently implemented
 *      using seqan::string.
 * \par Overview:
 *    This file is part of prolific, a library of useful C++ classes for
 *    working with genomic sequence data and Profile HMMs.  Please see the
 *    document CITING, which should have been included with this file.  You may
 *    use at will, subject to the license (Apache v2.0), but *please cite the
 *    relevant papers* in your documentation and publications associated with
 *    uses of this library.  Thank you!
 * \copyright
 *    Copyright &copy; 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
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

#ifndef __GALOSH_SEQUENCE_HPP__
#define __GALOSH_SEQUENCE_HPP__

#include "Prolific.hpp"

// For conversion between muscle's Seq class and our Sequence class
#ifdef __HAVE_MUSCLE
#include "muscle/muscle.h"
#include "muscle/seq.h"
#endif //__HAVE_MUSCLE

#include <string>
//using std::string; // TODO: Why is uncommenting this line causing a 'parse
//error' in g++?
#include <iostream>
#include <sstream>
using std::istringstream;
using std::ostringstream;

#include <seqan/sequence.h>
#include <seqan/file.h>

// TODO: REMOVE.  TESTING.
namespace galosh {

  template <typename ResidueType, typename AllocType>
  class Sequence;
}

namespace seqan {
  template <typename ResidueType, typename AllocType >
  inline typename seqan::Size<seqan::String<ResidueType, AllocType> const>::Type
  length ( galosh::Sequence<ResidueType, AllocType> const & seq )
  //length ( seqan::String<ResidueType, AllocType> const & seq )
  {
    //std::cout << "CALLING MY NEW LENGTH FUNCTION" << std::endl;
    return seqan::length( *( dynamic_cast<const seqan::String<ResidueType, AllocType> * const>( &seq ) ) );
  }
}

namespace galosh {

  template <typename ResidueType, typename AllocType=seqan::Alloc<> >
  class Sequence : public seqan::String<ResidueType, AllocType>
  {
    public:
  
    Sequence () :
      seqan::String<ResidueType, AllocType>()
    {
      // Do nothing else.
      assert( this->length() == 0 );
    } // <init>()
  
    Sequence ( uint32_t const length ) :
      seqan::String<ResidueType, AllocType>()
    {
      seqan::resize( *this, length );
      assert( this->length() == length );
    } // <init>( uint32_t )
  
    Sequence ( std::string str ) :
      seqan::String<ResidueType, AllocType>( str )
    {
      // Do nothing else
      // TODO: REMOVE
      //std::cout << "Sequence( " << str << " ): *this is " << *this << std::endl;
      //std::cout << "\tlength( *this ) is " << seqan::length( *this ) << std::endl;
      assert( this->length() == str.length() );
    } // <init>( string )

    Sequence ( char const * const str ) :
      seqan::String<ResidueType, AllocType>( str )
    {
      // Do nothing else
      // TODO: REMOVE
      //std::cout << "Sequence( " << str << " ): *this is " << *this << std::endl;
      //std::cout << "\tlength( *this ) is " << seqan::length( *this ) << std::endl;
      assert( this->length() == std::string( str ).length() );
    } // <init>( string )

    void
    reinitialize ()
    {
      seqan::clear( *this );
    } // reinitialize()

    void
    reinitialize ( uint32_t const length )
    {
      seqan::clear( *this );
      seqan::resize( *this, length );
    } // reinitialize( uint32_t )

    void
    reinitialize ( std::string const & str )
    {
      *this = str;
    } // reinitialize ( string const& )

    uint32_t
    length () const
    {
      return seqan::length( *( dynamic_cast<const seqan::String<ResidueType, AllocType> * const>( this ) ) );
    } // length()

    void
    clear ()
    {
      seqan::clear( *this );
    } // clear()

    /**
     * String reader.  Clobbers existing sequence, replacing it with what it
     * gleans from the given string.
     */
    void
    fromString ( std::string const & str )
    {
      *this = str;
    } // fromString( string const& )

    /**
     * String reader.  Appends to existing sequence what it gleans from the
     * given string.
     */
    void
    appendFromString ( std::string const& str )
    {
      seqan::append( *this, str, ( this->length() + str.length() ), seqan::Exact() );
    } // appendFromString( const& string )

#ifdef __HAVE_MUSCLE
    /**
     * For conversion from the Muscle Seq class
     */
    Sequence const &
    operator= ( Seq const & muscle_seq )
    {
      uint32_t len = muscle_seq.Length();
      reinitialize( len );
      for( uint32_t n = 0; n < len; n++ ) {
        ( *this )[ n ] = muscle_seq.GetChar( n );
      }
      return *this;
    } // operator= ( Seq const & )

    /**
     * For conversion to the Muscle Seq class
     */
    Seq &
    copyIntoSeq ( Seq & muscle_seq, galosh::Sequence<ResidueType> const & seq )
    {
      muscle_seq.FromString( seq.toString(), "" );

      return muscle_seq;
    } // copyIntoSeq ( Seq const &, Sequence<ResidueType> const & )
#endif //__HAVE_MUSCLE

  }; // End class Sequence

} // End namespace galosh

#endif // __GALOSH_SEQUENCE_HPP__
