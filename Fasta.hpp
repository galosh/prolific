/**
 * \file Fasta.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 *      galosh::prolific
 * \brief
 *      Class definition for the Galosh Fasta class, representing collections
 *      of sequences with associated descriptions, as is often represented in
 *      .fasta files.
 * \par Overview:
 *    This file is part of prolific, a library of useful C++ classes for
 *    working with genomic sequence data and Profile HMMs.  Please see the
 *    document CITING, which should have been included with this file.  You may
 *    use at will, subject to the license (Apache v2.0), but *please cite the
 *    relevant papers* in your documentation and publications associated with
 *    uses of this library.  Thank you!
 *
 *  \copyright &copy; 2006, 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson
 *    Cancer Research Center.
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

#ifndef __GALOSH_FASTA_HPP__
#define __GALOSH_FASTA_HPP__

#include "Prolific.hpp"

#include "Sequence.hpp"
using galosh::Sequence;

// For conversion between muscle's SeqVect class and our Fasta class
#ifdef __HAVE_MUSCLE
#include "muscle/muscle.h"
#include "muscle/seq.h"
#include "muscle/seqvect.h"
#endif // __HAVE_MUSCLE

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>

namespace galosh {

  template <typename ResidueType>
  class Fasta : public vector<Sequence<ResidueType> >
  {
  public:
    vector<string> m_descriptions;

    Fasta () :
      vector<Sequence<ResidueType> >(),
      m_descriptions()
    {
      // Do nothing else.
    } // <init>()

    Fasta (
      uint32_t const sequence_count
    ) :
      vector<Sequence<ResidueType> >( sequence_count ),
      m_descriptions( sequence_count )
    {
      // Do nothing else.
    } // <init>( uint32_t const )

    Fasta (
      string const& str
    ) :
      vector<Sequence<ResidueType> >(),
      m_descriptions()
    {
      fromString( str );
    } // <init>( string const& )

    void
    reinitialize ( 
      uint32_t const sequence_count
    )
    {
      if( this->size() != sequence_count ) {
        this->resize( sequence_count );
        m_descriptions.resize( sequence_count );
      }
      // Also initialize the contained sequences.
      for( uint32_t seq_i = 0; seq_i < sequence_count; seq_i++ ) {
        ( *this )[ seq_i ].reinitialize();
      }
    } // reinitialize( uint32_t const )


    uint32_t
    numSequences ()
    {
      return this->size();
    } // numSequences()

    /**
     * Append the sequences and descriptions contained in the other Fasta to
     * the end of this one.
     */
    Fasta &
    operator+= ( Fasta const& other_fasta )
    {
      // TODO: Use a std::vector append operation instead?
      for( uint32_t i = 0; i < other_fasta.size(); i++ ) {
        push_back( other_fasta[ i ] );
        m_descriptions.push_back( other_fasta.m_descriptions[ i ] );
      } // End foreach other fasta sequence

      return *this;
    } // operator+=( Fasta const& )

    /**
     * String reader.  Clobbers existing fasta, replacing it with what it
     * gleans from the given string.
     */
    void
    fromString ( string const& str )
    {
      this->clear();

      istringstream strm(( str ));

      operator>>( strm, *this );
    } // fromString( string const& )

    /**
     * File reader.  Clobbers existing fasta, replacing it with what it
     * gleans from the file with the given filename.
     */
    void
    fromFile ( string const & filename )
    {
      fromFile( filename.c_str() );
    } // fromFile( string )

    /**
     * File reader.  Clobbers existing fasta, replacing it with what it
     * gleans from the file with the given filename.
     */
    bool
    fromFile ( const char * filename )
    {
      std::ifstream fs ( filename );

      // TODO: Why can't I call fromFile( ifstream ) here?
      this->clear();

      if( !fs.is_open() ) {
        // TODO: REMOVE?
        std::cerr << "The fasta file '" << filename << "' could not be opened." << std::endl;
        return false;
      } else {
        operator>>( fs, *this );
        fs.close();
      }
      return true;
    } // fromFile( const char * )

    /**
     * File reader.  Clobbers existing fasta, replacing it with what it
     * gleans from the given file.
     */
    void
    fromFile ( std::ifstream fs )
    {
      this->clear();

      if( !fs.is_open() ) {
        // TODO: Complain.
        std::cerr << "The fasta file is not open." << endl;
      } else {
        operator>>( fs, *this );
      }

    } // fromFile( ifstream )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      //Fasta const& fasta )
      Fasta & fasta
    ) // TODO: Put back const.  Compiler doesn't like it.
    {
      fasta.writeFasta( os );
      return os;
    } // operator<<( basic_ostream&, Fasta const& )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    void
    writeFasta (
      std::basic_ostream<CharT,Traits>& os
    )
    {
      writeFasta( os, this->size() );
    } // writeFasta( basic_ostream& )

    /**
     * Stream writer.
     */
    template<class CharT, class Traits>
    void
    writeFasta (
      std::basic_ostream<CharT,Traits>& os,
      uint32_t num_sequences_to_write
    )
    {
      uint32_t last_seq = num_sequences_to_write - 1;
      uint32_t seq_i;
      for( seq_i = 0; seq_i <= last_seq; seq_i++ ) {
        os << "> " << m_descriptions[ seq_i ] << endl;
        // TODO: char wrap
        os << ( *this )[ seq_i ] << endl;
        os << endl; // Extra line at the end
      }
    } // writeFasta( basic_ostream&, uint32_t )

    /**
     * Stream reader.
     */
    friend std::istream&
    operator>> (
      std::istream& is,
      Fasta & fasta
    )
    {
      static const bool be_extra_verbose = false;//true;

      // TODO: MAKE MORE ROBUST
      string::size_type line_i;
      string line;

      string description;
      Sequence<ResidueType> sequence;

      while( !is.eof() ) {
        getline( is, line );
        if( be_extra_verbose ) {
          cout << "[Fasta]: " << line << endl;
        }
        if( line.size() == 0 ) {
          if( be_extra_verbose ) {
            cout << "[Fasta]: skipping blank line." << endl;
          }
          continue;
        }
        // description
        line_i = line.find( '>' );
        if( be_extra_verbose ) {
          cout << "[Fasta]: found '>' at pos " << line_i << endl;
        }
        if( line_i != string::npos ) {
          if( be_extra_verbose ) {
            cout << "[Fasta]: It's a desc line." << endl;
          }
          if( sequence.length() != 0 ) {
            if( be_extra_verbose ) {
              cout << "[Fasta]: We have a sequence already: " << sequence << endl;
            }
            fasta.m_descriptions.push_back( description );
            fasta.push_back( sequence );
            // Get ready for the next one.
            sequence.clear();
          }
          description =
            line.substr( line_i + 1, ( line.size() - ( line_i + 1 ) ) );
          if( be_extra_verbose ) {
            cout << "[Fasta]: The new desc is: " << description << endl;
          }
          continue;
        }
        sequence.appendFromString( line );
        if( be_extra_verbose ) {
          cout << "[Fasta]: Got some more sequence: " << sequence << endl;
        }
      } // End getting lines 'til eof()
      if( be_extra_verbose ) {
        cout << "[Fasta]: Reached eof." << endl;
      }
      if( sequence.length() != 0 ) {
        if( be_extra_verbose ) {
          cout << "[Fasta]: We have a sequence already: " << sequence << endl;
        }
        fasta.m_descriptions.push_back( description );
        fasta.push_back( sequence );
      }
      
      return is;
    } // operator>>( basic_istream&, Fasta const& )

#ifdef __HAVE_MUSCLE
    /**
     * For conversion from the Muscle SeqVect class.
     */
    Fasta const &
    operator= ( SeqVect const & muscle_seq_vect )
    {
      uint32_t len = muscle_seq_vect.size();

      // TODO: Assert it's the same alphabet type.
      //ALPHA alpha = muscle_seq_vect.GuessAlpha();
      // TODO: REMOVE
      //cout << "alpha is " << alpha << endl;

      reinitialize( len );
      for( uint32_t n = 0; n < len; n++ ) {
        ( *this )[ n ] = muscle_seq_vect.GetSeq( n );
        // TODO: REMOVE
        //cout << "Fasta.operator=(SeqVect): Converted \"" << muscle_seq_vect.GetSeqName( n ) << "\": " << ( *this )[ n ] << endl;
        m_descriptions[ n ] = muscle_seq_vect.GetSeqName( n );
      }
      return *this;
    } // operator= ( SeqVect const & )

    /**
     * For conversion to the Muscle SeqVect class
     */
    SeqVect &
    copyIntoSeqVect (
      SeqVect & muscle_seq_vect,
      Fasta const & fasta
    )
    {
      muscle_seq_vect.Clear();
      Seq muscle_seq;
      uint32_t last_seq = this->size() - 1;
      for( uint32_t seq_i = 0; seq_i <= last_seq; seq_i++ ) {
        muscle_seq.FromString(
          fasta[ seq_i ].toString().c_str(),
          fasta.m_descriptions[ seq_i ].c_str()
        );
        muscle_seq_vect.AppendSeq( muscle_seq );
      } // End foreach seq_i

      return muscle_seq_vect;
    } // copyIntoSeqVect ( SeqVect const &, Fasta const & )

#endif // __HAVE_MUSCLE

  }; // End class Fasta

} // End namespace galosh

#endif // __GALOSH_FASTA_HPP__
