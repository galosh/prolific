/**
 * \file Parameters.hpp
 * \author  D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * \par Library:
 * galosh::prolific
 * \brief command-line options and parameters for all kinds of galosh related programs.
 * Base class definition.
 * \copyright &copy; 2008, 2011, 2012, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *  All rights reserved.
 *****************************************************************************
##  Library:
##      galosh::prolific
##  File:
##      Parameters.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the Galosh Parameters class.
##  History:
##      Modified, Ted Holzman, 2013
##
#******************************************************************************
#*
#*    This file is part of prolific, a library of useful C++ classes for
#*    working with genomic sequence data and Profile HMMs.  Please see the
#*    document CITING, which should have been included with this file.  You may
#*    use at will, subject to the license (Apache v2.0), but *please cite the
#*    relevant papers* in your documentation and publications associated with
#*    uses of this library.  Thank you!
#*
#*    Copyright (C) 2006, 2008, 2011, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
#*    Research Center.
#*
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
#*****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PARAMETERS_HPP__
#define __GALOSH_PARAMETERS_HPP__

#include "Prolific.hpp"

#include <string>
using std::string;
#include <iostream>
using namespace std;

#include <iterator>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/version.hpp>
#include <boost/program_options.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/lexical_cast.hpp>

using boost::lexical_cast;
using boost::bad_lexical_cast;

#include "CommandlineParameters.hpp"
namespace po = boost::program_options;

namespace std {
  // Print a vector to an ostream.  Put this in namespace std so boost can use it.
  template<typename T>
  ostream& operator<< ( ostream& out, const vector<T>& v ) {
      out << "[";
      size_t last = v.size() - 1;
      for(size_t i = 0; i < v.size(); ++i) {
          out << v[i];
          if (i != last) 
              out << ", ";
      }
      out << "]";
      return out;
  }

} // namespace std


namespace galosh {

  enum DebugLevel {
    DEBUG_None = 0,
    DEBUG_Special = 2,
    DEBUG_Medium = 500,
    DEBUG_All = 1000
  };

  enum VerbosityLevel {
    VERBOSITY_None =     0,
    VERBOSITY_Meta =     5, // lower than low: show only info exterior to training
    VERBOSITY_Low =     10,
    VERBOSITY_Medium =  50,
    VERBOSITY_High =   100,
    VERBOSITY_All =   1000
  };

  class Parameters {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      /**
       * Members to be serialized
       *   TAH 9/13
       */
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
      ar & BOOST_SERIALIZATION_NVP( NAME )
      #include "GaloshOptions.hpp"  /// serialize local parameters
    } // serialize( Archive &, const unsigned int )

  public:
      po::variables_map m_galosh_options_map;
      po::options_description m_galosh_options_description;
      /**
       * Define Parameters "members".  These are tightly tied to the options.
       *   TAH 9/13
       */
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
      TYPE NAME
      #include "GaloshOptions.hpp"  /// declare basic Parameters members

    Parameters ()
    {
      #ifdef DEBUG
      cout << "[debug] Parameters::<init>()" << endl;
      cout << "[debug] using GaloshOptions.hpp" << endl;
      #endif
      /**
       *  Describe all options/parameters to the options_description object.  In the main
       *  routine (whatever it may be) these descriptions will be parsed from the commandline
       *  and possibly other sources.
       *    TAH 9/13
       **/
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
      m_galosh_options_description.add_options()(#NAME,po::value<TYPE>(&NAME)->default_value(DEFAULTVAL) TMP_EXTRA_STUFF,HELP)
      #include "GaloshOptions.hpp"  /// define all the commandline options for this module

      //resetToDefaults();  //This happens when the options map is populated in the main routine
    } // <init>()

    virtual ~Parameters () {};

    // Copy constructor; NOTE THAT THIS DOES NOT COPY THE "DEFAULTS" (THE KEY/VALUE PAIRS IN m_galosh_options_map.)
    Parameters ( const Parameters & copy_from )
    {
      if( copy_from.debug >= DEBUG_All ) {
        cout << "[debug] Parameters::<init>( copy_from )" << endl;
      }
      copyFrom( copy_from );
    } // <init>( Parameters & )

    // Copy constructor/operator
    Parameters & operator= (
      const Parameters & copy_from
    )
    {
      if( copy_from.debug >= DEBUG_All ) {
        cout << "[debug] Parameters::operator=( copy_from )" << endl;
      }
      copyFrom( copy_from );
      return *this;
    } // operator=( const Parameters & )

    void
    copyFromNonVirtual ( const Parameters & copy_from )
    {
      if( copy_from.debug >= DEBUG_All ) {
        cout << "[debug] Parameters::copyFromNonVirtual( copy_from )" << endl;
      }
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
      NAME = copy_from. NAME
      #include "GaloshOptions.hpp"  /// copy all local Parameters members

    } // copyFromNonVirtual( Parameters const & )

    virtual void
    copyFrom ( const Parameters & copy_from )
    {
      copyFromNonVirtual( copy_from );
    } // copyFrom( Parameters & )

    virtual void
    resetToDefaults ()
    {
       #ifdef DEBUG
       cout << "[debug] galosh::Parameters::resetToDefaults()" << endl;
       #endif
       #undef GALOSH_DEF_OPT
       #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) NAME = m_galosh_options_map[#NAME].as<TYPE>()
       #include "GaloshOptions.hpp"  /// copy all local Parameters members
    } // resetToDefaults()

    friend std::ostream &
    operator<< (
      std::ostream& os,
      Parameters const& parameters
    )
    {
      //TODO: REMOVE
      //cout << "in Parameters operator<<" << endl;
      parameters.writeParameters( os );
      return os;
    } // friend operator<< ( ostream, Parameters const& )

    void
    writeParameters (
       std::ostream &os
    ) const
    { // TAH 9/13  [Parameters] is commented out because tokens in square brackets have
      //TODO: REMOVE
      //cout << "in Parameters writeParameters()" << endl;
      // have a special meaning in "configuration files" sensu program_options
      //os << "#[Parameters]" << endl;

      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
      os << #NAME << " = " << lexical_cast<string>(NAME) << endl
      #include "GaloshOptions.hpp"  /// write all Parameters members to os

    } // writeParameters( ostream &) const

  }; // End class Parameters

  template <class ParametersType>
  class ParametersModifierTemplate {
    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      // Serialize the isModified_ stuff  TAH 9/13
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) ar & BOOST_SERIALIZATION_NVP(isModified_##NAME)
      #include "GaloshOptions.hpp"
      // Serialize the parameters
      ar & BOOST_SERIALIZATION_NVP( parameters );
    } // serialize( Archive &, const unsigned int )

  public:
    ParametersType parameters;
    /**
     * Declare isModified_<member> for all members of parameters
     * TAH 9/13
     */
    #undef GALOSH_DEF_OPT
    #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) bool isModified_##NAME
    #include "GaloshOptions.hpp"

    ParametersModifierTemplate () :
      parameters()
    {
      if( parameters.debug) {
        cout << "[debug] GALOSH::ParametersModifierTemplate::<init>()" << endl;
      }
      isModified_reset();
    } // <init>()

    // Copy constructor
    template <class AnyParametersModifierTemplate>
    ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from ) :
      parameters()
    {
      if( copy_from.parameters.debug >= DEBUG_All ) {
        cout << "[debug] GALOSH::ParametersModifierTemplate::<init>( copy_from )" << endl;
      }
      copyFromNonVirtual( copy_from );
    } // <init>( AnyParametersModifierTemplate const & )

    // Copy constructor/operator
    template <class AnyParametersModifierTemplate>
    ParametersModifierTemplate & operator= (
      const AnyParametersModifierTemplate & copy_from
    )
    {
      if( copy_from.parameters.debug >= DEBUG_All ) {
        cout << "[debug] GALOSH::ParametersModifierTemplate::operator=( copy_from )" << endl;
      }
      copyFromNonVirtual( copy_from );
      return *this;
    } // operator=( AnyParametersModifierTemplate const & )

    template <class AnyParametersModifierTemplate>
    void
    copyFromNonVirtual ( const AnyParametersModifierTemplate & copy_from )
    {
      if( copy_from.parameters.debug >= DEBUG_All ) {
        cout << "[debug] GALOSH::ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
      }
      isModified_copyFromNonVirtual( copy_from );

      parameters.copyFromNonVirtual( copy_from.parameters );
    } // copyFromNonVirtual( AnyParametersModifierTemplate const & )

    template <class AnyParametersModifierTemplate>
    void
    isModified_copyFromNonVirtual ( const AnyParametersModifierTemplate & copy_from )
    {
      /// TAH 9/13 copy all isModified members
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) isModified_##NAME = copy_from.isModified_##NAME
      #include "GaloshOptions.hpp"

    } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

    void reset ()
    {
      isModified_reset();
      parameters.resetToDefaults();
    } // reset()

    void
    isModified_reset ()
    {
      if( parameters.debug) {
        cout << "[debug] GALOSH::ParametersModifierTemplate::isModified_reset()" << endl;
      }
      /// TAH 9/13 set all isModified_<member> to false
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) isModified_##NAME = false;
      #include "GaloshOptions.hpp"

    } // isModified_reset()

    friend std::ostream&
    operator<< (
      std::ostream& os,
      ParametersModifierTemplate const& parameters_modifier
    )
    {
      parameters_modifier.writeParametersModifier( os );
      return os;
    } // friend operator<< ( ostream &, ParametersModifierTemplate const& )

    void
    writeParametersModifier (
      std::ostream os
    ) const
    { /// TAH 9/13 Parameter files in boost::program_options have special meanings for bracketed tokens
      os << "#[ParametersModifierTemplate]" << endl;
      /**
       * write out parameters iff they've been modified
       *   TAH 9/13
       */
       #undef GALOSH_DEF_OPT
       #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) if( isModified_##NAME ) os << #NAME << " = " << lexical_cast<string>(parameters. NAME)
       #include "GaloshOptions.hpp" // write out changed parameters

    } // writeParametersModifier ( ostream ) const

    template<class AnyParameters>
    void
    applyModifications ( AnyParameters & target_parameters ) const
    {
      /**
       * Set the parameters of a foreign Parameters object to this Parameter object's values
       * only if they have changed.
       *    TAH 9/13
       */
       #undef GALOSH_DEF_OPT
       #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) if( isModified_##NAME ) target_parameters. NAME = parameters. NAME
       #include "GaloshOptions.hpp" // copy changed parameters

    } // applyModifications( AnyParameters & ) const

  }; // End class ParametersModifierTemplate

  typedef ParametersModifierTemplate<galosh::Parameters> ParametersModifier;

} // End namespace galosh

#endif // __GALOSH_PARAMETERS_HPP__
