/**
 * \file ProlificParameters.hpp
 * \author D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
 * with some mods by Ted Holzman
 *  \par Library:
 *      galosh::prolific
 *  \brief
 * The galosh::Parameters descendant for programs using the
 * prolific library (inherits from galosh::DynamicProgramming::Parameters).  
 * Adds parameters relating to priors and initial values for Profile HMMs.
 *
 * \par Overview:
 *    This file is part of prolific, a library of useful C++ classes for
 *    working with genomic sequence data and Profile HMMs.  Please see the
 *    document CITING, which should have been included with this file.  You may
 *    use at will, subject to the license (Apache v2.0), but *please cite the
 *    relevant papers* in your documentation and publications associated with
 *    uses of this library.  Thank you!
 *
 *    \copyright &copy; 2008, 2011, 2013 by Paul T. Edlefsen, Fred Hutchinson Cancer
 *    Research Center.
 *
 *    \par License:
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
 *
 *    \par History
 *    9/13 TAH Modifications for boost::program_options to initialize 
 *    parameters and centralize parameter/option manipulations
 *****************************************************************************/

#if     _MSC_VER > 1000
#pragma once
#endif

#ifndef __GALOSH_PROLIFICPARAMETERS_HPP__
#define __GALOSH_PROLIFICPARAMETERS_HPP__

#include "Prolific.hpp"

#include "Parameters.hpp"
using galosh::Parameters;
using galosh::DebugLevel;
using galosh::VerbosityLevel;

#include "DynamicProgramming.hpp"
using galosh::DynamicProgramming;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

namespace galosh {

template <class ResidueType,
          class ProbabilityType,
          class ScoreType,
          class MatrixValueType>
  class ProlificParameters {
    public:

    class Parameters :
       public DynamicProgramming<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters
    {
    private: 
      typedef typename DynamicProgramming<ResidueType, ProbabilityType,ScoreType,MatrixValueType>::Parameters dynamic_programming_parameters_t;
      // Boost serialization
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( dynamic_programming_parameters_t );
        /**
         * ProlificParameters Members to be serialized
         *   TAH 9/13
         **/
         #undef GALOSH_DEF_OPT
         #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) ar & BOOST_SERIALIZATION_NVP( NAME )
         #include "ProlificParametersOptions.hpp"  /// serialize ProlificParameters parameters

      } // serialize( Archive &, const unsigned int )

   public: 
      /// PARAMETERS are now defined in ProlificParametersOptions.hpp
      /**
       * Define ProlificParametersProgramming::Parameters "members".  These are tightly tied to the options.
       *    TAH 9/13
       **/
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) TYPE NAME
      #include "ProlificParametersOptions.hpp"  /// declare Parameters members specific to ProlificParameters
        
      Parameters(); 

      virtual ~Parameters () {};
    
      // Copy constructor
      template <class AnyParameters>
      Parameters ( const AnyParameters & copy_from );
    
      // Copy constructor/operator
      template <class AnyParameters>
      Parameters &
      operator= (
        const AnyParameters & copy_from
      );
    
      template <class AnyParameters>
      void
      copyFromNonVirtual (
        AnyParameters const & copy_from
      );

      template <class AnyParameters>
      void
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      );

      virtual
      void
      copyFrom ( const Parameters & copy_from );
    
      virtual
      void
      resetToDefaults ();

      template<class CharT, class Traits>
      friend std::ostream&
      operator<< (
        std::ostream& os,
        Parameters const& parameters
      )
      {
        parameters.writeParameters( os );
        return os;
      } // friend operator<< ( ostream &, Parameters const& )

      template<class CharT, class Traits>
      void
      writeParameters (
        std::ostream& os
      ) const;

    }; // End inner class Parameters

    template <class ParametersType>
    class ParametersModifierTemplate :
      public DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::template ParametersModifierTemplate<ParametersType>
    {
      typedef typename DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::template ParametersModifierTemplate<ParametersType> base_parameters_modifier_t; 

    private:
      // Boost serialization
      friend class boost::serialization::access;
      template<class Archive>
      void serialize ( Archive & ar, const unsigned int /* file_version */ )
      {
        // save/load base class information.  This will serialize the
        // parameters too.
        ar & BOOST_SERIALIZATION_BASE_OBJECT_NVP( base_parameters_modifier_t );
        // Serialize the ProlificParameters::ParameterModifierTemplate specific isModified_ stuff  TAH 9/13
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) ar & BOOST_SERIALIZATION_NVP(isModified_##NAME)
        #include "ProlificParametersOptions.hpp"  //archive prolificparameters::ParameterModifierTemplate isModified_<member>s
      }

    public:
      /**
       * Declare isModified_<member>s for isModified_<member>s of ProlificParameters::ParameterModifierTemplate
       * TAH 9/13
       **/
      #undef GALOSH_DEF_OPT
      #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) bool isModified_##NAME
      #include "ProlificParametersOptions.hpp"  // Declare isModified<member>s for ProlificParameters::ParametersModifierTemplate

      // Base constructor
      ParametersModifierTemplate();

      // Copy constructor
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from );
    
      // Copy constructor/operator
      template <class AnyParametersModifierTemplate>
      ParametersModifierTemplate & operator= (
        const AnyParametersModifierTemplate & copy_from
      );
    
      template <class AnyParametersModifierTemplate>
      void
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );

      template <class AnyParametersModifierTemplate>
      void
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      );

      void
      reset ();

      void
      isModified_reset ();

      template<class CharT, class Traits>
      friend std::ostream&
      operator<< (
        std::ostream& os,
        ParametersModifierTemplate const& parameters_modifier
      )
      {
        parameters_modifier.writeParametersModifier( os );

        return os;
      } // friend operator<< ( ostream &, ParametersModifierTemplate const& )

      template<class CharT, class Traits>
      void
      writeParametersModifier (
        std::ostream& os
      ) const;

      template <class AnyParameters>
      void
      applyModifications ( AnyParameters & target_parameters );

    }; // End inner class ParametersModifierTemplate

    typedef ParametersModifierTemplate<typename ProlificParameters::Parameters> ParametersModifier;
  }; // End class ProlificParameters

  //======//// potentially non-inline implementations ////========//

  ////// Class galosh::ProlificParameters::Parameters ////
  // Base Constructor
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      Parameters ()
      {
        #ifdef DEBUG
        cout << "[debug] ProlificParameters::Parameters::<init>()" << endl;
        cout << "[debug] using ProlificParametersOptions.hpp" << endl; 
        #endif
        /**
         *  Describe all options/parameters to the options_description object.  In the main
         *  routine (whatever it may be) these descriptions will be parsed from the commandline
         *  and possibly other sources.
         *    TAH 9/13 
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP)          \
        this->galosh::Parameters::m_galosh_options_description.add_options()(#NAME,po::value<TYPE>(&NAME)->default_value(DEFAULTVAL) TMP_EXTRA_STUFF,HELP)
        #include "ProlificParametersOptions.hpp"  /// define all the commandline options for ProlificParameters::Parameters

        // TAH 9/13 No need for resetToDefaults here. 
        // This happens when m_galosh_options_map is populated when the main program starts up
        // resetToDefaults();
      } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      // Copy constructor
      Parameters ( const AnyParameters & copy_from )
      {
        //if( static_cast<galosh::Parameters>( copy_from ).debug >= DEBUG_All ) {
        //  cout << "[debug] ProlificParameters::Parameters::<init>( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_TRIVIAL
  typename ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters &
      // Copy constructor/operator
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      operator= (
        const AnyParameters & copy_from
      )
      {
        if( copy_from.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::Parameters::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_TRIVIAL
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFromNonVirtual (
        AnyParameters const & copy_from
      )
      {
        DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::Parameters::copyFromNonVirtual( copy_from );
        //if( copy_from.debug >= DEBUG_All ) {
        //  cout << "[debug] ProlificParameters::Parameters::copyFromNonVirtual( copy_from )" << endl;
        //} // End if DEBUG_All
        copyFromNonVirtualDontDelegate( copy_from );
      } // copyFromNonVirtual( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class AnyParameters>
  GALOSH_INLINE_COPY
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFromNonVirtualDontDelegate (
        AnyParameters const & copy_from
      )
      {
        /// TAH 9/13
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) NAME = copy_from. NAME
        #include "ProlificParametersOptions.hpp"  /// copy all ProlficParameters::Parameters members

      } // copyFromNonVirtualDontDelegate( AnyParameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  GALOSH_INLINE_TRIVIAL
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      copyFrom ( const Parameters & copy_from )
      {
        copyFromNonVirtual( copy_from );
      } // copyFrom( Parameters const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      resetToDefaults ()
      {
        DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::Parameters::resetToDefaults();
        if( this->galosh::Parameters::debug >= DEBUG_All ) {
           cout << "[debug] ProlificParameters::Parameters::resetToDefaults()" << endl;
        } // End if DEBUG_All
        /// TAH 9/13
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) NAME = this->galosh::Parameters::m_galosh_options_map[#NAME].template as<TYPE>()
        #include "ProlificParametersOptions.hpp"  /// reset all Parameters members 
                                                  /// (ProlificParameters and through inheritance tree)

      } // resetToDefaults()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::Parameters::
      writeParameters (
        std::ostream& os
      ) const
      {
        DynamicProgramming<ResidueType,ProbabilityType,ScoreType,MatrixValueType>::Parameters::writeParameters( os );
        os << endl;
        //Note: we must comment out [ProlificParameters] because it means something special to 
        //in configuration files sensu program_options
        os << "#[ProlificParameters]" << endl;

        /**
         * write out all ProlificParameters specific parameters in the style of a configuration
         * file, so that program_options parsers can read it back in
         *   TAH 9/13
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) os << #NAME << " = " << lexical_cast<string>(NAME) << endl 
        #include "ProlificParametersOptions.hpp"  /// write all ProlificParameters::Parameters members to os

      } // writeParameters ( ostream & ) const

  ////// Class galosh::ProlificParameters::ParametersModifierTemplate ////

  /// TAH 9/13 \todo figure out if we have to initialize the parameters member to ProlificParameters::Parameters  
  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      ParametersModifierTemplate ()
      {
        #ifdef DEBUG
        cout << "[debug] ProlificParameters::ParametersModifierTemplate::<init>()" << endl;
        #endif   
        isModified_reset();
      } // <init>()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_INIT
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor
      ParametersModifierTemplate ( const AnyParametersModifierTemplate & copy_from )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::ParametersModifierTemplate::<init>( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
      } // <init>( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_TRIVIAL
  typename ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::template ParametersModifierTemplate<ParametersType> &
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      // Copy constructor/operator
  operator= (
        const AnyParametersModifierTemplate & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::ParametersModifierTemplate::operator=( copy_from )" << endl;
        } // End if DEBUG_All
        copyFromNonVirtual( copy_from );
        return *this;
      } // operator=( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        if( copy_from.parameters.debug >= DEBUG_All ) {
          cout << "[debug] ProlificParameters::ParametersModifierTemplate::copyFromNonVirtual( copy_from )" << endl;
        } // End if DEBUG_All

        isModified_copyFromNonVirtual( copy_from );

        base_parameters_modifier_t::parameters.copyFromNonVirtual( copy_from.parameters );
      } // copyFromNonVirtual( AnyParametersModifierTemplate const & )


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParametersModifierTemplate>
  GALOSH_INLINE_COPY
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      isModified_copyFromNonVirtual (
        AnyParametersModifierTemplate const & copy_from
      )
      {
        /// Copy all parent isModified_<member>s
        base_parameters_modifier_t::isModified_copyFromNonVirtual( copy_from );
        /// TAH 9/13 copy all ProlficParameters::isModified_<member>s
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) isModified_##NAME = copy_from.isModified_##NAME
        #include "ProlificParametersOptions.hpp"  // Copy ProlificParameters::ParametersModifierTemplate::isModified_<member> members

      } // isModified_copyFromNonVirtual( AnyParametersModifierTemplate const & )

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_TRIVIAL
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      reset ()
      {
        isModified_reset();
        base_parameters_modifier_t::parameters.resetToDefaults();
      } // reset()

  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  GALOSH_INLINE_REINITIALIZE
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      isModified_reset ()
      {
        /// Set all parent ParametersModifierTemplate::isModified_<member>s to false
        base_parameters_modifier_t::isModified_reset();
        /// TAH 9/13 set all ProlificParameters::ParametersModifierTemplate::isModified_<member>s to false
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) isModified_##NAME = false;
        #include "ProlificParametersOptions.hpp" // Reset ProlificParameters::ParametersModifierTemplate::isModified_<member>s to false
  
      } // isModified_reset()

  template <typename ResidueType,
            typename ProbabilityType,
            typename ScoreType,
            typename MatrixValueType>
  template <class ParametersType>
  template<class CharT, class Traits>
  GALOSH_INLINE_OSTREAM
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      writeParametersModifier (
        std::ostream& os
      ) const
      {
        //base_parameters_modifier_t::operator<<( os, parameters_modifier );
        // Write out parent parameters if they've been changed
        base_parameters_modifier_t::writeParametersModifier( os );
        os << endl;
        /// TAH 9/13 must comment out tags in square braces for program_options config file parser 
        os << "#[ProlificParameters]" << endl;
        /**
         * write out ProlificParameters::parameters iff they've been modified
         *   TAH 9/13
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) if( isModified_##NAME ) os << #NAME << " = " << lexical_cast<string>(galosh::ParametersModifierTemplate<ParametersType>::parameters. NAME)
        #include "ProlificParametersOptions.hpp" // write out changed ProlificParameters::ParametersModifierTemplate parameters

      } // writeParametersModifier ( ostream & ) const


  template <class ResidueType,
            class ProbabilityType,
            class ScoreType,
            class MatrixValueType>
  template <class ParametersType>
  template <class AnyParameters>
  GALOSH_INLINE_PARAMETERSMODIFIER_APPLY_MODIFICATIONS
  void
  ProlificParameters<ResidueType, ProbabilityType, ScoreType, MatrixValueType>::ParametersModifierTemplate<ParametersType>::
      applyModifications ( AnyParameters & target_parameters )
      {
        /// Set the parameters of another object iff they've been changed in this one
        base_parameters_modifier_t::applyModifications( target_parameters );
        /**
         * Set the parameters of a foreign Parameters object to this Parameter object's values
         * iff they have changed.
         *    TAH 9/13
         **/
        #undef GALOSH_DEF_OPT
        #define GALOSH_DEF_OPT(NAME,TYPE,DEFAULTVAL,HELP) if( isModified_##NAME ) target_parameters. NAME = this->parameters. NAME
        #include "ProlificParametersOptions.hpp" // copy changed parameters

      } // applyModifications( Parameters & )
   
} // End namespace galosh

#endif // __GALOSH_PROLIFICPARAMETERS_HPP__
