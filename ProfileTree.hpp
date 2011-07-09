/*---------------------------------------------------------------------------##
##  Library:
##      galosh::prolific
##  File:
##      ProfileTree.hpp
##  Author:
##      D'Oleris Paul Thatcher Edlefsen   paul@galosh.org
##  Description:
##      Class definition for the Galosh ProfileTree class.  A ProfileTree is a
##      graph (a tree, actually) view of a collection of
##      ProfileTreeInternalNodes and a ProfileTreeRoot.
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
#*    Copyright (C) 2007, 2008, 2011 by Paul T. Edlefsen, Fred Hutchinson Cancer
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

#ifndef __GALOSH_PROFILETREE_HPP__
#define __GALOSH_PROFILETREE_HPP__

#include "Prolific.hpp"

#include "Profile.hpp"
using galosh::ProfileTreeRoot;
using galosh::ProfileTreeInternalNode;

#include "Random.hpp"
using galosh::Random;

#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/version.hpp>

#include <boost/tuple/tuple.hpp> // Boost version 1.44.0: This is supplies boost::tie, required by graph_as_tree, but for some reason not imported thereby.
#include <boost/graph/graph_as_tree.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/adj_list_serialize.hpp>
#include <boost/cstdlib.hpp>

//#include <boost/config.hpp>
//#include <iostream>                      // for std::cout
//#include <utility>                       // for std::pair
//#include <algorithm>                     // for std::for_each
//#include <boost/utility.hpp>             // for boost::tie
//#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
//#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/graphviz.hpp>

//#include <iostream>
using std::cout;
using std::endl;

#include <vector>
using std::vector;

using namespace boost;

namespace galosh {

template <typename ResidueType,
          typename ProbabilityType,
          typename InternalNodeType>// = ProfileTreeRoot<ResidueType, ProbabilityType> >
  class ProfileTree
  {
  protected:
    typedef adjacency_list<vecS, vecS, bidirectionalS, property<vertex_name_t, std::string> > graph_t;
    typedef graph_traits<graph_t>::vertex_descriptor vertex_t;
    typedef graph_traits<graph_t>::edge_descriptor edge_t;
    typedef property_map<graph_t, vertex_name_t>::type vertex_name_map_t;

    graph_t m_graph;
    ProfileTreeRoot<ResidueType, ProbabilityType> m_root;
    vector<InternalNodeType > m_internalNodes;

  public:
    struct ProfileTreeInternalNodePromise
    {
      ProfileTree * m_profileTree;
      vertex_t m_internalNodeVertex;

      ProfileTreeInternalNodePromise () :
        m_profileTree( NULL ),
        m_internalNodeVertex( 0 )
      {
        // Do nothing else.
      } // ProfileTreeInternalNodePromise()

      /**
       * This class works around the problem that when you add internal nodes,
       * the vector resizes and thus pointers to ProfileTreeInternalNodes
       * become invalid.  Instead we return an object of type
       * ProfileTreeInternalNodePromise, which will return the reference to the
       * internal node on operator().
       */
      ProfileTreeInternalNodePromise (
        ProfileTree * profile_tree,
        vertex_t internal_node_vertex
      ) :
        m_profileTree( profile_tree ),
        m_internalNodeVertex( internal_node_vertex )
      {
        // Do nothing else.
      } // ProfileTreeInternalNodePromise( ProfileTree *, vertex_t )

      InternalNodeType &
      operator() ()
      {
        return m_profileTree->m_internalNodes[ m_internalNodeVertex - 1 ];
      } // operator()()

      InternalNodeType const &
      operator() () const
      {
        return m_profileTree->m_internalNodes[ m_internalNodeVertex - 1 ];
      } // operator()() const

    }; // End inner class ProfileTreeInternalNodePromise

    // Boost serialization
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize ( Archive & ar, const unsigned int /* file_version */ )
    {
      ar & BOOST_SERIALIZATION_NVP( m_graph );
      ar & BOOST_SERIALIZATION_NVP( m_root );
      ar & BOOST_SERIALIZATION_NVP( m_internalNodes );
    } // serialize( Archive &, const unsigned int )

  public:
    ProfileTree () :
      m_graph(),
      m_root(),
      m_internalNodes()
    {
      addRootToGraph();
    } // <init>()

    ProfileTree (
      uint32_t length
    ) :
      m_graph(),
      m_root( length ),
      m_internalNodes()
    {
      addRootToGraph();
    } // <init>( uint32_t )

    // Copy constructor
    ProfileTree (
      ProfileTree<ResidueType, ProbabilityType,InternalNodeType> const& other_profile_tree
    ) :
      m_graph( other_profile_tree.m_graph ),
      m_root( other_profile_tree.m_root ),
      m_internalNodes( other_profile_tree.m_internalNodes )
    {
      // Do nothing else
      // TODO: REMOVE
      //m_root.copyFrom( other_profile_tree.m_root );
    } // <init>( ProfileTree const& )

    ProfileTree &
    operator= ( ProfileTree const& other_tree )
    {
      reinitialize( other_tree );

      return *this;
    } // operator=( AnyInternalNodeOrRoot const& )

    void
    reinitialize (
      uint32_t length
    )
    {
      m_graph.clear();
      m_root.reinitialize( length );
      m_internalNodes.resize( 0 );

      addRootToGraph();
    } // reinitialize( uint32_t )

    void
    reinitialize (
      ProfileTree<ResidueType, ProbabilityType, InternalNodeType> const& other_profile_tree
    )
    {
      m_graph = other_profile_tree.m_graph;
      //m_root = other_profile_tree.m_root;
      // TODO: REMOVE.  Workaround (see above).
      m_root.copyFrom( other_profile_tree.m_root );
      if( m_internalNodes.size() != other_profile_tree.m_internalNodes.size() ) {
        m_internalNodes.resize( other_profile_tree.m_internalNodes.size() );
      }

      //m_internalNodes = other_profile_tree.m_internalNodes;
      // TODO: REMOVE.  Workaround..  Not sure if this one is necessary.
      for( uint32_t internal_node_i = 0; internal_node_i < m_internalNodes.size(); internal_node_i++ ) {
        m_internalNodes[ internal_node_i ].reinitialize( other_profile_tree.m_internalNodes[ internal_node_i ] );
        m_internalNodes[ internal_node_i ].copyPositions( other_profile_tree.m_internalNodes[ internal_node_i ] );
      }

      // TODO: REMOVE
      //m_root.copyFrom( other_profile_tree.m_root );
    } // reinitialize( ProfileTree const& )

    ProfileTreeRoot<ResidueType, ProbabilityType> *
    getProfileTreeRoot ()
    {
      return &m_root;
    } // getProfileTreeRoot()

    ProfileTreeRoot<ResidueType, ProbabilityType> const *
    getProfileTreeRoot () const
    {
      return &m_root;
    } // getProfileTreeRoot() const

    /**
     * Return the size of this ProfileTree.
     */
    uint32_t
    nodeCount () const
    {
      return ( m_internalNodes.size() + 1 );
    } // nodeCount() const

    /**
     * Return a reference to the internal node with the given vertex.
     *
     * WARNING: holding onto a reference to the child is dangerous if the tree
     * is to be modified.  Use getProfileTreeInternalNodePromise instead if you
     * plan on changing the topology of the tree while holding a reference to a
     * child.
     */
    InternalNodeType &
    getProfileTreeInternalNode ( vertex_t const & vertex )
    {
      // TODO: This won't work if InternalNodeType isn't the same as the Root
      // type.
      if( vertex == 0 ) {
        return *getProfileTreeRoot();
      }
      if( ( vertex == 0 ) || ( vertex > m_internalNodes.size() ) ) {
        // TODO: REMOVE?
        cout << "WARNING: vertex out of range: " << vertex << "." << endl;
        assert( false );
        // TODO: Throw an exception!
      }
      // TODO: REMOVE
      //cout << "getProfileTreeInternalNode( " << vertex << " )" << endl;
      // vertex indices are offset by 1 from the internalNodes vector indices
      return m_internalNodes[ vertex - 1 ];
    } // getProfileTreeInternalNode( vertex_t const & )

    /**
     * Return a reference to the internal node with the given vertex.
     *
     * WARNING: holding onto a reference to the child is dangerous if the tree
     * is to be modified.  Use getProfileTreeInternalNodePromise instead if you
     * plan on changing the topology of the tree while holding a reference to a
     * child.
     */
    InternalNodeType const&
    getProfileTreeInternalNode ( vertex_t const & vertex ) const
    {
      // TODO: This won't work if InternalNodeType isn't the same as the Root
      // type.
      if( vertex == 0 ) {
        return *getProfileTreeRoot();
      }
      if( ( vertex == 0 ) || ( vertex > m_internalNodes.size() ) ) {
        // TODO: REMOVE?
        cout << "WARNING: vertex out of range: " << vertex << "." << endl;
        assert( false );
        // TODO: Throw an exception!
      }
      // TODO: REMOVE
      //cout << "getProfileTreeInternalNode( " << vertex << " )" << endl;
      // vertex indices are offset by 1 from the internalNodes vector indices
      return m_internalNodes[ vertex - 1 ];
    } // getProfileTreeInternalNode( vertex_t const & ) const

    ProfileTreeInternalNodePromise
    getProfileTreeInternalNodePromise ( vertex_t const & vertex )
    {
      // TODO: This won't work if InternalNodeType isn't the same as the Root
      // type.
      if( vertex == 0 ) {
        return ProfileTreeInternalNodePromise( this, 0 );
      }
      if( ( vertex == 0 ) || ( vertex > m_internalNodes.size() ) ) {
        // TODO: REMOVE?
        cout << "WARNING: vertex out of range: " << vertex << "." << endl;
        assert( false );
        // TODO: Throw an exception!
      }
      // TODO: REMOVE
      //cout << "getProfileTreeInternalNodePromise( " << vertex << " )" << endl;
      // vertex indices are offset by 1 from the internalNodes vector indices
      return ProfileTreeInternalNodePromise( this, vertex );
    } // getProfileTreeInternalNodePromise( vertex_t const & )

    vertex_t
    addChild (
      vertex_t const& parent_vertex
    )
    {
      if( parent_vertex == 0 ) {
        return addChildToRoot();
      } else {
        return addChildToInternalNode( parent_vertex );
      }
    } // addChild( vertex_t const& )

    vertex_t
    addChildToRoot ()
    {
      vector<int> variations =
        createDefaultParentPositionVariations( m_root.length() );
      // TODO: REMOVE
      //cout << "from createDefault..: variations are: [ ";
      //for ( int i = 0; i < variations.size(); i++ ) {
      //  if( i > 0 ) {
      //    cout << ", ";
      //  }
      //  cout << variations[ i ];
      //}
      //cout << " ]" << endl;

      return addChildToRoot( variations );
    } // addChildToRoot()

    vertex_t
    addChildToRoot (
      vector<int> const& parent_position_variations
    )
    {
      // This is the index into our vector, not into the graph.
      uint32_t new_child_index =
        m_internalNodes.size();
      // TODO: REMOVE
      //cout << "addChildToRoot: Parent Position variations are: [ ";
      //for ( int i = 0; i < parent_position_variations.size(); i++ ) {
      //  if( i > 0 ) {
      //    cout << ", ";
      //  }
      //  cout << parent_position_variations[ i ];
      //}
      //cout << " ]" << endl;

      InternalNodeType tmp( &m_root, parent_position_variations );
      m_internalNodes.push_back( tmp );
      //m_internalNodes.resize( new_child_index + 1 );
      //m_internalNodes[ new_child_index ].reinitialize(
      //   &m_root,
      //   parent_position_variations
      //);

      vertex_t new_child_vertex = add_vertex( m_graph );

      // assert( new_child_vertex == ( new_child_index + 1 ) )
      if( new_child_vertex != ( new_child_index + 1 ) ) {
        cout << "WARNING: new child vertex is " << new_child_vertex << ", NOT " << ( new_child_index + 1 ) << "!" << endl;
        assert( false );
      }

      // "ProfileTreeVertex" is the graph vertex (one greater than the
      // internalNodes vector index).
      m_internalNodes[ new_child_index ].setProfileTreeVertex(
        new_child_vertex // new_child_index + 1
      );
      // TODO: REMOVE.  TESTING.
      for( uint32_t i = 0; i <= new_child_index; i++ ) {
        assert( m_internalNodes[ i ].getProfileTreeVertex() == ( i + 1 ) );
      }

      // Get the root vertex.  How?
      //add_edge( root_vertex, new_child_vertex, m_graph );
      add_edge( 0, new_child_vertex, m_graph );

      // Name it.
      vertex_name_map_t name = get( vertex_name, m_graph );
      name[ new_child_vertex ] =
 	boost::lexical_cast<std::string>( new_child_vertex );

      return new_child_vertex;
    } // addChildToRoot( vector<int> const& )

    vertex_t
    addChildToInternalNode (
      vertex_t const& parent_vertex
    )
    {
      assert( parent_vertex > 0 );
      assert( parent_vertex <= m_internalNodes.size() );
      // TODO: REMOVE.  TESTING.
      for( uint32_t i = 0; i < m_internalNodes.size(); i++ ) {
        assert( m_internalNodes[ i ].getProfileTreeVertex() == ( i + 1 ) );
      }
      // TODO: REMOVE
      assert(
        parent_vertex ==
        m_internalNodes[ parent_vertex - 1 ].getProfileTreeVertex()
      );
      return
        addChildToInternalNode(
          m_internalNodes[ parent_vertex - 1 ]
        );
    } // addChildToInternalNode( InternalNodeType const& )

    vertex_t
    addChildToInternalNode (
      InternalNodeType const& parent
    )
    {
      vector<int> variations =
        createDefaultParentPositionVariations( parent.length() );
      return addChildToInternalNode( parent.getProfileTreeVertex(), variations );
    } // addChildToInternalNode( InternalNodeType const& )

    vertex_t
    addChildToInternalNode (
      vertex_t const& parent_vertex,
      vector<int> const& parent_position_variations
    )
    {
      assert( parent_vertex > 0 );
      assert( parent_vertex <= m_internalNodes.size() );
      // This is the index into our vector, not into the graph.
      uint32_t new_child_index =
        m_internalNodes.size();
      // TODO: REMOVE
      //cout << "Parent index is " << parent_vertex << endl;
      //cout << "Parent Position variations are: { " << endl;
      //for ( int i = 0; i < parent_position_variations.size(); i++ ) {
      //  if( i > 0 ) {
      //    cout << ", ";
      //  }
      //  cout << parent_position_variations[ i ];
      //}
      //cout << " }" << endl;

      InternalNodeType tmp(
        // graph/ProfileTree indices of internalNodes are offset by 1 (since 0
        // is root).  In the internalNodes array the indices are decremented
        // by 1.
        &m_internalNodes[ parent_vertex - 1 ],
        parent_position_variations
      );
      m_internalNodes.push_back( tmp );

      //m_internalNodes.resize( new_child_index + 1 );
      //m_internalNodes[ new_child_index ].reinitialize(
      //   // graph/ProfileTree indices of internalNodes are offset by 1 (since 0
      //   // is root).  In the internalNodes array the indices are decremented
      //   // by 1.
      //   &m_internalNodes[ parent_vertex - 1 ],
      //   parent_position_variations
      //);
      vertex_t new_child_vertex = add_vertex( m_graph );

      // assert( new_child_vertex == ( new_child_index + 1 ) )
      if( new_child_vertex != ( new_child_index + 1 ) ) {
        cout << "WARNING: new child vertex is " << new_child_vertex << ", NOT " << ( new_child_index + 1 ) << "!" << endl;
        assert( false );
      }

      // "ProfileTreeVertex" is the graph vertex (one greater than the
      // internalNodes vector index).
      m_internalNodes[ new_child_index ].setProfileTreeVertex(
        new_child_vertex // new_child_index + 1
      );
      // TODO: REMOVE.  TESTING.
      for( uint32_t i = 0; i <= new_child_index; i++ ) {
        assert( m_internalNodes[ i ].getProfileTreeVertex() == ( i + 1 ) );
      }

      add_edge( parent_vertex, new_child_vertex, m_graph );

      // Name it.
      vertex_name_map_t name = get( vertex_name, m_graph );
      name[ new_child_vertex ] =
 	boost::lexical_cast<std::string>( new_child_vertex );

      return new_child_vertex;
    } // addChildToInternalNode( vertex_t const &, vector<int> const& )

    /**
     * Return a Promise for one of the children of the Node with the given
     * vertex.  The child returned is the Nth, where N is the which_child
     * argument, and the order is the order in which the children were added.
     * Note that this means that children are indexed from 1, not 0.
     */
    ProfileTreeInternalNodePromise
    getChildPromise (
      vertex_t const & parent_vertex,
      uint32_t const which_child
    )
    {
      // TODO: Make sure that there are enough children?
      return getProfileTreeInternalNodePromise( getChildVertex( parent_vertex, which_child ) );
    } // getChildPromise( vertex_t const &, uint32_t const )

    /**
     * Return a reference to one of the children of the Node with the given
     * vertex.  The child returned is the Nth, where N is the which_child
     * argument, and the order is the order in which the children were added.
     * Note that this means that children are indexed from 1, not 0.
     *
     * WARNING: holding onto a reference to the child is dangerous if the tree
     * is to be modified.  Use getChildPromise instead if you plan on changing
     * the topology of the tree while holding a reference to a child.
     */
    InternalNodeType &
    getChild (
      vertex_t const & parent_vertex,
      uint32_t const which_child
    )
    {
      // TODO: Make sure that the vertices are valid.
      return getProfileTreeInternalNode( getChildVertex( parent_vertex, which_child ) );
    } // getChild( vertex_t const &, uint32_t const )

    /**
     * Return the vertex of one of the children of the Node with the given
     * vertex.  The child returned is the Nth, where N is the which_child
     * argument, and the order is the order in which the children were added.
     * Note that this means that children are indexed from 1, not 0.  If the
     * node with the given vertex does not have enough children to satisfy the
     * request, the return value will be 0.
     */
    vertex_t
    getChildVertex (
      vertex_t const & parent_vertex,
      uint32_t const which_child
    ) const
    {
      // TODO: REMOVE
      //cout << "getChildVertex( " << parent_vertex << ", " << which_child << " )" << endl;

      //typename property_map<Graph, vertex_index_t>::type 
      //  vertex_id = get( vertex_index, m_graph );
      typename graph_traits<graph_t>::adjacency_iterator ai, ai_end;
      uint32_t child_i = 1;
      // TODO: Is there a more efficient way?  The underlying type is vecS (a
      // vector), so shouldn't I be able to get constant-time access?
      for( tie( ai, ai_end ) =
             adjacent_vertices( parent_vertex, m_graph );
           ai != ai_end;
           ++ai, ++child_i
      ) {
        // TODO: REMOVE
        //cout << "child_i is " << child_i << endl;
        if( child_i == which_child ) {
          return *ai;
        }
      } // End foreach child of vertex v
      return 0;
    } // getChildVertex( vertex_t const &, uint32_t const ) const

    /**
     * This is the opposite of getChildVertex: given a vertex, return which
     * child it is.
     */
    uint32_t
    getChildIndexInParent (
      vertex_t const & parent_vertex,
      vertex_t const & child_vertex
    ) const
    {
      // TODO: REMOVE
      //cout << "getChildVertex( " << parent_vertex << ", " << which_child << " )" << endl;

      //typename property_map<Graph, vertex_index_t>::type 
      //  vertex_id = get( vertex_index, m_graph );
      typename graph_traits<graph_t>::adjacency_iterator ai, ai_end;
      uint32_t child_i = 1;
      // TODO: Is there a more efficient way?  The underlying type is vecS (a
      // vector), so shouldn't I be able to get constant-time access?
      for( tie( ai, ai_end ) =
             adjacent_vertices( parent_vertex, m_graph );
           ai != ai_end;
           ++ai, ++child_i
      ) {
        // TODO: REMOVE
        //cout << "child_i is " << child_i << endl;
        if( *ai == child_vertex ) {
          return child_i;
        }
      } // End foreach child of vertex v
      return 0;
    } // getChildIndexInParent( vertex_t const &, vertex_t const & ) const

    /**
     * Return a reference to the given node's parent.
     */
    // TODO: Make this work when InternalNodeType is not ProfileTreeRoot.
    InternalNodeType const &
    getParent (
      vertex_t const & child_vertex
    ) const
    {
      // TODO: Make sure that the vertices are valid.
      return getProfileTreeInternalNode( getParentVertex( child_vertex ) );
    } // getParent( vertex_t const & ) const

    /**
     * Return a reference to the given node's parent.
     */
    // TODO: Make this work when InternalNodeType is not ProfileTreeRoot.
    InternalNodeType &
    getParent (
      vertex_t const & child_vertex
    )
    {
      // TODO: Make sure that the vertices are valid.
      return getProfileTreeInternalNode( getParentVertex( child_vertex ) );
    } // getParent( vertex_t const & )

    /**
     * Return the vertex of the node's parent.  If the node is the root, it has
     * no parent, but 0 will be returned (which is the index of the root, so
     * it's like the root is it's own parent.  Be warned.
     */
    vertex_t
    getParentVertex (
      vertex_t const & node_vertex
    ) const
    {
      if( node_vertex == 0 ) {
        return node_vertex;
      }
      typename graph_traits<graph_t>::in_edge_iterator in_i, in_end;
      tie(in_i, in_end) = in_edges( node_vertex, m_graph );
      return source( *in_i, m_graph );
    } // getParentVertex( vertex_t const & ) const

    /**
     * Return the number of children of the Node with the given vertex.
     */
    uint32_t
    childCount ( vertex_t const & parent_vertex ) const
    {
      return out_degree( parent_vertex, m_graph );
    } // childCount( vertex_t const & ) const

    /**
     * Return a randomly-drawn child of the given vertex.  If the given parent
     * vertex has no children, the result will be 0.  This method uses a
     * uniform distribution over the children of the given parent_vertex.
     */
    vertex_t
    drawChild ( vertex_t const & parent_vertex, Random & random ) const
    {
      return
        drawChild( parent_vertex, NULL, random );
    } // drawChild( vertex_t const &, Random & ) const

    void
    setName (
      vertex_t const & node_vertex,
      string const & name
    ) {
      vertex_name_map_t name_map = get( vertex_name, m_graph );
      name_map[ node_vertex ] = name;
    } // setName( vertex_t const &, string const & )

    /**
     * Return a randomly-drawn child of the given vertex.  If the given parent
     * vertex has no children, the result will be 0.  This method uses a
     * uniform distribution over the children of the given parent_vertex when
     * its vector<ProbabilityType> * argument is NULL.  If non-null, that
     * vector must be of length childCount( parent_vertex ), and it is assumed
     * normalized (that is, all values are probabilities and they sum to 1).
     */
    vertex_t
    drawChild (
      vertex_t const & parent_vertex,
      vector<ProbabilityType> const * const probabilities,
      Random & random
    ) const
    {
      uint32_t child_count = childCount( parent_vertex );
      if( child_count == 0 ) {
        return 0;
      }
      
      // TODO: Check the assertion that probabilities (if non-NULL) is of
      // length child_count.
      ProbabilityType u = static_cast<ProbabilityType>( random.nextUniform() );
      ProbabilityType next_bin_boundary = ProbabilityType( 0 );
      for( uint32_t i = 0; i < child_count; i++ ) {
        if( probabilities != NULL ) {
          next_bin_boundary += ( *probabilities )[ i ];
        } else {
          next_bin_boundary += ( 1.0 / child_count );
        }
        if( u < next_bin_boundary ) {
          // We add 1 because which_child is indexed from 1 in getChildVertex().
          return getChildVertex( parent_vertex, ( i + 1 ) );
        }
      } // End foreach potential which_child index i

    } // drawChild( vertex_t const &, vector<ProbabilityType> const * const, Random & ) const

    /**
     * Return a randomly-drawn child of the given vertex, or the given vertex.
     * If the given parent vertex has no children, the result will always be
     * the parent.  This method uses a uniform distribution over the parent and
     * its children.
     */
    vertex_t
    drawChildOrParent ( vertex_t const & parent_vertex, Random & random ) const
    {
      return
        drawChildOrParent( parent_vertex, NULL, random );
    } // drawChildOrParent( vertex_t const &, Random & ) const

    /**
     * Return a randomly-drawn child of the given vertex, or the given vertex.
     * If the given parent vertex has no children, the result will always be
     * the parent.  This method uses a uniform distribution over the parent and
     * children when its vector<ProbabilityType> argument is NULL.  If
     * non-null, that vector must be of length ( 1 + childCount( parent_vertex
     * ) ), and it is assumed normalized (that is, all values are probabilities
     * and they sum to 1).  The 0th index of that vector will be used for the
     * probability of the parent.
     */
    vertex_t
    drawChildOrParent (
      vertex_t const & parent_vertex,
      vector<ProbabilityType> const * const probabilities,
      Random & random
    ) const
    {
      uint32_t child_count = childCount( parent_vertex );
      if( child_count == 0 ) {
        return parent_vertex;
      }

      // TODO: Check the assertion that probabilities (if non-NULL) is of
      // length (child_count + 1 ).
      ProbabilityType u = static_cast<ProbabilityType>( random.nextUniform() );
      ProbabilityType next_bin_boundary = ProbabilityType( 0 );
      for( uint32_t i = 0; i <= child_count; i++ ) {
        if( probabilities != NULL ) {
          next_bin_boundary += ( *probabilities )[ i ];
        } else {
          next_bin_boundary += ProbabilityType( 1.0 / ( child_count + 1 ) );
        }
        if( u < next_bin_boundary ) {
          if( i == 0 ) {
            return parent_vertex;
          } else {
            return getChildVertex( parent_vertex, i );
          }
        }
      } // End foreach potential which_child index i

    } // drawChildOrParent( vertex_t const &, vector<ProbabilityType> const * const, Random & ) const

    template<class CharT, class Traits>
    friend std::basic_ostream<CharT,Traits>&
    operator<< (
      std::basic_ostream<CharT,Traits>& os,
      ProfileTree & profile_tree
    )
    {
      profile_tree.showGraph( os );
      os << "Root:" << endl;
      os << ( *profile_tree.getProfileTreeRoot() ) << endl;
      uint32_t tree_size = profile_tree.nodeCount();
      for( uint32_t internal_node_i = 1; internal_node_i < tree_size; internal_node_i++ ) {
        os << "Node " << boost::lexical_cast<string>( internal_node_i ) << ":" << endl;
        os << profile_tree.getProfileTreeInternalNode( internal_node_i ) << endl;
      }

      return os;
    } // friend operator<< ( basic_ostream, ProfileTreeInternalNode const&)

    template <typename StreamType>
    void
    showGraph ( StreamType & stream )
    {
      boost::property_map<graph_t, vertex_index_t>::type 
        vertex_id = get(vertex_index, m_graph);
    
      vertex_name_map_t
        vertex_name_map = get(vertex_name, m_graph);
      typedef graph_traits<graph_t>::vertex_iterator vertex_iter;
      std::pair<vertex_iter, vertex_iter> vp;

      vertex_t root_vertex = 0;
      //stream << "vertices(m_graph) = ";
      for (vp = vertices(m_graph); vp.first != vp.second; ++vp.first) {
        //stream << vertex_name_map[get(vertex_id, *vp.first)] <<  " ";
        if( get(vertex_id, *vp.first) == 0 ) {
          root_vertex = *vp.first;
        }
      }
      //stream << std::endl;
      
      //stream << "edges(m_graph) = ";
      //graph_traits<graph_t>::edge_iterator ei, ei_end;
      //for (tie(ei,ei_end) = edges(m_graph); ei != ei_end; ++ei) {
      //  stream << "(" << vertex_name_map[get(vertex_id, source(*ei, m_graph))]
      //            << "," << vertex_name_map[get(vertex_id, target(*ei, m_graph))]
      //            << ") ";
      //}
      //stream << std::endl;
      
      //std::for_each(vertices(m_graph).first, vertices(m_graph).second,
      //              exercise_vertex<graph_t, StreamType>(m_graph, stream));
    
    
    
      typedef iterator_property_map<std::vector<vertex_t>::iterator,
        property_map<graph_t, vertex_index_t>::type> parent_map_t;
      std::vector<vertex_t> parent(num_vertices(m_graph));
      typedef graph_as_tree<graph_t, parent_map_t> tree_t;

      tree_t t(m_graph, root_vertex, make_iterator_property_map(parent.begin(), 
                                                get(vertex_index, m_graph)));
      
      tree_printer<StreamType> vis( stream );
      traverse_tree( root_vertex, t, vis );
      stream << endl;
    } // showGraph( StreamType stream )

    void
    showGraph ()
    {
      showGraph( std::cout );
    } // showGraph()

  protected:
    void addRootToGraph ()
    {
      vertex_t root = add_vertex( m_graph );
      // Name it.
      vertex_name_map_t name = get( vertex_name, m_graph );
      name[ root ] = "R";
    } // addRootToGraph()

    /**
     * Create a vector of ints to pass to the ProfileTreeInternalNode
     * constructor indicating that the child will have the same length as the
     * parent and will replace every parent position with its own.
     */
    vector<int>
    createDefaultParentPositionVariations ( uint32_t parent_length )
    {
      // Replace every position of the parent profile with a different one in
      // the child profile.
      vector<int> variations = vector<int>( parent_length + 1, -1 );
      // The first index is special: 0 means don't add new positions before
      // the first position.
      variations[ 0 ] = 0;

      // TODO: REMOVE
      //cout << "IN createDefault..( " << parent_length << " ): variations are: [ ";
      //for ( int i = 0; i < variations.size(); i++ ) {
      //  if( i > 0 ) {
      //    cout << ", ";
      //  }
      //  cout << variations[ i ];
      //}
      //cout << " ]" << endl;

      return variations;
    } // createDefaultParentPositionVariations()

    template <class Graph, typename StreamType>
    struct exercise_vertex {
      StreamType & m_stream;
      exercise_vertex(Graph& g_, StreamType & stream) : g(g_), m_stream(stream) { }
      typedef typename graph_traits<Graph>::vertex_descriptor Vertex;
      void operator()(const Vertex& v) const
      {
        using namespace boost;
        typename property_map<Graph, vertex_index_t>::type 
          vertex_id = get(vertex_index, g);
        m_stream << "vertex: " << get(vertex_id, v) << std::endl;
    
        // Write out the outgoing edges
        m_stream << "\tout-edges: ";
        typename graph_traits<Graph>::out_edge_iterator out_i, out_end;
        typename graph_traits<Graph>::edge_descriptor e;
        for (tie(out_i, out_end) = out_edges(v, g); 
             out_i != out_end; ++out_i)
        {
          e = *out_i;
          Vertex src = source(e, g), targ = target(e, g);
          m_stream << "(" << get(vertex_id, src)
                    << "," << get(vertex_id, targ)
                    << ") ";
        }
        m_stream << std::endl;
    
        // Write out the incoming edges    
        m_stream << "\tin-edges: ";
        typename graph_traits<Graph>::in_edge_iterator in_i, in_end;
        for (tie(in_i, in_end) = in_edges(v, g); in_i != in_end; ++in_i)
        {
          e = *in_i;
          Vertex src = source(e, g), targ = target(e, g);
          m_stream << "(" << get(vertex_id, src)
                    << "," << get(vertex_id, targ) << ") ";
        }
        m_stream << std::endl;
    
        // Write out all adjacent vertices    
        //m_stream << "\tadjacent vertices: ";
        //typename graph_traits<Graph>::adjacency_iterator ai, ai_end;
        //for (tie(ai,ai_end) = adjacent_vertices(v, g);  ai != ai_end; ++ai)
        //  m_stream << get(vertex_id, *ai) <<  " ";
        //m_stream << std::endl;
      }
      Graph& g;
    }; // End inner struct exercise_vertex

    template <typename StreamType>
    class tree_printer {
      StreamType & m_stream;

    public:
      tree_printer ( StreamType & stream ) :
        m_stream( stream )
      {
        // Do nothing else
      } // <init>( StreamType )

      template <typename Node, typename Tree> 
      void preorder(Node, Tree&) {
        m_stream << "(";
      }
      template <typename Node, typename Tree> 
      void inorder(Node n, Tree& t)
      {
        m_stream << get(boost::vertex_name, t)[n];
      }
      template <typename Node, typename Tree> 
      void postorder(Node, Tree&) {
        m_stream << ")";
      }
      
    }; // End inner class tree_printer

  }; // End class ProfileTree

} // End namespace galosh

#endif // __GALOSH_PROFILETREE_HPP__
