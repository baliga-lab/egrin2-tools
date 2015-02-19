// ====================================================================================================
//
// compute_tanimoto.cpp
//
// author : Antoine Allard (antoine.allard.1@gmail.com)
// affil. : Université Laval, Québec, Canada
// www    : dynamica.phy.ulaval.ca (or is.gd/allard)
// created: 2011/07/13
// modif. : 2011/08/08
//
// This program computes the Tanimoto coefficients between contiguous edges.
//
// This code has been greatly inspired by "calcAndWrite_Jaccard.ccp" from Jim Bagrow (see below).
//
// ====================================================================================================
// Original disclaimer by Jim Bagrow

// calcAndWrite_Jaccards.cpp
// Jim Bagrow
// Last Modified: 2008-12-30

/*
Copyright 2008,2009,2010 James Bagrow


This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


// USAGE:
//      g++ -O5 -o calc calcAndWrite_Jaccards.cpp
//      ./calc network.pairs network.jaccs
//
//  -- network.pairs is an integer edgelist (one edge, two nodes
//  per line)
//  -- network.jaccs will contain the jaccard coefficient for each
//  pair of edges compared, of the form:
//      i_0 i_1 j_0 j_1 jaccard<newline>
//      ...
//  for edges (i_0,i_1) and (j_0,j_1)


// same as calcAndWrite_Jaccards.cpp, but without the graph code...
// but not much quicker (:-()

// all this does is calculate the jaccard for "each" edge pair and
// write it to a file.  Two make this into real code will take some
// more work (the next step can be the hierarchical clustering over
// the output jacc file...)

// CORRECTNESS:
//  Returns same jaccard file as calcAndWrite_Jaccards.cpp, as shown
//  by compareTwoJaccs.py...

// ====================================================================================================

#include <fstream>
#include <iostream>
#include <set>
#include <utility>   // for pairs
#include <string>
#include <cstdlib>   // for exit
#include <vector>
//#include <cmath>     // for pow
#include <boost/math/special_functions/binomial.hpp> // for binomial coefficient
#include <ctime>     // to keep track of the time elaspsed
#include "boost/multi_array.hpp"
#include <boost/lexical_cast.hpp>


// ====================================================================================================
// function computing the inner product between two adjacency lists
double inner_product( const std::set<std::pair<int,double> > &A, const std::set<std::pair<int,double> > &B )
{
  // only get number of elements, don't build a new set
  // this assumes the sets are ordered, which std::sets are!
  double value = 0.0;
  std::set<std::pair<int,double> >::const_iterator As = A.begin(), Af = A.end(),
                                                   Bs = B.begin(), Bf = B.end();
  while( As != Af && Bs != Bf )
  {
    if( As->first < Bs->first )
      ++As;

    else if( Bs->first < As->first )
      ++Bs;

    else // happens if the two nodes have a common neighbors
    {
      value += As->second * Bs->second;
      ++As;
      ++Bs;
    }
  }

  return value;
}



// ====================================================================================================
// main function
int main (int argc, char const *argv[])
{

  // names of the files that are used as input/ouput
  std::string wpairs_filename(argv[1]); wpairs_filename += ".wpairs";       // input file
  std::string tanimoto_filename(argv[1]); tanimoto_filename += ".tanimoto"; // output file


  // objects to compute the duration of the calculation
  time_t start_time, stop_time;
  time(&start_time);


  // open the input file to load the network
  std::ifstream wpairs_file(wpairs_filename.c_str());


  // check if the input file as been openned correctly
  if (!wpairs_file)
  {
    std::cout << "ERROR: unable to open input file" << std::endl;
    exit(1); // terminate with error
  }


  // choice of definition for the Tanimoto Coefficient
  const int tanimoto_definition(atoi(argv[2]));
  if(tanimoto_definition == 0)
    std::cout << "Using Ahn et al. definition for the Tanimoto Coefficient" << std::endl;
  else if(tanimoto_definition == 1)
    std::cout << "Using Kalinka et al. definition for the Tanimoto Coefficient" << std::endl;
  else
    std::cout << "The definition of the Tanimoto Coefficients to be computed has not been entered properly (must be either 0 or 1)." << std::endl;
  


  // scan the input file once to compute the number of nodes
  std::cout << std::endl << "Counting the number of nodes..." << std::endl;
  std::cout.flush();
  int node1, node2, num_nodes, max_node(-1);
  double weight;
  while (wpairs_file >> node1 >> node2 >> weight)
  {
    if (node1 > max_node){  max_node = node1;  }
    if (node2 > max_node){  max_node = node2;  }
  }
  num_nodes = max_node + 1;
  wpairs_file.close();


  // set the keystone range covered during this calculation
  int first_node(atoi(argv[3])), last_node(atoi(argv[4]));
  if(last_node>num_nodes)
    last_node = num_nodes;
  tanimoto_filename += "_";
  tanimoto_filename += argv[3];
  tanimoto_filename += "_";
  tanimoto_filename += boost::lexical_cast<std::string>(last_node-1);

  // ===========================================================================
  // builds the network

  // displays the filename of the file containing the network
  std::cout << "Building the network from file: " << wpairs_filename << std::endl;
  std::cout.flush();

  // initiate the array that stores the network itself
  std::set<std::pair<int,double> > *nodes = NULL;
  nodes = new std::set<std::pair<int,double> >[num_nodes];


  // fills the network structure from the input file
  wpairs_file.open( wpairs_filename.c_str() );
  while(wpairs_file >> node1 >> node2 >> weight)
  {
    nodes[node1].insert(std::make_pair(node2,weight));
    nodes[node2].insert(std::make_pair(node1,weight));
  }
  wpairs_file.close();
  

  // makes the nodes' adjacency list inclusive
  std::set<std::pair<int,double> >::const_iterator it,end;
  if(tanimoto_definition == 0) // using the definition of Ahn et al.
  {  
    int k;
    double avg_weight;
    for(int i(0); i<num_nodes; ++i)
    {
      // computes the average weights of the node's edges
      avg_weight = 0;
      it = nodes[i].begin();
      end = nodes[i].end();
      for(; it!=end; ++it)
        avg_weight += it->second;
      avg_weight /= nodes[i].size();

      // adds the nodes itself to the adjacency list
      nodes[i].insert( std::make_pair(i,avg_weight) );
    }
  }

  if(tanimoto_definition == 1) // using the definition of Kalinka et al.
    for(int i(0); i<num_nodes; ++i)
      nodes[i].insert( std::make_pair(i,1.0) );


  // computes the number of contiguous edge pairs to be considered (number of lines in the output file)
  unsigned long int num_edges(0), num_edge_pairs(0), deg, progress_step, counter(0);
  double avg_degree(0);
  for(int i(0); i<num_nodes; ++i)
  {
    deg = nodes[i].size()-1; // the neighborhood of nodes is inclusive
    num_edges += deg;
    if(first_node <= i && i < last_node)
    {
      if(deg>1)
        num_edge_pairs += boost::math::binomial_coefficient<double>(deg, 2);
    }
  }
  avg_degree = num_edges / static_cast<double>(num_nodes); // needs to be calculated before "num_edges" gets divided by two
  num_edges /= 2;
  progress_step = num_edge_pairs / 100 + 1; // will be used to display the progress of the calculation


  // displays data about the network and the calculation
  std::cout << "  - number of nodes:\t\t\t" << num_nodes << std::endl;
  std::cout << "  - number of edges:\t\t\t" << num_edges << std::endl;
  std::cout << "  - average degree of nodes:\t\t" << std::fixed << std::setprecision(2) << avg_degree << std::endl;
  std::cout << "Considering keystones in [" << first_node << "," << last_node-1 << "]" << std::endl;
  std::cout << "  - number of contiguous edge pairs:\t" << num_edge_pairs << std::endl;


//// Version A: the inner_product between two different nodes is computed for each contiguous pairs of edges
//  // ===========================================================================
//  // computes the squared norm of the vector containing each nodes' neighbors
//  std::vector<double> squared_norm(num_nodes);
//  for(int i(0); i<num_nodes; ++i)
//    squared_norm[i] = inner_product(nodes[i],nodes[i]);

// Version B: the inner_product is computed for each edge pair beforehand
  // ===========================================================================
  // computes the inner_product for each pair of nodes (faster computation time for medium-size networks)
  // ATTENTION: ASSUMES UNDIRECTED NETWORK
  std::cout << "Computing the 'inner product matrix'...";
  std::cout.flush();
  int inner_prod_mat_numel = boost::math::binomial_coefficient<double>(num_nodes+1,2);
  progress_step = inner_prod_mat_numel / 100 + 1;
  boost::multi_array<double,2> inner_prod_mat(boost::extents[num_nodes][num_nodes]);
  for(int i(0); i<num_nodes; ++i)
  {
    for(int j(i); j<num_nodes; ++j)
    {
      if(counter % progress_step == 0)
      {
        std::cout << "\rComputing the 'inner product matrix'... ";
        std::cout << std::fixed << std::setprecision(0) << 100*static_cast<double>(counter)/static_cast<double>(inner_prod_mat_numel) << "%";
        std::cout.flush();
      }
      inner_prod_mat[i][j] = inner_product(nodes[i],nodes[j]);
      ++counter;
    }
  }
  std::cout << "\rComputing the 'inner product matrix'... 100%" << std::endl;

/*
  // displays the network's adjacency list
  for(int i(0); i<num_nodes; ++i)
  {
    it = nodes[i].begin();
    end = nodes[i].end();
    std::cout<<"Node "<<i<<std::endl;
    for(; it!=end; ++it)
      std::cout  <<  it->first  <<  "  "  <<  it->second  <<  std::endl;
    std::cout<<"Squared norm: "<<squared_norm[i]<<std::endl;
    std::cout<<std::endl;
  }
*/


  // ===========================================================================
  // the calculation itself

  std::cout << "Computing edge similarities...0%";
  std::cout.flush();
  progress_step = num_edge_pairs / 100 + 1; // will be used to display the progress of the calculation
  counter = 0;

  FILE * jaccFile = fopen(tanimoto_filename.c_str(),"w");
  double tanimoto_coeff, inner_prod;
  std::set<std::pair<int,double> >::const_iterator node1_it, node2_it;

  for (int keystone(first_node); keystone<last_node; ++keystone)
  { // loop over keystones 

    end = nodes[keystone].end();

    // loop over all neighbors of keystone
    node1_it = nodes[keystone].begin();
    for(; node1_it!=end; ++node1_it)
    {
      node1 = node1_it->first;
      if(node1 == keystone)
        continue;

      // loop over all neighbors of keystone
      node2_it = nodes[keystone].begin();
      for(; node2_it!=end; ++node2_it)
      {
        node2 = node2_it->first;
        if(node2 == keystone || node1 >= node2)
          continue;

        // updates the displayed counter is necessary
        if(counter % progress_step == 0)
        {
          std::cout << "\rComputing edge similarities...  ";
          std::cout << std::fixed << std::setprecision(0) << 100*static_cast<double>(counter)/static_cast<double>(num_edge_pairs-1) << "%";
          std::cout.flush();
        }
        ++counter;

// Version A: the inner_product between two different nodes is computed for each contiguous pairs of edges
//        // this point being reached means that "keystone-node1" and "keystone-node2" are a proper contiguous edge pair
//        // computes the Tanimoto coefficient for this pair of edges
//        inner_prod = inner_product(nodes[node1],nodes[node2]);
//        tanimoto_coeff = inner_prod / (squared_norm[node1] + squared_norm[node2] - inner_prod);

// Version B: the inner_product is computed for each edge pair beforehand
        // this point being reached means that "keystone-node1" and "keystone-node2" are a proper contiguous edge pair
        // computes the Tanimoto coefficient for this pair of edges
        tanimoto_coeff = inner_prod_mat[node1][node2] / (inner_prod_mat[node1][node1] + inner_prod_mat[node2][node2] - inner_prod_mat[node1][node2]);


        // writes the tanimoto coefficient in the output file
        if(keystone < node1 && keystone < node2)
        {
          fprintf( jaccFile, "%i\t%i\t%i\t%i\t%.5f\n", keystone, node1, keystone, node2, tanimoto_coeff );
        } 
        else if (keystone < node1 && keystone > node2)
        {
          fprintf( jaccFile, "%i\t%i\t%i\t%i\t%.5f\n", keystone, node1, node2, keystone, tanimoto_coeff );
        }
        else if (keystone > node1 && keystone < node2)
        {
          fprintf( jaccFile, "%i\t%i\t%i\t%i\t%.5f\n", node1, keystone, keystone, node2, tanimoto_coeff );
        }
        else
        {
          fprintf( jaccFile, "%i\t%i\t%i\t%i\t%.5f\n", node1, keystone, node2, keystone, tanimoto_coeff );
        }
      }
    }
  } // done loop over keystones

  fclose(jaccFile);

  delete [] nodes;


  // updates the information displayed
  std::cout << "\rComputing edge similarities... 100%" << std::endl;
  std::cout << "  - results written in: " << tanimoto_filename << std::endl;

  // gets the stop time and displays the duration of the calculation
  time(&stop_time);
  std::cout << "Elapsed time:\t" << difftime(stop_time,start_time) << " seconds (";
  std::cout << std::fixed << std::setprecision(2) << difftime(stop_time,start_time)/3600 << " hrs.)" << std::endl << std::endl;

  return 0;
}
