// ====================================================================================================
//
// compute_tanimoto.cpp
//
// author : Antoine Allard (antoine.allard.1@gmail.com)
// affil. : Université Laval, Québec, Canada
// www    : dynamica.phy.ulaval.ca (or is.gd/allard)
// created: 2011/07/??
// modif. : 2011/10/31
//
// This program computes the Tanimoto coefficients between contiguous edges.
//
// This code is mostly identical to the original "clusterJaccsFile.cpp" from Jim Bagrow (see below).
// New density definitions have been added and inputs/ouputs have been modified.
// Most of the modifications I've made are identified by "AA"
//
// ====================================================================================================
// Original disclaimer by Jim Bagrow

// clusterJaccsFile.cpp
// Jim Bagrow
// Last Modified: 2009-03-10

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
//
//  g++ -O5 -o cluster clusterJaccsFile.cpp
//  ./cluster network.pairs network.jaccs network.clusters threshold
//   
//     -- network.pairs is an integer edgelist (one edge, two nodes
//     per line)
//     
//     -- network.jaccs contains network.jaccs the jaccard 
//     coefficient for each pair of edges compared, of the form:
//         i_0 i_1 j_0 j_1 jaccard<newline>
//         ...
//     for edges (i_0,i_1) and (j_0,j_1), etc.
//     
//     -- network.clusters will contain one cluster of edges per line 
//     (edge nodes are comma-separated and edges are space-separated)
//     
//     --threshold is the [0,1] threshold for the clustering

// CORRECTNESS:
//  ??????

// OPTIMIZATION:
//  * Probably pretty slow with all the O(log n) .find map lookups
//  * Better to store edge2iterator? (edge2pointer really)  But do 
//  iterators get invalidated after an insert?
//  * According to http://en.wikipedia.org/wiki/Map_(C%2B%2B_container):
//  "Iterators are not invalidated by insert and erase operations
//  which don't remove the object to which the iterator points."
//  implemented edge2iter but there's no speed up!  Should get fancy...

// ====================================================================================================

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <set>
#include <map>
#include <utility>   // for pairs
#include <algorithm> // for swap
#include <ctime>     // to keep track of the time elaspsed
using namespace std;


int main (int argc, char const *argv[])
{

  // value of the threshold
  double threshold = atof( argv[2] );

  // names of the files that are used as input/ouput
  std::string wpairs_filename(argv[1]); wpairs_filename += ".wpairs";       // input file
  std::string tanimoto_filename(argv[1]); tanimoto_filename += ".tanimoto"; // input file
  std::string clusters_filename(argv[1]); clusters_filename += ".clusters_"; clusters_filename += argv[2]; // output file
  std::string cluster_stats_filename(argv[1]); cluster_stats_filename += ".cluster_stats_"; cluster_stats_filename += argv[2]; // output file
  std::string density_filename(argv[1]); density_filename += ".density_"; density_filename += argv[2]; // output file

  // makes sure that the threshold is in [0,1]
  if (threshold < 0.0 || threshold > 1.0)
  {
    cout << "ERROR: specified threshold not in [0,1]" << endl;
    exit(1);
  }
  

  // objects to compute the duration of the calculation
  time_t start_time, stop_time;
  time(&start_time);    
    
  
  // starts loads edgelist
  ifstream inFile;
  inFile.open( wpairs_filename.c_str() );
  if (!inFile)
  {
    cout << "ERROR: unable to open input file" << endl;
    exit(1); // terminate with error
  }

  // index should be iterator not integer????
  map< int,           set<pair<int,int> > > index2cluster; // O(log n) access too slow?
  map< pair<int,int>, map< int,set<pair<int,int> > >::iterator > edge2iter;
  std::map< pair<int,int>, double > edge2weight; // AA

  int ni, nj, index = 0;
  double weight, Mw(0);

  while (inFile >> ni >> nj >> weight)
  { // scan edgelist to populate
    if(ni >= nj)
      swap(ni,nj); // undirected!
        
    index2cluster[ index ].insert( make_pair(ni,nj) );         // build cluster index to set of edge-pairs map
    edge2iter[ make_pair(ni,nj) ] = index2cluster.find(index); // build edge pair to cluster iter map ******????
    index++;

    edge2weight[ make_pair(ni,nj) ] = weight; // AA
    Mw += weight; // AA
  }
  inFile.close();
  inFile.clear();


  //************* loop over jaccards file and do the clustering
  ifstream jaccFile;  jaccFile.open( tanimoto_filename.c_str() );
  if (!jaccFile)
  {
    cout << "ERROR: unable to open tanimoto file" << endl;
    exit(1); // terminate with error
  }

  int i0,i1,j0,j1; double jacc;
  int idx_i, idx_j;
  map< int, set<pair<int,int> > >::iterator iter_i,iter_j;
  set<pair<int,int> >::iterator iterS;
  while ( jaccFile >> i0 >> i1 >> j0 >> j1 >> jacc )
  {
    if ( jacc >= threshold )
    {
      if (i0 >= i1)
        swap(i0,i1); // undirected!
      if (j0 >= j1)
        swap(j0,j1); // undirected!
            
      iter_i = edge2iter[ make_pair(i0,i1) ];
      iter_j = edge2iter[ make_pair(j0,j1) ];
      if ( iter_i != iter_j )
      {
        // always merge smaller cluster into bigger:
        if ( (*iter_j).second.size() > (*iter_i).second.size() )
          swap(iter_i, iter_j);

        // merge cluster j into i and update index for all elements in j:
        for (iterS = iter_j->second.begin(); iterS != iter_j->second.end(); iterS++)
        {
          iter_i->second.insert( *iterS );
          edge2iter[ *iterS ] = iter_i;
        }
                
        // delete cluster j:
        index2cluster.erase(iter_j);
      } 
    } // done merging clusters i and j
  }
  jaccFile.close();
  //************* done looping over jaccards file

    
    //************* write the clusters to file:
    //cout << "There were " << index2cluster.size() << " clusters at threshold " << threshold << "." << endl;
    
    // all done clustering, write to file (and calculated partition density):
    FILE * clustersFile     = fopen( clusters_filename.c_str(), "w" );
    FILE * clusterStatsFile = fopen( cluster_stats_filename.c_str(), "w" );
    
    set<int> clusterNodes;
    int mc, nc;
    int M = 0, Mns = 0;
    double wSum = 0.0;
    double m_cw_c, average_weighted_density(0.0), average_weighted_density_var(0.0); // AA
    
    set< pair<int,int> >::iterator S;
    map< int,set<pair<int,int> > >::iterator it;
    
    for ( it = index2cluster.begin(); it != index2cluster.end(); it++ ){

        m_cw_c = 0.0; // AA
        clusterNodes.clear();
        for (S = it->second.begin(); S != it->second.end(); S++ ){
            fprintf( clustersFile, "%i,%i ", S->first, S->second ); // this leaves a trailing space...!
            
            clusterNodes.insert(S->first);
            clusterNodes.insert(S->second);

            m_cw_c += edge2weight[ make_pair(S->first,S->second) ]; // AA
        }

        // counts the number of nodes and edges in cluster c, and updates the total number of edges
        mc = it->second.size();
        nc = clusterNodes.size();
        M += mc;
        
        // computes the density according to Ahn et al. definition
        if (nc != 2) {
            Mns  += mc;
            wSum += mc * (mc - (nc-1.0)) / ((nc-2.0)*(nc-1.0));
        }

        // computes the weighted density // AA
        if(nc != 2)
          average_weighted_density += m_cw_c * m_cw_c * (mc - (nc-1.0)) / ((nc-2.0)*(nc-1.0)) / mc;

        // computes the weighted density // AA
        if(nc != 2)
          average_weighted_density_var += m_cw_c * (mc - (nc-1.0)) / ((nc-2.0)*(nc-1.0));
        
       
        fprintf( clustersFile, "\n" );
        fprintf( clusterStatsFile, "%i %i %.5f\n", mc, nc, m_cw_c/mc);
    }
    fclose(clustersFile);
    fclose(clusterStatsFile);
    //*************

  // writes the average cluster density in the different file (for plotting)
  std::ofstream density_file( density_filename.c_str() );
  density_file << std::fixed << std::setprecision(7);
  density_file << threshold << "\t";                                   // value of the threshold used
  density_file << M << "\t";                                           // number of edges found in communities
  density_file << 2.0 * wSum / M << "\t";                              // density def. #1 (see Progress Report II)
  density_file << Mns << "\t";                                         // number of edges found in communities of size greater than 2
  density_file << 2.0 * wSum / Mns << "\t";                            // density def. #2 (see Progress Report II)
  density_file << index2cluster.size() << "\t";                        // number of communities detected
  density_file << 2.0 * average_weighted_density / Mw << "\t";         // density def. #5 (see Progress Report II)
  density_file << 2.0 * average_weighted_density_var / Mw << "\t";     // density def. #3 (see Progress Report II)
  density_file << 2.0 * average_weighted_density_var / M << std::endl; // density def. #4 (see Progress Report II)
  density_file.close();

  // gets the stop time and displays the duration of the calculation
  time(&stop_time); 
  std::cout << "Elapsed time:\t" << difftime(stop_time,start_time) << " seconds (";
  std::cout << std::fixed << std::setprecision(2) << difftime(stop_time,start_time)/3600 << " hrs.)" << std::endl;
  std::cout.flush();  

  return 0;
}
