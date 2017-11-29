// ====================================================================================================
//
// adjmat2wpairs.cpp
//
// author : Antoine Allard (antoine.allard.1@gmail.com)
// affil. : Université Laval, Québec, Canada
// www    : dynamica.phy.ulaval.ca (or is.gd/allard)
// created: 2011/07/??
// modif. : 2011/10/31
//
// This program computes a "wpairs" file from an weighted adjacency matrix (or edge list).
//
// ====================================================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <cstdlib>


int main(int argc, char* argv[])
{

  std::cout << "Formating adjacency matrix file into the weighted pairs (.wpairs) format...";
  std::cout.flush();

  /// ====================================================================================================
  /// definition of useful objets/variables

  // output streams
  std::string pairs_weighted_filename(argv[1]); pairs_weighted_filename += ".wpairs";
  std::string id2name_filename(argv[1]); id2name_filename += ".numid2name";

  // indicates whether the format of the input file is an adjacency matrix (0, default) or an edge list (1)
  int input_format(0);
  if(argc >= 3)
    input_format = atoi(argv[2]);

  std::ifstream in_file(argv[1]);
  std::ofstream pairs_weighted_file(pairs_weighted_filename.c_str());
  std::ofstream id2name_file(id2name_filename.c_str());
  pairs_weighted_file.precision(15);

  // variables counting the number of nodes and edges
  int nb_nodes(0), nb_edges(0);

  // vector containing the nodes' name
  std::vector<std::string> names;
  names.reserve(2400); // has been added because we know the maximum size of the data set

  int node1, node2;
  std::string name1, name2;
  double weight, num_zero(1e-12); // default threshold value
  bool new_node;
  time_t start_time, stop_time;
  time(&start_time);

  if(argc == 4)
    if(atof(argv[3]) > num_zero)
      num_zero = atof(argv[3]);

  /// ====================================================================================================
  /// reading files and checking if everything's okay

  if (!in_file)
  {
    std::cout << "ERROR: unable to open input file" << std::endl;
    abort();
  }

  if (!pairs_weighted_file)
  {
    std::cout << "ERROR: unable to open output file" << std::endl;
    abort();
  }

  if (!id2name_file)
  {
    std::cout << "ERROR: unable to open output file" << std::endl;
    abort();
  }


  // reading stream
  if(in_file.is_open() && pairs_weighted_file.is_open() && id2name_file.is_open())
  {

    // for every nonempty line of the adjacency matrix file
    while (in_file >> name1 >> name2 >> weight)
    {

      if(weight > num_zero)
      {
        // check if "node1" has already been registered (i.e. received a numerical node id)
        new_node = true;
        for(int i(0); i<nb_nodes; ++i)
        {
          if(name1 == names[i])
          {
            new_node = false;
            node1 = i;
            continue;
          }
        }

        // if "node1" is a new node, adds it to the list
        if(new_node)
        {
          node1 = nb_nodes;
          ++nb_nodes;
          names.push_back(name1);
          id2name_file << node1 << "\t" << name1 << std::endl;
        }

        // check if "node2" has already been registered (i.e. received a numerical node id)
        new_node = true;
        for(int i(0); i<nb_nodes; ++i)
        {
          if(name2 == names[i])
          {
            new_node = false;
            node2 = i;
            continue;
          }
        }

        // if "node2" is a new node, adds it to the list
        if(new_node)
        {
          node2 = nb_nodes;
          ++nb_nodes;
          names.push_back(name2);
          id2name_file << node2 << "\t" << name2 << std::endl;
        }

        // writes the pair and weight
        // CAUTION: assumes that the adjacency matrix is symmetrical and that edges are undirected (i.e. a same edge appears twice in the file)!

        if(input_format == 0) // if an adjacency matrix is provided
        {
          if(node1 <= node2)
          {
            ++nb_edges;
            pairs_weighted_file << std::fixed << node1 << "\t\t" << node2 << "\t\t" << weight << std::endl;
          }
        }

        if(input_format == 1) // if an edge list is provided
        {
          if(node1 < node2)
          {
            ++nb_edges;
            pairs_weighted_file << std::fixed << node1 << "\t\t" << node2 << "\t\t" << weight << std::endl;
          }
          else if(node1 > node2)
          {
            ++nb_edges;
            pairs_weighted_file << std::fixed << node2 << "\t\t" << node1 << "\t\t" << weight << std::endl;
          }
        }
      }

    }

    in_file.close();
    pairs_weighted_file.close();
    id2name_file.close();

    time(&stop_time);

    std::cout << "done!" << std::endl;
    std::cout << std::endl;
    std::cout << "Number of nodes:\t" << nb_nodes << std::endl;
    std::cout << "Number of edges:\t" << nb_edges << std::endl;
    std::cout << std::endl;
    std::cout << "Elapsed time:\t" << difftime(stop_time,start_time) << " seconds (" << difftime(stop_time,start_time)/3600 << " hrs.)" << std::endl;

  }

    return 0;
}
