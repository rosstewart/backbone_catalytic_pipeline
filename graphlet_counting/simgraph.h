/**
 * Simple Graph - auxiliary data structure. 
 *
 * Modified from Vladimir Vacic graphlet kernel
 *
 * Modified by:
 * Jose Lugo-Martinez, jlugomar@indiana.edu
 * School of Informatics and Computing
 * Indiana University-Bloomington
 *
 * Jun-24-2013
 *
 * Copyright (c) 2014 Jose Lugo-Martinez,
 * Vladimir Vacic, and Predrag Radivojac.
 *
 */


#ifndef __SIMPLE_GRAPH_H__
#define __SIMPLE_GRAPH_H__

#include <iostream>
#include <vector>
using namespace std;


class SimpleGraph  {
public:
    SimpleGraph()  {}
    ~SimpleGraph()  {}    

    /** Read an adjacency list file and create a graph. */
    static SimpleGraph read_graph(const char*, const char*);

    /** Generates a GraphViz file. */
    void print_dot(ostream&);

    /** Generates a GraphViz file, but includes only */
    void print_dot(ostream &out, const vector<unsigned>&, unsigned);

    /** Breadth-first assignment of distances from a source node. */ 
    vector<unsigned> breadth_first_sort(unsigned) const;

    string nodes;                  // Vertex labels.
    vector<vector<unsigned> > adj; // Adjacency lists.
};

#endif

