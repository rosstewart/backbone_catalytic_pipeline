#include "simgraph.h"
#include "string.h"
#include <fstream>
#include <set>
#include <queue>


SimpleGraph SimpleGraph::read_graph(const char *nlabels_file, const char *graph_file)  {
    ifstream lin(nlabels_file, ios::in);
    ifstream gin(graph_file, ios::in);

    if (lin.fail())  {
        cerr << "ERROR: Node labels file " << nlabels_file << " could not be opened." << endl; exit(1);
    }

    if (gin.fail())  {
        cerr << "ERROR: Graph file " << graph_file << " could not be opened." << endl; exit(1);
    }

    SimpleGraph g;

    // Read graph labels
    string line;
    if (getline(lin,line))
        g.nodes = line;    

    // Size of the graph is determined by the number of node labels
    g.adj.resize(g.nodes.size());

    // Read graph
    // Parse adjacency list
    while(getline(gin, line))  {
        vector<string> tokens = split(line, '\t');
        unsigned i = to_i(tokens[0]);

        if (i >= g.adj.size())  {
            cerr << "ERROR: Node index " << i << " >= graph size " << g.adj.size() << " in graph file " << graph_file << "." << endl; exit(1);
        }

        for (unsigned j=1; j<tokens.size(); j++)  {
            unsigned k = to_i(tokens[j]);
            bool found(false);

            if (k >= g.adj.size())  {
                cerr << "ERROR: Node index " << k << " >= graph size " << g.adj.size() << " in graph file " << graph_file << "." << endl; exit(1);
            }

            for (vector<unsigned>::iterator it = g.adj[i].begin(); it < g.adj[i].end(); it++)  {
                if (*it == k) {
                    found = true;
                }
            }
            if (!found)  {
                g.adj[i].push_back(k);
            }
        }
    }

    return g;
}


void SimpleGraph::print_dot(ostream &out)  {
    out << "graph PCG {\n";
    out << "    node [shape=circle,style=filled,color=lightgray]; {node [label=\"" << nodes[0] << "\"] n0; }\n";
    out << "\n";

    for (unsigned i=0; i<adj.size(); i++)
        for (unsigned j=0; j<adj[i].size(); j++)
            if (i<adj[i][j])
                out << "n" << i << " -- n" << adj[i][j] << ";\n";

    out << "}";
}


void SimpleGraph::print_dot(ostream &out, const vector<unsigned> &dist, unsigned bound)  {
    out << "graph PCG {\n";
    out <<  "   graph [splines=true overlap=false]\n";
    out << "    node [shape=doublecircle]; {node [label=\"" << nodes[0] << "\"] n0; }\n";
    out << "    node [shape=circle];\n";
    out << "\n";

    for (unsigned i=0; i<nodes.size(); i++)
        if (dist[i]<bound)
            out << "    n" << i << " [label=\"" << nodes[i] << "\"];\n";
    out << "\n";

    for (unsigned i=0; i<bound; i++)  {
        out << "{rank=same;"; 
        for (unsigned j=0; j<dist.size(); j++)
            if (i == dist[j])
                out << " n" << j; 
        out << ";}\n";  
    }
    out << "\n";

    for (unsigned i=0; i<adj.size(); i++)
        for (unsigned j=0; j<adj[i].size(); j++)
            if (i<adj[i][j] && dist[i]<bound && dist[adj[i][j]]<bound)  {
                out << "    n" << i << " -- n" << adj[i][j] << ";\n";
            }

    out << "}";
}


vector<unsigned> SimpleGraph::breadth_first_sort(unsigned root) const  {
    vector<unsigned> dist(nodes.size(), UINT_MAX);    
    dist[root] = 0;
   
    queue<unsigned> Q;
    Q.push(root);
    
    while (!Q.empty())  {
        unsigned i = Q.front();
        
        for (unsigned j=0; j<adj[i].size(); j++)  {
            unsigned k = adj[i][j];
            
            if (UINT_MAX==dist[k])  {
                dist[k] = dist[i] + 1;
                Q.push(k);
            }
        }
        Q.pop();
    }
    return dist;
}

