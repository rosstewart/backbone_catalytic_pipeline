/**
 * Generates a kernel matrix according to the
 * user-specified graph-based kernel method
 * over set of vertices of interest in a graph.
 *
 * Modified from Vladimir Vacic graphlet kernel
 *
 * Jose Lugo-Martinez, jlugomar@indiana.edu
 * School of Informatics and Computing
 * Indiana University-Bloomington 
 *
 * Aug-6-2012
 *
 * Copyright (c) 2014 Jose Lugo-Martinez,
 * Vladimir Vacic, and Predrag Radivojac.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "gkernel.h"
#include "string.h"
#include <iostream>
#include <fstream>
using namespace std;


void print_help()  {
    cout << "Usage: run_kernel -p FILE -n FILE -g G_FILE -l L_FILE -t TYPE -[k|s] OUTPUT [...]\n";
    cout << "Options:\n\n";

    cout << "  -h         Displays this message.\n\n";

    cout << "  -t TYPE    Kernel type (0-Cumulative Random Walk, 1-Random Walk, 2-Standard Graphlet, 3-Label Substitutions Graphlet, 4-Edge Indels Graphlet, 5-Edit Distance Graphlet).\n";
    cout << "             Defaults to standard graphlet.\n\n";

    cout << "  -p FILE    List of positive (vertices) examples.\n";
    cout << "  -n FILE    List of negative (vertices) examples.\n";
    cout << "  -g G_FILE  Input graph file.\n";
    cout << "  -l L_FILE  Vertex labels file for input graph.\n\n";

    cout << "  -N         Normalize the kernel matrix.\n";
    cout << "             Defaults to false.\n\n";

    cout << "  -k KERNEL  Output file for the kernel matrix in standard output.\n";
    cout << "   or\n";
    cout << "  -s SPARSE  Output file for the sparse attribute matrix (SVML).\n";
    cout << "             Defaults to KERNEL.\n\n";

    cout << "  -I STEPS   Number of steps. (Needed for Random Walk Kernels)\n";
    cout << "             Defaults to 100,000 steps.\n\n";

    cout << "  -R RESTART Restart probability. (Needed for Random Walk Kernels)\n\n";
    cout << "             Defaults to 0.3\n\n";

    cout << "  -S SIMMAT  Similarity matrix for weighting each possible vertex label substitution. (Needed for Label Substitutions and Edit Distance Kernels)\n";
    cout << "             Defaults to uniform weights (i.e. all label mismatches are equally weighted to 1).\n\n";

    cout << "  -M LABMIS  Fraction of nodes in the n-graphlet allowed to have vertex label mismatches. (Needed for Label Substitutions and Edit Distance Kernels)\n";
    cout << "             For example, M=0.34 allows 1-vertex label mismatch for n=3,4 and 5 graphlets, whereas M=0.26 allows 1-vertex label mismatch for n=4 and 5 graphlets.\n";
    cout << "             Defaults to 0.0 (i.e. no vertex label mismatches allowed).\n\n";

    cout << "  -E EDGMIS  Total number of edge insertions and/or deletions allowed between graphlets. (Needed for Edge Indels Kernel and Edit Distance Kernels)\n";
    cout << "             Defaults to 0.\n\n";

    cout << "  -A ALPHA   Vertex labels alphabet over problem statement is defined (i.e. all possible labels for a vertex). (Needed for Label Substitutions and Edit Distance Kernels)\n";
    cout << "             Defaults to 20 amino acid residues alphabet.\n\n";

    cout << "  -c LABELS  Output file for each example class label.\n\n";

    cout << "  -v         Verbose (prints progress messages).\n\n"; 
}

int main(int argc, char* argv[])  {
    typedef enum kerneltype  {
        RANDOM_WALK_CUMULATIVE,
        RANDOM_WALK,
        STANDARD_GRAPHLET,
        LABEL_MISMATCH,
        EDGE_MISMATCH,
        EDIT_DISTANCE
     } KernelType;

    typedef enum outformat  {
        KERNEL,
        SPARSE_SVML
     } OutputFormat;

    string pos_file;
    string neg_file;
    string l_file, g_file;
    string output_file;
    OutputFormat format(KERNEL);
    KernelType k_type(STANDARD_GRAPHLET);    
    string labels_file;
    bool normalize(false);
    bool verbose(false);

    // Random Walk Kernels Parameters
    int steps(100000); 
    double restart(0.3);

    // Label Substitutions Graphlet Kernel Parameters
    float mismatches(0.0);
    string alphabet; 
    string root_alphabet; 
    string sim_matrix_file;

    // Edge Indels Graphlet Kernel Parameter
    unsigned edgmis(0);

    // Parse command line arguments.
    for (int i=1; i<argc && (argv[i])[0] == '-'; i++)  {
        switch ((argv[i])[1])  {
            case 'h': print_help(); exit(0);
            case 't': 
                i++; 
                switch (to_i(argv[i]))  {
                    case 0:
                        k_type=RANDOM_WALK_CUMULATIVE;
                        break;
                    case 1:
                        k_type=RANDOM_WALK;
                        break;
                    case 2:
                        k_type=STANDARD_GRAPHLET;
                        break;
                    case 3:
                        k_type=LABEL_MISMATCH;
                        break;
                    case 4:
                        k_type=EDGE_MISMATCH;
                        break;
                    case 5:
                        k_type=EDIT_DISTANCE;
                        break;
                    default:
                        k_type=STANDARD_GRAPHLET;
                }
                break;
            case 'p': i++; pos_file=argv[i]; break;
            case 'n': i++; neg_file=argv[i]; break;
            case 'g': i++; g_file=argv[i]; break;
            case 'l': i++; l_file=argv[i]; break;
            case 'N': normalize=true; break;
            case 'k': i++; format=KERNEL; output_file=argv[i]; break;
            case 's': i++; format=SPARSE_SVML; output_file=argv[i]; break;
            // Kernel-specific parameters                    
            case 'I': i++; steps=to_i(argv[i]); break;
            case 'R': i++; restart=to_f(argv[i]); break;

            case 'S': i++; sim_matrix_file=argv[i]; break;
            case 'M': 
                i++; 
                mismatches=to_f(argv[i]);
                if(mismatches < 0.0 || mismatches > 1.0)  {
                    cerr << "ERROR: Fraction of nodes allowed to have vertex label mismatches, M, must be 0<=M<=1, but you entered " << argv[i] << endl;
                    print_help();  exit(1);
                }
                break;
            case 'E':
                i++;
                edgmis=to_i(argv[i]);
                if(edgmis < 0 || edgmis > 2)  {
                    cerr << "ERROR: Total number of edge insertions and deletions must be either 1 or 2, but you entered " << argv[i] << endl;
                    print_help();  exit(1);
                }
                break;
            case 'A': i++; alphabet=argv[i]; root_alphabet=argv[i]; break;
            case 'c': i++; labels_file=argv[i]; break;
            case 'v': verbose=true; break;
            default: 
                cerr << "ERROR: Unknown option " << argv[i] << endl;
                print_help();  exit(1);
        }
    }

    if (0 == output_file.size())  {
        cerr << "ERROR: Output file name not specified." << endl;  print_help();  exit(1);
    }

    if (0 == sim_matrix_file.size() && (k_type == LABEL_MISMATCH || k_type == EDIT_DISTANCE))  {
        // User-defined probability similarity matrix
        // Use default matrix
        sim_matrix_file = "user_defined.matrix";
    }

    if ((k_type == LABEL_MISMATCH || k_type == EDIT_DISTANCE) && mismatches > 0.0 && (0 == alphabet.size() || 0 == root_alphabet.size()))  {
        cerr << "ERROR: Alphabet for the vertex labels not specified. It is required for selected kernel type." << endl;  print_help();  exit(1);
    }

    GraphKernel gk;

    string line;
    vector<unsigned> examples;
    vector<int> labels;

    // Read list of positive examples.
    ifstream p(pos_file.c_str(), ios::in);
    if (p.fail()) {
        cerr << "WARNING: Positives file " << pos_file << " cannot be opened." << endl;
    }
    else  {
        while(getline(p, line))  {
            vector<string> tokens = split(line, '\t');            
            examples.push_back(to_i(strip(tokens[0])));
            labels.push_back(1);
        }
    }
    p.close();

	// Read list of negative examples.
    ifstream n(neg_file.c_str(), ios::in);
    if (n.fail())  {
        cerr << "WARNING: Negatives file " << neg_file << " cannot be opened." << endl;
    }
    else  {
        while(getline(n, line))  {
            vector<string> tokens = split(line, '\t');        
            examples.push_back(to_i(strip(tokens[0])));
            labels.push_back(-1);
        }
    }
    n.close();

    if (examples.size() < 1)  {
        cerr << "ERROR: Too few examples." << endl << endl; print_help(); exit(1);
    }

    if (normalize)  gk.set_normalize();
    if (verbose)  gk.set_verbose();

    gk.read_graphs(l_file, g_file, examples);
    gk.set_labels(labels);

    switch (format)  {
        case KERNEL:
			switch (k_type)  {
				case RANDOM_WALK_CUMULATIVE:
					gk.compute_random_walk_cumulative_matrix(steps, restart);
					break;
				case RANDOM_WALK:
					gk.compute_random_walk_matrix(steps, restart);
					break;
				case STANDARD_GRAPHLET:
					gk.set_number_label_mismatches(0.0);
					gk.compute_label_mismatch_matrix();
					break;
				case LABEL_MISMATCH:
					gk.set_number_label_mismatches(mismatches);
					gk.set_label_mismatches_alphabet(alphabet);
					gk.set_label_mismatches_root_alphabet(root_alphabet);
					gk.read_sim_matrix(sim_matrix_file);
					gk.compute_label_mismatch_matrix();
					break;
				case EDGE_MISMATCH:
					gk.set_number_edges_mismatches(edgmis);
                    gk.compute_edge_mismatch_matrix();
					break;
                case EDIT_DISTANCE:
                    gk.set_number_label_mismatches(mismatches);
                    gk.set_label_mismatches_alphabet(alphabet);
                    gk.set_label_mismatches_root_alphabet(root_alphabet);
                    gk.read_sim_matrix(sim_matrix_file);
                    gk.set_number_edges_mismatches(edgmis);
                    if (edgmis == 2)
                        gk.compute_edit_distance2_matrix();
                    else
                        gk.compute_edit_distance_matrix();
                    break;
			}
			gk.write_matrix(output_file.c_str());
            break;
        case SPARSE_SVML:
            switch (k_type)  {
				case RANDOM_WALK_CUMULATIVE:
                    gk.compute_random_walk_cumulative_matrix(steps, restart);
                    break;
                case RANDOM_WALK:
                    gk.compute_random_walk_matrix(steps, restart);
					break;
				case STANDARD_GRAPHLET:
					gk.set_number_label_mismatches(0.0);
					gk.write_sparse_svml_lm(output_file.c_str());
					break;
				case LABEL_MISMATCH:
					gk.set_number_label_mismatches(mismatches);
					gk.set_label_mismatches_alphabet(alphabet);
					gk.set_label_mismatches_root_alphabet(root_alphabet);
					gk.read_sim_matrix(sim_matrix_file);
					gk.write_sparse_svml_lm(output_file.c_str());
					break;
				case EDGE_MISMATCH:
                    gk.set_number_edges_mismatches(edgmis);
					gk.write_sparse_svml_em(output_file.c_str());
                    break;
                case EDIT_DISTANCE:
                    gk.set_number_label_mismatches(mismatches);
                    gk.set_label_mismatches_alphabet(alphabet);
                    gk.set_label_mismatches_root_alphabet(root_alphabet);
                    gk.read_sim_matrix(sim_matrix_file);
                    gk.set_number_edges_mismatches(edgmis);
                    if (edgmis == 2)
                        gk.write_sparse_svml_ed2(output_file.c_str()); 
                    else
                        gk.write_sparse_svml_ed(output_file.c_str());
            }
    }

    if (labels_file.size() > 0)
        gk.write_labels(labels_file.c_str());

    exit(0);
}

