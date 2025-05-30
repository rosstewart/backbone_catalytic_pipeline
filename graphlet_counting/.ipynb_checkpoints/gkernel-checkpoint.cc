#include "gkernel.h"
#include "string.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <iomanip>


/*********************** GraphKernel methods ***********************/
void GraphKernel::read_graphs(string nlabels_file, string graph_file, const vector<unsigned> &vertices)  {
    if (VERBOSE)  cerr << "Reading input data ... ";

    graph = SimpleGraph::read_graph(nlabels_file.c_str(), graph_file.c_str());

    for (unsigned i=0; i<vertices.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;
                
        roots.push_back(vertices[i]);
    }

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::read_sim_matrix(string sim_matrix_file)  {
    if (VERBOSE)  cerr << "Reading probability similarity matrix for vertex labels file ... ";
    string column, row, key;
	
    // Read vertex labels similarity matrix.
    ifstream p(sim_matrix_file.c_str(), ios::in);
    if (p.fail()) {
        cerr << "ERROR: Vertex labels similarity matrix file " << sim_matrix_file << " cannot be opened." << endl; exit(1);
    }
    else  {
        if (VERBOSE)  cerr << sim_matrix_file.c_str();
        getline(p, column);
        while(getline(p, row))  {
            vector<string> tokens = split(row, '\t'); 
            for(unsigned i=0 ; i<column.size() ; i++)  {
                key = tokens[0] + column[i];
                map<string, float>::iterator it = sim_vlm_matrix.find(key);
                if(it == sim_vlm_matrix.end())  {
                    sim_vlm_matrix[key] = to_f(tokens[i+1]);
                }
            }
        }
    }
    p.close();

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::compute_random_walk_cumulative_matrix(int steps, double restart)  {
    if (VERBOSE)  cerr << "Computing Cumulative Random Walk Graph Kernel for #steps = " << steps << " restart prob = " << restart << " ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);

    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;
        
        kernel[i][i] = random_walk_cumulative(graph, roots[i], roots[i], steps, restart);
        for (unsigned j=0; j<i; j++)  {			
            kernel[i][j] = random_walk_cumulative(graph, roots[i], roots[j], steps, restart);
        }        
    }
    if (VERBOSE)  cerr << endl;
}

void GraphKernel::compute_random_walk_matrix(int steps, double restart)  {
    if (VERBOSE)  cerr << "Computing Random Walk Graph Kernel for #steps = " << steps << " restart prob = " << restart << " ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);

    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;

        kernel[i][i] = random_walk(graph, roots[i], roots[i], steps, restart);
        for (unsigned j=0; j<i; j++)  {
            kernel[i][j] = random_walk(graph, roots[i], roots[j], steps, restart);
        }
    }
    if (VERBOSE)  cerr << endl;
}

void GraphKernel::compute_label_mismatch_matrix()  {
    if (VERBOSE)  {
        if (set_k((GRAPHLET_TYPES-1), SF) > 0)
            cerr << "Computing Label Substitutions Graphlet Kernel ... ";
        else
            cerr << "Computing Standard Graphlet Kernel ... ";
    }

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);
    
	for (unsigned i=0; i<roots.size(); i++)  {
		for (unsigned j=0; j<=i; j++)
			kernel[i][j] = 0.0;
	}

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(graph, roots[i]));
    }

	unsigned long g_type = 0;
	while (g_type < GRAPHLET_TYPES)  {
        if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;            
            
		    for (unsigned i=0; i<roots.size(); i++)  {
			    if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;
                
                add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, VLM, false);

                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, false, VLM, false);

			    kernel[i][i] = kernel[i][i] + distance_hash_join(hashes[i][g_type], hashes[i][g_type], g_type);
			    for (unsigned j=0; j<i; j++)  {
				    kernel[i][j] = kernel[i][j] + distance_hash_join(hashes[i][g_type], hashes[j][g_type], g_type);
			    }
		    }
			vl_mismatch_neighborhood.clear();

			for (unsigned i=0; i<roots.size(); i++)  {
				hashes[i][g_type].clear();
            }            
        }
        g_type = g_type + 1; 
	}

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::compute_edge_mismatch_matrix()  {
    if (VERBOSE)   cerr << "Computing Edge Indels Graphlet Kernel ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(graph, roots[i]));
        add_edge_mismatch_counts(hashes[i]);
    }

    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;
        
        kernel[i][i] = distance_hash_join(hashes[i], hashes[i]);
        for (unsigned j=0; j<i; j++)  {
            kernel[i][j] = distance_hash_join(hashes[i], hashes[j]);
        }
    }

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::compute_edit_distance_matrix()  {
    if (VERBOSE)  cerr << "Computing Edit Distance Graphlet Kernel (d=1) ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);
    
	for (unsigned i=0; i<roots.size(); i++)  {
		for (unsigned j=0; j<=i; j++)
			kernel[i][j] = 0.0;
	}

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(graph, roots[i]));
        add_edge_mismatch_counts(hashes[i]);
    }

	unsigned long g_type = 0;
	while (g_type < GRAPHLET_TYPES)  {
        if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;

		    for (unsigned i=0; i<roots.size(); i++)  {
			    if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;
                
                if (VLM >= 1)  {
                    add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, 1, false);
                }

                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, false, 1, true);

			    kernel[i][i] = kernel[i][i] + distance_hash_join(hashes[i][g_type], hashes[i][g_type], g_type);
			    for (unsigned j=0; j<i; j++)  {
				    kernel[i][j] = kernel[i][j] + distance_hash_join(hashes[i][g_type], hashes[j][g_type], g_type);
			    }
		    }
			vl_mismatch_neighborhood.clear();

            for (unsigned i=0; i<roots.size(); i++)  {
                hashes[i][g_type].clear();
            }
        }
        g_type = g_type + 1; 
	}

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::compute_edit_distance2_matrix()  {
    if (VERBOSE)  cerr << "Computing Edit Distance Graphlet Kernel (d=2) ... ";

    kernel.resize(roots.size());
    for (unsigned i=0; i<roots.size(); i++)
        kernel[i].resize(i+1);
    
	for (unsigned i=0; i<roots.size(); i++)  {
		for (unsigned j=0; j<=i; j++)
			kernel[i][j] = 0.0;
	}

    for (unsigned i=0; i<roots.size(); i++)  {
        hashes.push_back(get_graphlets_counts(graph, roots[i]));
        add_1_edge_mismatch_counts(hashes[i]);
    }

	unsigned long g_type = 0;
	while (g_type < GRAPHLET_TYPES)  {
        if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;

		    for (unsigned i=0; i<roots.size(); i++)  {
                if (VLM >= 1)  {
                    add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, 1, true);
                }
               
                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, true, 1, true);
            }
        }
		vl_mismatch_neighborhood.clear();
        g_type = g_type + 1;
    }
    for (unsigned i=0; i<roots.size(); i++)  {
        add_2_edge_mismatch_counts(hashes[i]);
    }

	g_type = 0;
	while (g_type < GRAPHLET_TYPES)  {
        if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
            int VLM = set_k(g_type, SF);
            map<Key,MismatchInfo> mismatch_hash;

            for (unsigned i=0; i<roots.size(); i++)  {
                if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;

                if (VLM >= 2)  {
                    add_vertex_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, 2, false);
                }

                update_label_mismatch_counts(hashes[i][g_type], mismatch_hash, g_type, false, 2, true);

                kernel[i][i] = kernel[i][i] + distance_hash_join(hashes[i][g_type], hashes[i][g_type], g_type);
			    for (unsigned j=0; j<i; j++)  {
				    kernel[i][j] = kernel[i][j] + distance_hash_join(hashes[i][g_type], hashes[j][g_type], g_type);
			    }
		    }
			vl_mismatch_neighborhood.clear();

            for (unsigned i=0; i<roots.size(); i++)  {
                hashes[i][g_type].clear();
            }
        }
        g_type = g_type + 1; 
	}
    if (VERBOSE)  cerr << endl;
}

#if OUTPUT_FORMAT == 0
void GraphKernel::write_matrix(const char *file)  {
    ofstream out(file, ios::out | ios::binary);

    unsigned g_size = roots.size();
    out.write((char*) &g_size, sizeof(unsigned)); 

    for (unsigned i=0; i<g_size; i++)
        for (unsigned j=0; j<=i; j++) 
            out.write((char*) &kernel[i][j], sizeof(float));

    out.close();
}
#elif OUTPUT_FORMAT == 1
void GraphKernel::write_matrix(const char *file)  {
    ofstream out(file, ios::out );

    unsigned g_size = roots.size();

    for (unsigned i=0; i<g_size; i++)  {
        for (unsigned j=0; j<=i; j++)   {
            out  << kernel[i][j] << "\t";
        }
        out  <<  endl;
    }
    out.close();
}
#else
void GraphKernel::write_matrix(const char *file)  {
	ofstream out(file, ios::out );
	
	unsigned g_size = roots.size();
	
	for (unsigned i=0; i<g_size; i++)  {
		for (unsigned j=0; j<g_size; j++)   {
		    if (j>i) {
                out  << setprecision (10) << kernel[j][i] <<  "\t"  ;
		    } else  {
                out  << setprecision (10) << kernel[i][j]  <<  "\t"  ;
            }
        }
        out  <<  endl;
    }
	out.close();
}
#endif

void GraphKernel::write_sparse_svml_lm(const char *file)  {
    if (VERBOSE)   { 
        if (set_k((GRAPHLET_TYPES-1), SF) > 0)
            cerr << "Computing attributes for Label Substitutions Graphlet Kernel ... ";
        else
            cerr << "Computing attributes for Standard Graphlet Kernel ... ";
    }

	ofstream out(file, ios::out);	
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(graph, roots[i]);
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
			if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
                int VLM= set_k(g_type, SF);
                map<Key,MismatchInfo> mismatch_hash;

                add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, VLM, false);

                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, false, VLM, false);                                
				for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)  {
                    out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first);
                }
			}
			vl_mismatch_neighborhood.clear();
        }
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::write_sparse_svml_em(const char *file)  {
    if (VERBOSE)  cerr << "Computing attributes for Edge Indels Graphlet Kernel ... ";

	ofstream out(file, ios::out);
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(graph, roots[i]);
        add_edge_mismatch_counts(g_hash);        
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
                for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)  {
                    out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_edge_mismatch_count(g_hash[g_type], it->first);
                }
            }
        }
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::write_sparse_svml_ed(const char *file)  {
    if (VERBOSE)  cerr << "Computing attributes for Edit Distance Graphlet Kernel (d=1) ... ";

    ofstream out(file, ios::out);
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(graph, roots[i]);
        add_edge_mismatch_counts(g_hash);
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {            
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
                int VLM = set_k(g_type, SF);
                map<Key,MismatchInfo> mismatch_hash;
                if(VLM >= 1)                    
                    add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, 1, false);
                    
                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, false, 1, true);
                
                for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)  {                    
                    if (retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first) > 0.0)  {
                        out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first);
                    }
                }
            }
			vl_mismatch_neighborhood.clear();
        }       
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::write_sparse_svml_ed2(const char *file)  {
    if (VERBOSE)  cerr << "Computing attributes for Edit Distance Graphlet Kernel (d=2) ... ";

    ofstream out(file, ios::out);
    for (unsigned i=0; i<roots.size(); i++)  {
        if (VERBOSE && (i+1) % 100 == 0)  cerr << " " << i+1;

        out << labels[i];

        vector<map<Key,MismatchInfo> > g_hash = get_graphlets_counts(graph, roots[i]);
        add_1_edge_mismatch_counts(g_hash);
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {            
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
                int VLM = set_k(g_type, SF);
                map<Key,MismatchInfo> mismatch_hash;

                if(VLM >= 1)                    
                    add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, 1, true);
                                    
                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, true, 1, true);
                vl_mismatch_neighborhood.clear();

                add_2_edge_mismatch_counts(g_hash);

                if (VLM >= 2)
                    add_vertex_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, 2, false);

                update_label_mismatch_counts(g_hash[g_type], mismatch_hash, g_type, false, 2, true);
                vl_mismatch_neighborhood.clear();

                for (map<Key,MismatchInfo>::iterator it = g_hash[g_type].begin(); it != g_hash[g_type].end(); it++)  {                    
                    if (retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first) > 0.0)  {
                        out << " " << get_feature_id(it->first, g_type) << ":" << retrieve_label_mismatch_count(g_type, g_hash[g_type], it->first);
                    }
                }
            }
        }
        out << " #" << i << endl;
    }
    out.close();

    if (VERBOSE)  cerr << endl;
}

void GraphKernel::write_labels(const char *file)  {
    unsigned num_pos(0), num_neg(0);

    ofstream out(file, ios::out);

    out << labels[0];
    labels[0]==1 ? num_pos++ : num_neg++;

    for (unsigned i=1; i<labels.size(); i++)  {
        out << "\n" << labels[i];
        labels[i]==1 ? num_pos++ : num_neg++;
    }
    out.close();

    if (VERBOSE)  {
        cerr << "Positive examples: " << num_pos << endl;
        cerr << "Negative examples: " << num_neg << endl;
    }
}

float GraphKernel::random_walk_cumulative(SimpleGraph &g, unsigned g1_root, unsigned g2_root, int STEPS, double RESTART)  {
    unsigned i1, i2, j1, j2;
    unsigned g1_next_node, g2_next_node; 
	float value = compare_labels(g.nodes[g1_root], g.nodes[g2_root]);

    if(g.adj[g1_root].size() == 0 || g.adj[g2_root].size() == 0)  {
        return value;
    }

    i1 = g1_root;
    i2 = g2_root;
    int step(1);
    while(step < STEPS)  {
        if(g.adj[i1].size() == 0 || g.adj[i2].size() == 0)  {
            return value;
        }
        j1 = randint(g.adj[i1].size());
        j2 = randint(g.adj[i2].size());
        g1_next_node = g.adj[i1][j1];
        g2_next_node = g.adj[i2][j2];
        char node1_label = g.nodes[g1_next_node];
        char node2_label = g.nodes[g2_next_node];
        value = value + compare_labels(node1_label, node2_label);
        double prob = randdouble();
        // Restart random walk
        if(prob < RESTART)  {
            i1 = g1_root;
            i2 = g2_root;
        }
        // Continue random walk
        else  {
            i1 = g1_next_node;
            i2 = g2_next_node;
        }
        step += 1;
    }

    return value;
}

float GraphKernel::random_walk(SimpleGraph &g, unsigned g1_root, unsigned g2_root, int STEPS, double RESTART)  {
    unsigned i1, i2, j1, j2;
    unsigned g1_next_node, g2_next_node; 
	float value(0.0);
    string g1_walk_seq, g2_walk_seq;

    if (g.nodes[g1_root] != g.nodes[g2_root] || g.adj[g1_root].size() == 0 || g.adj[g2_root].size() == 0)  {
        return value;
    }

    g1_walk_seq.push_back(g.nodes[g1_root]);
    g2_walk_seq.push_back(g.nodes[g2_root]);
    i1 = g1_root;
    i2 = g2_root;
    int step(1);
    while (step < STEPS)  {
        if (g.adj[i1].size() == 0 || g.adj[i2].size() == 0)  {
            return value;
        }
        j1 = randint(g.adj[i1].size());
        j2 = randint(g.adj[i2].size());
        g1_next_node = g.adj[i1][j1];
        g2_next_node = g.adj[i2][j2];
        g1_walk_seq.push_back(g.nodes[g1_next_node]);
        g2_walk_seq.push_back(g.nodes[g2_next_node]);
        double prob = randdouble();
        // Restart random walk
        if(prob < RESTART)  {
            i1 = g1_root;
            i2 = g2_root;
            if (g1_walk_seq.compare(g2_walk_seq) == 0) {
                value = value + 1.0;
            }
            g1_walk_seq.clear();
            g2_walk_seq.clear();
            g1_walk_seq.push_back(g.nodes[g1_root]);
            g2_walk_seq.push_back(g.nodes[g2_root]);
        }
        // Continue random walk
        else  {
            i1 = g1_next_node;
            i2 = g2_next_node;
        }
        step += 1;
    }

    return value;
}

// Count graphlets starting from root.
vector<map<Key,MismatchInfo> > GraphKernel::get_graphlets_counts(SimpleGraph &g, unsigned g_root)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > hash(GRAPHLET_TYPES, T);
    vector<unsigned> dist = g.breadth_first_sort(g_root);
	vector<Key> mismatches;
	Key key;

	unsigned i, j, k, l;
    char root, a, b, c, d;
    
    root = g.nodes[g_root];

	// 1-graphlets, case 0
    #if GRAPHLETS_1	
	key = create_permutations_subset(mismatches, root, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, 0);	
	increment_match_hash(hash[0], key, mismatches);

    #endif

    for (unsigned i_=0; i_<g.adj[g_root].size(); i_++)  {
        i = g.adj[g_root][i_];
        a = g.nodes[i];

        // 2-graphlets, case 01
		#if GRAPHLETS_2
		key = create_permutations_subset(mismatches, root, a, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, 1);	
		increment_match_hash(hash[1], key, mismatches);        
		#endif
		
        for (unsigned j_=0; j_<i_; j_++)  {
            // 3-graphlets, case 011
            j = g.adj[g_root][j_];
            b = g.nodes[j];

			#if GRAPHLETS_3
            bool found(false);
            unsigned t(0);
            while (!found && t<g.adj[i].size())  {
                found = (j == g.adj[i][t++]);
            }
            if (found)  {
				key = create_permutations_subset(mismatches, root, a, b, ZERO_CHAR, ZERO_CHAR, 4);
				increment_match_hash(hash[4], key, mismatches);                
            } 
            else  {
				key = create_permutations_subset(mismatches, root, a, b, ZERO_CHAR, ZERO_CHAR, 3);
				increment_match_hash(hash[3], key, mismatches);                
            }
			#endif

            // 4-graphlets, case 0111
            for (unsigned k_=0; k_<j_; k_++)  {
                k = g.adj[g_root][k_];
                c = g.nodes[k];
				
				#if GRAPHLETS_4
                bool found_ij(false), found_ik(false), found_jk(false);
                unsigned t(0);
                while(!found_ij && t<g.adj[i].size())
                    found_ij = (j == g.adj[i][t++]);

                t = 0;
                while (!found_ik && t<g.adj[i].size())
                    found_ik = (k == g.adj[i][t++]);

                t = 0;
                while (!found_jk && t<g.adj[j].size())
                    found_jk = (k == g.adj[j][t++]);

                if (found_ij)  {
                    if (found_ik)  {
                        if (found_jk)  {
							key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 15);
							increment_match_hash(hash[15], key, mismatches);                            
                        } else  {
							key = create_permutations_subset(mismatches, root, b, c, a, ZERO_CHAR, 14);
							increment_match_hash(hash[14], key, mismatches);
                        }
                    }
                    else  {
                        if (found_jk)  {
							key = create_permutations_subset(mismatches, root, a, c, b, ZERO_CHAR, 14);
							increment_match_hash(hash[14], key, mismatches);
                        } else  {
							key = create_permutations_subset(mismatches, root, c, a, b, ZERO_CHAR, 10);
							increment_match_hash(hash[10], key, mismatches);
                        }
                    }
                }
                else  {
                    if (found_ik)  {
                        if (found_jk)  {
							key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 14);
							increment_match_hash(hash[14], key, mismatches);
                        } else  {
							key = create_permutations_subset(mismatches, root, b, a, c, ZERO_CHAR, 10);
							increment_match_hash(hash[10], key, mismatches);                            
                        }
                    }
                    else  {
                        if (found_jk)  {
							key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 10);
							increment_match_hash(hash[10], key, mismatches);
                        } else  {
							key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 8);
							increment_match_hash(hash[8], key, mismatches);
                        }
                    }
                }
				#endif
				
				// 5-graphlets, case 01111
				for (unsigned l_=0; l_<k_; l_++)  {
					l = g.adj[g_root][l_];
					d = g.nodes[l];                    
					
					#if GRAPHLETS_5
                    bool found_ij(false), found_ik(false), found_il(false), found_jk(false), found_jl(false), found_kl(false);
					unsigned t(0);
					while(!found_ij && t<g.adj[i].size())  {
						found_ij = (j == g.adj[i][t++]);
                    }
					
					t = 0;
					while (!found_ik && t<g.adj[i].size())
						found_ik = (k == g.adj[i][t++]);
					
					t = 0;
					while(!found_il && t<g.adj[i].size())
						found_il = (l == g.adj[i][t++]);
					
					t = 0;
					while (!found_jk && t<g.adj[j].size())
						found_jk = (k == g.adj[j][t++]);
					
					t = 0;
					while (!found_jl && t<g.adj[j].size())
						found_jl = (l == g.adj[j][t++]);
					
					t = 0;
					while (!found_kl && t<g.adj[k].size())  {
						found_kl = (l == g.adj[k][t++]);
                    }

					if (found_ij)  {
						if (found_ik)  {
							if (found_il)  {
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 26);
											increment_match_hash(hash[26], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 25);
											increment_match_hash(hash[25], key, mismatches);
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, c, b, d, 25);
											increment_match_hash(hash[25], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, d, b, c, a, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
									}
								}
								else  {                                    
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, d, b, c, 25);
											increment_match_hash(hash[25], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, c, b, d, a, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, c, d, a, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, b, c, d, a, 22);
											increment_match_hash(hash[22], key, mismatches);
										}
									}
								}
							}
							else  {                                    
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, c, a, d, 25);
											increment_match_hash(hash[25], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, d, a, c, b, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, d, a, b, c, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, d, a, b, c, 21);
											increment_match_hash(hash[21], key, mismatches);
										}
									}
								}
								else  {                                    
									if (found_jl)  {
										if (found_kl)  {                                            
                                            key = create_permutations_subset(mismatches, root, a, d, b, c, 23);
											increment_match_hash(hash[23], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 20);
											increment_match_hash(hash[20], key, mismatches);
										}
									}
									else  {                                            
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, a, c, b, d, 20);
											increment_match_hash(hash[20], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, d, b, c, a, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
								}
							}
						}
						else  {                                
							if (found_il)  {
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, d, a, c, 25);
											increment_match_hash(hash[25], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, c, a, d, b, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
									}
									else  {
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, a, c, b, d, 23);
											increment_match_hash(hash[23], key, mismatches);											
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, a, b, d, c, 20);
											increment_match_hash(hash[20], key, mismatches);
										}
									}
								}
								else  {                                    
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, c, a, b, d, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, c, a, b, d, 21);
											increment_match_hash(hash[21], key, mismatches);
										}
									}
									else  {
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, a, d, b, c, 20);
											increment_match_hash(hash[20], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, c, b, d, a, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
								}
							}
							else  {                                    
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, c, d, b, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, a, c, d, b, 22);
											increment_match_hash(hash[22], key, mismatches);
										}
									}
									else  {                                            
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, b, c, a, d, 20);
											increment_match_hash(hash[20], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, d, a, c, b, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
								}
								else  {                                    
									if (found_jl)  {
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, b, d, a, c, 20);
											increment_match_hash(hash[20], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, c, a, d, b, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
									else  { 
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, a, b, c, d, 18);
											increment_match_hash(hash[18], key, mismatches);											
										}
										else  {
											key = create_permutations_subset(mismatches, root, c, d, a, b, 17);
											increment_match_hash(hash[17], key, mismatches);
										}
									}
								}
							}
						}
					}
					else  {
						if (found_ik)  {
							if (found_il)  {
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, c, d, a, b, 25);
											increment_match_hash(hash[25], key, mismatches);
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, a, b, c, d, 23);
											increment_match_hash(hash[23], key, mismatches);											
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, a, d, c, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, a, c, d, b, 20);
											increment_match_hash(hash[20], key, mismatches);                                            
										}
									}
								}
								else  {                                    
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, a, c, d, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, a, d, c, b, 20);
											increment_match_hash(hash[20], key, mismatches);                                            
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, a, c, d, 21);
											increment_match_hash(hash[21], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, b, c, d, a, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
								}
							}
							else  {                                    
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, d, c, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, b, c, d, a, 20);
											increment_match_hash(hash[20], key, mismatches);                                            
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, d, c, 22);
											increment_match_hash(hash[22], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, d, a, b, c, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
								}
								else  {
									if (found_jl)  {
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, c, d, a, b, 20);
											increment_match_hash(hash[20], key, mismatches);                                            
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, a, c, b, d, 18);
											increment_match_hash(hash[18], key, mismatches);											
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, a, d, c, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, b, d, a, c, 17);
											increment_match_hash(hash[17], key, mismatches);
										}
									}
								}
							}
						}
						else  {                                
							if (found_il)  {
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 24);
											increment_match_hash(hash[24], key, mismatches);
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, b, d, c, a, 20);
											increment_match_hash(hash[20], key, mismatches);                                            
										}
									}
									else  {
										if (found_kl)  {
                                            key = create_permutations_subset(mismatches, root, c, d, b, a, 20);
											increment_match_hash(hash[20], key, mismatches);                                            
										}
										else  {
                                            key = create_permutations_subset(mismatches, root, a, d, b, c, 18);
											increment_match_hash(hash[18], key, mismatches);											
										}
									}
								}
								else  {                                  
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 22);
											increment_match_hash(hash[22], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, c, a, b, d, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
									else  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, b, a, c, d, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, b, c, a, d, 17);
											increment_match_hash(hash[17], key, mismatches);
										}
									}
								}
							}
							else  {
								if (found_jk)  {
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 21);
											increment_match_hash(hash[21], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, a, c, d, b, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
									}
									else  {                                            
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, d, c, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, a, d, b, c, 17);
											increment_match_hash(hash[17], key, mismatches);
										}
									}
								}
								else  {                                    
									if (found_jl)  {
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 19);
											increment_match_hash(hash[19], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, a, c, b, d, 17);
											increment_match_hash(hash[17], key, mismatches);
										}
									}
									else  {             
										if (found_kl)  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 17);
											increment_match_hash(hash[17], key, mismatches);
										}
										else  {
											key = create_permutations_subset(mismatches, root, a, b, c, d, 16);
											increment_match_hash(hash[16], key, mismatches);
										}
									}
								}
							}
						}		
					}
					#endif
				}
            }
        }
		
		for (unsigned j_=0; j_<g.adj[g_root].size(); j_++)  {
			if (i_ == j_)  continue;
			
			j = g.adj[g_root][j_];
			b = g.nodes[j];
			
			for (unsigned k_=0; k_<g.adj[g_root].size(); k_++)  {
				if (k_ == i_ || k_ == j_)  continue;
				k = g.adj[g_root][k_];
				c = g.nodes[k];
				
                // 5-graphlets, case 01112
				for (unsigned l_=0; l_<g.adj[i].size(); l_++)  {
					l = g.adj[i][l_];
					if (2 != dist[l])  continue;
					d = g.nodes[l];
					
					#if GRAPHLETS_5
                    bool found_ij(false), found_ik(false), found_jk(false), found_jl(false), found_kl(false);
					unsigned t(0);
					while(!found_ij && t<g.adj[i].size())
						found_ij = (j == g.adj[i][t++]);
					
					t = 0;
					while (!found_ik && t<g.adj[i].size())
						found_ik = (k == g.adj[i][t++]);
					
					t = 0;
					while (!found_jk && t<g.adj[j].size())
						found_jk = (k == g.adj[j][t++]);
					
					t = 0;
					while (!found_jl && t<g.adj[j].size())
						found_jl = (l == g.adj[j][t++]);
					
					t = 0;
					while(!found_kl && t<g.adj[k].size())
						found_kl = (l == g.adj[k][t++]);

					if (found_ij)  {
						if (found_ik)  {
							if (found_jk)  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, a, b, c, d, 58);
											increment_match_hash(hash[58], key, mismatches);
									    }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 57);
											increment_match_hash(hash[57], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 57);
											increment_match_hash(hash[57], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 54);
											increment_match_hash(hash[54], key, mismatches);
									    }
                                    }
								}
							}
							else  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, a, b, c, d, 56);
											increment_match_hash(hash[56], key, mismatches);
									    }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 53);
											increment_match_hash(hash[53], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 53);
											increment_match_hash(hash[53], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 49);
											increment_match_hash(hash[49], key, mismatches);
									    }
                                    }
								}
							}
						}
						else  {
							if (found_jk)  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 56);
											increment_match_hash(hash[56], key, mismatches);
									    }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, b, a, d, 53);
											increment_match_hash(hash[53], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 55);
											increment_match_hash(hash[55], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, b, a, d, 47);
											increment_match_hash(hash[47], key, mismatches);
									    }
                                    }
								}
							}
							else  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 52);
											increment_match_hash(hash[52], key, mismatches);
									    }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 50);
											increment_match_hash(hash[50], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 48);
											increment_match_hash(hash[48], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, b, a, d, 45);
											increment_match_hash(hash[45], key, mismatches);
									    }
                                    }
								}
							}
						}
					}
					else  {
						if (found_ik)  {
							if (found_jk)  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 56);
											increment_match_hash(hash[56], key, mismatches);
									    }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 55);
											increment_match_hash(hash[55], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 53);
											increment_match_hash(hash[53], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 47);
											increment_match_hash(hash[47], key, mismatches);
									    }
                                    }
								}
							}
							else  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 52);
											increment_match_hash(hash[52], key, mismatches);
									    }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 48);
											increment_match_hash(hash[48], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 50);
											increment_match_hash(hash[50], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 45);
											increment_match_hash(hash[45], key, mismatches);
									    }
                                    }
								}
							}
						}
						else  {
							if (found_jk)  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, a, b, c, d, 52);
											increment_match_hash(hash[52], key, mismatches);
                                        }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, b, a, d, 48);
											increment_match_hash(hash[48], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 48);
											increment_match_hash(hash[48], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 44);
											increment_match_hash(hash[44], key, mismatches);
									    }
                                    }
								}
							}
							else  {
								if (found_jl)  {
									if (found_kl)  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, a, b, c, d, 51);
											increment_match_hash(hash[51], key, mismatches);
									    }
                                    }
									else  {
                                        if (i_ < j_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, c, a, b, d, 46);
											increment_match_hash(hash[46], key, mismatches);
									    }
                                    }
								}
								else  {
									if (found_kl)  {
                                        if (i_ < k_ && j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, a, c, d, 46);
											increment_match_hash(hash[46], key, mismatches);
									    }
                                    }
									else  {
                                        if (j_ < k_)  {
										    key = create_permutations_subset(mismatches, root, b, c, a, d, 43);
											increment_match_hash(hash[43], key, mismatches);
									    }
                                    }
								}
							}
						}
					}
					#endif
				}
			}
			
			// 4-graphlets, case 0112
			for (unsigned k_=0; k_<g.adj[i].size(); k_++)  {
				k = g.adj[i][k_];
				if (2 != dist[k])  continue;
				c = g.nodes[k];
				
				#if GRAPHLETS_4
                bool found_ij(false), found_jk(false);
				unsigned t(0);
				while(!found_ij && t<g.adj[i].size())
					found_ij = (j == g.adj[i][t++]);
				
				t = 0;
				while (!found_jk && t<g.adj[j].size())
					found_jk = (k == g.adj[j][t++]);

				if (found_ij)  {                    
					if (found_jk)  {						
						if(i_ < j_)  {
							key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 13);
							increment_match_hash(hash[13], key, mismatches);
						}						
                    } else  {
						key = create_permutations_subset(mismatches, root, b, a, c, ZERO_CHAR, 11);
						increment_match_hash(hash[11], key, mismatches);
					}
				} else  {
					if (found_jk)  {						
						if(i_ < j_)  {
							key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 12);
							increment_match_hash(hash[12], key, mismatches);
						}						
					} else  {
						key = create_permutations_subset(mismatches, root, b, a, c, ZERO_CHAR, 6);
						increment_match_hash(hash[6], key, mismatches);
					}
				}				
				#endif
				
				// 5-graphlets, case 01122, Type 1
				for (unsigned l_=0; l_<g.adj[i].size(); l_++)  {
					l = g.adj[i][l_];
					if (2 != dist[l] || l_ == k_)  continue;
					d = g.nodes[l];
					
					#if GRAPHLETS_5					
                    bool found_ij(false), found_jk(false), found_jl(false), found_kl(false);
					unsigned t(0);
					while(!found_ij && t<g.adj[i].size())
						found_ij = (j == g.adj[i][t++]);
					
					t = 0;
					while (!found_jk && t<g.adj[j].size())
						found_jk = (k == g.adj[j][t++]);
					
					t = 0;
					while (!found_jl && t<g.adj[j].size())
						found_jl = (l == g.adj[j][t++]);
					
					t = 0;
					while (!found_kl && t<g.adj[k].size())
						found_kl = (l == g.adj[k][t++]);

					if (!found_jl && !found_jk)  {
						if (found_ij)  {
							if (found_kl)  {
								if(k_ < l_)  {
									key = create_permutations_subset(mismatches, root, b, a, c, d, 42);
									increment_match_hash(hash[42], key, mismatches);
								}
							}
							else  {
								if(k_ < l_)  {
									key = create_permutations_subset(mismatches, root, b, a, c, d, 41);
									increment_match_hash(hash[41], key, mismatches);
								}
							}
						}
						else  {
							if (found_kl)  {
								if(k_ < l_)  {
									key = create_permutations_subset(mismatches, root, b, a, c, d, 40);
									increment_match_hash(hash[40], key, mismatches);
								}
							}
							else  {
								if(k_ < l_)  {
									key = create_permutations_subset(mismatches, root, b, a, c, d, 39);
									increment_match_hash(hash[39], key, mismatches);
								}
							}
						}
					}
					#endif
				}
				
				// 5-graphlets, case 01122, Type 2
				for (unsigned l_=0; l_<g.adj[j].size(); l_++)  {
					l = g.adj[j][l_];
					if (2 != dist[l] || l == k)  continue;
					d = g.nodes[l];
					
					#if GRAPHLETS_5
                    bool found_ij(false), found_il(false), found_jk(false), found_kl(false);
					unsigned t(0);
					while(!found_ij && t<g.adj[i].size())
						found_ij = (j == g.adj[i][t++]);
					
					t = 0;
					while(!found_il && t<g.adj[i].size())
						found_il = (l == g.adj[i][t++]);
					
					t = 0;
					while (!found_jk && t<g.adj[j].size())
						found_jk = (k == g.adj[j][t++]);
					
					t = 0;
					while (!found_kl && t<g.adj[k].size())
						found_kl = (l == g.adj[k][t++]);

					if (found_ij)  {
						if (found_il)  {
							if (found_jk)  {
								if (found_kl)  {
									if(i_ < j_ && k < l)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 38);
										increment_match_hash(hash[38], key, mismatches);
									}
								}
								else  {
									if(i_ < j_ && k < l)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 36);
										increment_match_hash(hash[36], key, mismatches);
									}
								}
							}
							else {
								if (found_kl)  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, b, a, c, d, 37);
										increment_match_hash(hash[37], key, mismatches);
									}
								}
								else  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, b, a, c, d, 32);
										increment_match_hash(hash[32], key, mismatches);
									}
								}	
							}
						} 
						else  {
							if (found_jk)  {
								if (found_kl)  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, d, c, 37);
										increment_match_hash(hash[37], key, mismatches);
									}
								}
								else  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, d, c, 32);
										increment_match_hash(hash[32], key, mismatches);
									}
								}
							}
							else  {
								if (found_kl)  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 31);
										increment_match_hash(hash[31], key, mismatches);                                        
									}
								}
								else  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 28);
										increment_match_hash(hash[28], key, mismatches);                                        
									}
								}
							}							
						}
					} 
					else  {
						if (found_il)  {
							if (found_jk)  {
								if (found_kl)  {
									if(i_ < j_ && k < l)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 35);
										increment_match_hash(hash[35], key, mismatches);
									}
								}
								else  {
									if(i_ < j_ && k < l)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 34);
										increment_match_hash(hash[34], key, mismatches);
									}
								}
							}
							else {
								if (found_kl)  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, b, a, c, d, 33);
										increment_match_hash(hash[33], key, mismatches);
									}
								}
								else  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, b, a, c, d, 30);
										increment_match_hash(hash[30], key, mismatches);
									}
								}
							}
						} 
						else  {
							if (found_jk)  {
								if (found_kl)  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, d, c, 33);
										increment_match_hash(hash[33], key, mismatches);
									}
								}
								else  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, d, c, 30);
										increment_match_hash(hash[30], key, mismatches);
									}
								}
							}
							else  {
								if (found_kl)  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 29);
										increment_match_hash(hash[29], key, mismatches);                                        
									}
								}
								else  {
									if(i_ < j_)  {
										key = create_permutations_subset(mismatches, root, a, b, c, d, 27);
										increment_match_hash(hash[27], key, mismatches);                                        
									}
								}
							}							
						}						
					}
					#endif
				}
				
				// 5-graphlets, case 01123 
				for (unsigned l_=0; l_<g.adj[k].size(); l_++)  {
					l = g.adj[k][l_];
					//Generalize to include all "local paths"
					if(dist[l] <= 1 || l == i || l == j || l == k)  continue;
					d = g.nodes[l];
					
					#if GRAPHLETS_5
                    bool found_il(false), found_jl(false), found_ij(false), found_jk(false);
					unsigned t(0);
					while (!found_il && t<g.adj[i].size())
						found_il = (l == g.adj[i][t++]);
					
					t = 0;
					while (!found_jl && t<g.adj[j].size())
						found_jl = (l == g.adj[j][t++]);
					
					t = 0;
					while(!found_ij && t<g.adj[i].size())
						found_ij = (j == g.adj[i][t++]);
					
					t = 0;
					while (!found_jk && t<g.adj[j].size())
						found_jk = (k == g.adj[j][t++]);

					if (!found_il && !found_jl)  {
						if (found_ij)  {
							if (found_jk)  {
								if (i_ < j_)  {
									key = create_permutations_subset(mismatches, root, a, b, c, d, 66);
									increment_match_hash(hash[66], key, mismatches);
								}
							}
							else  {
								key = create_permutations_subset(mismatches, root, b, a, c, d, 64);
								increment_match_hash(hash[64], key, mismatches);
							}
						}
						else  {
							if (found_jk)  {
								if (i_ < j_)  {
									key = create_permutations_subset(mismatches, root, a, b, c, d, 65);
									increment_match_hash(hash[65], key, mismatches);
								}
							}
							else  {
								key = create_permutations_subset(mismatches, root, b, a, c, d, 63);
								increment_match_hash(hash[63], key, mismatches);
							}
						}
					}
					#endif
				}					
			}
		}
		
		// 3-graphlets, case 012
		for (unsigned j_=0; j_<g.adj[i].size(); j_++)  {
			j = g.adj[i][j_];
			if (2 != dist[j])  continue;
			b = g.nodes[j];
			
			#if GRAPHLETS_3
			key = create_permutations_subset(mismatches, root, a, b, ZERO_CHAR, ZERO_CHAR, 2);
			increment_match_hash(hash[2], key, mismatches);
			#endif
			
			// 4-graphlets, case 0122
			for (unsigned k_=0; k_<j_; k_++)  {
				k = g.adj[i][k_];
				if (2 != dist[k])  continue;
				c = g.nodes[k];
				
				#if GRAPHLETS_4
				bool found(false);
				unsigned t(0);
				while (!found && t<g.adj[j].size())
					found = (k == g.adj[j][t++]);
				
				if (found)  {
					key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 9);
					increment_match_hash(hash[9], key, mismatches);
				} else  {
					key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 7);
					increment_match_hash(hash[7], key, mismatches);
				}
				#endif
				
				// 5-graphlets, case 01222
				for (unsigned l_=0; l_<k_; l_++)  {
					l = g.adj[i][l_];
					if (2 != dist[l] || l == i)  continue;
					d = g.nodes[l];
					
					#if GRAPHLETS_5
					bool found_jk(false), found_jl(false), found_kl(false);
					unsigned t(0);
					while (!found_jk && t<g.adj[j].size())
						found_jk = (k == g.adj[j][t++]);
					
					t = 0;
					while (!found_jl && t<g.adj[j].size())
						found_jl = (l == g.adj[j][t++]);
					
					t = 0;
					while (!found_kl && t<g.adj[k].size())
						found_kl = (l == g.adj[k][t++]);
					if (found_jk)  {
						if (found_jl)  {
							if (found_kl)  {
								key = create_permutations_subset(mismatches, root, a, b, c, d, 62);
								increment_match_hash(hash[62], key, mismatches);
							}
							else  {
								key = create_permutations_subset(mismatches, root, a, b, c, d, 61);
								increment_match_hash(hash[61], key, mismatches);
							}
						}
						else  {
							if (found_kl)  {
								key = create_permutations_subset(mismatches, root, a, c, b, d, 61);
								increment_match_hash(hash[61], key, mismatches);
							}
							else  {
								key = create_permutations_subset(mismatches, root, a, d, b, c, 60);
								increment_match_hash(hash[60], key, mismatches);
							}							
						}
					}
					else  {
						if (found_jl)  {
							if (found_kl)  {
								key = create_permutations_subset(mismatches, root, a, d, b, c, 61);
								increment_match_hash(hash[61], key, mismatches);
							}
							else  {
								key = create_permutations_subset(mismatches, root, a, c, b, d, 60);
								increment_match_hash(hash[60], key, mismatches);
							}
						}
						else  {
							if (found_kl)  {
								key = create_permutations_subset(mismatches, root, a, b, c, d, 60);
								increment_match_hash(hash[60], key, mismatches);
							}
							else  {
								key = create_permutations_subset(mismatches, root, a, b, c, d, 59);
								increment_match_hash(hash[59], key, mismatches);
							}							
						}						
					}
					#endif	
				}	
			}
			
			// 5-graphlets, case 01223
            for (unsigned k_=0; k_<g.adj[i].size(); k_++)  {
				k = g.adj[i][k_];
				if (2 != dist[k] || j_ == k_)  continue;
				c = g.nodes[k];
				
				for (unsigned l_=0; l_<g.adj[j].size(); l_++)  {
					l = g.adj[j][l_];					
                    // Generalize to include all "local paths"
					if(dist[l] <= 1 || l == i || l == j || l == k)  continue;                    
					d = g.nodes[l];
					
					#if GRAPHLETS_5
					bool found_il(false), found_jk(false), found_kl(false);
					unsigned t(0);
					while (!found_il && t<g.adj[i].size())
						found_il = (l == g.adj[i][t++]);
					
					t = 0;
					while(!found_jk && t<g.adj[j].size())
						found_jk = (k == g.adj[j][t++]);
					
					t = 0;
					while (!found_kl && t<g.adj[k].size())
						found_kl = (l == g.adj[k][t++]);

                    if (!found_il)  {
					    if (found_jk)  {
						    if (found_kl)  {
                                if (j_ < k_)  {
							        key = create_permutations_subset(mismatches, root, a, b, c, d, 70);
									increment_match_hash(hash[70], key, mismatches);
						        }
                            }
						    else  {
							    key = create_permutations_subset(mismatches, root, a, c, b, d, 69);
								increment_match_hash(hash[69], key, mismatches);
                            }
					    }
					    else  {
						    if (found_kl)  {
                                if (j_ < k_)  {
							        key = create_permutations_subset(mismatches, root, a, b, c, d, 68);
									increment_match_hash(hash[68], key, mismatches);
						        }
                            }
						    else  {
							    key = create_permutations_subset(mismatches, root, a, c, b, d, 67);
								increment_match_hash(hash[67], key, mismatches);
                            }
					    }
                    }
					#endif
				}
            }
			
            // 4-graphlets, case 0123
            for (unsigned k_=0; k_<g.adj[j].size(); k_++)  {
				k = g.adj[j][k_];
                // Expanded Vacic et al orbit 5 case to include all "local" graphlets whose path satisfies the 0123 property.
				if(dist[k] <= 1 || k == i || k == j)  continue;
				c = g.nodes[k];
				
				#if GRAPHLETS_4
				bool found_ik(false);
				unsigned t(0);
				while (!found_ik && t<g.adj[i].size())  {
					found_ik = (k == g.adj[i][t++]);
				}
				
                // Verify that local path does not classifies under previous graphlet types.
				if (!found_ik)  {
                    key = create_permutations_subset(mismatches, root, a, b, c, ZERO_CHAR, 5);
					increment_match_hash(hash[5], key, mismatches);
                }
				#endif
				
				// 5-graphlets, case 01233 
				for (unsigned l_=0; l_<k_; l_++)  {
					l = g.adj[j][l_];
                    // Expanded to include all "local" graphlets whose path satisfies the 0123 property.
                    if(dist[l] <= 1 || l == i || l == j)  continue;										
					d = g.nodes[l];
		
					#if GRAPHLETS_5
					bool found_il(false), found_ik(false), found_kl(false);
					unsigned t(0);
					while (!found_il && t<g.adj[i].size())
						found_il = (l == g.adj[i][t++]);
					
                    t = 0;
                    while (!found_ik && t<g.adj[i].size())
                        found_ik = (k == g.adj[i][t++]);
					
					t = 0;
					while (!found_kl && t<g.adj[k].size())
						found_kl = (l == g.adj[k][t++]);

                    // Verify that local path does not classifies under previous graphlet types.
                    if(!found_il && !found_ik)  {
					    if(found_kl)  {
						    key = create_permutations_subset(mismatches, root, a, b, c, d, 72);
							increment_match_hash(hash[72], key, mismatches);
					    }
					    else  {
						    key = create_permutations_subset(mismatches, root, a, b, c, d, 71);
							increment_match_hash(hash[71], key, mismatches);
					    }
                    }
					#endif
				}
				
				// 5-graphlets, case 01234 
				for (unsigned l_=0; l_<g.adj[k].size(); l_++)  {
					l = g.adj[k][l_];
                    // Orbit 73 case includes all "local" graphlets whose path satisfies the 01234 property.
					if(dist[l] <= 1 || l == i || l == j || l == k)  continue;					
					d = g.nodes[l];
					
					#if GRAPHLETS_5
					bool found_ik(false), found_il(false), found_jl(false);
					unsigned t(0);
					while (!found_ik && t<g.adj[i].size())
						found_ik = (k == g.adj[i][t++]);

                    t = 0;
					while (!found_il && t<g.adj[i].size())
						found_il = (l == g.adj[i][t++]);
					
					t = 0;
					while (!found_jl && t<g.adj[j].size())
						found_jl = (l == g.adj[j][t++]);
					
                    // Verify that local path does not classifies under previous graphlet types.
                    if(!found_ik && !found_il && !found_jl)  {
					    key = create_permutations_subset(mismatches, root, a, b, c, d, 73);
						increment_match_hash(hash[73], key, mismatches);
                    }
					#endif
				}
			}
		}
    } 
    
    if (NORMALIZE)
        normalize_spectral(hash);

    return hash;
}

// Add inexact graphlets by allowing vertex and edge label mismatches upto VLM.
void GraphKernel::add_vertex_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long g_type, int VLM, bool option)  {
	// Update counts to include vertex label mismacthes.
	if(VLM > 0 && ((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5)))  {
        // For each exact graphlet, generate all mismatch graphlets upto vertex label distance VLM.
        for (map<Key,MismatchInfo>::iterator git = hash.begin(); git != hash.end(); git++)  {
            if (hash[git->first].matches > 0)
                generate_vertex_label_mismatch_graphlets(vl_mismatch_neighborhood, hash, mismatch_hash, git->first, g_type, ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, VLM);
            else  {
                if (option)
                    generate_vertex_label_mismatch_graphlets(vl_mismatch_neighborhood, hash, mismatch_hash, git->first, g_type, ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, VLM);
            }
        }
    }
}

void GraphKernel::update_label_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, unsigned long g_type, bool option, int VLM, bool eq)  { 
    if((VLM > 0) && ((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5)))  {
        // Update mismatch counts for newly created graphlets.
        for (map<Key,MismatchInfo>::iterator git = hash.begin(); git != hash.end(); git++)  {
            if (option)  {
                update_mismatch_count(hash, git->first, (hash[git->first].matches + hash[git->first].mismatches), g_type, sim_vlm_matrix, VLM, eq);
                update_mismatch_count(mismatch_hash, git->first, (hash[git->first].matches + hash[git->first].mismatches), g_type, sim_vlm_matrix, VLM, eq);
            }
            else  {
                update_mismatch_count(hash, git->first, hash[git->first].matches, g_type, sim_vlm_matrix, VLM, eq);
                update_mismatch_count(mismatch_hash, git->first, hash[git->first].matches, g_type, sim_vlm_matrix, VLM, eq);
            }
        }
        // Merge all non-zero graphlets into one hashable list. 
        for (map<Key,MismatchInfo>::iterator git = mismatch_hash.begin(); git != mismatch_hash.end(); git++)  {
            insert_mismatch_counts(hash, mismatch_hash, git->first);
        }
        mismatch_hash.clear();
    }
}

// Add inexact graphlets by allowing edge insertions and deletions upto EM.
void GraphKernel::add_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > mismatch_hash(GRAPHLET_TYPES, T);
    unsigned vindex = 0;
    list<pair <unsigned long, Key> > L;
    float mult_factor;

    if (EM > 0)  {
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  { 
                for (map<Key,MismatchInfo>::iterator it = hash[g_type].begin(); it != hash[g_type].end(); it++)  { 
                    vindex = 0;
                    vector<list<pair <unsigned long, Key> > > EM_set(EM, L);
                    pair <unsigned long, Key> p (g_type, it->first); 
                    mult_factor = retrieve_exact_matches_count(hash[g_type], it->first);

                    // Generate edge mismacth neighborhood for this graphlet.
                    update_edge_mismatch_count(EM_set, it->first, g_type, EM, vindex);

                    list<pair <unsigned long, Key> > :: iterator list_it, search_it;
                    bool already_added(false);

                    // Add graphlets within edge mismacth neighborhood of current graphlet.
                    for (unsigned em = 0; em < EM; em++)  { 
                        EM_set[em].sort();
                        for (list_it = EM_set[em].begin(); list_it != EM_set[em].end(); list_it++)  {
                            already_added = false;
                            for (int previous = em-1; previous >= 0 ; previous--)  {
                                for (search_it = EM_set[previous].begin(); search_it != EM_set[previous].end(); search_it++)  {
                                    if ((list_it->first == p.first && list_it->second == p.second) || (list_it->first == search_it->first && list_it->second == search_it->second))  {
                                        already_added = true;
                                        break;
                                    }
                                }
                            }
                            if (!already_added)  { 
                                vector<Key> mismatches;
                                char root, a, b, c, d;
			                    initialize_vertices_labels(list_it->second, root, a, b, c, d);
			                    Key k = create_permutations_subset(mismatches, root, a, b, c, d, list_it->first);
                                if (k != list_it->second)  {
                                    cerr << "ERROR: Graphlets keys do not match " << print_key(k) << " vs " << print_key(list_it->second) << " ; " << list_it->first << endl; exit(1);
                                }
                                increment_edge_mismatch_hash(hash[list_it->first], mismatch_hash[list_it->first], list_it->second, mult_factor, mismatches);
                            }
                        }
                    }
                }
            }
        }
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
                for (map<Key,MismatchInfo>::iterator mit = mismatch_hash[g_type].begin(); mit != mismatch_hash[g_type].end(); mit++)  {
					insert_mismatch_counts(hash[g_type], mismatch_hash[g_type], mit->first);
                }
            }
        }
        mismatch_hash.clear();
    }
}

// Add inexact graphlets by allowing 1-edge insertion and deletion.
void GraphKernel::add_1_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > mismatch_hash(GRAPHLET_TYPES, T);
    unsigned vindex = 0;
    list<pair <unsigned long, Key> > L;
    float mult_factor;

    if (EM > 0)  {
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  { 
                for (map<Key,MismatchInfo>::iterator it = hash[g_type].begin(); it != hash[g_type].end(); it++)  { 
                    vindex = 0;
                    vector<list<pair <unsigned long, Key> > > EM_set(EM, L);
                    pair <unsigned long, Key> p (g_type, it->first); 
                    mult_factor = retrieve_exact_matches_count(hash[g_type], it->first);

                    // Generate edge mismacth neighborhood for this graphlet.
                    update_edge_mismatch_count(EM_set, it->first, g_type, 1, vindex);

                    list<pair <unsigned long, Key> > :: iterator list_it;

                    // Add graphlets within edge mismacth neighborhood of current graphlet.
                    EM_set[0].sort();
					for (list_it = EM_set[0].begin(); list_it != EM_set[0].end(); list_it++)  {
                        vector<Key> mismatches;
                        char root, a, b, c, d;
			            initialize_vertices_labels(list_it->second, root, a, b, c, d);
			            Key k = create_permutations_subset(mismatches, root, a, b, c, d, list_it->first);
                        if (k != list_it->second)  {
                            cerr << "ERROR: Graphlets keys do not match " << print_key(k) << " vs " << print_key(list_it->second) << endl; exit(1);
                        }
                        increment_edge_mismatch_hash(hash[list_it->first], mismatch_hash[list_it->first], list_it->second, mult_factor, mismatches);
                    }
                }
            }
        }
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
                for (map<Key,MismatchInfo>::iterator mit = mismatch_hash[g_type].begin(); mit != mismatch_hash[g_type].end(); mit++)  {                    
					insert_mismatch_counts(hash[g_type], mismatch_hash[g_type], mit->first);
                }
            }
        }
        mismatch_hash.clear();
    }
}

// Add inexact graphlets by allowing 2-edge insertions and deletions.
void GraphKernel::add_2_edge_mismatch_counts(vector<map<Key,MismatchInfo> > &hash)  {
    map<Key,MismatchInfo> T;
    vector<map<Key,MismatchInfo> > mismatch_hash(GRAPHLET_TYPES, T);
    unsigned vindex = 0;
    list<pair <unsigned long, Key> > L;
    float mult_factor;

    if (EM > 0)  {
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  { 
                for (map<Key,MismatchInfo>::iterator it = hash[g_type].begin(); it != hash[g_type].end(); it++)  { 
                    vindex = 0;
                    vector<list<pair <unsigned long, Key> > > EM_set(EM, L);
                    pair <unsigned long, Key> p (g_type, it->first); 
                    mult_factor = retrieve_exact_matches_count(hash[g_type], it->first);

                    // Generate edge mismacth neighborhood for this graphlet.
                    update_edge_mismatch_count(EM_set, it->first, g_type, 2, vindex);

                    list<pair <unsigned long, Key> > :: iterator list_it, search_it;
                    bool already_added(false);

                    // Add graphlets within edge mismacth neighborhood of current graphlet.
                    for (unsigned em = 1; em < EM; em++)  { 
                        EM_set[em].sort();
                        for (list_it = EM_set[em].begin(); list_it != EM_set[em].end(); list_it++)  {
                            already_added = false;
                            for (int previous = em-1; previous >= 0 ; previous--)  {
                                for (search_it = EM_set[previous].begin(); search_it != EM_set[previous].end(); search_it++)  {
                                    if ((list_it->first == p.first && list_it->second == p.second) || (list_it->first == search_it->first && list_it->second == search_it->second))  {
                                        already_added = true;
                                        break;
                                    }
                                }
                            }
                            if (!already_added)  {
                                vector<Key> mismatches;
                                char root, a, b, c, d;
			                    initialize_vertices_labels(list_it->second, root, a, b, c, d);
			                    Key k = create_permutations_subset(mismatches, root, a, b, c, d, list_it->first);
                                if (k != list_it->second)  {
                                    cerr << "ERROR: Graphlets keys do not match " << print_key(k) << " vs " << print_key(list_it->second) << endl; exit(1);
                                }
                                increment_edge_mismatch_hash(hash[list_it->first], mismatch_hash[list_it->first], list_it->second, mult_factor, mismatches);
                            }
                        }
                    }
                }
            }
        }
        for (unsigned g_type=0; g_type<GRAPHLET_TYPES; g_type++)  {
            if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
                for (map<Key,MismatchInfo>::iterator mit = mismatch_hash[g_type].begin(); mit != mismatch_hash[g_type].end(); mit++)  {
					insert_mismatch_counts(hash[g_type], mismatch_hash[g_type], mit->first);
                }
            }
        }
        mismatch_hash.clear();
    }
}

void GraphKernel::normalize_spectral(map<Key,MismatchInfo> &hash, unsigned long g_type)  {
    float norm = sqrt(distance_hash_join(hash, hash, g_type));    
    if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
		for (map<Key,MismatchInfo>::iterator it = hash.begin(); it != hash.end(); it++)  {
			it->second.matches = retrieve_label_mismatch_count(g_type, hash, it->first)/norm;
		}
	}
}

float GraphKernel::distance_hash_join(map<Key,MismatchInfo> g_hash, map<Key,MismatchInfo> h_hash, unsigned long g_type)  {
    float sum(0);
	
	if((g_type == 0 && (GRAPHLETS_1)) || (g_type == 1 && GRAPHLETS_2) || ((g_type >= 2 && g_type <= 4) && GRAPHLETS_3) || ((g_type >= 5 && g_type <= 15) && GRAPHLETS_4) || ((g_type >= 16 && g_type <= 73) && GRAPHLETS_5))  {
		for (map<Key,MismatchInfo>::iterator git = g_hash.begin(); git != g_hash.end(); git++)  {
			map<Key,MismatchInfo>::iterator hit = h_hash.find(git->first);
			if (hit != h_hash.end())  {
				if (NORMALIZE)  {
					sum += git->second.matches * hit->second.matches;
                }
				else  {                
					sum += retrieve_label_mismatch_count(g_type, g_hash, git->first) * retrieve_label_mismatch_count(g_type, h_hash, hit->first);
                }
			}
		}
	}
    return sum;
}

void GraphKernel::normalize_spectral(vector<map<Key,MismatchInfo> > &hash)  {
    float norm = sqrt(distance_hash_join(hash, hash));
    
    for (unsigned i=0; i<GRAPHLET_TYPES; i++)  {
        for (map<Key,MismatchInfo>::iterator it = hash[i].begin(); it != hash[i].end(); it++)
            it->second.matches = retrieve_edge_mismatch_count(hash[i], it->first)/norm;
    }
}

float GraphKernel::distance_hash_join(vector<map<Key,MismatchInfo> > g_hash, vector<map<Key,MismatchInfo> > h_hash)  {
    float sum(0);
    
    for (unsigned i=0; i<GRAPHLET_TYPES; i++)  {
        if((i == 0 && (GRAPHLETS_1)) || (i == 1 && GRAPHLETS_2) || ((i >= 2 && i <= 4) && GRAPHLETS_3) || ((i >= 5 && i <= 15) && GRAPHLETS_4) || ((i >= 16 && i <= 73) && GRAPHLETS_5))  {
            for (map<Key,MismatchInfo>::iterator git = g_hash[i].begin(); git != g_hash[i].end(); git++)  {
                map<Key,MismatchInfo>::iterator hit = h_hash[i].find(git->first);
            
                if (hit != h_hash[i].end())  {
                    if (NORMALIZE)
                        sum += git->second.matches * hit->second.matches;
                    else
                        sum += retrieve_edge_mismatch_count(g_hash[i], git->first) * retrieve_edge_mismatch_count(h_hash[i], hit->first);
                }
            }
        }
    }
    return sum;
}
