#include "mismatches.h"
#include "string.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <queue>
#include <iomanip>


/************************ Auxiliary functions for Graph Kernel class *************************/
string get_key(Key k)  {	
	string temp(GRAPHLET_SIZE,ZERO_CHAR);
    for (unsigned i=GRAPHLET_SIZE-1; i > 0; i--)  {
        temp[i] = (k & ALPHABET_SIZE) + ZERO_CHAR; k = k >> LOG_ALPHABET_SIZE;
    }
    temp[0] = (k & ALPHABET_SIZE) + ZERO_CHAR; //root  
	
    return temp;
}

string print_key(Key k)  {
	string temp(GRAPHLET_SIZE,ZERO_CHAR);
    for (unsigned i=GRAPHLET_SIZE-1; i > 0; i--)  {
        temp[i] = (k & ALPHABET_SIZE) + ZERO_CHAR; k = k >> LOG_ALPHABET_SIZE;
    }
    temp[0] = (k & ALPHABET_SIZE) + ZERO_CHAR; //root  
	
    ostringstream s;
    s << temp ;
    return s.str();
}

Key make_key(char root, char a, char b, char c, char d, unsigned long g_type)  {
    Key curr_key(0);

	curr_key = ((unsigned long) (root-ZERO_CHAR)) << LOG_ALPHABET_SIZE;
    curr_key = (curr_key+a-ZERO_CHAR) << LOG_ALPHABET_SIZE;
	curr_key = (curr_key+b-ZERO_CHAR) << LOG_ALPHABET_SIZE;
    curr_key = (curr_key+c-ZERO_CHAR) << LOG_ALPHABET_SIZE;
    curr_key = (curr_key+d-ZERO_CHAR);

    return curr_key;
}

void initialize_vertices_labels(Key key, char &root, char &a, char &b, char &c, char &d)  {
	d = (key & ALPHABET_SIZE) + ZERO_CHAR;  key = key >> LOG_ALPHABET_SIZE;
    c = (key & ALPHABET_SIZE) + ZERO_CHAR;  key = key >> LOG_ALPHABET_SIZE;
    b = (key & ALPHABET_SIZE) + ZERO_CHAR;  key = key >> LOG_ALPHABET_SIZE;
    a = (key & ALPHABET_SIZE) + ZERO_CHAR;  key = key >> LOG_ALPHABET_SIZE;
    root = (key & ALPHABET_SIZE) + ZERO_CHAR;
}

Key get_feature_id(Key k, unsigned long g_type)  {
    Key feature_id(0);
    feature_id = (k << LOG_GRAPHLET_TYPES_SIZE) + g_type;
    return feature_id;
}

void increment_match_hash(map<Key,MismatchInfo> &hash, const Key &k, vector<Key> &mismatches_list)  {
    map<Key,MismatchInfo>::iterator it; 
    if ((it = hash.find(k)) == hash.end())  {
        hash[k].matches = 1.0;
        hash[k].mismatches = 0.0;
        for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
            hash[k].mismatchesGraph[*mismatch_key] = 0.0;
        }
    }
    else  {
        hash[k].matches = it->second.matches + 1.0;
    }
	mismatches_list.clear();
}

float retrieve_exact_matches_count(map<Key,MismatchInfo> &hash, const Key &k)  {
    float perfect_matches = 0.0;
    map<Key,MismatchInfo>::iterator it = hash.find(k);
    if (it != hash.end())  {
        perfect_matches = hash[k].matches ;
    }
    return perfect_matches;
}

float retrieve_edge_mismatch_count(map<Key,MismatchInfo> &hash, const Key &k)  {
    float counts = 0.0;
    map<Key,MismatchInfo>::iterator it = hash.find(k);
    if (it != hash.end())  {
        counts = hash[k].matches + hash[k].mismatches ;
    }
    return counts;
}

float retrieve_label_mismatch_count(unsigned long g_type, map<Key,MismatchInfo> &hash, const Key key)  {
    map<Key,MismatchInfo>::iterator it = hash.find(key);
    float counts = 0.0;
    if (it != hash.end())  {
	    counts = hash[key].matches + hash[key].mismatches;
	    for (map<Key,float>::iterator mismatches = hash[key].mismatchesGraph.begin(); mismatches != hash[key].mismatchesGraph.end(); mismatches++)  {
		    counts += mismatches->second;
	    }
    }
	return counts;
}

void insert_mismatch_counts(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k)  {
    map<Key,MismatchInfo>::iterator it;    
    if ((it = hash.find(k)) == hash.end())  {
        hash[k].matches = mismatch_hash[k].matches;
        hash[k].mismatches = mismatch_hash[k].mismatches;
        for (map<Key,float>::iterator mismatch_key = mismatch_hash[k].mismatchesGraph.begin(); mismatch_key != mismatch_hash[k].mismatchesGraph.end(); mismatch_key++)  {
            hash[k].mismatchesGraph[mismatch_key->first] = mismatch_hash[k].mismatchesGraph[mismatch_key->first];
        }
    }
    else  {
        cerr << "Warning: While merging found a double counted graphlet " << print_key(k) << endl;
    }
}

void insert_mismatches_hash(map<Key,MismatchInfo> &mismatch_hash, const Key &k, vector<Key> mismatches_list)  {
    map<Key,MismatchInfo>::iterator it;
	mismatch_hash[k].matches = 0.0;
    mismatch_hash[k].mismatches = 0.0;
	for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
		mismatch_hash[k].mismatchesGraph[*mismatch_key] = 0.0;
    }
}

// Create all vertex-labeled graphlet permutations for each graphlet found.
Key create_permutations_subset(vector<Key> &mismatches, char root, char a, char b, char c, char d, unsigned long g_type)  {
	Key k(0), curr_key(0); //For label mismatches labels.
	char a1, b1, c1, d1; 

	mismatches.clear();

	switch (g_type)  {
		case 0:
			// Make key label
			// P (Pivot) 
			k = make_key(root, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);			
			break;
			
		case 1:
			// Make key label
			// P-A 
			k = make_key(root, a, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			break;
			
		case 3:
		case 4:
			// Make key label
			// P-A1-A2
			compare_two(a, b, a1, b1);
			k = make_key(root, a1, b1, ZERO_CHAR, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A2-A1
			curr_key = make_key(root, b1, a1, ZERO_CHAR, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 2:
			// Make key label
			// P-A-B                 
			k = make_key(root, a, b, ZERO_CHAR, ZERO_CHAR, g_type);			
			//Insert into permutations subset
			mismatches.push_back(k);
			break;
			
		case 5:
		case 6:
		case 11:
			// Make key label
			// P-A-B-C 
			k = make_key(root, a, b, c, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			break;
			
		case 12:
		case 13:
		case 14:
			// Make key label
			// P-A1-A2-B 
			compare_two(a, b, a1, b1);
			k = make_key(root, a1, b1, c, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A2-A1-B 
			curr_key = make_key(root, b1, a1, c, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 7:
		case 9:
		case 10:
			// Make key label
			// P-A-B1-B2 
			compare_two(b, c, b1, c1);
			k = make_key(root, a, b1, c1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A-B2-B1 
			curr_key = make_key(root, a, c1, b1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 8:
		case 15:
			// Make key label
			// P-A1-A2-A3 
			compare_three(a, b, c, a1, b1, c1);
			k = make_key(root, a1, b1, c1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A1-A3-A2 
			curr_key = make_key(root, a1, c1, b1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A1-A3 
			curr_key = make_key(root, b1, a1, c1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A3-A1 
			curr_key = make_key(root, b1, c1, a1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A1-A2 
			curr_key = make_key(root, c1, a1, b1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A2-A1 
			curr_key = make_key(root, c1, b1, a1, ZERO_CHAR, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 30:
		case 32:
		case 33:
		case 37:
		case 45:
		case 47:
		case 48:
		case 53:
		case 63:
		case 64:
		case 67:
		case 69:
		case 73:
			// Make key label
			// P-A-B-C-D 
			k = make_key(root, a, b, c, d, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			break;

		case 20:
		case 27:
		case 28:
		case 29:
		case 31:
			// Make key label
			if (a < b)  {
				// P-A1-A2-B1-B2 
				k = make_key(root, a, b, c, d, g_type);
				//Insert into permutations subset
				mismatches.push_back(k);
				// Generate additional allowed permutations for this graphlet instance.
				// P-A2-A1-B2-B1 
				curr_key = make_key(root, b, a, d, c, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
			}
			else if (a > b)  {
				// P-A1-A2-B1-B2 
				k = make_key(root, b, a, d, c, g_type);
				//Insert into permutations subset
				mismatches.push_back(k);
				// Generate additional allowed permutations for this graphlet instance.
				// P-A2-A1-B2-B1 
				curr_key = make_key(root, a, b, c, d, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
			}
			else  {
				if (c <= d)  {
					// P-A1-A2-B1-B2 
					k = make_key(root, a, b, c, d, g_type);
					//Insert into permutations subset
					mismatches.push_back(k);
					// Generate additional allowed permutations for this graphlet instance.
					// P-A2-A1-B2-B1 
					curr_key = make_key(root, b, a, d, c, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
				}
				else  {
					// P-A1-A2-B1-B2 
					k = make_key(root, b, a, d, c, g_type);
					//Insert into permutations subset
					mismatches.push_back(k);
					// Generate additional allowed permutations for this graphlet instance.
					// P-A2-A1-B2-B1 
					curr_key = make_key(root, a, b, c, d, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
				}
			}
			break;
			
		case 17:
		case 25:
		case 34:
		case 35:
		case 36:
		case 38:
			// Make key label
			// P-A1-A2-B1-B2 
			compare_two(a, b, a1, b1);
			compare_two(c, d, c1, d1);
			k = make_key(root, a1, b1, c1, d1, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A1-A2-B2-B1 
			curr_key = make_key(root, a1, b1, d1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A1-B1-B2 
			curr_key = make_key(root, b1, a1, c1, d1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A1-B2-B1 
			curr_key = make_key(root, b1, a1, d1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 43:
		case 44:
		case 49:
		case 54:
		case 65:
		case 66:
			// Make key label
			// P-A1-A2-B-C 
			compare_two(a, b, a1, b1);
			k = make_key(root, a1, b1, c, d, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A2-A1-B-C 
			curr_key = make_key(root, b1, a1, c, d, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 19:
		case 24:
		case 46:
		case 50:
		case 52:
		case 55:
		case 56:
		case 57:
		case 68:
		case 70:
			// Make key label
			// P-A-B1-B2-C 
			compare_two(b, c, b1, c1);
			k = make_key(root, a, b1, c1, d, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A-B2-B1-C 
			curr_key = make_key(root, a, c1, b1, d, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 39:
		case 40:
		case 41:
		case 42:
		case 60:
		case 61:
		case 71:
		case 72:
			// Make key label
			// P-A-B-C1-C2 
			compare_two(c, d, c1, d1);
			k = make_key(root, a, b, c1, d1, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A-B-C2-C1 
			curr_key = make_key(root, a, b, d1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 18:
        case 23:
			// Make key label
            compare_two(a, b, a1, b1);
            compare_two(c, d, c1, d1);
            if (a1 < c1)  {			
				// P-A1-A2-a1-a2                         
				k = make_key(root, a1, b1, c1, d1, g_type);
				//Insert into permutations subset
				mismatches.push_back(k);
				// Generate additional allowed permutations for this graphlet instance.
				// P-A1-A2-a2-a1
				curr_key = make_key(root, a1, b1, d1, c1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-A2-A1-a1-a2 
				curr_key = make_key(root, b1, a1, c1, d1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-A2-A1-a2-a1 
				curr_key = make_key(root, b1, a1, d1, c1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a1-a2-A1-A2 
				curr_key = make_key(root, c1, d1, a1, b1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a1-a2-A2-A1 
				curr_key = make_key(root, c1, d1, b1, a1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a2-a1-A1-A2
				curr_key = make_key(root, d1, c1, a1, b1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a2-a1-A2-A1
				curr_key = make_key(root, d1, c1, b1, a1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
            }
            else if (a1 > c1)  {
				// P-A1-A2-a1-a2                         
				k = make_key(root, c1, d1, a1, b1, g_type);
				//Insert into permutations subset
				mismatches.push_back(k);
				// Generate additional allowed permutations for this graphlet instance.
				// P-A1-A2-a2-a1
				curr_key = make_key(root, a1, b1, d1, c1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-A2-A1-a1-a2 
				curr_key = make_key(root, b1, a1, c1, d1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-A2-A1-a2-a1 
				curr_key = make_key(root, b1, a1, d1, c1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a1-a2-A1-A2 
				curr_key = make_key(root, a1, b1, c1, d1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a1-a2-A2-A1 
				curr_key = make_key(root, c1, d1, b1, a1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a2-a1-A1-A2
				curr_key = make_key(root, d1, c1, a1, b1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
				// P-a2-a1-A2-A1
				curr_key = make_key(root, d1, c1, b1, a1, g_type);
				//Insert into permutations subset
				insert_permutation(curr_key, mismatches);
            }
            else  {
                if (b1 <= d1)  {
					// P-A1-A2-a1-a2                         
					k = make_key(root, a1, b1, c1, d1, g_type);
					//Insert into permutations subset
					mismatches.push_back(k);
					// Generate additional allowed permutations for this graphlet instance.
					// P-A1-A2-a2-a1
					curr_key = make_key(root, a1, b1, d1, c1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-A2-A1-a1-a2 
					curr_key = make_key(root, b1, a1, c1, d1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-A2-A1-a2-a1 
					curr_key = make_key(root, b1, a1, d1, c1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a1-a2-A1-A2 
					curr_key = make_key(root, c1, d1, a1, b1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a1-a2-A2-A1 
					curr_key = make_key(root, c1, d1, b1, a1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a2-a1-A1-A2
					curr_key = make_key(root, d1, c1, a1, b1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a2-a1-A2-A1
					curr_key = make_key(root, d1, c1, b1, a1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
                }
                else  {
					// P-A1-A2-a1-a2                         
					k = make_key(root, c1, d1, a1, b1, g_type);
					//Insert into permutations subset
					mismatches.push_back(k);
					// Generate additional allowed permutations for this graphlet instance.
					// P-A1-A2-a2-a1
					curr_key = make_key(root, a1, b1, d1, c1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-A2-A1-a1-a2 
					curr_key = make_key(root, b1, a1, c1, d1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-A2-A1-a2-a1 
					curr_key = make_key(root, b1, a1, d1, c1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a1-a2-A1-A2 
					curr_key = make_key(root, a1, b1, c1, d1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a1-a2-A2-A1 
					curr_key = make_key(root, c1, d1, b1, a1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a2-a1-A1-A2
					curr_key = make_key(root, d1, c1, a1, b1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
					// P-a2-a1-A2-A1
					curr_key = make_key(root, d1, c1, b1, a1, g_type);
					//Insert into permutations subset
					insert_permutation(curr_key, mismatches);
                }
            }
            break;

		case 22:
		case 51:
		case 58:
			// Make key label
			// P-A1-A2-A3-B 
			compare_three(a, b, c, a1, b1, c1);
			k = make_key(root, a1, b1, c1, d, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A1-A3-A2-B 
			curr_key = make_key(root, a1, c1, b1, d, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A1-A3-B 
			curr_key = make_key(root, b1, a1, c1, d, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A3-A1-B 
			curr_key = make_key(root, b1, c1, a1, d, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A1-A2-B 
			curr_key = make_key(root, c1, a1, b1, d, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A2-A1-B 
			curr_key = make_key(root, c1, b1, a1, d, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;				
			
		case 21:
		case 59:
		case 62:
			// Make key label
			// P-A-B1-B2-B3 
			compare_three(b, c, d, b1, c1, d1);
			k = make_key(root, a, b1, c1, d1, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A1-A3-A2-B 
			curr_key = make_key(root, a, b1, d1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A1-A3-B 
			curr_key = make_key(root, a, c1, b1, d1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A3-A1-B 
			curr_key = make_key(root, a, c1, d1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A1-A2-B 
			curr_key = make_key(root, a, d1, b1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A2-A1-B 
			curr_key = make_key(root, a, d1, c1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;
			
		case 16:
		case 26:
			// Make key label
			// P-A1-A2-A3-A4 
			compare_four(a, b, c, d, a1, b1, c1, d1);
			k = make_key(root, a1, b1, c1, d1, g_type);
			//Insert into permutations subset
			mismatches.push_back(k);
			// Generate additional allowed permutations for this graphlet instance.
			// P-A1-A2-A4-A3 
			curr_key = make_key(root, a1, b1, d1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A1-A3-A2-A4 
			curr_key = make_key(root, a1, c1, b1, d1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A1-A3-A4-A2 
			curr_key = make_key(root, a1, c1, d1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A1-A4-A2-A3 
			curr_key = make_key(root, a1, d1, b1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A1-A4-A3-A2 
			curr_key = make_key(root, a1, d1, c1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A1-A3-A4 
			curr_key = make_key(root, b1, a1, c1, d1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A1-A3-A4 
			curr_key = make_key(root, b1, a1, d1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A3-A1-A4 
			curr_key = make_key(root, b1, c1, a1, d1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A3-A4-A1 
			curr_key = make_key(root, b1, c1, d1, a1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A4-A1-A3 
			curr_key = make_key(root, b1, d1, a1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A2-A4-A3-A1 
			curr_key = make_key(root, b1, d1, c1, a1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A1-A2-A4 
			curr_key = make_key(root, c1, a1, b1, d1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A1-A4-A2 
			curr_key = make_key(root, c1, a1, d1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A2-A1-A4 
			curr_key = make_key(root, c1, b1, a1, d1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A2-A4-A1 
			curr_key = make_key(root, c1, b1, d1, a1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A4-A1-A2 
			curr_key = make_key(root, c1, d1, a1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A3-A4-A2-A1 
			curr_key = make_key(root, c1, d1, b1, a1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A4-A1-A2-A3 
			curr_key = make_key(root, d1, a1, b1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A4-A1-A3-A2 
			curr_key = make_key(root, d1, a1, c1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A4-A2-A1-A3 
			curr_key = make_key(root, d1, b1, a1, c1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A4-A2-A3-A1 
			curr_key = make_key(root, d1, b1, c1, a1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A4-A3-A1-A2 
			curr_key = make_key(root, d1, c1, a1, b1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			// P-A4-A3-A2-A1 
			curr_key = make_key(root, d1, c1, b1, a1, g_type);
			//Insert into permutations subset
			insert_permutation(curr_key, mismatches);
			break;			
	} //End of switch

	return k;
}

void insert_graphlet_mismatch_neighborhood(list<Key> &neighborhood, Key k)  {
    list<Key>::iterator list_it;
    bool found(false);
    for (list_it = neighborhood.begin(); list_it != neighborhood.end(); list_it++)  {
        if (k == *list_it)  {
            found = true;
            return;
        }
    }
    neighborhood.push_back(k);
}

inline float get_sim_score(char vlabel1, char vlabel2, map<string, float> &sim_vlm_matrix)  {    
    string lookup;
    lookup += toupper(vlabel1);
    lookup += toupper(vlabel2);
    map<string,float>::iterator git = sim_vlm_matrix.find(lookup);
    if(git != sim_vlm_matrix.end())  {
        return git->second;
    }
    return 0.0;
}

void generate_graphlet_mismatch_neighborhood_m1(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, unsigned long g_type, Key key)  {
    Key k;
    char root, a, b, c, d;
	vector<Key> mismatches;
    float sim_score(0.0);

	initialize_vertices_labels(key, root, a, b, c, d);
    
    for (unsigned i=0; i<ALPHABET_ROOT.length(); i++)  {
        sim_score = get_sim_score(root, ALPHABET_ROOT[i], sim_vlm_matrix);
        if (ALPHABET_ROOT[i] != root && sim_score >= SIMILARITY_THRESHOLD)  {
            k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], a, b, c, d, g_type);
            insert_graphlet_mismatch_neighborhood(neighborhood, k);
        }
    }

    for (unsigned i=0; i<ALPHABET.length(); i++)  {
        sim_score = get_sim_score(a, ALPHABET[i], sim_vlm_matrix);
        if (ALPHABET[i] != a && sim_score >= SIMILARITY_THRESHOLD)  {
            k = create_permutations_subset(mismatches, root, ALPHABET[i], b, c, d, g_type);
            insert_graphlet_mismatch_neighborhood(neighborhood, k);
        }

        if (get_graphlet_length(g_type) > 2 && ALPHABET[i] != b)  {
            sim_score = get_sim_score(b, ALPHABET[i], sim_vlm_matrix);
            if (sim_score >= SIMILARITY_THRESHOLD)  {
                k = create_permutations_subset(mismatches, root, a, ALPHABET[i], c, d, g_type);
                insert_graphlet_mismatch_neighborhood(neighborhood, k);
            }
        }

        if (get_graphlet_length(g_type) > 3 && ALPHABET[i] != c)  {
            sim_score = get_sim_score(c, ALPHABET[i], sim_vlm_matrix);
            if (sim_score >= SIMILARITY_THRESHOLD)  {
                k = create_permutations_subset(mismatches, root, a, b, ALPHABET[i], d, g_type);
                insert_graphlet_mismatch_neighborhood(neighborhood, k);
            }
        }

        if (get_graphlet_length(g_type) > 4 && ALPHABET[i] != d)  {
            sim_score = get_sim_score(d, ALPHABET[i], sim_vlm_matrix);        
            if (sim_score >= SIMILARITY_THRESHOLD)  {
                k = create_permutations_subset(mismatches, root, a, b, c, ALPHABET[i], g_type);
                insert_graphlet_mismatch_neighborhood(neighborhood, k);
            }
        }
    }
}

void generate_graphlet_mismatch_neighborhood_m2(list<Key> &neighborhood, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, unsigned long g_type, Key key)  {
    Key k;
    char root, a, b, c, d;
	vector<Key> mismatches;
    float sim_score1(0.0), sim_score2(0.0);

	initialize_vertices_labels(key, root, a, b, c, d);

	for (unsigned i=0; i<ALPHABET_ROOT.length(); i++)  {
        sim_score1 = get_sim_score(root, ALPHABET_ROOT[i], sim_vlm_matrix);
        if (sim_score1 >= SIMILARITY_THRESHOLD)  {
            for (unsigned j=0; j<ALPHABET.length(); j++)  {
                sim_score2 = get_sim_score(a, ALPHABET[j], sim_vlm_matrix);
                if (sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], ALPHABET[j], b, c, d, g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }

                if (get_graphlet_length(g_type) > 2)  {
                    sim_score2 = get_sim_score(b, ALPHABET[j], sim_vlm_matrix);
                    if (sim_score2 >= SIMILARITY_THRESHOLD)  {
                        k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], a, ALPHABET[j], c, d, g_type);
                        insert_graphlet_mismatch_neighborhood(neighborhood, k);
                    }
                }

                if (get_graphlet_length(g_type) > 3)  {
                    sim_score2 = get_sim_score(c, ALPHABET[j], sim_vlm_matrix);
                    if (sim_score2 >= SIMILARITY_THRESHOLD)  {
                        k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], a, b, ALPHABET[j], d, g_type);
                        insert_graphlet_mismatch_neighborhood(neighborhood, k);
                    }
                }

                if (get_graphlet_length(g_type) > 4)  {
                    sim_score2 = get_sim_score(d, ALPHABET[j], sim_vlm_matrix);
                    if (sim_score2 >= SIMILARITY_THRESHOLD)  {
                        k = create_permutations_subset(mismatches, ALPHABET_ROOT[i], a, b, c, ALPHABET[j], g_type);
                        insert_graphlet_mismatch_neighborhood(neighborhood, k);
                    }
                }
            }
        }
    }

    for (unsigned i=0; i<ALPHABET.length(); i++)  { 
        for (unsigned j=0; j<ALPHABET.length(); j++)  {
            if (get_graphlet_length(g_type) > 2)  {
                sim_score1 = get_sim_score(a, ALPHABET[i], sim_vlm_matrix);
                sim_score2 = get_sim_score(b, ALPHABET[j], sim_vlm_matrix);
                if (sim_score1 >= SIMILARITY_THRESHOLD && sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, root, ALPHABET[i], ALPHABET[j], c, d, g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }
            }
            if (get_graphlet_length(g_type) > 3)  {
                sim_score1 = get_sim_score(a, ALPHABET[i], sim_vlm_matrix);
                sim_score2 = get_sim_score(c, ALPHABET[j], sim_vlm_matrix);
                if (sim_score1 >= SIMILARITY_THRESHOLD && sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, root, ALPHABET[i], b, ALPHABET[j], d, g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }

                sim_score1 = get_sim_score(b, ALPHABET[i], sim_vlm_matrix);
                sim_score2 = get_sim_score(c, ALPHABET[j], sim_vlm_matrix);
                if (sim_score1 >= SIMILARITY_THRESHOLD && sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, root, a, ALPHABET[i], ALPHABET[j], d, g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }
            }
            if (get_graphlet_length(g_type) > 4)  {
                sim_score1 = get_sim_score(a, ALPHABET[i], sim_vlm_matrix);
                sim_score2 = get_sim_score(d, ALPHABET[j], sim_vlm_matrix);
                if (sim_score1 >= SIMILARITY_THRESHOLD && sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, root, ALPHABET[i], b, c, ALPHABET[j], g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }

                sim_score1 = get_sim_score(b, ALPHABET[i], sim_vlm_matrix);
                sim_score2 = get_sim_score(d, ALPHABET[j], sim_vlm_matrix);
                if (sim_score1 >= SIMILARITY_THRESHOLD && sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, root, a, ALPHABET[i], c, ALPHABET[j], g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }

                sim_score1 = get_sim_score(c, ALPHABET[i], sim_vlm_matrix);
                sim_score2 = get_sim_score(d, ALPHABET[j], sim_vlm_matrix);
                if (sim_score1 >= SIMILARITY_THRESHOLD && sim_score2 >= SIMILARITY_THRESHOLD)  {
                    k = create_permutations_subset(mismatches, root, a, b, ALPHABET[i], ALPHABET[j], g_type);
                    insert_graphlet_mismatch_neighborhood(neighborhood, k);
                }
            }
        }
    }
}

// Generate corresponding vertex label mismatch graphlets for each graphlet found.	 
void generate_vertex_label_mismatch_graphlets(map<Key, list<Key> > &vl_mismatch_neighborhood, map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, Key key, unsigned long g_type, string ALPHABET_ROOT, string ALPHABET, map<string, float> &sim_vlm_matrix, int VLM)  {
    char root, a, b, c, d;
	list<Key>::iterator list_it;
	map<Key, list<Key> >::iterator lit;

	if ((lit = vl_mismatch_neighborhood.find(key)) == vl_mismatch_neighborhood.end())  {
		if (VLM == 1)  {
			generate_graphlet_mismatch_neighborhood_m1(vl_mismatch_neighborhood[key], ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, g_type, key);
		}
		if (VLM == 2)  {
			generate_graphlet_mismatch_neighborhood_m2(vl_mismatch_neighborhood[key], ALPHABET_ROOT, ALPHABET, sim_vlm_matrix, g_type, key);
		}
    }

    for (list_it = vl_mismatch_neighborhood[key].begin(); list_it != vl_mismatch_neighborhood[key].end(); list_it++)  {
		map<Key,MismatchInfo>::iterator it;
		map<Key,MismatchInfo>::iterator mit;
		if (((it = hash.find(*list_it)) == hash.end()) && ((mit = mismatch_hash.find(*list_it)) == mismatch_hash.end()))  {
            vector<Key> mismatches;
			initialize_vertices_labels(*list_it, root, a, b, c, d);
			Key k = create_permutations_subset(mismatches, root, a, b, c, d, g_type);
			//Insert mismatch graphlets into hash
			insert_mismatches_hash(mismatch_hash, k, mismatches);
		}
	}		
}

void increment_mismatch_count(map<Key,MismatchInfo> &hash, const Key &key, const Key &mismatch_key, float sim_score, float mult_factor)  {
	map<Key,MismatchInfo>::iterator it = hash.find(key);
    if (sim_score >= SIMILARITY_THRESHOLD)  {
        if (it != hash.end())  {
		    map<Key,float>::iterator mismatch = hash[it->first].mismatchesGraph.find(mismatch_key);
		    if (mismatch != hash[key].mismatchesGraph.end())  {
			    hash[it->first].mismatchesGraph[mismatch_key] = hash[it->first].mismatchesGraph[mismatch_key] + (mult_factor * sim_score);
		    }
			else  {
				cerr << "ERROR: Mismatch graphlet permutation " << print_key(mismatch_key) << " cannot be found among subset of allowed permutations for mismatch graphlet " << print_key(key) << "." << endl; exit(1);
			}
	    }
		else  {
            cerr << "ERROR: Mismatch derived graphlet " << print_key(key) << " cannot be found." << endl; exit(1);
        }
    }
}

float compare_graphlets(Key key1, Key key2, unsigned long g_type, map<string, float> &sim_vlm_matrix, float &sim_score)  {
	string key1_str = get_key(key1), key2_str = get_key(key2);
	float distance(0.0);
    sim_score = 1.0;

	for (unsigned i=0; i<get_graphlet_length(g_type); i++)  {
        if(key1_str[i] != key2_str[i])  {
            distance = distance + 1.0;
            string lookup;
		    lookup += toupper(key2_str[i]);
		    lookup += toupper(key1_str[i]);
            map<string,float>::iterator git = sim_vlm_matrix.find(lookup);
		    if(git != sim_vlm_matrix.end())  {
                sim_score *= git->second;
		    }
		    else {
		        cerr << "ERROR: Substitution pair  " << lookup << " cannot be found in vertex-labels similarity matrix." << endl; exit(1);
		    }
        }
	}
	return distance;
}

void update_mismatch_count(map<Key,MismatchInfo> &hash, Key k, float mult_factor, unsigned long g_type, map<string, float> sim_vlm_matrix, int VLM, bool eq)  {
	float curr_dist, min_dist;
    float min_score, sim_score;
	Key min_key, min_graphlet;
	
	if (mult_factor > 0.0 && VLM > 0)  {
        for (map<Key,MismatchInfo>::iterator it = hash.begin(); it != hash.end(); it++)  {
			if (it->first != k) {
                min_dist = float(get_graphlet_length(g_type)) + 1.0;
                min_score = 0.0;
				min_graphlet = make_key(ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, g_type);
				min_key = make_key(ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, ZERO_CHAR, g_type);
				for (map<Key,float>::iterator mismatches_hash = hash[(it->first)].mismatchesGraph.begin(); mismatches_hash != hash[(it->first)].mismatchesGraph.end(); mismatches_hash++)  {
					curr_dist = compare_graphlets(mismatches_hash->first, k, g_type, sim_vlm_matrix, sim_score);
					if (curr_dist < min_dist || (curr_dist == min_dist && sim_score > min_score))  {
						min_dist = curr_dist;
						min_key = it->first;
						min_graphlet = mismatches_hash->first;
						min_score = sim_score;
					}
				}
                if (eq)  {
			        if (min_dist == VLM)  {
				        increment_mismatch_count(hash, min_key, min_graphlet, min_score, mult_factor);
			        }
                }
                else  {
			        if (min_dist <= VLM)  {
				        increment_mismatch_count(hash, min_key, min_graphlet, min_score, mult_factor);
			        }
                }
            }
		}
    }
}

void increment_edge_mismatch_hash(map<Key,MismatchInfo> &hash, map<Key,MismatchInfo> &mismatch_hash, const Key &k, float mult_factor, vector<Key> &mismatches_list)  {
    map<Key,MismatchInfo>::iterator it = hash.find(k);
    map<Key,MismatchInfo>::iterator mit = mismatch_hash.find(k);
    if (it != hash.end())  {
        hash[k].mismatches = it->second.mismatches + mult_factor;
    }
    else  {
        if (mit == mismatch_hash.end())  {
            mismatch_hash[k].matches = 0.0;
            mismatch_hash[k].mismatches = mult_factor;
            for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
                mismatch_hash[k].mismatchesGraph[*mismatch_key] = 0.0;
            }
            mismatches_list.clear();
        }
        else  {
            mismatch_hash[k].mismatches = mit->second.mismatches + mult_factor;
        }
    }
}

void insert_edge_mismatch_graphlet(list<pair <unsigned long, Key> > &EM_set, pair <unsigned long, Key> graphlet, unsigned index)  {
    list<pair <unsigned long, Key> > :: iterator list_it;
    bool found(false);
    for (list_it = EM_set.begin(); list_it != EM_set.end(); list_it++)  {
        if (list_it->first == graphlet.first && list_it->second == graphlet.second)  {
            found = true; 
            break;
        }
    }
    if (!found)
        EM_set.push_back(graphlet);
}

// Update mismatch graphlet counter.
void update_edge_mismatch_count(vector<list<pair <unsigned long, Key> > > &EM_set, Key k, unsigned g_type, unsigned EDGE_MISMATCHES_ALLOWED, unsigned vindex)  {
	Key key(k);
	char root, a, b, c, d;
    pair <unsigned long, Key> p (0,0);
	vector<Key> mismatches;

    if (EDGE_MISMATCHES_ALLOWED <= 0)
        return;    

	initialize_vertices_labels(key, root, a, b, c, d);
	
	switch (g_type)  {
		case 0: // P 
			// There are no edges addition/removal for this case.
			break;
			
		case 1: // P-A 
			// There are no edges addition/removal for this case.
			break;
			
		case 2: // P-A-B 
			// P-A-A (Edge Mismatch by adding edge from P to B)
            p.first = 4;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);        
            update_edge_mismatch_count(EM_set, p.second, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 3: // P-A1-A2 
			// P-A-A (Edge mismatch by adding edge from A1 to A2)
            p.first = 4;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 4, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 4: // P-A1-A2
			// P-A-B (Edge Mismatch by removing edge from P to A1)
            p.first = 2;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);            
            update_edge_mismatch_count(EM_set, p.second, 2, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A (Edge Mismatch by removing edge from A1 to A2)
            p.first = 3;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);         
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 3, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B (Edge Mismatch by removing edge from P to A2)
            p.first = 2;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);        
            update_edge_mismatch_count(EM_set, p.second, 2, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 5: // P-A-B-C 
			// P-A-B-C (Edge Mismatch by adding edge P to B)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
			update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B (Edge Mismatch by adding edge P to C)
            p.first = 12;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-B (Edge Mismatch by adding edge A to C)
            p.first = 9;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);			
            update_edge_mismatch_count(EM_set, p.second, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 6: // P-A-B-C 
			// P-A-B-B (Edge Mismatch by adding edge P to C)
            p.first = 10;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			        
			// P-A-B-C (Edge Mismatch by adding edge A to B)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-A-B (Edge Mismatch by adding edge A to C)
            p.first = 12;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);			
            update_edge_mismatch_count(EM_set, p.second, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 7: // P-A-B1-B2 
			// P-A-B-C (Edge Mismatch by adding edge P to B1)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B (Edge Mismatch by adding edge B1 to B2)
            p.first = 9;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by adding edge P to B2)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 8: // P-A1-A2-A3
			// P-A-B-B (Edge Mismatch by adding edge A1 to A2)
            p.first = 10;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B (Edge Mismatch by adding edge A2 to A3)
            p.first = 10;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B (Edge Mismatch by adding edge A1 to A3)
            p.first = 10;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 9: // P-A-B1-B2 
			// P-A-A-B (Edge Mismatch by adding edge P to B1)
            p.first = 13;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
		
			// P-A-B-C (Edge Mismatch by removing edge A to B1)
            p.first = 5;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B (Edge Mismatch by removing edge B1 to B2)
            p.first = 7;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B (Edge Mismatch by adding edge P to B2)
            p.first = 13;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by removing edge A to B2)
            p.first = 5;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 10: // P-A-B1-B2 
			// P-A-A-B (Edge Mismatch by adding edge A to B1)
            p.first = 14;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by removing edge P to B1)
            p.first = 6;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-A (Edge Mismatch by removing edge B1 to B2)
            p.first = 8;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 8, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B (Edge Mismatch by adding edge A to B2)
            p.first = 14;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by removing edge P to B2)
            p.first = 6;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 11: // P-A-B-C 
			// P-A-A-B (Edge Mismatch by adding edge P to C)
            p.first = 14;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-A-B (Edge Mismatch by adding edge A to C)
            p.first = 13;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B (Edge Mismatch by removing edge P to A)
            p.first = 7;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 7, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-C (Edge Mismatch by removing edge P to B)
            p.first = 5;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by removing edge A to B)
            p.first = 6;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 12: // P-A1-A2-B
            // P-A-A-B (Edge Mismatch by adding edge P to B)
            p.first = 14;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B (Edge Mismatch by adding edge A1 to A2)
            p.first = 13;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by removing edge P to A1)
            p.first = 5;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-C (Edge Mismatch by removing edge A1 to B)
            p.first = 6;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-C (Edge Mismatch by removing edge P to A2)
            p.first = 5;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 5, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
			// P-A-B-C (Edge Mismatch by removing edge A2 to B)
            p.first = 6;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 6, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 13: // P-A1-A2-B 
			// P-A-A-A (Edge Mismatch by adding edge P to B)
            p.first = 15;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B (Edge Mismatch by removing edge P to A1)
            p.first = 9;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
            // P-A-A-B (Edge Mismatch by removing edge A1 to A2)
            p.first = 12;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by removing edge A1 to B)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-B (Edge Mismatch by removing edge P to A2)
            p.first = 9;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 9, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
			// P-A-B-C (Edge Mismatch by removing edge A2 to B)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 14: // P-A1-A2-B 
			// P-A-A-A (Edge Mismatch by adding edge A1 to A2)
            p.first = 15;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 15, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C (Edge Mismatch by removing edge P to A1)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-A-B (Edge Mismatch by removing edge P to B)
            p.first = 12;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 12, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-B (Edge Mismatch by removing edge A1 to B)
            p.first = 10;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-C (Edge Mismatch by removing edge P to A2)
            p.first = 11;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 11, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B (Edge Mismatch by removing edge A2 to B)
            p.first = 10;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 10, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
			
		case 15: // P-A1-A2-A3 
			// P-A-A-B (Edge Mismatch by removing edge P to A1)
            p.first = 13;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-A-B (Edge Mismatch by removing edge A1 to A2)
            p.first = 14;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B (Edge Mismatch by removing edge P to A2)
            p.first = 13;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B (Edge Mismatch by removing edge A2 to A3)
            p.first = 14;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-A-B (Edge Mismatch by removing edge P to A3)
            p.first = 13;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 13, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B (Edge Mismatch by removing edge A1 to A3)
            p.first = 14;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 14, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;
        
		case 16: // P-A1-A2-A3-A4 
			// P-A-A-B-B (Edge Mismatch by adding edge A1 to A2)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, c, d, a, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-A-B-B (Edge Mismatch by adding edge A2 to A3)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-A-B-B (Edge Mismatch by adding edge A2 to A4)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
	        // P-A-A-B-B (Edge Mismatch by adding edge A1 to A3)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge A3 to A4)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
            // P-A-A-B-B (Edge Mismatch by adding edge A1 to A4)	
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);		
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;                

		case 17: // P-A1-A2-B1-B2 
			// P-A-B-C-D (Edge Mismatch by adding edge A1 to A2)
            p.first = 18;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 18, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-B-C (Edge Mismatch by adding edge A1 to B1)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to B1)
            p.first = 43;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 43, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-A-A (Edge Mismatch by removing edge B1 to B2)
            p.first = 16;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 16, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by adding edge A1 to B2)            
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-A-B-C (Edge Mismatch by removing edge P to B2)
            p.first = 43;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 43, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));            
			
			// P-A-B-B-C (Edge Mismatch by adding edge A2 to B1)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by adding edge A2 to B2)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;        

		case 18: // P-A-B-C-D 
            // P-A-B-C-D (Edge Mismatch by adding edge A to C)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-A-B-C (Edge Mismatch by removing edge P to A)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, c, d, b, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge A to B)            
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge C to D)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, c, d, a, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to C)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge A to D)                        
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-A-B-C (Edge Mismatch by removing edge P to D)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-C-D (Edge Mismatch by adding edge B to C)                        
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to B)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, c, d, a, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge B to D)                        
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;        

		case 19: // P-A-B1-B2-C
            // P-A-B-C-D (Edge Mismatch by adding edge A to B1)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-A-B (Edge Mismatch by adding edge A to C)
            p.first = 22;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-B (Edge Mismatch by adding edge B1 to B2)
            p.first = 21;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge P to B1)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to C)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge B1 to C)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge A to B2)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, c, d, a, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-C-D (Edge Mismatch by removing edge P to B2)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-A-B-B (Edge Mismatch by removing edge B2 to C)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
    
		case 20: // P-A-B-C-D
            // P-A-B-B-C (Edge Mismatch by adding edge A to D)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, c, b, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by adding edge B to C)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, d, a, c, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-C (Edge Mismatch by adding edge C to D)
            p.first = 23;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, d, b, c, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, c, a, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge P to C)            
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, d, b, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge P to D)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge A to B)
            p.first = 18;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 18, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge A to C)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, c, a, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge B to D)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, d, b, c, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 21: // P-A-B1-B2-B3
            // P-A-B-B-C (Edge Mismatch by adding edge A to B1)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to B1)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge B1 to B2)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-B-C (Edge Mismatch by adding edge A to B2) 
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-B-B-C (Edge Mismatch by removing edge P to B2)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-B-B-C (Edge Mismatch by removing edge B2 to B3)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-B-C (Edge Mismatch by adding edge A to B3)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to B3)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-B-B-C (Edge Mismatch by removing edge B1 to B3)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 22: // P-A1-A2-A3-B
            // P-A-B-B-C (Edge Mismatch by adding edge A1 to A2)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to A1)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-A-B (Edge Mismatch by removing edge P to B)
            p.first = 51;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 51, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge A1 to B)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by adding edge A2 to A3)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-A-B-C (Edge Mismatch by removing edge P to A2)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-B-B-C (Edge Mismatch by removing edge A2 to B)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-B-C (Edge Mismatch by adding edge A1 to A3)            
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to A3)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge A3 to B)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 23: // P-A-B-C1-C2
            // P-A-A-B-B (Edge Mismatch by adding edge A to B)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge C1 to C2)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, c, d, a, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to A)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to C1)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge A to C1)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, d, c, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to B)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge B to C1)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, a, d, c, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-B-C-D (Edge Mismatch by removing edge B to C2)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to C2)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge A to C2)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 24: // P-A-B1-B2-C
            // P-A-A-B-B (Edge Mismatch by adding edge A to B1)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to A)
            p.first = 54;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 54, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge P to B1)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, d, c, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-B-C (Edge Mismatch by removing edge P to C)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-B (Edge Mismatch by removing edge A to C)
            p.first = 21;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-A-B (Edge Mismatch by removing edge B1 to B2)
            p.first = 22;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge B1 to C)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, c, d, b, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-A-B-B (Edge Mismatch by adding edge A to B2)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, c, d, a, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-C-D (Edge Mismatch by removing edge P to B2)            
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));            

            // P-A-B-C-D (Edge Mismatch by removing edge B2 to C)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, d, c, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 25: // P-A1-A2-B1-B2
            // P-A-A-A-A (Edge Mismatch by adding edge B1 to B2)
            p.first = 26;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to A1)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge P to B1)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-C (Edge Mismatch by removing edge A1 to A2)
            p.first = 23;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge A1 to B1)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, c, a, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-B-C (Edge Mismatch by removing edge P to A2)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-B-B-C (Edge Mismatch by removing edge A2 to B1)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, c, b, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-B-B-C (Edge Mismatch by removing edge A2 to B2)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, d, b, c, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
        
            // P-A-B-B-C (Edge Mismatch by removing edge P to B2)            
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-B-C (Edge Mismatch by removing edge A1 to B2)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, d, a, c, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 26: // P-A1-A2-A3-A4
            // P-A-A-A-B (Edge Mismatch by removing edge P to A1)
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge A1 to A2)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, c, d, a, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-A-A-B (Edge Mismatch by removing edge P to A2)
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge A2 to A3)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-A-B-B (Edge Mismatch by removing edge A2 to A4)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-A-B (Edge Mismatch by removing edge P to A3)            
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge A1 to A3)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
            // P-A-A-B-B (Edge Mismatch by removing edge A3 to A4)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-A-B (Edge Mismatch by removing edge P to A4)            
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by removing edge A1 to A4)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 27: // P-A-B-C-D
			// P-A-A-B-C (Edge Mismatch by adding edge P to C)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge P to D)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to B)
            p.first = 28;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-C-D (Edge Mismatch by adding edge A to D)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
            // P-A-B-C-D (Edge Mismatch by adding edge B to C)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge C to D)
            p.first = 29;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
            
		case 28: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to C)            
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to D)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge C to D)
            p.first = 31;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B)
            p.first = 27;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)            
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, d, b, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));          

            // P-A-B-C-D (Edge Mismatch by adding edge B to C)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 29: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to C)            
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)            
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge A to B)
            p.first = 31;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to D)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
            // P-A-B-C-D (Edge Mismatch by adding edge B to C)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, b, d, c, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to C)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to D)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge C to D)
            p.first = 27;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 30: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to C) 
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, c, b, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)            
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge A to B)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by adding edge A to C)
            p.first = 34;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 34, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
            // P-A-B-C-D (Edge Mismatch by adding edge C to D)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge A to D)
            p.first = 39;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 39, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to D)
            p.first = 27;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 31: // P-A-B-C-D
			// P-A-B-B-C (Edge Mismatch by adding edge P to C)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge P to D)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge A to D)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge B to C)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-B-C (Edge Mismatch by removing edge P to A)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge P to B)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B)
            p.first = 29;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to C)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to D)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge C to D)
            p.first = 28;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 32: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to C)            
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, c, b, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge P to D)
            p.first = 54;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 54, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge A to C)
            p.first = 36;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 36, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge C to D)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-C-C (Edge Mismatch by removing edge P to A)
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to D)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to D)
            p.first = 28;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
            
		case 33: // P-A-B-C-D
			// P-A-B-B-C (Edge Mismatch by adding edge P to C)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)            
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge A to B)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge A to C)
            p.first = 35;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 35, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A)        
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to B)
            p.first = 72;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 72, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge A to D)
            p.first = 40;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 40, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B to C)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to D)
            p.first = 29;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge C to D)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 34: // P-A1-A2-B1-B2
			// P-A-B-B-C (Edge Mismatch by adding edge P to B1)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge A1 to A2)
            p.first = 36;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 36, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by adding edge B1 to B2)
            p.first = 35;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 35, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			
			// P-A-B-B-C (Edge Mismatch by removing edge P to A1)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B1)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-B-C (Edge Mismatch by adding edge P to B2)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B2)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-B-C (Edge Mismatch by removing edge P to A2)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B1)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));          
                
			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B2)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 35: // P-A1-A2-B1-B2
			// P-A-B-B-C (Edge Mismatch by adding edge P to B1)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge A1 to A2)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge P to A1)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B1)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by removing edge B1 to B2)
            p.first = 34;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 34, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
        
			// P-A-B-B-C (Edge Mismatch by adding edge P to B2)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B2)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-B-C (Edge Mismatch by removing edge P to A2)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B1)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B2)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 36: // P-A1-A2-B1-B2
			// P-A-B-B-C (Edge Mismatch by adding edge P to B1)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge B1 to B2)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to A1)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by removing edge A1 to A2)
            p.first = 34;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 34, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B1)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-B-C (Edge Mismatch by adding edge P to B2)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
              
			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B2)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
        
			// P-A-B-C-C (Edge Mismatch by removing edge P to A2)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B1)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B2)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 37: // P-A-B-C-D
			// P-A-B-B-C (Edge Mismatch by adding edge P to C)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge P to D)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-B (Edge Mismatch by adding edge A to C)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to A)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge P to B)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge A to D)
            p.first = 42;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 42, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B to C)
            p.first = 66;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 66, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to D)
            p.first = 31;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge C to D)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 38: // P-A1-A2-B1-B2
			// P-A-A-A-B (Edge Mismatch by adding edge P to B1)
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-B (Edge Mismatch by removing edge P to A1)
            p.first = 62;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 62, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by removing edge A1 to A2)
            p.first = 35;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 35, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B1)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by removing edge B1 to B2)
            p.first = 36;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 36, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-A-B (Edge Mismatch by adding edge P to B2)
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B2)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
        
			// P-A-B-B-B (Edge Mismatch by removing edge P to A2)
            p.first = 62;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 62, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B1)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B2)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 39: // P-A-B-C1-C2
			// P-A-B-C-D (Edge Mismatch by adding edge P to C1)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge A to B)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C1)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge C1 to C2)
            p.first = 40;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 40, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge P to C2)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge A to C2)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 40: // P-A-B-C1-C2
			// P-A-B-B-C (Edge Mismatch by adding edge P to C1)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge A to B)
            p.first = 42;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 42, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C1)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge B to C1)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge C1 to C2)
            p.first = 39;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 39, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge P to C2)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C2)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                        
			// P-A-B-C-C (Edge Mismatch by removing edge B to C2)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 41: // P-A-B-C1-C2
			// P-A-A-B-C (Edge Mismatch by adding edge P to C1)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge C1 to C2)
            p.first = 42;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 42, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C1)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to A)
            p.first = 59;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 59, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to B)
            p.first = 71;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 71, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge A to B)
            p.first = 39;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 39, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-A-B-C (Edge Mismatch by adding edge P to C2)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C2)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 42: // P-A-B-C1-C2
			// P-A-B-C-D (Edge Mismatch by adding edge P to C1)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C1)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to A)            
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to B)            
            p.first = 72;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 72, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge A to B)
            p.first = 40;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 40, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge B to C1)            
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-C (Edge Mismatch by removing edge C1 to C2)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge P to C2)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C2)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge B to C2)            
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 43: // P-A1-A2-B-C
			// P-A-A-B-B (Edge Mismatch by adding edge P to C)
            p.first = 17;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 17, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge A1 to A2)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A1 to B)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A1 to C)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                   
			// P-A-B-C-D (Edge Mismatch by adding edge A2 to B)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
         
            // P-A-B-B-C (Edge Mismatch by adding edge A2 to C)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 44: // P-A1-A2-B-C
			// P-A-B-C-D (Edge Mismatch by adding edge P to C)
            p.first = 18;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 18, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A1 to B)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A1 to C)            
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A1)
            p.first = 27;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge A1 to A2)
            p.first = 43;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 43, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge A2 to B)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
         
            // P-A-B-C-D (Edge Mismatch by adding edge A2 to C)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A2)
            p.first = 27;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 27, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 45: // P-A-B-C-D
			// P-A-B-B-C (Edge Mismatch by adding edge P to D)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to B)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge A to C)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge A to D)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
         
            // P-A-B-B-C (Edge Mismatch by adding edge B to D)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to B)
            p.first = 39;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 39, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to C)			
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B to C)
            p.first = 43;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 43, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 46: // P-A-B1-B2-C
			// P-A-B-B-C (Edge Mismatch by adding edge P to D)
            p.first = 19;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 19, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to B1)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-A-B (Edge Mismatch by adding edge A to D)
            p.first = 51;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 51, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
         
            // P-A-B-B-C (Edge Mismatch by adding edge B1 to B2)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B1)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B1 to C)
            p.first = 43;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 43, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to B2)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by removing edge P to B2)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B2 to C)
            p.first = 43;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 43, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 47: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to D)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge A to C)
            p.first = 54;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 54, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to D)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge B to D)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
         
            // P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 28;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to C)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B to C)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
            
		case 48: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to D)
            p.first = 20;
            p.second = create_permutations_subset(mismatches, root, b, d, a, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 20, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to C)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to D)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge B to C)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
         
            // P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, c, b, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 29;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to C)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to B)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B to D)
            p.first = 44;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 44, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge C to D)			
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 49: // P-A1-A2-B-C
			// P-A-A-A-B (Edge Mismatch by adding edge P to C)
            p.first = 22;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge A1 to A2)
            p.first = 54;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 54, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A1 to C)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
           
            // P-A-B-C-C (Edge Mismatch by removing edge P to A1)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge P to B)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A2 to C)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                     
			// P-A-B-C-C (Edge Mismatch by removing edge P to A2)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 50: // P-A-B1-B2-C
			// P-A-B-B-B (Edge Mismatch by adding edge P to D)
            p.first = 21;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 21, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to B1)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to D)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
         
            // P-A-B-C-C (Edge Mismatch by removing edge P to B1)
            p.first = 40;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 40, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge B1 to B2)            
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B1 to D)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to B2)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-C (Edge Mismatch by removing edge P to B2)
            p.first = 40;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 40, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B2 to D)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 51: // P-A1-A2-A3-B
			// P-A-A-A-B (Edge Mismatch by adding edge P to B)
            p.first = 22;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 22, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A1 to A2)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to A1)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A1 to B)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by adding edge A2 to A3)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge P to A2)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A2 to B)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
            // P-A-B-B-C (Edge Mismatch by adding edge A1 to A3)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                    
			// P-A-A-B-C (Edge Mismatch by removing edge P to A3)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A3 to B)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 52: // P-A-B1-B2-C
			// P-A-B-B-C (Edge Mismatch by adding edge P to C)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to B1)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-A-B-C (Edge Mismatch by removing edge P to A)
            p.first = 66;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 66, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B1)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to C)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-A-B (Edge Mismatch by removing edge B1 to B2)
            p.first = 51;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 51, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B1 to C)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-B-C (Edge Mismatch by adding edge A to B2)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            			        
			// P-A-B-C-D (Edge Mismatch by removing edge P to B2)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B2 to C)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, c, b, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 53: // P-A-B-C-D
			// P-A-B-B-C (Edge Mismatch by adding edge P to D)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to C)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to D)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            			        
            // P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, c, b, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to C)
            p.first = 42;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 42, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to B)
            p.first = 50;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 50, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to C)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to D)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge C to D)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 54: // P-A1-A2-B-C
			// P-A-B-B-C (Edge Mismatch by adding edge P to C)
            p.first = 24;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 24, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A1 to C)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge P to A1)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, b, c, d, a, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge P to B)
            p.first = 66;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 66, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge A1 to A2)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A2 to C)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            			        
			// P-A-B-C-D (Edge Mismatch by removing edge P to A2)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 55: // P-A-B1-B2-C
			// P-A-B-C-C (Edge Mismatch by adding edge P to C)
            p.first = 23;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 23, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to C)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge B1 to B2)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            			        
            // P-A-A-B-B (Edge Mismatch by removing edge P to A)
            p.first = 34;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 34, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B1)
            p.first = 31;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B1)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B1 to C)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B2)
            p.first = 31;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B2)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B2 to C)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 56: // P-A-B1-B2-C
			// P-A-A-B-B (Edge Mismatch by adding edge P to C)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-A-B (Edge Mismatch by adding edge B1 to B2)
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
           			        
            // P-A-A-B-B (Edge Mismatch by removing edge P to A)
            p.first = 35;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 35, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B1)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to B1)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to C)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B1 to C)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B2)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to B2)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B2 to C)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
            
		case 57: // P-A-B1-B2-C
			// P-A-A-B-B (Edge Mismatch by adding edge P to C)
            p.first = 25;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 25, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-A-B (Edge Mismatch by adding edge A to C)
            p.first = 58;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 58, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
           			        
            // P-A-A-B-B (Edge Mismatch by removing edge P to A)
            p.first = 36;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 36, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B1)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B1)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge B1 to B2)
            p.first = 55;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 55, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B1 to C)
            p.first = 54;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 54, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by removing edge P to B2)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B2)
            p.first = 53;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 53, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge B2 to C)
            p.first = 54;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 54, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));            
            break;

		case 58: // P-A1-A2-A3-B
			// P-A-A-A-B (Edge Mismatch by adding edge P to C)
            p.first = 26;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 26, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			          			        
            // P-A-A-B-B (Edge Mismatch by removing edge P to A1)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge A1 to A2)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A1 to B)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-A-B-B (Edge Mismatch by removing edge P to A2)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A2 to A3)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
                
			// P-A-B-B-C (Edge Mismatch by removing edge A2 to B)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by removing edge P to A3)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A1 to A3)
            p.first = 56;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 56, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A3 to B)
            p.first = 57;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 57, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 59: // P-A-B1-B2-B3
            // P-A-B-C-C (Edge Mismatch by adding edge P to B1)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge B1 to B2)
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge P to B2)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge B2 to B3)
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge P to B3)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));            
			          			                    
			// P-A-B-C-C (Edge Mismatch by adding edge B1 to B3)
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;

		case 60: // P-A-B-C1-C2
			// P-A-B-C-C (Edge Mismatch by adding edge P to B)
            p.first = 42;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 42, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to C1)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, c, a, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-C (Edge Mismatch by adding edge B to C1)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to C1)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-B (Edge Mismatch by removing edge C1 to C2)
            p.first = 59;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 59, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge P to C2)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			          			                    
			// P-A-B-C-C (Edge Mismatch by adding edge B to C2)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to C2)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 61: // P-A-B-C1-C2
			// P-A-B-C-C (Edge Mismatch by adding edge P to B)
            p.first = 36;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 36, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to C1)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, c, a, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-B (Edge Mismatch by adding edge C1 to C2)
            p.first = 62;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 62, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-B-C (Edge Mismatch by removing edge A to B)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by removing edge A to C1)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-C (Edge Mismatch by removing edge B to C1)
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to C2)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, d, a, c, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			          			                    
			// P-A-B-C-D (Edge Mismatch by removing edge A to C2)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge B to C2)
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 62: // P-A-B1-B2-B3
			// P-A-A-B-B (Edge Mismatch by adding edge P to B1)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to B1)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge B1 to B2)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-A-B-B (Edge Mismatch by adding edge P to B2)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to B2)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge B2 to B3)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by adding edge P to B3)
            p.first = 38;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 38, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge A to B3)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge B1 to B3)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 63: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to C)
            p.first = 45;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 45, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge P to D)
            p.first = 46;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 46, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to B)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge A to C)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to D)
            p.first = 29;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge B to D)
            p.first = 40;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 40, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;

		case 64: // P-A-B-C-D 
			// P-A-B-C-D (Edge Mismatch by adding edge P to C)
            p.first = 47;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 47, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)
            p.first = 48;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 48, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge A to C)
            p.first = 66;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 66, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to D)
            p.first = 31;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge B to D)
            p.first = 42;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 42, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to B)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 65: // P-A1-A2-B-C
			// P-A-A-B-C (Edge Mismatch by adding edge P to B)
            p.first = 49;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 49, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-A-B (Edge Mismatch by adding edge P to C)
            p.first = 51;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 51, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge A1 to A2)
            p.first = 66;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 66, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A1 to C)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to A1)
            p.first = 71;
            p.second = create_permutations_subset(mismatches, root, b, c, a, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 71, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A2 to C)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge P to A2)
            p.first = 71;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 71, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B)
            p.first = 63;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 63, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));            
            break;

		case 66: // P-A1-A2-B-C
			// P-A-A-B-C (Edge Mismatch by adding edge P to B)
            p.first = 54;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 54, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge P to C)
            p.first = 52;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 52, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A1 to C)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, b, a, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A1)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by removing edge A1 to A2)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A1 to B)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A2 to C)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge P to A2)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A2 to B)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;
            
		case 67: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to B)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, b, a, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to C)
            p.first = 28;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 28, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge A to D)
            p.first = 60;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 60, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge B to C)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge B to D)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 68: // P-A-B1-B2-C
			// P-A-B-C-D (Edge Mismatch by adding edge P to B1)
            p.first = 31;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by adding edge P to C)
            p.first = 34;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 34, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge A to C)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge B1 to B2)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B1)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B1 to C)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge P to B2)
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 31, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to B2)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B2 to C)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 69: // P-A-B-C-D
			// P-A-A-B-C (Edge Mismatch by adding edge P to B)
            p.first = 66;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 66, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to C)
            p.first = 32;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 32, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, d, a, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge A to D)
            p.first = 61;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 61, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge B to D)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge A to B)
            p.first = 71;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 71, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge A to C)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to C)
            p.first = 67;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 67, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
			break;

		case 70: // P-A-B1-B2-C
			// P-A-B-C-D (Edge Mismatch by adding edge P to B1)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-B (Edge Mismatch by adding edge P to C)
            p.first = 35;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 35, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-B (Edge Mismatch by adding edge A to C)
            p.first = 62;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 62, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge A to B1)
            p.first = 72;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 72, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by removing edge B1 to B2)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge B1 to C)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to B2)
            p.first = 37;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 37, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-C (Edge Mismatch by removing edge A to B2)
            p.first = 72;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 72, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by adding edge B2 to C)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 71: // P-A-B-C1-C2
			// P-A-B-C-C (Edge Mismatch by adding edge P to B)
            p.first = 41;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 41, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge P to C1)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge A to C1)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, c, b, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge C1 to C2)
            p.first = 72;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 72, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-A-B-C (Edge Mismatch by adding edge P to C2)
            p.first = 65;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 65, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C2)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 72: // P-A-B-C1-C2
			// P-A-B-C-C (Edge Mismatch by adding edge P to B)
            p.first = 42;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 42, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to C1)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to C1)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by removing edge B to C1)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by removing edge C1 to C2)
            p.first = 71;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 71, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

            // P-A-B-C-D (Edge Mismatch by adding edge P to C2)
            p.first = 33;
            p.second = create_permutations_subset(mismatches, root, a, d, c, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 33, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to C2)
            p.first = 70;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 70, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            
			// P-A-B-C-D (Edge Mismatch by removing edge B to C2)
            p.first = 73;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 73, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;

		case 73: // P-A-B-C-D
			// P-A-B-C-D (Edge Mismatch by adding edge P to B)
            p.first = 64;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 64, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to C)
            p.first = 30;
            p.second = create_permutations_subset(mismatches, root, a, c, d, b, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 30, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge P to D)
            p.first = 29;
            p.second = create_permutations_subset(mismatches, root, a, d, b, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 29, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-D (Edge Mismatch by adding edge A to C)
            p.first = 69;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 69, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-B-C (Edge Mismatch by adding edge A to D)
            p.first = 68;
            p.second = create_permutations_subset(mismatches, root, a, b, d, c, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 68, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));

			// P-A-B-C-C (Edge Mismatch by adding edge B to D)
            p.first = 72;
            p.second = create_permutations_subset(mismatches, root, a, b, c, d, p.first);
            insert_edge_mismatch_graphlet(EM_set[vindex], p, vindex);
            update_edge_mismatch_count(EM_set, p.second, 72, (EDGE_MISMATCHES_ALLOWED-1), (vindex+1));
            break;        
	}
}

