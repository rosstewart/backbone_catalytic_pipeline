#include "utils.h"
#include <stdlib.h>
#include <iostream>
using namespace std;


float compare_labels(char label1, char label2)  {
    if (label1 == label2)
        return 1.0;
    return 0.0;
}


int randint(int max)  {
    if( max == 0)  {
        cerr << "Potential error in randint()" << endl;
        return -1;
    }
    return rand() % max;
}


double randdouble()  {
    return randint(101) / 100.0;
}


unsigned long get_graphlet_length(unsigned long g_type)  {
    unsigned long g_length;
    if (g_type == 0)  {
        g_length = 1;
    }
    else if (g_type == 1)  {
        g_length = 2;
    }
    else if (g_type >= 2 && g_type <= 4)  {
        g_length = 3; 
    }
    else if (g_type >= 5 && g_type <= 15)  {
        g_length = 4;
    }
    else if (g_type >= 16 && g_type <= 73)  {
        g_length = 5;
    }
    else  {
        cerr << "ERROR: Orbit: " << g_type << " is unsupported." << endl; exit(1);
    }
    return g_length;
}


void compare_two(char &a, char &b, char &a1, char &b1)  {
	if (a < b)  { // a,b
		a1 = a; b1 = b;
	} 
	else  { // b,a
		a1 = b; b1 = a;
	}
}


void compare_three(char &a, char &b, char &c, char &a1, char &b1, char &c1)  {
	if (a < b)  {
		if (b < c)  { // a,b,c
			a1 = a; b1 = b; c1 = c;
		} else if (a < c)  { // a,c,b
			a1 = a; b1 = c; c1 = b;
		} else  { //  c,a,b
			a1 = c; b1 = a; c1 = b;
		}
	} 
	else  {
		if (a < c)  { // b,a,c
			a1 = b; b1 = a; c1 = c;
		} else if (b < c)  { // b,c,a
			a1 = b; b1 = c; c1 = a;
		} else  { // c,b,a
			a1 = c; b1 = b; c1 = a;
		}
	}
}


void compare_four(char &a, char &b, char &c, char &d, char &a1, char &b1, char &c1, char &d1)  {
	if (a < b)  {
		if (b < c)  {
			if(c < d)  { // a,b,c,d
				a1 = a; b1 = b; c1 = c; d1 = d;
			}
			else if (b < d) { // a,b,d,c
				a1 = a; b1 = b; c1 = d; d1 = c;
			}
			else if (a < d) { // a,d,b,c
				a1 = a; b1 = d; c1 = b; d1 = c;
			}
			else { // d,a,b,c
				a1 = d; b1 = a; c1 = b; d1 = c;
			}
		}
		else if (a < c) {
			if(b < d)  { // a,c,b,d
				a1 = a; b1 = c; c1 = b; d1 = d;
			}
			else if (c < d) { // a,c,d,b
				a1 = a; b1 = c; c1 = d; d1 = b;
			}
			else if (a < d) { // a,d,c,b
				a1 = a; b1 = d; c1 = c; d1 = b;
			}
			else { //d,a,c,b
				a1 = d; b1 = a; c1 = c; d1 = b;
			}
		}
		else  {
			if(d < c)  { // d,c,a,b
				a1 = d; b1 = c; c1 = a; d1 = b;
			}
			else if (d < a) { // c,d,a,b
				a1 = c; b1 = d; c1 = a; d1 = b;
			}
			else if (d < b) { // c,a,d,b
				a1 = c; b1 = a; c1 = d; d1 = b;
			}
			else { // c,a,b,d
				a1 = c; b1 = a; c1 = b; d1 = d;
			}	
		}
	}
	else {
		if (a < c)  {
			if(d < b)  { // d,b,a,c
				a1 = d; b1 = b; c1 = a; d1 = c;
			}
			else if (d < a) { // b,d,a,c
				a1 = b; b1 = d; c1 = a; d1 = c;
			}
			else if (d < c) { // b,a,d,c
				a1 = b; b1 = a; c1 = d; d1 = c;
			}
			else { // b,a,c,d
				a1 = b; b1 = a; c1 = c; d1 = d;
			}
		}
		else if (b < c) {
			if(a < d)  { // b,c,a,d
				a1 = b; b1 = c; c1 = a; d1 = d;
			}
			else if (d < b) { // d,b,c,a
				a1 = d; b1 = b; c1 = c; d1 = a;
			}
			else if (d < c) { // b,d,c,a
				a1 = b; b1 = d; c1 = c; d1 = a;
			}
			else { // b,c,d,a
				a1 = b; b1 = c; c1 = d; d1 = a;
			}
		}
		else  {
			if(a < d)  { // c,b,a,d
				a1 = c; b1 = b; c1 = a; d1 = d;
			}
			else if (b < d) { // c,b,d,a
				a1 = c; b1 = b; c1 = d; d1 = a;
			}
			else if (c < d) { // c,d,b,a
				a1 = c; b1 = d; c1 = b; d1 = a;
			}
			else { // d,c,b,a
				a1 = d; b1 = c; c1 = b; d1 = a;
			}	
		}
	}
}


unsigned int set_k(unsigned long g_type, float sf)  {
    unsigned int k_mismatches = 0;
   
    if ((g_type >= 0 && g_type < GRAPHLET_TYPES) && (sf >= 0.0 && sf <= 1.0))  {
        k_mismatches = unsigned(float(get_graphlet_length(g_type)) * sf);
    }
    else  {
        if (g_type < 0 || g_type >= GRAPHLET_TYPES)  {
            cerr << "ERROR: Orbit: " << g_type << " is unsupported." << endl; exit(1);
        }
        else  {
            cerr << "ERROR: Scaling fraction: " << sf << " is invalid." << endl; exit(1);
        }
    }

    return k_mismatches;
}


/* This function first checks if graphlet permutation is already listed.
 * If not, it inserts it into a set of allowed permutations.
 */  
void insert_permutation(const Key &target, vector<Key> &mismatches_list)  {
	bool found(false);

	for (vector<Key>::iterator mismatch_key = mismatches_list.begin(); mismatch_key < mismatches_list.end(); mismatch_key++)  {
		if (*mismatch_key == target)  {
			found = true;
			break;
		}
	}
	if (!found)  {
	    mismatches_list.push_back(target);
	}
}

