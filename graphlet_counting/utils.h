/**
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

#ifndef UTILS_H
#define	UTILS_H

#include "config.h"
#include <vector>
using namespace std;


float compare_labels(char label1, char label2);

int randint( int max); //[0.00-1.00)

double randdouble();

unsigned long get_graphlet_length(unsigned long g_type);

void compare_two(char &a, char &b, char &a1, char &b1);

void compare_three(char &a, char &b, char &c, char &a1, char &b1, char &c1);

void compare_four(char &a, char &b, char &c, char &d, char &a1, char &b1, char &c1, char &d1);

unsigned int set_k(unsigned long g_type, float sf);

void insert_permutation(const Key &target, vector<Key> &mismatches_list);

#endif

