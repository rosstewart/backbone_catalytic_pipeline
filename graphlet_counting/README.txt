Graphlet Kernel Framework version 1.0 (Adapted 
from Vladimir Vacic's graphlet kernel code)

Jose Lugo-Martinez, jlugomar@indiana.edu
Department of Computer Science
School of Informatics and Computing
Indiana University, Bloomington, IN, USA

This is version 1.0 of the graphlet kernels 
software for vertex classification. This 
version provides a framework for users to 
run a variety of graph-based kernel methods 
on (sparse) vertex-labeled undirected graphs. 
Kernels implemented are: cumulative random 
walk, standard random walk, standard graphlet 
kernel, edit distance graphlet kernel, label 
substitutions graphlet kernel and edge indels 
graphlet kernel.

Please direct all comments and bug reports
of this version to jlugomar@indiana.edu

A description of each method and discussion 
of parameters can be found in the paper. 

If you use this framework, please cite our paper:

Lugo-Martinez J. and Radivojac P. "Generalized
Graphlet Kernels for Probabilistic Inference 
in Sparse Graphs", Network Science, 2014.


--------------------------------------------
COMPILATION
--------------------------------------------

To compile the graphlet kernel framework, 
type "make" on the command prompt. Program 
file "run_kernel" will be generated. There 
is a configuration file where the user can 
set up some meta-parameters (graphlet 
combinations of interest, vertex labels 
alphabet size, etc.)

Running the binary with the -h switch will 
list all the command line options.


--------------------------------------------
PROGRAM OPTIONS
--------------------------------------------

Usage: run_kernel -p FILE -n FILE -g G_FILE -l L_FILE -t TYPE -[k|s] OUTPUT [...]
Options:

  -h         Displays this message.

  -t TYPE    Kernel type (0-Cumulative Random Walk, 1-Random Walk, 2-Standard Graphlet, 3-Label Substitutions Graphlet, 4-Edge Indels Graphlet, 5-Edit Distance Graphlet)
             Defaults to standard graphlet.

  -p FILE    List of positive (vertices) examples.
  -n FILE    List of negative (vertices) examples.
  -g G_FILE  Input graph file.
  -l L_FILE  Vertex labels file for input graph.

  -N         Normalize the kernel matrix.
             Defaults to false.

  -k KERNEL  Output file for the kernel matrix in standard output.
  -s SPARSE  Output file for the sparse attribute matrix (SVML).
             Defaults to KERNEL.

  -I STEPS   Number of steps. (Needed for Random Walk Kernels)
             Defaults to 100,000 steps.

  -R RESTART Restart probability. (Needed for Random Walk Kernels)
             Defaults to 0.3

  -S SIMMAT  Similarity matrix for weighting each possible vertex label substitution. (Needed for Label Substitutions and Edit Distance Kernels)
             Defaults to uniform weights (i.e. all label mismatches are equally weighted to 1).

  -M LABMIS  Fraction of nodes in the n-graphlet allowed to have vertex label mismatches. (Needed for Label Substitutions and Edit Distance Kernels)
             For example, M=0.34 allows 1-vertex label mismatch for n=3,4 and 5 graphlets, whereas M=0.26 allows 1-vertex label mismatch for n=4 and 5 graphlets.
             Defaults to 0.0 (i.e. no vertex label mismatches allowed).

  -E EDGMIS  Total number of edge insertions and/or deletions allowed between graphlets. (Needed for Edge Indels and Edit Distance Kernels)
             Defaults to 1.

  -A ALPHA   Vertex labels alphabet over problem statement is defined (i.e. all possible labels for a vertex). (Needed for Label Substitutions and Edit Distance Kernels)

  -c LABELS  Output file for each example class label.

  -v         Verbose (prints progress messages).


--------------------------------------------
EXAMPLE
--------------------------------------------

The example subdirectory contains two files with 
positive and negative vertices of interest within 
a vertex-labeled undirected graph, namely 'example.pos' 
and 'example.neg'. These files contain only one field, 
which indicates the vertex of interest within the 
input labeled undirected graph.

Preparing a kernel matrix from graph file is a one 
step process:

run_kernel is used to generate the user-specified 
graph-based kernel method. There are several output 
options, out of which SVM^Light format is the easiest 
to use, because it can readily be used with SVM^Light
(see http://svmlight.joachims.org). In a nutshell, 
this is a space-separated file with the first field 
equal to 1 (for positives) or -1 (for negatives), and 
all other non-zero entries are feature_key:feature_value
pairs.  

A sample "run_example.sh" shell script is provided 
as an illustrative example of the process and command 
line parameters required to run each kernel method.


--------------------------------------------
COMMENT REGARDING LARGE DATASETS
--------------------------------------------

Using the SVM^Light format directly on very large 
datasets is not efficient, because the dot-product 
between feature vectors has to be unnecessarily 
recomputed over and over again. There is an option 
to precompute the kernel matrix and save it as a 
binary file (or standard file), however, using it 
with SVM^Light is slightly more complicated, because 
the reader for this file format has to be custom 
coded and under the SVM^Light copyright agreement 
we are not allowed to distribute the modified version
of the SVM^Light code. 
