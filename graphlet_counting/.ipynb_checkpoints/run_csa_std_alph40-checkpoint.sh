#!/bin/bash
DATASET=../../mat_graphs
#GRAPH_FILE=../../mat_graphs/1arz_A.graph
#LABELS_FILE=../../mat_graphs/1arz_A.labels
RESULTS_DIR=../../graphlet_std_results_alph40

echo
echo Generating kernel matrices for each kernel method ...

for file in ${DATASET}/*.graph
do
    filename=$(basename -- "$file")
    id="${filename%.*}"

    GRAPH_FILE=${DATASET}/${id}.graph
    LABELS_FILE=${DATASET}/${id}.cons_labels

    ./run_kernel -p ${DATASET}/${id}.pos -n ${DATASET}/${id}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -A AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy -s ${RESULTS_DIR}/${id}_sgk.svml -N true -v

    #echo "File: $file, ID: $id, done"
done

#Running Cumulative Random Walk Kernel
#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 0 -I 100000 -R 0.4 -k ${RESULTS_DIR}/${DATASET}_crw.dat -v 

#Running Random Walk Kernel
#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 1 -I 100000 -R 0.2 -k ${RESULTS_DIR}/${DATASET}_rw.dat -v

#Running Standard Graphlet Kernel Kernel
#./run_kernel -p ../../${DATASET}/1arz_A.pos -n ../../${DATASET}/1arz_A.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -k ${RESULTS_DIR}/${DATASET}_sgk.dat -v

#./run_kernel -p ../../${DATASET}/1arz_A.pos -n ../../${DATASET}/1arz_A.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${DATASET}_sgk.svml -v

#Running Label Substitutions Graphlet Kernel
#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 3 -M 0.34 -A ACDEFGHIKLMNPQRSTVWY -k ${RESULTS_DIR}/${DATASET}_lsk.dat -v

#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 3 -M 0.34 -A ACDEFGHIKLMNPQRSTVWY -s ${RESULTS_DIR}/${DATASET}_lsk.svml -v

#Running Edge Indels Graphlet Kernel
#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 4 -E 1 -k ${RESULTS_DIR}/${DATASET}_eik.dat -v

#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 4 -E 1 -s ${RESULTS_DIR}/${DATASET}_eik.svml -v

#Running Edit Distance Graphlet Kernel
#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 5 -E 1 -M 0.34 -A ACDEFGHIKLMNPQRSTVWY -k ${RESULTS_DIR}/${DATASET}_edk.dat -v

#./run_kernel -p ${DATASET}/${DATASET}.pos -n ${DATASET}/${DATASET}.neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 5 -E 1 -M 0.34 -A ACDEFGHIKLMNPQRSTVWY -s ${RESULTS_DIR}/${DATASET}_edk.svml -v

echo

date

