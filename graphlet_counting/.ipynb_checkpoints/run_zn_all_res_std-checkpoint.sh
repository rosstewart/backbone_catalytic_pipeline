#!/bin/bash
#cd 
DATASET=../../xin_lit_backed/xin_cpps_supplementary/all_graphs_single_chain
RESULTS_DIR=/data/ross/phos_dna_zn/zn_all_res_graphlet_std_results

echo
echo Generating kernel matrices for each kernel method ...

for file in ${DATASET}/*.zn_pos
do
    filename=$(basename -- "$file")
    id="${filename%.*}"
    #id=$(basename "$file" .phos_pos)
    #echo ${id}

    GRAPH_FILE=${DATASET}/${id}.graph
    LABELS_FILE=${DATASET}/${id}.labels

    #echo ${GRAPH_FILE}
    #echo "./run_kernel -p ${DATASET}/${id}.phos_pos -n ${DATASET}/${id}.phos_neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${id}_sgk.svml -N true -v"
    #./run_kernel -p ${DATASET}/${id}.zn_pos -n ${DATASET}/${id}.zn_neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${id}_sgk.svml -N true
    ./run_kernel -p ${DATASET}/${id}.all_res -n ${DATASET}/${id}.all_res_neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${id}_sgk.svml -N true

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

