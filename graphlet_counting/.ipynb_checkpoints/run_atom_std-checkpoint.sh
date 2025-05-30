#!/bin/bash
#cd 
DATASET=../../xin_lit_backed/xin_cpps_supplementary/atom_graphlets/atom_graphs

method=$1
echo $method
# for method in cat phos dna zn; do
RESULTS_DIR=/data/ross/atom_graphlets/${method}_graphlet_std_results

echo
echo Generating kernel matrices for $method kernel method ...

for file in ${DATASET}/*.${method}_atom_pos
do
    filename=$(basename -- "$file")
    id="${filename%.*}"
    #id=$(basename "$file" .phos_pos)
    #echo ${id}

    GRAPH_FILE=${DATASET}/${id}.graph
    LABELS_FILE=${DATASET}/${id}.atom_labels

    #echo ${GRAPH_FILE}
    #echo "./run_kernel -p ${DATASET}/${id}.phos_pos -n ${DATASET}/${id}.phos_neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${id}_sgk.svml -N true -v"
    # echo ${RESULTS_DIR}/${id}_sgk.svml
    ./run_kernel -p ${DATASET}/${id}.${method}_atom_pos -n ${DATASET}/${id}.${method}_atom_neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 2 -A "ABCDEGHMNOPQSTXYZh?" -s ${RESULTS_DIR}/${id}_sgk.svml -N true -v
    # ./run_kernel -p ${DATASET}/${id}.${method}_atom_pos -n ${DATASET}/${id}.${method}_atom_neg -g ${GRAPH_FILE} -l ${LABELS_FILE} -t 3 -M 0.00 -A "ABCDEGHMNOPQSTXYZh?" -s ${RESULTS_DIR}/${id}_sgk.svml -N true -v

    #echo "File: $file, ID: $id, done"
done
# done


echo

date

