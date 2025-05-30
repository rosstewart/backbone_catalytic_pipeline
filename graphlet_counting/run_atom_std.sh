#!/bin/bash

pdb_id=$1
DATASET=$2
RESULTS_DIR=$2

run_job() {
    id=$1

    GRAPH_FILE="${DATASET}/${id}.graph"
    LABELS_FILE="${DATASET}/${id}.atom_labels"
    OUTPUT_FILE="${RESULTS_DIR}/${id}_sgk.svml"

    ./run_kernel -p "${DATASET}/${id}.atom_pos" \
                 -n "${DATASET}/${id}.atom_neg" \
                 -g "${GRAPH_FILE}" \
                 -l "${LABELS_FILE}" \
                 -t 2 \
                 -A "NAC" \
                 -s "${OUTPUT_FILE}" \
                 -N true \
                 -v
}

run_job "${pdb_id}"
