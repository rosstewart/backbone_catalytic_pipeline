#!/bin/bash
#cd 

method=interaction_loss

DATASET=../../ppi_lossgain/af3_graphs
RESULTS_DIR=/data/ross/ppi_lossgain/interaction_loss/${method}_4graphlet_sep_lsk_results

mkdir $RESULTS_DIR

echo
echo Generating kernel matrices for each kernel method ...

for file in ${DATASET}/*.${method}_{pos,neg}
do
    filename=$(basename -- "$file")
    id="${filename%.*}"
    #id=$(basename "$file" .phos_pos)
    #echo ${id}

    GRAPH_FILE=${DATASET}/${id}.graph
    WT_LABELS_FILE=${DATASET}/${id}.labels_separated

    # wild-type
    if [[ "$file" == *_pos* ]]; then
        if [[ -e "${file/_pos/_neg}" ]]; then
            ./run_kernel -p ${DATASET}/${id}.${method}_pos -n ${DATASET}/${id}.${method}_neg -g ${GRAPH_FILE} -l ${WT_LABELS_FILE} -t 3 -M 0.34 -A AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy -s ${RESULTS_DIR}/${id}_sgk.svml
        else
            touch ${DATASET}/${id}.${method}_filler
            #echo "./run_kernel -p ${DATASET}/${id}.${method}_pos -n ${DATASET}/${id}.${method}_filler -g ${GRAPH_FILE} -l ${WT_LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${id}_sgk.svml"
            ./run_kernel -p ${DATASET}/${id}.${method}_pos -n ${DATASET}/${id}.${method}_filler -g ${GRAPH_FILE} -l ${WT_LABELS_FILE} -t 3 -M 0.34 -A AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy -s ${RESULTS_DIR}/${id}_sgk.svml
        fi
    else
        if [[ -e "${file/_neg/_pos}" ]]; then
            ./run_kernel -p ${DATASET}/${id}.${method}_pos -n ${DATASET}/${id}.${method}_neg -g ${GRAPH_FILE} -l ${WT_LABELS_FILE} -t 3 -M 0.34 -A AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy -s ${RESULTS_DIR}/${id}_sgk.svml
        else
            touch ${DATASET}/${id}.${method}_filler
            ./run_kernel -p ${DATASET}/${id}.${method}_filler -n ${DATASET}/${id}.${method}_neg -g ${GRAPH_FILE} -l ${WT_LABELS_FILE} -t 3 -M 0.34 -A AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy -s ${RESULTS_DIR}/${id}_sgk.svml
        fi
    fi

    #echo "wt done"

    # read all variants
    while read -r mt_idx variant; do
        VT_LABELS_FILE="${DATASET}/${id}_${method}_variant_${variant}.labels_separated"
        #echo $VT_LABELS_FILE

        if [[ ! -f "$VT_LABELS_FILE" ]]; then
            echo "Error: File $VT_LABELS_FILE does not exist."
            exit 1
        fi

        temp_site_file=$(mktemp)
        trap 'rm -f "$temp_site_file"' EXIT

        echo "${mt_idx}" > "$temp_site_file"

        if [[ "$file" == *_pos* ]]; then
            touch ${DATASET}/${id}.${method}_filler
            ./run_kernel -p $temp_site_file -n ${DATASET}/${id}.${method}_filler -g ${GRAPH_FILE} -l ${VT_LABELS_FILE} -t 3 -M 0.34 -A AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy -s ${RESULTS_DIR}/${id}_variant_${variant}_sgk.svml
        else
            touch ${DATASET}/${id}.${method}_filler
            ./run_kernel -p ${DATASET}/${id}.${method}_filler -n $temp_site_file -g ${GRAPH_FILE} -l ${VT_LABELS_FILE} -t 3 -M 0.34 -A AaCcDdEeFfGgHhIiKkLlMmNnPpQqRrSsTtVvWwYy -s ${RESULTS_DIR}/${id}_variant_${variant}_sgk.svml
        fi

    done < "$file"

    # variants
#    for VT_LABELS_FILE in ${DATASET}/${id}_${method}_variant_*.labels
#    do
#        variant=$(basename "$VT_LABELS_FILE" | cut -d'.' -f1 | awk -F'_' '{print $NF}')
#
#        if [[ "$file" == *_pos* ]]; then
#            touch ${DATASET}/${id}.${method}_filler
#            ./run_kernel -p ${DATASET}/${id}.${method}_pos -n ${DATASET}/${id}.${method}_filler -g ${GRAPH_FILE} -l ${VT_LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${id}_variant_${variant}_sgk.svml -N false
#        else
#            touch ${DATASET}/${id}.${method}_filler
#            ./run_kernel -p ${DATASET}/${id}.${method}_filler -n ${DATASET}/${id}.${method}_neg -g ${GRAPH_FILE} -l ${VT_LABELS_FILE} -t 2 -s ${RESULTS_DIR}/${id}_variant_${variant}_sgk.svml -N false 
#        fi
# 
#
#    done

done

echo

date

