import os
import glob
import pandas as pd
import numpy as np
import sys
from scipy.sparse import coo_matrix, csr_matrix, save_npz

import numpy as np
from scipy.sparse import csr_matrix, load_npz, coo_matrix
from sklearn.svm import LinearSVC
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_curve, roc_auc_score, precision_recall_curve, average_precision_score
import matplotlib.pyplot as plt
import sys
import os, glob
from joblib import Parallel, delayed
import pickle

def svml_to_sparse(save_dir):

    svml_results_dir = save_dir
    results_dir = save_dir
    
    pdb_ids_out_f = f"{results_dir}/pdb_ids.npy" # pdb id corresponding to every data point, useful if you are stratifying by protein. leave empty if you don't wish to save
    features_out_f = f"{results_dir}/features.npz" # feature matrix file name
    labels_out_f = f"{results_dir}/labels.npy" # label file name
    os.makedirs(results_dir,exist_ok=True)
    
    
    '''
    read graphlet features from each `.svml` file
    '''
    
    df = pd.DataFrame()
    
    residue_idx = 0 # the row of the sparse matrix
    labels = []
    pdb_ids = []
    sparse_feats = {}
    for svml in glob.glob(f'{svml_results_dir}/*.svml'):
        pdb_id = svml.split('/')[-1].split('.')[0]
        with open(svml,'r') as f:
            for line in f:
                vals = line.strip().split()
                label = int(vals[0])
                idx = int(vals[-1].replace('#',''))
                feats_and_counts  = vals[1:-1]
                for feat_and_count in feats_and_counts:
                    feat,count = feat_and_count.split(':')
                    # assert int(feat) not in sparse_feats
                    sparse_feats[(residue_idx,int(feat))] = float(count)
                residue_idx += 1
                labels.append(label)
                pdb_ids.append(pdb_id)
    
    if len(pdb_ids_out_f) != 0:
        np.save(pdb_ids_out_f, pdb_ids)
        # print(f'pdb ids saved in {pdb_ids_out_f} with shape {np.array(pdb_ids).shape}')
    
    np.save(labels_out_f, labels)
    # print(f'labels saved in {labels_out_f} with shape {np.array(labels).shape}')
    
    
    '''
    convert into `scipy.sparse` format
    '''
    
    rows,cols,vals = [],[],[]
    for (row,col), val in sparse_feats.items():
        rows.append(row)
        cols.append(col)
        vals.append(val)
    
    '''
    map feature indices to start at zero as the integers overflow (are too large)
    '''
    
    unique_str_feat = set()
    for col in cols:
        unique_str_feat.add(str(col))
    
    feature_mapping = {old_index: new_index for new_index, old_index in enumerate(unique_str_feat)}
    new_indices = np.array([feature_mapping[str(old_index)] for old_index in cols])
    
    n_rows = max(rows) + 1
    n_cols = len(unique_str_feat)
    
    sparse_mat = coo_matrix((vals, (rows, new_indices)), shape=(n_rows, n_cols)).tocsr()
    
    
    '''
    save feature file
    '''
    
    save_npz(features_out_f.replace('.npz','_formatted.npz'), sparse_mat)
    # print(f'sparse matrix saved as {features_out_f} with shape {sparse_mat.shape}')



def run_svm_inference(pdb_id, save_dir, model_dir):
    sparse_dir = save_dir
    graph_dir = save_dir

    X = load_npz(f'{save_dir}/features_formatted.npz')
    y = np.load(f'{save_dir}/labels.npy')
    
    # per-atom pdb_id
    pdb_ids = np.array([pdb_id.replace('_sgk','') for pdb_id in np.load(f'{save_dir}/pdb_ids.npy')])
    
    
    pos_atom_indices = []
    neg_atom_indices = []
    with open(f'{graph_dir}/{pdb_id}.atom_pos','r') as f:
        for line in f:
            pos_atom_indices.append(int(line.strip()))
    # THIS BLOCK IS FOR LABELED DATA (TRAINING)
    # with open(f'{graph_dir}/{pdb_id}.atom_neg','r') as f:
    #     for line in f:
    #         neg_atom_indices.append(int(line.strip()))
    with open(f'{graph_dir}/{pdb_id}.atom_res_indices','r') as f:
        for line in f:
            atom_res_indices = np.array(line.strip().split(',')).astype(int)
            break
    
    splits_res_indices = []
    for i in range(len(y)):
        if pdb_id == pdb_ids[i]:
            # get the residue idx corresponding to the atom
            # all the positives are first, then all the negatives
            if len(pos_atom_indices) != 0:
                assert y[i] == 1
                atom_idx = pos_atom_indices.pop(0)
            else:
                assert len(neg_atom_indices) != 0
                assert y[i] == -1
                atom_idx = neg_atom_indices.pop(0)
                
            res_idx = atom_res_indices[atom_idx]
            splits_res_indices.append(f'{pdb_id}_{res_idx}')
    
    assert len(pos_atom_indices) == 0 and len(neg_atom_indices) == 0
    
    # with open(f'{model_dir}/svm_model.pkl','rb') as f:
    #     svm = pickle.load(f)
    with open(f'{model_dir}/svm_model.pkl', 'rb') as f:
        calibrated_svm = pickle.load(f)

    # posteriors of catalytic residues
    y_pred = calibrated_svm.predict_proba(X)[:, 1]
    
    assert len(splits_res_indices) == len(y_pred)
    
    prev_res_idx = ''
    res_avg_preds = []
    all_res_preds, all_res_labels = [],[]
    for i in range(len(splits_res_indices)):
        res_idx = splits_res_indices[i]
        pred = y_pred[i]
        label = y[i]
    
        if res_idx != prev_res_idx:
            # avg predictions for prev residue
            if len(res_avg_preds) != 0:
                avg_pred = np.mean(res_avg_preds)
                all_res_preds.append(avg_pred)
                all_res_labels.append(y[i-1])
    
            # start new residue
            res_avg_preds = [pred]
        else:
            res_avg_preds.append(pred)
            assert label == y[i-1]
        
        prev_res_idx = res_idx
    
    if len(res_avg_preds) != 0:
        avg_pred = np.mean(res_avg_preds)
        all_res_preds.append(avg_pred)
        all_res_labels.append(label) # label will hold last value of y
     
    
    all_res_preds = np.array(all_res_preds)
    all_res_labels = np.array(all_res_labels)
    
    assert np.all(all_res_labels == 1)
    
    return all_res_preds











