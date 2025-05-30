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

assert len(sys.argv) == 7, f"Usage: python {sys.argv[0]} <sparse_dir> <split_file> <method> <save_dir> <is_4graphlet> <graph_dir>\n<sparse_dir> must end in <edge_threshold>"

is_4graphlet = int(sys.argv[5])
graph_dir = sys.argv[6]

edge_threshold = sys.argv[1].split('sparse_')[-1].replace('/','')
# Load data
file_name = sys.argv[1]+'/features_formatted.npz'
if '.npz' in file_name:
    X = load_npz(file_name)
else:
    X = np.load(file_name)
y = np.load(sys.argv[1]+'/labels.npy')#'pca_test_labels.npy')
pdb_ids = np.array([pdb_id.replace('_sgk','') for pdb_id in np.load(sys.argv[1]+'/pdb_ids.npy')])

split_f = '/home/rcstewart/gnn/cv_splits/cat_splits.tsv'
# split_f = sys.argv[2]

# graph_dir = f'xin_lit_backed/xin_cpps_supplementary/atom_graphlets/atom_graphs_{edge_threshold}'
method = sys.argv[3]
save_dir = sys.argv[4]
os.makedirs(save_dir,exist_ok=True)

# splits by group (chain)
num_weird = 0
with open(split_f,'r') as f:
    splits = [([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[])]
    splits_res_indices = [([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[]),([],[])]
    # seen_1hio_D = False
    for line in f:
        # print(line)
        fold, train_test, pdb_id = line.strip().split('\t')
        
        pos_atom_indices = []
        neg_atom_indices = []
        with open(f'{graph_dir}/{pdb_id}.{method}_atom_pos','r') as f:
            for line in f:
                pos_atom_indices.append(int(line.strip()))
        with open(f'{graph_dir}/{pdb_id}.{method}_atom_neg','r') as f:
            for line in f:
                neg_atom_indices.append(int(line.strip()))
        with open(f'{graph_dir}/{pdb_id}.atom_res_indices','r') as f:
            for line in f:
                atom_res_indices = np.array(line.strip().split(',')).astype(int)
                break
        
        for i in range(len(y)):
            # print(pdb_id,pdb_ids[i])
            if pdb_id == pdb_ids[i]:# or pdb_ids[i] in pdb_id or pdb_ids[i] in pdb_id: # old: some cat ones don't have chain identifier, smh
                # if pdb_id == '1hio_A':
                #     seen_1hio_D = True
                # found entry with pdb id in X
                # print(fold, train_test, pdb_id,pdb_ids[i])
                splits[int(fold[-1])][0 if train_test == 'train' else 1].append(i)

                ''' get the residue idx corresponding to the atom '''
                # all the positives are first, then all the negatives
                if len(pos_atom_indices) != 0:
                    assert y[i] == 1
                    atom_idx = pos_atom_indices.pop(0)
                else:
                    assert len(neg_atom_indices) != 0
                    assert y[i] == -1
                    atom_idx = neg_atom_indices.pop(0)
                    
                res_idx = atom_res_indices[atom_idx]
                splits_res_indices[int(fold[-1])][0 if train_test == 'train' else 1].append(f'{pdb_id}_{res_idx}')

        # if train_test == 'test':
        if len(pos_atom_indices) == 0 and len(neg_atom_indices) == 0:#, f"fold {fold} {train_test} {pdb_id}: pos {pos_atom_indices} len {len(pos_atom_indices)}, neg {neg_atom_indices} len {len(neg_atom_indices)}, seen_1hio_D {seen_1hio_D}"
            # num_weird += 1
            pass
        else:
            num_weird += 1

print('num with no pos/neg site labels',num_weird)

for i in range(len(splits)):
    splits[i] = (np.array(splits[i][0]).astype(np.int32),np.array(splits[i][1]).astype('int32'))
    splits_res_indices[i] = (np.array(splits_res_indices[i][0]).astype(str),np.array(splits_res_indices[i][1]).astype(str))


# SVM with a linear kernel
svm = LinearSVC(max_iter=500000000, C=1, penalty='l2', dual=True)

# Function to process each fold
def process_fold(train_index, test_index, fold_splits_res_indices):
    X_train, X_test = X[train_index], X[test_index]
    y_train, y_test = y[train_index], y[test_index]
    
    print('fitting fold...',flush=True)
    svm.fit(X_train, y_train)
    print('getting fold preds...',flush=True)
    
    # Get decision function scores for the test set
    y_pred = svm.decision_function(X_test)
    test_fold_splits_res_indices = fold_splits_res_indices[1]
    assert len(test_fold_splits_res_indices) == len(y_pred)

    prev_res_idx = ''
    res_avg_preds = []
    all_res_preds, all_res_labels = [],[]
    for i in range(len(test_fold_splits_res_indices)):
        res_idx = test_fold_splits_res_indices[i]
        pred = y_pred[i]
        label = y_test[i]

        if res_idx != prev_res_idx:
            # avg predictions for prev residue
            if len(res_avg_preds) != 0:
                avg_pred = np.mean(res_avg_preds)
                all_res_preds.append(avg_pred)
                all_res_labels.append(y_test[i-1])

            # start new residue
            res_avg_preds = [pred]
        else:
            res_avg_preds.append(pred)
            assert label == y_test[i-1]
        
        prev_res_idx = res_idx
    
    if len(res_avg_preds) != 0:
        avg_pred = np.mean(res_avg_preds)
        all_res_preds.append(avg_pred)
        all_res_labels.append(label) # label will hold last value of y_test
    
    print('finished fold',flush=True)    

    all_res_preds = np.array(all_res_preds)
    all_res_labels = np.array(all_res_labels)

    return all_res_labels, all_res_preds

# Parallelize the cross-validation process
results = Parallel(n_jobs=-1)(
    delayed(process_fold)(train_index, test_index, splits_res_indices[fold]) for fold, (train_index, test_index) in enumerate(splits) #cv.split(X, y)
)


# Collect the true and predicted labels
y_true, y_pred = zip(*results)
y_true = np.concatenate(y_true)
y_pred = np.concatenate(y_pred)

# Calculate ROC and PR AUC
fpr, tpr, _ = roc_curve(y_true, y_pred)
roc_auc = roc_auc_score(y_true, y_pred)

precision, recall, _ = precision_recall_curve(y_true, y_pred)
pr_auc = average_precision_score(y_true, y_pred)

print(f"ROC AUC: {roc_auc:.4f}",flush=True)
print(f"PR AUC: {pr_auc:.4f}",flush=True)

graphlet_count = '4' if is_4graphlet else '5'
np.save(f'{save_dir}/{method}_{edge_threshold}_atom_avg_{graphlet_count}graphlet_svm_preds.npy',y_pred)
np.save(f'{save_dir}/{method}_{edge_threshold}_atom_avg_{graphlet_count}graphlet_svm_labels.npy',y_true)

with open(f'{save_dir}/{method}_{edge_threshold}_{graphlet_count}graphlet_svm.pkl','wb') as f:
    pickle.dump(svm, f)
