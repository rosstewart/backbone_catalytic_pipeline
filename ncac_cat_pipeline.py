'''
Ross Stewart 05/30/25

ncac_cat_pipeline

predicts catalytic residue probabilities from an input pdb file

args:
figure it out

output:
posteriors
'''


import numpy as np
from Bio import PDB
import tempfile
from contact_graph_utils import make_graph, write_graph_and_labels, make_pos_file
from svm_utils import svml_to_sparse, run_svm_inference
import subprocess
import os
import sys

if len(sys.argv) != 3:
    print(f'Usage: python {sys.argv[0]} <input.pdb> <output.npy>')
    sys.exit(1)


pdb_f = sys.argv[1]
assert os.path.exists(pdb_f)
preds_f = sys.argv[2]

three_letter_to_one = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
    "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
    "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
    "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O", "Asx": "B", "Glx": "Z",
    "Xaa": "X", "Ter": "*"
}

pdb_atom_mapping = {"N": "N", "CA": "A", "C": "C"}


parser = PDB.PDBParser(QUIET=True)
wd = os.getcwd()
graphlet_wd = './graphlet_counting'
model_dir = './svm'
edge_dist_threshold = 7.5

with tempfile.TemporaryDirectory() as save_dir:
    
    # generate contact graph
    pdb_id, mat_data = make_graph(pdb_f, edge_dist_threshold, parser, three_letter_to_one, pdb_atom_mapping)
    
    # write files needed downstream
    write_graph_and_labels(pdb_id, mat_data, save_dir)
    make_pos_file(pdb_id, save_dir)

    # generate graphlet features (7.5 Å edge dist threshold, only up to 4-graphlet for speed)
    os.chdir(graphlet_wd) # need to be in graphlet dir unless i want to refactor more stuff
    subprocess.run([f'./run_atom_std.sh', pdb_id, save_dir], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.chdir(wd)
    
    svml_to_sparse(save_dir)
    residue_cat_preds = run_svm_inference(pdb_id, save_dir, model_dir)

    np.save(preds_f, residue_cat_preds)

print(f'\ncatalytic residue predictions saved to {preds_f}\n')
