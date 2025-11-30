import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import torch
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem

import torch
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from torch_geometric.utils import to_networkx
from chem_repr.graph_utils import smile_to_graph

def compare_molecule_to_graph(smiles, title="Molecule"):
    """
    Plots the RDKit image (Left) vs the NetworkX Graph (Right).
    """
    # 1. Get the Chemistry Representation (RDKit)
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print(f"Invalid SMILES: {smiles}")
        return

    # Add index labels to the RDKit image so we can match them to the graph
    for idx, atom in enumerate(mol.GetAtoms()):
        atom.SetProp('molAtomMapNumber', str(idx))
    
    # 2. Get the Graph Representation (PyG -> NetworkX)
    data = smile_to_graph(smiles)
    G = to_networkx(data, to_undirected=True)

    # --- PLOTTING ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # LEFT: The Chemist's View
    img = Draw.MolToImage(mol)
    axes[0].imshow(img)
    axes[0].set_title(f"{title}\n(Chemist View)", fontsize=14)
    axes[0].axis('off')

    # RIGHT: The Computer's View (Graph)
    # We will color nodes by atomic number to make it readable
    atomic_nums = [mol.GetAtomWithIdx(i).GetAtomicNum() for i in range(mol.GetNumAtoms())]
    # Simple color map: C=grey, O=red, F=green, Li=purple, others=blue
    colors = []
    labels = {}
    for i, atomic_num in enumerate(atomic_nums):
        labels[i] = f"{i}\n{Chem.PeriodicTable.GetElementSymbol(Chem.GetPeriodicTable(), atomic_num)}"
        if atomic_num == 6: colors.append('#808080')   # Carbon
        elif atomic_num == 8: colors.append('#FF0000') # Oxygen
        elif atomic_num == 9: colors.append('#00FF00') # Fluorine
        elif atomic_num == 3: colors.append('#800080') # Lithium
        elif atomic_num == 15: colors.append('#FFA500')# Phosphorus
        else: colors.append('#0000FF')

    pos = nx.kamada_kawai_layout(G) # Force-directed layout
    
    nx.draw(G, pos, ax=axes[1], 
            node_color=colors, 
            with_labels=True, 
            labels=labels,
            node_size=800, 
            font_color="white",
            font_weight="bold")
    
    axes[1].set_title(f"Graph Topology\n(Model Input: edge_index)", fontsize=14)
    axes[1].axis('off')

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    # Example 1: DME (Solvent)
    print("Visualizing DME...")
    compare_molecule_to_graph("COCCOC", "DME (Dimethoxyethane)")

    # Example 2: LiFSI (Salt)
    # Note: RDKit treats salts usually as disconnected, but let's see how the graph handles it.
    print("Visualizing LiFSI...")
    compare_molecule_to_graph("[Li+].F[S](=O)(=O)[N-][S](=O)(=O)F", "LiFSI")

    print("Visualizing TTE...")
    compare_molecule_to_graph("C(C(C(F)F)(F)F)OC(C(F)F)(F)F", "TTE")

    print("Visualizing perfluoro-THF...")
    compare_molecule_to_graph("C(COCC(F)F)OCC(F)F", "perfluoro-THF")