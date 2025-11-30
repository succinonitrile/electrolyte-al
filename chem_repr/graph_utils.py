import torch
from rdkit import Chem
from torch_geometric.data import Data

# SAFE TO MODIFY: You can add more atomic features here.
def smile_to_graph(smile: str) -> Data:
    """
    Converts SMILES to PyG Data object.
    """
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smile}")
        
    # Atoms
    # Feature vector: [AtomicNum, Degree, IsAromatic]
    # Keep features minimal for desktop memory efficiency
    x = []
    for atom in mol.GetAtoms():
        x.append([
            atom.GetAtomicNum(),
            atom.GetTotalDegree(),
            1 if atom.GetIsAromatic() else 0
        ])
    x = torch.tensor(x, dtype=torch.float)

    # Bonds
    edge_index = []
    edge_attr = []
    for bond in mol.GetBonds():
        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        edge_index.append([start, end])
        edge_index.append([end, start])
        
        # Feature: [BondTypeAsDouble]
        b_type = bond.GetBondTypeAsDouble()
        edge_attr.append([b_type])
        edge_attr.append([b_type])

    if len(edge_index) == 0:
        edge_index = torch.empty((2, 0), dtype=torch.long)
        edge_attr = torch.empty((0, 1), dtype=torch.float)
    else:
        edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
        edge_attr = torch.tensor(edge_attr, dtype=torch.float)

    return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)