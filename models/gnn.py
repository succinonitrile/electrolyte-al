import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool

class MoleculeEncoder(nn.Module):
    """
    GCN to embed single molecules.
    SAFE TO MODIFY: Hidden dims and layers.
    """
    def __init__(self, node_in_dim: int, hidden_dim: int, out_dim: int):
        super().__init__()
        self.conv1 = GCNConv(node_in_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.conv3 = GCNConv(hidden_dim, out_dim)
        self.relu = nn.ReLU()

    def forward(self, x, edge_index, batch):
        x = self.relu(self.conv1(x, edge_index))
        x = self.relu(self.conv2(x, edge_index))
        x = self.conv3(x, edge_index)
        # Global pooling to get one vector per molecule
        return global_mean_pool(x, batch)