import torch
import torch.nn as nn

class FormulationNet(nn.Module):
    """
    Aggregates component embeddings into a cycle life prediction.
    Uses a DeepSets-like approach (sum/mean pooling + MLP).
    """
    def __init__(self, mol_embed_dim: int, hidden_dim: int):
        super().__init__()
        # Input: Mol embedding + Concentration (1 float) + Role (1 float encoding)
        self.component_process = nn.Linear(mol_embed_dim + 2, hidden_dim)
        
        self.predictor = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1) # Predicts log(cycles) or raw cycles
        )

    def forward(self, mol_embeddings, concentrations, roles, batch_indices):
        """
        mol_embeddings: (Total Components in Batch, Dim)
        concentrations: (Total Components in Batch, 1)
        roles: (Total Components in Batch, 1) -> 0 for solvent, 1 for salt
        batch_indices: (Total Components in Batch) -> indicates which formulation
        """
        # Concat component features
        features = torch.cat([mol_embeddings, concentrations, roles], dim=1)
        
        # Process individual components
        h = torch.relu(self.component_process(features))
        
        # Aggregation (Sum pooling per formulation)
        # In a loop for clarity, but strictly should be scatter_add for GPU speed
        # Using a simple scatter add implementation:
        out = torch.zeros(batch_indices.max() + 1, h.size(1), device=h.device)
        out.index_add_(0, batch_indices, h)
        
        # Predict
        prediction = self.predictor(out)
        return prediction.squeeze(-1)