import torch
import numpy as np
from typing import List, Tuple
from .gnn import MoleculeEncoder
from .aggregator import FormulationNet

class FormulationSurrogateModel:
    """
    Wrapper for an Ensemble of models.
    DO NOT MODIFY: The predict interface logic.
    """
    def __init__(self, config: dict, device: torch.device):
        self.config = config
        self.device = device
        self.models = []
        self._build_ensemble()

    def _build_ensemble(self):
        # Create N identical models
        for _ in range(self.config['model']['ensemble_size']):
            encoder = MoleculeEncoder(
                node_in_dim=3, # Matches smile_to_graph
                hidden_dim=self.config['model']['hidden_dim'],
                out_dim=self.config['model']['embedding_dim']
            )
            aggregator = FormulationNet(
                mol_embed_dim=self.config['model']['embedding_dim'],
                hidden_dim=self.config['model']['hidden_dim']
            )
            # Combine into a container
            model = torch.nn.ModuleDict({'encoder': encoder, 'agg': aggregator})
            model.to(self.device)
            self.models.append(model)

    def train_ensemble(self, train_loader):
        """
        Trains models sequentially to save VRAM on desktop.
        """
        for idx, model in enumerate(self.models):
            print(f"Training ensemble member {idx+1}/{len(self.models)}...")
            optimizer = torch.optim.Adam(model.parameters(), lr=self.config['training']['learning_rate'])
            criterion = torch.nn.MSELoss() # Or Survival Loss (TODO)
            
            model.train()
            for epoch in range(self.config['training']['epochs']):
                for batch in train_loader:
                    # Move batch to device
                    batch = batch.to(self.device)
                    optimizer.zero_grad()
                    
                    # Forward pass
                    mols_embed = model['encoder'](batch.x, batch.edge_index, batch.batch)
                    preds = model['agg'](mols_embed, batch.concs, batch.roles, batch.formulation_idx)
                    
                    # Simple Regression Mode (ignoring censoring for MVP)
                    # TODO: Implement Survival Loss using batch.censored
                    loss = criterion(preds, batch.y)
                    loss.backward()
                    optimizer.step()

    def predict(self, batch) -> Tuple[np.ndarray, np.ndarray]:
        """
        Returns (mean, std)
        """
        batch = batch.to(self.device)
        preds_list = []
        
        with torch.no_grad():
            for model in self.models:
                model.eval()
                mols_embed = model['encoder'](batch.x, batch.edge_index, batch.batch)
                pred = model['agg'](mols_embed, batch.concs, batch.roles, batch.formulation_idx)
                preds_list.append(pred.cpu().numpy())
        
        preds_stack = np.stack(preds_list) # (n_models, batch_size)
        return np.mean(preds_stack, axis=0), np.std(preds_stack, axis=0)