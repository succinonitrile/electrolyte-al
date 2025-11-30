import torch
import numpy as np
import pandas as pd
from .acquisition import select_candidates
from .generator import CandidateGenerator
from models.surrogate import FormulationSurrogateModel
from constraints.filters import SafetyFilter
from constraints.purchasability import PurchasabilityManager

class ActiveLearner:
    """
    Orchestrates the loop: Train -> Generate -> Filter -> Score -> Recommend.
    """
    def __init__(self, config, device, data_path, purchasability_path):
        self.config = config
        self.device = device
        self.data_path = data_path
        
        # Constraints
        self.safety = SafetyFilter(config['constraints']['forbidden_smarts'])
        self.purchasing = PurchasabilityManager(purchasability_path)
        
        # Model
        self.surrogate = FormulationSurrogateModel(config, device)

    def load_data_and_train(self):
        # 1. Load CSV
        df = pd.read_csv(self.data_path)
        # TODO: Convert DF to PyG DataLoader (omitted for brevity, requires custom collate)
        # train_loader = create_loader(df, batch_size=self.config['training']['batch_size'])
        
        # 2. Train
        # self.surrogate.train_ensemble(train_loader)
        print("Simulated training complete.")

    def run_cycle(self):
        # 1. Generate Candidates
        generator = CandidateGenerator(self.purchasing, self.safety)
        # Generate raw candidate objects
        candidates = generator.generate_random_pool(self.config['active_learning']['candidate_pool_size'])
        
        # 2. Score Candidates (Batched Inference)
        # Convert candidates to PyG batch (omitted)
        # means, stds = self.surrogate.predict(candidate_batch)
        
        # Mocking predictions for code structure validity
        means = np.random.uniform(50, 600, len(candidates))
        stds = np.random.uniform(10, 100, len(candidates))
        
        # 3. Select
        results = select_candidates(
            candidates, means, stds, 
            n_high=self.config['active_learning']['n_high_performance'],
            n_low=self.config['active_learning']['n_low_performance'],
            n_unc=self.config['active_learning']['n_uncertainty'],
            beta=self.config['active_learning']['ucb_beta']
        )
        
        return results