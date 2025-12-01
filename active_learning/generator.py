import pandas as pd
import numpy as np
import random
import uuid
from typing import List, Optional
from rdkit import Chem
from rdkit import RDLogger  # <--- Import Logger
from data_model.types import Formulation, Component

class CandidateGenerator:
    def __init__(self, purchasability_manager, safety_filter):
        self.purchasing = purchasability_manager
        self.safety = safety_filter
        
        # --- SILENCE RDKIT ERRORS ---
        # This stops the "SMILES Parse Error" blocks from printing to console
        RDLogger.DisableLog('rdApp.*') 
        
        # Load and categorize inventory
        self.salts = []
        self.solvents = []
        self.additives = []
        self._load_and_categorize_inventory()

    def _load_and_categorize_inventory(self):
        try:
            df = pd.read_csv(self.purchasing.csv_path)
            
            # Clean Headers
            df.columns = df.columns.str.strip()
            
            # Clean Boolean
            if df['purchasable'].dtype == 'object':
                df['purchasable'] = df['purchasable'].astype(str).str.strip().str.upper() == 'TRUE'
            
            # Filter
            df = df[df['purchasable'] == True]
            
            for _, row in df.iterrows():
                notes = str(row.get('notes', '')).lower()
                name = row['name']
                smiles = str(row['smiles']).strip()
                
                # Validation
                if "INVALID" in smiles or "FALSE" in smiles or not smiles:
                    continue

                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    # RDKit error logs are now silenced, so this will skip quietly
                    continue
                
                item = {'name': name, 'smiles': smiles}

                if 'salt' in notes:
                    self.salts.append(item)
                elif 'solvent' in notes or 'ether' in notes or 'carbonate' in notes:
                    if 'additive' in notes:
                        self.additives.append(item)
                    else:
                        self.solvents.append(item)
                elif 'additive' in notes:
                    self.additives.append(item)
            
            print(f"Debug: Loaded {len(self.salts)} salts, {len(self.solvents)} solvents, {len(self.additives)} additives.")

        except Exception as e:
            print(f"Error loading inventory: {e}")
            self.salts, self.solvents, self.additives = [], [], []

    def generate_random_pool(self, n_candidates: int) -> List[Formulation]:
        candidates = []
        attempts = 0
        max_attempts = n_candidates * 10 
        
        if not self.salts or not self.solvents:
            print("Error: No valid ingredients found.")
            return []
        
        while len(candidates) < n_candidates and attempts < max_attempts:
            attempts += 1
            try:
                salt = random.choice(self.salts)
                n_solvents = random.randint(1, 3)
                solvents = random.sample(self.solvents, k=min(n_solvents, len(self.solvents)))
                
                components = []
                components.append(Component('salt', salt['name'], salt['smiles'], round(random.uniform(0.5, 2.0), 2), 'molar'))
                
                raw_vols = [random.uniform(0.2, 1.0) for _ in solvents]
                total_vol = sum(raw_vols)
                for s, v in zip(solvents, raw_vols):
                    components.append(Component('solvent', s['name'], s['smiles'], round(v / total_vol, 2), 'volume_fraction'))
                
                formulation = Formulation(f"gen_{uuid.uuid4().hex[:8]}", components, 0.0, False, "Generated")
                
                if not any(self.safety.is_forbidden(c.smiles) for c in components):
                    candidates.append(formulation)
                    
            except Exception:
                continue

        print(f"Generated {len(candidates)} unique candidates.")
        return candidates