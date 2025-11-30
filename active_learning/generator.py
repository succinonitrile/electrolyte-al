import pandas as pd
import numpy as np
import random
import uuid
from typing import List, Optional
from data_model.types import Formulation, Component

class CandidateGenerator:
    def __init__(self, purchasability_manager, safety_filter):
        self.purchasing = purchasability_manager
        self.safety = safety_filter
        
        # Load and categorize inventory
        self.salts = []
        self.solvents = []
        self.additives = []
        self._load_and_categorize_inventory()

    def _load_and_categorize_inventory(self):
        """
        Reads the inventory CSV from the path stored in the purchasability manager
        and categorizes chemicals into salts, solvents, and additives based on notes.
        """
        try:
            # We access the CSV path stored in the manager
            df = pd.read_csv(self.purchasing.csv_path)
            
            # Filter for purchasable items only
            df = df[df['purchasable'] == True]
            
            for _, row in df.iterrows():
                # Basic heuristic to categorize components based on the 'notes' column
                notes = str(row.get('notes', '')).lower()
                name = row['name']
                smiles = row['smiles']
                
                item = {
                    'name': name, 
                    'smiles': smiles
                }

                if 'salt' in notes:
                    self.salts.append(item)
                elif 'solvent' in notes or 'ether' in notes or 'carbonate' in notes:
                    # Note: 'ether' or 'carbonate' often implies solvent if not explicitly additive
                    if 'additive' in notes:
                        self.additives.append(item)
                    else:
                        self.solvents.append(item)
                elif 'additive' in notes:
                    self.additives.append(item)
                else:
                    # Fallback/Default behavior (optional: skip or treat as solvent)
                    pass
                    
        except Exception as e:
            print(f"Error loading inventory for generation: {e}")
            self.salts, self.solvents, self.additives = [], [], []

    def generate_random_pool(self, n_candidates: int) -> List[Formulation]:
        """
        Generates N random, valid, and safe formulations.
        """
        candidates = []
        attempts = 0
        max_attempts = n_candidates * 10 # Avoid infinite loops
        
        while len(candidates) < n_candidates and attempts < max_attempts:
            attempts += 1
            
            try:
                # 1. Select Components
                # Strategy: 1 Salt + 1-3 Solvents + 0-2 Additives
                if not self.salts or not self.solvents:
                    break # Cannot generate without ingredients
                
                salt = random.choice(self.salts)
                
                n_solvents = random.randint(1, 3)
                solvents = random.sample(self.solvents, k=min(n_solvents, len(self.solvents)))
                
                n_additives = random.randint(0, 2)
                if self.additives:
                    additives = random.sample(self.additives, k=min(n_additives, len(self.additives)))
                else:
                    additives = []

                # 2. Assign Amounts
                components = []
                
                # Salt (Molar: 0.5 to 2.0 M)
                components.append(Component(
                    role='salt',
                    name=salt['name'],
                    smiles=salt['smiles'],
                    amount=round(random.uniform(0.5, 2.0), 2),
                    amount_type='molar'
                ))
                
                # Solvents (Volume Fraction: sum to 1.0)
                raw_vols = [random.uniform(0.2, 1.0) for _ in solvents]
                total_vol = sum(raw_vols)
                for s, v in zip(solvents, raw_vols):
                    components.append(Component(
                        role='solvent',
                        name=s['name'],
                        smiles=s['smiles'],
                        amount=round(v / total_vol, 2),
                        amount_type='volume_fraction'
                    ))
                    
                # Additives (Molar: 0.01 to 0.2 M)
                for a in additives:
                    components.append(Component(
                        role='additive',
                        name=a['name'],
                        smiles=a['smiles'],
                        amount=round(random.uniform(0.01, 0.2), 3),
                        amount_type='molar'
                    ))

                # 3. Create Formulation
                formulation = Formulation(
                    formulation_id=f"gen_{uuid.uuid4().hex[:8]}",
                    components=components,
                    cycles_to_80=0.0, # Placeholder for candidate
                    censored=False,
                    notes="Generated Candidate"
                )

                # 4. Safety Check
                # If any component is forbidden, discard the whole formulation
                is_safe = True
                for comp in components:
                    if self.safety.is_forbidden(comp.smiles):
                        is_safe = False
                        break
                
                if is_safe:
                    candidates.append(formulation)
                    
            except Exception as e:
                continue

        print(f"Generated {len(candidates)} unique candidates.")
        return candidates