import pandas as pd
import numpy as np
import random
import uuid
from typing import List, Optional
from rdkit import Chem # <--- Added to validate molecules
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
        Reads the inventory CSV and categorizes chemicals.
        Filters out invalid SMILES automatically to prevent crashes.
        """
        try:
            df = pd.read_csv(self.purchasing.csv_path)
            
            # 1. Clean Headers
            df.columns = df.columns.str.strip()
            
            # 2. Clean Boolean Column
            if df['purchasable'].dtype == 'object':
                df['purchasable'] = df['purchasable'].astype(str).str.strip().str.upper() == 'TRUE'
            
            # 3. Filter Purchasable
            df = df[df['purchasable'] == True]
            
            valid_count = 0
            skipped_count = 0

            for _, row in df.iterrows():
                notes = str(row.get('notes', '')).lower()
                name = row['name']
                smiles = str(row['smiles']).strip()
                
                # --- VALIDATION CHECK ---
                # Skip explicit placeholders
                if "INVALID" in smiles or "FALSE" in smiles or not smiles:
                    skipped_count += 1
                    continue

                # Check if RDKit can parse it
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    # Log only the first few errors to avoid spamming
                    if skipped_count < 3:
                        print(f"Warning: Skipping invalid SMILES: {smiles}")
                    elif skipped_count == 3:
                        print("Warning: ... and more invalid SMILES skipped.")
                    skipped_count += 1
                    continue
                # ------------------------
                
                item = {
                    'name': name, 
                    'smiles': smiles
                }

                if 'salt' in notes:
                    self.salts.append(item)
                elif 'solvent' in notes or 'ether' in notes or 'carbonate' in notes:
                    if 'additive' in notes:
                        self.additives.append(item)
                    else:
                        self.solvents.append(item)
                elif 'additive' in notes:
                    self.additives.append(item)
                else:
                    pass
            
            print(f"Debug: Loaded {len(self.salts)} salts, {len(self.solvents)} solvents, {len(self.additives)} additives. (Skipped {skipped_count} invalid entries)")

        except Exception as e:
            print(f"Error loading inventory for generation: {e}")
            self.salts, self.solvents, self.additives = [], [], []

    def generate_random_pool(self, n_candidates: int) -> List[Formulation]:
        """
        Generates N random, valid, and safe formulations.
        """
        candidates = []
        attempts = 0
        max_attempts = n_candidates * 10 
        
        # Check if we have enough ingredients to work with
        if not self.salts:
            print("Error: No valid salts found in inventory.")
            return []
        if not self.solvents:
            print("Error: No valid solvents found in inventory.")
            return []
        
        while len(candidates) < n_candidates and attempts < max_attempts:
            attempts += 1
            
            try:
                # 1. Select Components
                salt = random.choice(self.salts)
                
                n_solvents = random.randint(1, 3)
                # Ensure we don't try to sample more solvents than exist
                k_solv = min(n_solvents, len(self.solvents))
                solvents = random.sample(self.solvents, k=k_solv)
                
                n_additives = random.randint(0, 2)
                additives = []
                if self.additives:
                    k_add = min(n_additives, len(self.additives))
                    additives = random.sample(self.additives, k=k_add)

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
                    cycles_to_80=0.0, 
                    censored=False,
                    notes="Generated Candidate"
                )

                # 4. Safety Check
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