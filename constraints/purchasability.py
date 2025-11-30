import pandas as pd
from typing import Set

class PurchasabilityManager:
    def __init__(self, csv_path: str):
        self.csv_path = csv_path
        self.purchasable_smiles = self._load()

    def _load(self) -> Set[str]:
        try:
            df = pd.read_csv(self.csv_path)
            # Filter for strict purchasability
            available = df[df['purchasable'] == True]
            return set(available['smiles'].tolist())
        except FileNotFoundError:
            print("Warning: Purchasability DB not found. Assuming empty.")
            return set()

    def is_purchasable(self, smiles: str) -> bool:
        return smiles in self.purchasable_smiles