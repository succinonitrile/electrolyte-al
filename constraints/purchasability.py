import pandas as pd
from typing import Set

class PurchasabilityManager:
    def __init__(self, csv_path: str):
        self.csv_path = csv_path
        self.purchasable_smiles = self._load()

    def _load(self) -> Set[str]:
        try:
            df = pd.read_csv(self.csv_path)
            
            # Fix 1: Strip whitespace from column headers (e.g., " purchasable" -> "purchasable")
            df.columns = df.columns.str.strip()

            # Fix 2: Debugging check
            if 'purchasable' not in df.columns:
                print(f"ERROR: 'purchasable' column missing in {self.csv_path}")
                print(f"Found columns: {list(df.columns)}")
                
                # Check for the common 'entire line quoted' issue
                if len(df.columns) == 1:
                    print("HINT: It looks like your CSV has quotes wrapping the entire header row.")
                    print("      Please open the CSV in a text editor and remove outer quotes.")
                
                # Return empty to prevent immediate crash, or let it raise specific error
                raise KeyError(f"'purchasable' column not found. Columns are: {list(df.columns)}")

            # Filter for strict purchasability
            available = df[df['purchasable'] == True]
            return set(available['smiles'].tolist())
            
        except FileNotFoundError:
            print("Warning: Purchasability DB not found. Assuming empty.")
            return set()