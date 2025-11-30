from rdkit import Chem
from typing import List

class SafetyFilter:
    def __init__(self, forbidden_smarts: List[str]):
        self.patterns = [Chem.MolFromSmarts(s) for s in forbidden_smarts]
        
    def is_forbidden(self, smiles: str) -> bool:
        """Returns True if molecule contains a forbidden motif."""
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return True # Invalid molecules are unsafe
        
        for pat in self.patterns:
            if mol.HasSubstructMatch(pat):
                return True
        return False