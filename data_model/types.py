from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Any
import pandas as pd
import json

@dataclass
class Component:
    """
    Represents a single salt or solvent.
    """
    role: str # 'salt' or 'solvent'
    name: str
    smiles: str
    amount: float # concentration
    amount_type: str # 'molar', 'volume_fraction', etc.

@dataclass
class Formulation:
    """
    Represents a complete electrolyte mixture.
    CORE INTERNAL: Do not modify structure casually.
    """
    formulation_id: str
    components: List[Component]
    cycles_to_80: float
    censored: bool
    notes: Optional[str] = None

    @classmethod
    def from_row(cls, row: pd.Series) -> 'Formulation':
        """Parses a CSV row into a Formulation object."""
        try:
            comps_data = json.loads(row['components'])
            components = [Component(**c) for c in comps_data]
            
            # Handle potential missing/NaN values
            cycles = float(row['cycles_to_80'])
            censored = bool(row['censored'])
            
            return cls(
                formulation_id=str(row['formulation_id']),
                components=components,
                cycles_to_80=cycles,
                censored=censored,
                notes=row.get('notes', "")
            )
        except Exception as e:
            raise ValueError(f"Error parsing row {row.name}: {e}")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)