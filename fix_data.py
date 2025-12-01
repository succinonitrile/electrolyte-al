import os
import pandas as pd

# Mappings of Invalid SMILES -> Valid RDKit SMILES or "INVALID_ENTRY"
REPLACEMENTS = {
    # --- Salts ---
    "[Li+].[BF4-]": "[Li+].F[B-](F)(F)F",
    "[Li+].[NO3-]": "[Li+].[O-][N+](=O)[O-]",
    "[Li+].B(OC(CF3)F)2F2": "[Li+].F[B-]1(F)OC(=O)C(=O)O1",
    "[Li+].N(S(=O)2CF3)S(=O)2CF3": "[Li+].N(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F",
    
    # --- Specific Fixes from your logs ---
    # Fix: C6H3Cl2F is a formula, not SMILES. Structure: 2,3-Dichlorofluorobenzene
    "C6H3Cl2F": "Clc1c(F)c(Cl)ccc1",
    
    # Fix: Extra parenthesis in Diethyl tetrafluorosuccinate
    "CCOC(=O)C(C(=O)OCC)(F)F)(F)F": "CCOC(=O)C(F)(F)C(F)(F)C(=O)OCC",
    
    # Fix: Bad syntax for 1,1-Difluoro-2-iodoethylene
    "C(=C(I)H)(F)F": "FC(F)=CI",
    
    # Fix: Typo in Perfluoro alkyl phosphate (Formula provided instead of SMILES)
    "C27H25F34N2O8PS2": "INVALID_ENTRY", 
    
    # --- General Cleanup ---
    # Fix "nan" appearing in SMILES columns
    "nan": "INVALID_ENTRY",
    
    # Fix common shorthand typos in your file (CF3 -> C(F)(F)F)
    "C(F)CF3": "C(F)C(F)(F)F",
    "COC(CF3)(CF3)": "COC(C(F)(F)F)(C(F)(F)F)"
}

def patch_file(filename):
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return

    print(f"Patching {filename}...")
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    
    count = 0
    for bad, good in REPLACEMENTS.items():
        if bad in content:
            matches = content.count(bad)
            count += matches
            content = content.replace(bad, good)
    
    # Write back
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"  -> Applied {count} fixes.")

if __name__ == "__main__":
    patch_file("inventory.csv")
    patch_file("experiments.csv")
    print("Done. Data is repaired.")