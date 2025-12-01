import pandas as pd
import os

# --- Mappings of Bad SMILES -> Good SMILES ---
# We use simple string replacement to ensure we catch them inside JSON strings (experiments.csv)
# and standard columns (inventory.csv).
REPLACEMENTS = {
    # 1. Lithium Tetrafluoroborate (LiBF4)
    # Old: [Li+].[BF4-] -> New: Explicit connectivity
    "[Li+].[BF4-]": "[Li+].F[B-](F)(F)F",
    
    # 2. Lithium Difluoro(oxalato)borate (LiDFOB)
    # Old: [Li+].B(OC(CF3)F)2F2 (Invalid syntax) -> New: Correct structure
    "[Li+].B(OC(CF3)F)2F2": "[Li+].F[B-]1(F)OC(=O)C(=O)O1",
    
    # 3. Lithium Nitrate (LiNO3)
    # Old: [Li+].[NO3-] -> New: Explicit resonance structure
    "[Li+].[NO3-]": "[Li+].[O-][N+](=O)[O-]",
    
    # 4. Diphenyl Carbonate
    # Old: Broken ring closure -> New: Standard Canonical
    "c1ccc(Oc2ccccc2)OC(=O)c3ccccc3": "O=C(Oc1ccccc1)Oc2ccccc2",
    
    # 5. Perfluoro alkyl phosphate
    # Old: Formula instead of SMILES -> New: Set to False (exclude from purchasable)
    "C27H25F34N2O8PS2": "FALSE_SMILES_INVALID" 
}

def clean_file(filename):
    if not os.path.exists(filename):
        print(f"Skipping {filename} (not found)")
        return

    print(f"Processing {filename}...")
    
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Perform replacements
    count = 0
    for bad, good in REPLACEMENTS.items():
        if bad in content:
            matches = content.count(bad)
            count += matches
            content = content.replace(bad, good)
            print(f"  - Replaced {matches} occurance(s) of '{bad}'")
            
    # Write back to file
    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)
        
    print(f"  -> Saved {filename} with {count} corrections.\n")

if __name__ == "__main__":
    # Fix both files
    clean_file("inventory.csv")
    clean_file("experiments.csv")
    
    print("Done! You can now run your active learning script.")