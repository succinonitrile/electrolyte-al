import os

# Mappings of Invalid SMILES -> Valid RDKit SMILES
REPLACEMENTS = {
    # Salts
    "[Li+].[BF4-]": "[Li+].F[B-](F)(F)F",  # LiBF4
    "[Li+].[NO3-]": "[Li+].[O-][N+](=O)[O-]",  # LiNO3
    "[Li+].B(OC(CF3)F)2F2": "[Li+].F[B-]1(F)OC(=O)C(=O)O1",  # LiDFOB (corrected to standard)
    "[Li+].N(S(=O)2CF3)S(=O)2CF3": "[Li+].N(S(=O)(=O)C(F)(F)F)S(=O)(=O)C(F)(F)F", # LiTFSI (fixed CF3 group)
    
    # Common Bad Patterns in your file
    "C(F)CF3": "C(F)C(F)(F)F", # Fix compressed CF3
    "C(=C(I)H)(F)F": "FC(F)=CI", # 1,1-difluoro-2-iodoethylene
    "C6H3Cl2F": "Clc1c(F)c(Cl)ccc1", # 2,3-Dichlorofluorobenzene (Guessing structure from formula)
    "nan": "INVALID_ENTRY"
}

def patch_file(filename):
    if not os.path.exists(filename):
        print(f"File not found: {filename}")
        return

    print(f"Patching {filename}...")
    with open(filename, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Apply all replacements
    count = 0
    for bad, good in REPLACEMENTS.items():
        if bad in content:
            matches = content.count(bad)
            count += matches
            content = content.replace(bad, good)
            print(f"  - Fixed {matches} occurrences of {bad.split(']')[0]}...")

    with open(filename, 'w', encoding='utf-8') as f:
        f.write(content)
    print(f"  > Saved with {count} fixes.\n")

if __name__ == "__main__":
    patch_file("inventory.csv")
    patch_file("experiments.csv")
    print("Data repair complete. Now run your active learning script.")