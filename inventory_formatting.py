import os

input_path = "inventory.csv"
output_path = "inventory_clean.csv"

print(f"Processing {input_path}...")

with open(input_path, "r", encoding="utf-8") as f_in, open(output_path, "w", encoding="utf-8", newline="") as f_out:
    lines = f_in.readlines()
    
    # 1. Fix the Header
    # Replace the malformed first line with the standard clean header
    f_out.write("name,smiles,purchasable,vendor,notes\n")
    
    # 2. Fix the Rows
    # Skip the original bad header (index 0) and process the rest
    for i, line in enumerate(lines[1:]):
        # Replace triple quotes """ with single standard quotes "
        # This converts """[Li+].F...""" to "[Li+].F..." which is valid CSV
        clean_line = line.replace('"""', '"')
        
        # Remove any leading/trailing whitespace issues if present
        clean_line = clean_line.strip() + "\n"
        
        f_out.write(clean_line)

print("Done!")
print(f"Clean file saved to: {output_path}")
print("Please rename 'inventory_clean.csv' to 'inventory.csv' before running your code.")