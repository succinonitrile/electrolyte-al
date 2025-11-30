import argparse
import yaml
import pandas as pd
from utils.device import get_device
from active_learning.learner import ActiveLearner

def main():
    parser = argparse.ArgumentParser(description="Electrolyte Active Learning (Desktop)")
    parser.add_argument("--config", type=str, default="config/default.yaml", help="Path to config YAML")
    parser.add_argument("--data", type=str, required=True, help="Path to experimental data CSV")
    parser.add_argument("--purchase", type=str, required=True, help="Path to purchasability CSV")
    parser.add_argument("--output", type=str, default="recommendations.csv", help="Output file")
    
    args = parser.parse_args()

    # 1. Load Config
    with open(args.config, "r") as f:
        config = yaml.safe_load(f)

    # 2. Setup Device
    device = get_device(
        use_gpu=config['hardware']['use_gpu'], 
        device_str=config['hardware']['cuda_device']
    )
    print(f"Running on device: {device}")

    # 3. Initialize Learner
    learner = ActiveLearner(config, device, args.data, args.purchase)

    # 4. Train
    print("Starting Model Training...")
    learner.load_data_and_train()

    # 5. Active Learning Step
    print("Generating and Scoring Candidates...")
    recommendations = learner.run_cycle()

    # 6. Save
    rows = []
    for category, formulated_list in recommendations.items():
        for f in formulated_list:
            # Flatten for CSV
            rows.append({
                "category": category,
                "formulation": str(f) # Convert object to string/json representation
            })
    
    pd.DataFrame(rows).to_csv(args.output, index=False)
    print(f"Recommendations saved to {args.output}")

if __name__ == "__main__":
    main()