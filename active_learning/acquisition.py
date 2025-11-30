import numpy as np

def select_candidates(candidates, means, stds, n_high, n_low, n_unc, beta=1.0):
    """
    CPU-friendly sorting and selection.
    """
    # LCB for High Performance (Optimistic) or UCB (depending on def)
    # Usually for maximization: UCB = Mean + Beta*Std
    ucb = means + beta * stds
    lcb = means - beta * stds
    
    # 1. High Performance (High UCB)
    high_indices = np.argsort(ucb)[-n_high:][::-1]
    
    # 2. Low Performance (Low LCB - predicted to fail)
    low_indices = np.argsort(lcb)[:n_low]
    
    # 3. Uncertainty (High Std)
    unc_indices = np.argsort(stds)[-n_unc:][::-1]
    
    return {
        "high_performance": [candidates[i] for i in high_indices],
        "low_performance": [candidates[i] for i in low_indices],
        "uncertainty": [candidates[i] for i in unc_indices]
    }