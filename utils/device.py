import torch
import logging

# DO NOT MODIFY: This logic ensures single-GPU usage on desktop
def get_device(use_gpu: bool = True, device_str: str = "cuda:0") -> torch.device:
    """
    Returns the torch device. Enforces single-GPU or CPU fallback.
    """
    if use_gpu and torch.cuda.is_available():
        # Check if the specific device index is valid
        if "cuda" in device_str:
            try:
                # Test the device
                _ = torch.tensor([1]).to(device_str)
                return torch.device(device_str)
            except Exception as e:
                logging.warning(f"Requested {device_str} not available ({e}). Falling back to cpu.")
                return torch.device("cpu")
        return torch.device("cuda")
    
    return torch.device("cpu")