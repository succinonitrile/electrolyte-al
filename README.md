# Single-Desktop Electrolyte Active Learning System

A lightweight, GPU-accelerated active learning framework for designing Li-metal battery electrolytes. Designed to run on a standard consumer desktop (e.g., Intel i5 + Nvidia RTX 4060).

## Hardware Requirements
- **CPU**: Modern consumer CPU (Intel i5/AMD Ryzen 5 or better).
- **GPU**: Nvidia GPU (RTX 4060 class with ~8GB VRAM recommended).
- **RAM**: 16GB+.

## Installation
1. Install Python 3.10+.
2. Install PyTorch (ensure CUDA version matches your driver):
   ```bash
   pip install torch torchvision --index-url [https://download.pytorch.org/whl/cu118](https://download.pytorch.org/whl/cu118)