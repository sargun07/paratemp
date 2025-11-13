# parallel tempering

A lightweight, model-agnostic Python library for **Parallel Tempering (Replica Exchange MCMC)**.

## Features
- Generic engine — works for any model (energy-based, probabilistic, physics).
- Plug in your own energy function + local sampling kernel.
- Swaps between replicas follow the correct Metropolis acceptance rule.
- Clean API similar to modern samplers.
- Installable via `pip install -e .`

## Installation

```bash
pip install -e .
