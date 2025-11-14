# Parallel Tempering (Replica Exchange Monte Carlo)
A documentation overview for the `paratemp` Python library.

Parallel Tempering (PT), also known as Replica Exchange Monte Carlo, is a sampling technique designed to efficiently explore multimodal or rugged energy landscapes. PT is widely used in:

- statistical physics  
- computational chemistry  
- biomolecular simulations  
- Bayesian inference  
- optimization problems such as Max-Cut  

The key idea:

- Run **multiple replicas** of the same system at different temperatures.
- High-temperature replicas explore broadly.
- Low-temperature replicas explore low-energy basins.
- Periodically **swap** configurations between neighboring temperatures using a Metropolis-style acceptance rule.

This dramatically improves the ability to cross energy barriers and sample complex distributions.

---

## Features of the `paratemp` Library

- Generic PT implementation (any state representation)
- Pluggable energy function and proposal move
- Geometric temperature ladder
- Boltzmann and Tsallis distributions
- Replica-count auto-calibration
- Energy histogram utilities
- Full swap acceptance diagnostics

> [!TIP]
> PT is most useful when standard MCMC gets trapped in local minima.

---

## Documentation Structure

This folder contains the following reference pages:

- [`theory.md`](theory.md) — PT theory and equations  
- [`temperature_and_replicas.md`](temperature_and_replicas.md) — choosing temperatures & replica count  
- [`distributions.md`](distributions.md) — Boltzmann vs Tsallis  
- [`implementation.md`](implementation.md) — how the library is implemented  
- [`diagnostics.md`](diagnostics.md) — acceptance statistics & histograms  
- [`references.md`](references.md) — scientific references  
