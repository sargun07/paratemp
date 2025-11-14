# Diagnostics, Analysis & Debugging in PT

This page explains how to interpret acceptance rates, histograms, and trajectories.

---

## 1. Global Swap Acceptance Rate

Defined as:

$$
A_{\text{global}} = \frac{\text{accepted swaps}}{\text{attempted swaps}}.
$$

Use:

```python
sampler.swap_acceptance_rate()
```

---

## 2. Per-Pair Acceptance

Useful for diagnosing temperature spacing problems.

```python
sampler.pair_acceptance_rates()
```

---

## 3. Energy Histograms

Energy histograms indicate:

- sampling quality  
- overlap of distributions  
- multimodality at low temperature  

Our library includes a textual histogram helper.

---

## 4. Low-Temperature Trajectories

Tracking:

```python
replicas[0].state
```

gives insight into:

- barrier crossings  
- mixing behavior  
- sampling stability  

---

## 5. What Good PT Looks Like

- swap acceptance ≈ 0.20–0.30  
- smooth movement of configurations across ladder  
- energy histogram with appropriate overlap  
- low-T replica jumps between modes over time  

