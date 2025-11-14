# Probability Distributions Supported in PT

Parallel Tempering is not restricted to the Boltzmann distribution.
Our library supports:

- **Boltzmann Distribution** (standard)
- **Tsallis Distribution** (heavy-tailed)

---

## 1. Boltzmann Distribution (default)

$$
\pi(x) \propto e^{-\beta E(x)}.
$$

Used when sampling physical systems or Bayesian posteriors.

---

## 2. Tsallis Distribution

A generalized, heavy-tailed distribution:

$$
\pi(x) \propto [1-(1-q)\beta E(x)]^{\frac{1}{1-q}}.
$$

Useful when:

- the landscape is extremely rugged  
- you want broader exploration  
- Boltzmann declines too sharply  

---

## How to Select a Distribution

```python
pt = ParallelTempering(
    energy_fn=energy_fn,
    local_step_fn=local_step_fn,
    distribution="tsallis"
)
```

Default is:

```
distribution="boltzmann"
```

---

> [!CAUTION]
> Tsallis distribution requires choosing the parameter `q`.  
> Users should experiment or consult literature for appropriate values.
