# Implementation Details of the `paratemp` Library

This document explains the software architecture and design choices.

---

## 1. Replica Dataclass

```python
@dataclass
class Replica:
    state: Any
    energy: float
    beta: float
    index: int
```

Notes:

- `index` is the fixed **ladder position** for the temperature slot.
- Only `state` and `energy` move during swaps.
- Temperatures never move.

---

## 2. Key Components

### **local_step_fn**
User-defined MCMC proposal:

```python
x_new, E_new = local_step_fn(x, beta, rng)
```

### **step_local(n)**
Perform `n` MCMC updates per replica.

### **attempt_swaps()**
Try swapping every adjacent pair.

### **run()**
Main driver:

1. `step_local(n_local_steps)`
2. `attempt_swaps()`
3. optional callback
4. repeat for N iterations

---

## 3. Acceptance Rules

### Local Metropolis:
$$
A_{\text{local}} = \min(1, e^{-\beta(E' - E)}).
$$

### Swap Metropolis:
$$
A_{\text{swap}} = \min(1, e^{(\beta_i-\beta_j)(E_j - E_i)}).
$$

---

## 4. Initial States

Two options:

### A. Single initial state
```python
init_state = 0.0
```

Distributed to all replicas.

### B. Multiple initial states
```python
init_states = [0.0, 1.0, -1.0, ...]
```

> [!NOTE]
> Starting all replicas at the same state is common and correct —  
> PT rapidly decorrelates them.

