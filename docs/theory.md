# Parallel Tempering Theory

Parallel Tempering simulates a joint distribution over multiple **replicas**, each at a different temperature.

Let replicas be indexed:

$$
T_1 < T_2 < \cdots < T_M
$$

with inverse temperatures:

$$
\beta_i = \frac{1}{k_B T_i}.
$$

---

## 1. Single-Temperature MCMC

Standard Metropolis sampling uses the Boltzmann distribution:

$$
\pi(x) \propto e^{-\beta E(x)}.
$$

The local acceptance rule is:

$$
A_{\text{local}} = \min\left(1, e^{-\beta(E(x') - E(x))}\right).
$$

> [!NOTE]
> A random number $u \sim U(0,1)$ determines whether a proposal is accepted.

---

## 2. Joint Distribution in PT

For \(M\) replicas with states \( x_1, x_2, \ldots, x_M \):

$$
\Pi(x_1,\ldots,x_M) = \prod_{i=1}^M e^{-\beta_i E(x_i)}.
$$

Each replica independently samples from its own Boltzmann distribution.

---

## 3. Replica Exchange (Swap Move)

Swapping the states of replicas \(i\) and \(j=i+1\):

$$
(x_i, x_j) \rightarrow (x_j, x_i)
$$

is accepted with:

$$
A_{\text{swap}} = \min\left(
1,\,
e^{(\beta_i - \beta_j)\,\left(E(x_j) - E(x_i)\right)}
\right).
$$

> [!IMPORTANT]
> Swap acceptance depends on **energy histogram overlap** between adjacent temperatures.

If temperatures are too far apart → no overlap → swaps rarely accepted → PT fails.

---

## 4. Why PT Works

- High-temperature replicas can cross energy barriers.
- Low-temperature replicas refine low-energy configurations.
- Swaps let configurations **travel across the temperature ladder**.

PT allows a Markov chain to escape local minima and explore a rugged landscape efficiently.

---

## 5. Detailed Balance

Both:

- local moves  
- swap moves  

are constructed to satisfy detailed balance with respect to the joint distribution.

Thus PT is an **exact sampling method**, unlike simulated annealing.

