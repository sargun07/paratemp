# Advanced Topics in Parallel Tempering  
_A supplementary guide to advanced concepts mentioned in “Parallel Tempering: Theory, Applications, and New Perspectives” (Earl & Deem, 2005)._

This document summarizes **additional theoretical elements** of Parallel Tempering (PT) that were part of our study but not included in the main documentation pages.

These topics include:  
- analytical acceptance-rate predictions  
- iterative temperature-selection formulas  
- multidimensional PT  
- Hamiltonian parallel tempering  
- multicanonical PT  
- density-of-states (Wang–Landau) PT  
- partial-system swaps  
- scaling behavior with system size  
- Tsallis (non-Boltzmann) PT  
- virtual-move PT  
- free-energy PT  

They are provided as theoretical background, future-extension ideas, and references for advanced users.

---

# 1. Acceptance Probability Theory (Kofke, 2002; 2005)

One of the most important theoretical results is the Gaussian-overlap expression for swap acceptance rates.

If the energy distributions at adjacent temperatures are approximately Gaussian, with heat capacity \( C_v \) assumed constant across the interval, Kofke derived:

$$
\langle A \rangle = \operatorname{erfc} \!\left[
\frac{1}{2} \left( \frac{C_v}{2} \right)^{1/2}\,
\frac{1 - \beta_j / \beta_i}
{\sqrt{1 + (\beta_j / \beta_i)^2}}
\right].
$$

Where:

- \( \beta_i, \beta_j \) are inverse temperatures  
- \( C_v \) is heat capacity  
- `erfc` is the complementary error function  

> [!IMPORTANT]
> This formula explains why **temperature spacing** determines swap acceptance.

---

# 2. Iterative Temperature-Adjustment Schemes

Several iterative formulas exist for adjusting temperatures so swap acceptances reach a target.

## 2.1 Sanbonmatsu Equation

Given a target acceptance \( A_{\text{target}} \):

$$
A_{\text{target}} = \exp(\Delta\beta \, \Delta E)
$$

This allows solving for:

- the temperature spacing  
- or differences in mean energies across temperatures  

---

## 2.2 Rathore–Mannar–de Pablo σ-based System

They propose solving:

$$
\frac{\Delta E}{\sigma_m}\bigg|_{T_j}
=
\left( \frac{\Delta E}{\sigma} \right)_{\text{target}}
$$

where:

- \( \sigma_m = \frac{\sigma(T_i) + \sigma(T_j)}{2} \)

Ensures consistent swap acceptance across all pairs.

---

# 3. Multidimensional Parallel Tempering (Yan & de Pablo)

PT is not limited to a **1D temperature ladder**.  
It can operate in *multiple orthogonal dimensions*:

- temperature  
- chemical potential  
- Hamiltonian parameters  
- pulling forces  
- cutoffs in potentials  

A replica becomes an element of an \( n \)-dimensional grid.

### Swaps may occur:
- within dimensions  
- between dimensions  
- across arbitrary axes

> [!NOTE]
> This dramatically increases sampling power but also computational cost.

---

# 4. Hamiltonian Parallel Tempering (Fukunishi, Watanabe, Takada)

Instead of altering temperatures, PT can alter the **Hamiltonian**:

$$
H_i(x) = H_0(x) + \lambda_i \Delta H(x)
$$

Examples:
- scale van der Waals interactions  
- scale hydrophobic interactions  
- modify part of the potential energy term  

Acceptance rule:

$$
A = \min\Big(1,\,
e^{ -\beta [ (H_i(x') + H_j(x)) - (H_i(x) + H_j(x')) ] }
\Big)
$$

Used heavily in biomolecular simulations.

---

# 5. Multicanonical + PT (de Pablo Group)

An alternative ensemble:

$$
p_{\text{multi}}(x) = e^{-\beta E(x)} \, w(E)
$$

where \( w(E) \) is a weight that flattens the energy histogram.

### Benefits:
- lower effective energy barriers  
- fewer replicas required  
- broader overlap of histograms  

### PT usage:
- each replica has its own multicanonical weight  
- swaps prevent simulations from sticking in narrow energy windows  
- weights are refined during equilibration using histogram reweighting

---

# 6. Density-of-States / Wang–Landau PT

Wang–Landau sampling performs a **random walk in energy space**, estimating the density of states \( g(E) \).

In this hybrid:
- PT runs multiple overlapping energy windows  
- each window estimates a region of the density of states  
- swaps allow communication between windows  

Used for:
- protein folding  
- solid–liquid equilibria  
- polymer transitions

---

# 7. Partial-System Swaps (Swendsen–Wang, 1986)

In spin systems, you can exchange **only part of the configuration** (e.g., clusters of spins).  
This solves the problem:

$$
\text{Number of replicas} \propto \sqrt{N}
$$

But extending this to atomistic systems is difficult due to surface energies when splitting molecules.

---

# 8. Virtual-Move Parallel Tempering (Coluzza & Frenkel)

Instead of only exchanging configurations that were *actually attempted*:

- **all possible swaps** between all replica pairs contribute statistical weight  
- even unattempted moves are included virtually  
- a generalization of the “waste-recycling” Monte Carlo method

Reported to improve statistical efficiency by **up to 20×**.

---

# 9. Free-Energy PT (λ-PT)

Used when the Hamiltonian depends on a coupling parameter \( \lambda \):

$$
U_\lambda = (1-\lambda) U_0 + \lambda U_1
$$

Replicas are placed at different λ.  
Swaps allow fast movement through alchemical states.

Useful for:
- free energy differences  
- chemical transformations  
- ligand binding calculations  

---

# 10. Tsallis (Non-Boltzmann) PT (Whitfield et al.)

Tsallis distribution:

$$
\pi_q(x) \propto [1 - (1-q)\beta E(x)]^{\frac{1}{1-q}}
$$

Properties:
- heavy tails  
- weaker suppression of high-energy states  
- better exploration in extremely rugged landscapes  

Used in:
- peptide conformational search  
- generalized ensembles  

---

# 11. Scaling Behavior with System Size

For canonical PT:

- histogram width ∝ √N  
- mean energy ∝ N  
- temperature spacing must scale accordingly  
- number of replicas needed:

$$
M \propto \sqrt{N}
$$

This becomes a bottleneck for large biomolecular systems.

---

# Summary

This document consolidates advanced PT concepts:

- analytical acceptance formulas  
- optimal temperature selection  
- multidimensional PT  
- Hamiltonian PT  
- multicanonical and density-of-states hybrids  
- partial swaps  
- Tsallis statistics  
- free-energy PT  
- scaling and performance considerations  

These ideas expand on the main theory in the PT review paper and represent established extensions of the method.

They are intended for expert users and future development of the `paratemp` library.

