# Temperature Selection & Number of Replicas

Choosing temperatures and replica count is one of the most important steps in PT.

---

## 1. Requirements for Good Temperature Ladders

1. **High temperature** must be high enough to escape local minima.  
2. **Energy histograms of adjacent replicas must overlap.**  
3. **Swap acceptance should be approximately 20–30%.**  
4. **Number of replicas** should be large enough for smooth mixing,  
   but small enough to avoid wasted computation.

---

## 2. Geometric Temperature Ladder

A common and effective choice:

$$
T_i = T_{\min} \cdot r^i, \quad 
r = \left(\frac{T_{\max}}{T_{\min}}\right)^{1/(M-1)}.
$$

This often yields **nearly uniform swap acceptance** (Kofke 2002).

---

## 3. Target Acceptance Rate ≈ 0.25

Studies (Rathore et al., Kone & Kofke) show optimal global mixing when:

$$
A_{\text{swap}} \approx 0.23.
$$

> [!TIP]
> 0.20–0.30 is the “sweet spot.”  
> Below 0.15 → poor mixing.  
> Above 0.50 → unnecessary replicas.

---

## 4. Auto-Choosing the Number of Replicas

The library provides a helper function that:

- tries candidate values of `n_replicas`
- runs a short PT trial for each
- measures swap acceptance
- selects the value closest to the target (≈ 0.25)

This matches best practices from PT literature.

---

## 5. When to Use More Replicas

Use more replicas when:

- the system size is large  
- energy barriers are very high  
- energy histograms do not overlap  
- acceptance drops too low

Use fewer replicas when:

- acceptance is > 0.5  
- temperatures differ only slightly  
- computation cost is high  

