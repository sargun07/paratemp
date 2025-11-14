from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable, List, Optional, Sequence, Tuple, Dict
import math
import random


# -------------------------------------------------------------------
# Type aliases
# 1. EnergyFn: it takes a state and returns a float energy
# 2.  it takes state, beta, rng and returns new state and new energy (this is the single-temperature MCMC step)
# 3. function called after each PT iteration
# -------------------------------------------------------------------

EnergyFn = Callable[[Any], float]
LocalStepFn = Callable[[Any, float, random.Random], Tuple[Any, float]]
CallbackFn = Callable[[int, List["Replica"], "ParallelTempering"], None]


# -------------------------------------------------------------------
# Distributions (Boltzmann, Tsallis)
# -------------------------------------------------------------------

class BaseDistribution:

    def log_weight(self, energy: float, beta: float) -> float:
        """
        Return log π(x) up to an additive constant.
        Must be implemented by subclasses.
        """
        raise NotImplementedError


class BoltzmannDistribution(BaseDistribution):
    """
    Standard Boltzmann distribution: π(x) ∝ exp(-β E(x))
    """

    def log_weight(self, energy: float, beta: float) -> float:
        return -beta * energy


class TsallisDistribution(BaseDistribution):
    """
    Tsallis distribution with parameter q.

    π(x) ∝ [1 - (1 - q) β E(x)]^(1 / (1 - q))

    For q -> 1 this approaches the Boltzmann distribution.
    """

    def __init__(self, q: float = 1.2) -> None:
        if q == 1.0:
            raise ValueError(
                "For q = 1.0, Tsallis reduces to Boltzmann; "
                "use BoltzmannDistribution instead."
            )
        self.q = q

    def log_weight(self, energy: float, beta: float) -> float:
        base = 1.0 - (1.0 - self.q) * beta * energy
        if base <= 0.0:
            # Outside the support → zero probability → log π = -∞
            return float("-inf")
        return math.log(base) / (1.0 - self.q)


# -------------------------------------------------------------------
# Replica data structure
# -------------------------------------------------------------------

@dataclass
class Replica:
    """
    One replica in the parallel tempering ensemble.

    Attributes: state, energy, beta (inverse temperature), index (index of this replica in the ladder)
    """

    state: Any
    energy: float
    beta: float
    index: int


# -------------------------------------------------------------------
# Temperature ladder helpers
# -------------------------------------------------------------------

def geometric_temperatures(
    T_min: float,
    T_max: float,
    n_replicas: int,
) -> List[float]:
    """
    Generate a geometric temperature ladder between T_min and T_max.

    T_i = T_min * r^(i),  i = 0, ..., n_replicas-1
    where r = (T_max / T_min)^(1 / (n_replicas - 1))
    """
    if n_replicas < 2:
        raise ValueError("n_replicas must be at least 2.")
    if T_min <= 0.0 or T_max <= 0.0:
        raise ValueError("Temperatures must be positive.")
    if T_max <= T_min:
        raise ValueError("T_max must be > T_min.")

    r = (T_max / T_min) ** (1.0 / (n_replicas - 1))
    temps = [T_min * (r ** i) for i in range(n_replicas)]
    return temps


# -------------------------------------------------------------------
# Parallel Tempering engine
# -------------------------------------------------------------------

class ParallelTempering:
    """
    Generic Parallel Tempering sampler.

    It holds M replicas at different temperatures, run local moves in each replica, attempts swaps between neighbouring replicas, track swap statistics

    Parameters
    ----------
    energy_fn : callable
        energy_fn(state) -> float
    local_step_fn : callable
        local_step_fn(state, beta, rng) -> (new_state, new_energy)
        This function should implement a valid single-temperature MCMC step
        that leaves π(x) (e.g. Boltzmann at that β) invariant.
    temperatures : Sequence[float], optional
        Explicit list of temperatures (T). If not provided, T_min, T_max,
        and n_replicas must be given.
    T_min, T_max : float, optional
        Minimum and maximum temperatures for ladder construction if
        `temperatures` is not provided.
    n_replicas : int, optional
        Number of replicas if using ladder construction. Default: 8.
    ladder : {"geometric"}, optional
        Method to construct temperatures if `temperatures` is None.
        For now we support only "geometric".
    distribution : {"boltzmann", "tsallis"}, optional
        Probability distribution type. Default: "boltzmann".
    tsallis_q : float, optional
        q parameter for Tsallis distribution. Only used if distribution="tsallis".
    init_state : Any, optional
        Single initial state to copy across all replicas.
    init_states : Sequence[Any], optional
        List of initial states, one per replica. Overrides init_state if given.
    rng : random.Random, optional
        Random number generator. If None, a new RNG with default seed is used.
    """

    def __init__(
        self,
        energy_fn: EnergyFn,
        local_step_fn: LocalStepFn,
        temperatures: Optional[Sequence[float]] = None,
        T_min: Optional[float] = None,
        T_max: Optional[float] = None,
        n_replicas: int = 8,
        ladder: str = "geometric",
        distribution: str = "boltzmann",
        tsallis_q: float = 1.2,
        init_state: Optional[Any] = None,
        init_states: Optional[Sequence[Any]] = None,
        rng: Optional[random.Random] = None,
    ) -> None:
        # ----------------------------
        # Store core callables and RNG
        # ----------------------------
        self.energy_fn: EnergyFn = energy_fn
        self.local_step_fn: LocalStepFn = local_step_fn
        self.rng: random.Random = rng or random.Random()

        # ----------------------------
        # Choose distribution model
        # ----------------------------
        if distribution.lower() == "boltzmann":
            self.distribution: BaseDistribution = BoltzmannDistribution()
        elif distribution.lower() == "tsallis":
            self.distribution = TsallisDistribution(q=tsallis_q)
        else:
            raise ValueError(
                f"Unknown distribution '{distribution}'. "
                "Supported: 'boltzmann', 'tsallis'."
            )

        # ----------------------------
        # Build temperature ladder
        # ----------------------------
        if temperatures is not None:
            temps = list(temperatures)
            if len(temps) < 2:
                raise ValueError("At least two temperatures are required.")
        else:
            if T_min is None or T_max is None:
                raise ValueError(
                    "If `temperatures` is not given, T_min and T_max must be provided."
                )
            if ladder == "geometric":
                temps = geometric_temperatures(T_min, T_max, n_replicas)
            else:
                raise ValueError(f"Unsupported ladder type: '{ladder}'")

        # Sort temperatures ascending for clarity:
        temps = sorted(temps)
        self.temperatures: List[float] = temps
        self.betas: List[float] = [1.0 / T for T in temps]

        M = len(self.temperatures)

        # ----------------------------
        # Initialize replicas
        # ----------------------------
        if init_states is not None:
            if len(init_states) != M:
                raise ValueError(
                    "len(init_states) must equal number of replicas (len(temperatures))."
                )
            states = list(init_states)
        else:
            if init_state is None:
                raise ValueError(
                    "Either `init_state` or `init_states` must be provided."
                )
            states = [init_state for _ in range(M)]

        self.replicas: List[Replica] = []
        for idx, (state, beta) in enumerate(zip(states, self.betas)):
            energy = self.energy_fn(state)
            self.replicas.append(Replica(state=state, energy=energy, beta=beta, index=idx))

        # ----------------------------
        # Swap statistics
        # ----------------------------
        self.n_swap_attempts: int = 0
        self.n_swaps_accepted: int = 0
        # Per-pair statistics: key = (i, j), i<j
        self.pair_stats: Dict[Tuple[int, int], Dict[str, int]] = {}

    # ------------------------------------------------------------------
    # Properties & helpers
    # ------------------------------------------------------------------

    @property
    def n_replicas(self) -> int:
        return len(self.replicas)

    def swap_acceptance_rate(self) -> float:
        """
        Overall swap acceptance rate across all neighbor pairs.
        """
        if self.n_swap_attempts == 0:
            return 0.0
        return self.n_swaps_accepted / self.n_swap_attempts

    def pair_acceptance_rates(self) -> Dict[Tuple[int, int], float]:
        """
        Acceptance rate for each neighbor pair (i, j).
        Returns a dict: {(i, j): rate}
        """
        rates: Dict[Tuple[int, int], float] = {}
        for pair, stats in self.pair_stats.items():
            att = stats.get("attempts", 0)
            acc = stats.get("accepted", 0)
            rates[pair] = 0.0 if att == 0 else acc / att
        return rates

    # ------------------------------------------------------------------
    # Core PT operations
    # ------------------------------------------------------------------

    def step_local(self, n_steps: int = 1) -> None:
        """
        Perform `n_steps` local MCMC moves on each replica independently.

        """
        for _ in range(n_steps):
            for i, rep in enumerate(self.replicas):
                new_state, new_energy = self.local_step_fn(
                    rep.state, rep.beta, self.rng
                )
                # local_step_fn is responsible for accept/reject; we just store.
                self.replicas[i].state = new_state
                self.replicas[i].energy = new_energy

    def _log_weight(self, energy: float, beta: float) -> float:
        return self.distribution.log_weight(energy, beta)

    def _attempt_swap_pair(self, i: int, j: int) -> None:
        """
        Attempt a replica-exchange (swap) between replicas i and j.

        Generic Metropolis–Hastings acceptance:

            A = min(1, [ π_i(x_j) π_j(x_i) ] / [ π_i(x_i) π_j(x_j) ])

        where π_k uses the distribution model (Boltzmann, Tsallis, etc.)
        associated with replica k's β.
        """
        if i == j:
            return
        if not (0 <= i < self.n_replicas and 0 <= j < self.n_replicas):
            return

        rep_i = self.replicas[i]
        rep_j = self.replicas[j]

        beta_i, beta_j = rep_i.beta, rep_j.beta
        E_i, E_j = rep_i.energy, rep_j.energy

        log_pi_xi = self._log_weight(E_i, beta_i)
        log_pi_xj = self._log_weight(E_j, beta_i)
        log_pj_xj = self._log_weight(E_j, beta_j)
        log_pj_xi = self._log_weight(E_i, beta_j)

        log_num = log_pi_xj + log_pj_xi
        log_den = log_pi_xi + log_pj_xj
        log_A = log_num - log_den

        # A = exp(log_A), but clip at 1
        A = 1.0 if log_A >= 0.0 else math.exp(log_A)

        # Update global stats
        self.n_swap_attempts += 1
        pair = (min(i, j), max(i, j))
        stats = self.pair_stats.setdefault(pair, {"attempts": 0, "accepted": 0})
        stats["attempts"] += 1

        if self.rng.random() < A:
            # Swap states *and* energies; betas stay with replicas i & j
            new_rep_i = Replica(
                state=rep_j.state,
                energy=rep_j.energy,
                beta=beta_i,
                index=rep_i.index,
            )
            new_rep_j = Replica(
                state=rep_i.state,
                energy=rep_i.energy,
                beta=beta_j,
                index=rep_j.index,
            )
            self.replicas[i] = new_rep_i
            self.replicas[j] = new_rep_j

            self.n_swaps_accepted += 1
            stats["accepted"] += 1

    def attempt_swaps(self, scheme: str = "even-odd") -> None:
        """
        Attempt swaps between neighboring replicas.

        Parameters
        ----------
        scheme : {"even-odd", "all"}
            'all'      – try all neighbor pairs (0,1), (1,2), ..., (M-2, M-1).
            'even-odd' – on each call, randomly choose between:
                          (0,1), (2,3), ...
                       or (1,2), (3,4), ...
        """
        M = self.n_replicas
        if M < 2:
            return

        if scheme == "all":
            pairs = [(i, i + 1) for i in range(M - 1)]
        elif scheme == "even-odd":
            if self.rng.random() < 0.5:
                # even pairs: (0,1), (2,3), ...
                pairs = [(i, i + 1) for i in range(0, M - 1, 2)]
            else:
                # odd pairs: (1,2), (3,4), ...
                pairs = [(i, i + 1) for i in range(1, M - 1, 2)]
        else:
            raise ValueError(f"Unknown swap scheme: '{scheme}'")

        for i, j in pairs:
            self._attempt_swap_pair(i, j)

    def run(
        self,
        n_iterations: int,
        n_local_steps: int = 1,
        swap_scheme: str = "all",
        callback: Optional[CallbackFn] = None,
    ) -> None:
        """
        Run the parallel tempering simulation.

        Each iteration:
          1. runs `n_local_steps` local moves on all replicas
          2. attempts swaps between neighboring replicas

        Parameters
        ----------
        n_iterations : int
            Number of outer iterations.
        n_local_steps : int
            Number of local MCMC steps per iteration per replica.
        swap_scheme : {"even-odd", "all"}
            Swap pattern for neighbor pairs.
        callback : callable, optional
            Function called after each iteration:
                callback(iteration_index, replicas, sampler)
        """
        for it in range(n_iterations):
            self.step_local(n_local_steps)
            self.attempt_swaps(scheme=swap_scheme)
            if callback is not None:
                callback(it, self.replicas, self)
