from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, List, Optional, Sequence, Tuple
import math
import random


# Type aliases for clarity
EnergyFn = Callable[[Any], float]
LocalStepFn = Callable[[Any, float, random.Random], Tuple[Any, float]]


@dataclass
class Replica:
    """
    A single replica in the parallel tempering ensemble.

    Attributes
    ----------
    state : Any
        The configuration of the system (user-defined type).
    energy : float
        The energy E(state) at the current state.
    beta : float
        Inverse temperature, beta = 1 / (k_B * T).
    """
    state: Any
    energy: float
    beta: float


class ParallelTempering:
    """
    Generic Parallel Tempering (Replica Exchange) sampler.

    This class is model-agnostic. You provide:
        - an energy function E(x)
        - a local move function that performs one step at fixed beta

    The sampler:
        - holds M replicas at different betas
        - runs local moves in each replica
        - attempts swaps between neighboring replicas

    Parameters
    ----------
    init_states : Sequence[Any]
        Initial state for each replica.
    betas : Sequence[float]
        Inverse temperature values, sorted in increasing order.
    energy_fn : Callable[[Any], float]
        Function computing the energy E(x).
    local_step_fn : Callable[[Any, float, random.Random], (Any, float)]
        One local MCMC/MD step at fixed beta. Returns (new_state, new_energy).
        This function is responsible for internal accept/reject logic.
    rng : random.Random, optional
        Random number generator to use. If None, a new RNG is created.
    """

    def __init__(
        self,
        init_states: Sequence[Any],
        betas: Sequence[float],
        energy_fn: EnergyFn,
        local_step_fn: LocalStepFn,
        rng: Optional[random.Random] = None,
    ) -> None:
        if len(init_states) != len(betas):
            raise ValueError("init_states and betas must have the same length")

        if sorted(betas) != list(betas):
            raise ValueError("betas must be in increasing order")

        self.energy_fn = energy_fn
        self.local_step_fn = local_step_fn
        self.rng = rng or random.Random()

        # Initialize replicas
        self.replicas: List[Replica] = []
        for state, beta in zip(init_states, betas):
            energy = energy_fn(state)
            self.replicas.append(Replica(state=state, energy=energy, beta=beta))

        # Swap statistics
        self.n_swap_attempts: int = 0
        self.n_swaps_accepted: int = 0

    # ------------- Properties & helpers -------------

    @property
    def betas(self) -> List[float]:
        """Current list of inverse temperatures associated with replicas."""
        return [r.beta for r in self.replicas]

    def swap_acceptance_rate(self) -> float:
        """
        Overall swap acceptance rate (accepted / attempted), or 0 if no swaps attempted.
        """
        if self.n_swap_attempts == 0:
            return 0.0
        return self.n_swaps_accepted / self.n_swap_attempts

    # ------------- Local moves -------------

    def step_local(self, n_steps: int = 1) -> None:
        """
        Perform `n_steps` local moves on each replica independently.

        We assume `local_step_fn` already enforces detailed balance at fixed beta.
        """
        for _ in range(n_steps):
            for i, rep in enumerate(self.replicas):
                new_state, new_energy = self.local_step_fn(
                    rep.state, rep.beta, self.rng
                )
                # local_step_fn decides accept/reject; we just take its output.
                self.replicas[i].state = new_state
                self.replicas[i].energy = new_energy

    # ------------- Swap moves -------------

    def _attempt_swap_pair(self, i: int, j: int) -> None:
        """
        Attempt a replica-exchange (swap) between replicas i and j.

        Acceptance probability:
            A = min(1, exp[(β_i - β_j) * (E_j - E_i)])
        """
        rep_i = self.replicas[i]
        rep_j = self.replicas[j]

        beta_i, beta_j = rep_i.beta, rep_j.beta
        E_i, E_j = rep_i.energy, rep_j.energy

        delta = (beta_i - beta_j) * (E_j - E_i)
        accept = False
        if delta >= 0:
            accept = True
        else:
            u = self.rng.random()
            accept = (u < math.exp(delta))

        self.n_swap_attempts += 1
        if accept:
            # Swap just the states + energies; betas stay associated with positions i, j
            self.replicas[i], self.replicas[j] = (
                Replica(state=rep_j.state, energy=rep_j.energy, beta=beta_i),
                Replica(state=rep_i.state, energy=rep_i.energy, beta=beta_j),
            )
            self.n_swaps_accepted += 1

    def attempt_swaps(self, scheme: str = "even-odd") -> None:
        """
        Attempt swaps between neighboring replicas.

        Parameters
        ----------
        scheme : {'all', 'even-odd'}
            'all'      – try all neighbor pairs (0,1), (1,2), ..., (M-2, M-1).
            'even-odd' – alternate between even and odd neighbor pairs
                         to avoid correlations.
        """
        M = len(self.replicas)
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
            raise ValueError(f"Unknown swap scheme: {scheme}")

        for i, j in pairs:
            self._attempt_swap_pair(i, j)

    # ------------- Main driver -------------

    def run(
        self,
        n_iterations: int,
        n_local_steps: int = 1,
        swap_scheme: str = "even-odd",
        callback: Optional[
            Callable[[int, List[Replica], "ParallelTempering"], None]
        ] = None,
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
            Number of local steps per iteration per replica.
        swap_scheme : {'even-odd', 'all'}
            Swap pattern for neighbor pairs.
        callback : callable, optional
            Function called after each iteration:
                callback(iteration_index, replicas, sampler)
        """
        for it in range(n_iterations):
            self.step_local(n_local_steps)
            self.attempt_swaps(swap_scheme)
            if callback is not None:
                callback(it, self.replicas, self)
