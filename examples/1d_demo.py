import math
import random
import numpy as np

from paratemp import ParallelTempering

# defining 1 D energy function
def energy_fn(x: float) -> float:
    return float((x*x - 1.0) ** 2)

# single temperature local monte carlo move -> it returns new state and new energy
def local_step_fn(state: float, beta: float, rng: random.Random):
    
    x = state
    E_x = energy_fn(x)

    step_size = 0.5
    proposal = x + rng.normalvariate(0.0, step_size)
    E_prop = energy_fn(proposal)

    # acceptance rule
    delta = beta * (E_prop - E_x)
    A = min(1,math.exp(-beta * delta))

    if A==1.0 or rng.random() < A:
        return proposal, E_prop
    else:
        return x, E_x


def main():

    # taking min temp, max temp and number of replicas as input
    T_min = 0.5
    T_max = 5.0
    n_replicas = 8

    init_state = 0.0
    rng = random.Random(42)

    pt = ParallelTempering(
        energy_fn=energy_fn,
        local_step_fn=local_step_fn,
        T_min=T_min,
        T_max=T_max,
        n_replicas=n_replicas,
        ladder="geometric",
        distribution="boltzmann",  
        init_state=init_state,
        rng=rng,
    )

    # we will store the values of the lowest temperature 
    samples_lowT = []

    # it prints overall acceptance rate after each swapping period
    def cb(it, replicas, sampler):
        samples_lowT.append(float(replicas[0].state))
        if (it + 1) % 200 == 0:
            print(
                f"Iter {it+1:5d} | "
                f"swap acc = {sampler.swap_acceptance_rate():.3f}"
            )

    pt.run(n_iterations=2000, n_local_steps=10, callback=cb)


if __name__ == "__main__":
    main()
