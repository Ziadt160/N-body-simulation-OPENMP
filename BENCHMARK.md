# Benchmarks

Measured on 28 Oct 2026 from commit `ed7a735`.

**Machine:** AMD Ryzen 5 3600 (6 cores, 12 threads), WSL2 Ubuntu, g++ 15.2.0, CMake 4.2.3.
**Build:** `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release`, plus a serial build with `-DBARNES_HUT_ENABLE_OPENMP=OFF`.
**Method:** wall time printed by `barnes-hut` ("Simulation completed in … seconds"), best of 3 runs, `--seed 42`.

## Thread scaling

20,000 bodies, 50 steps, θ = 0.5:

| Build / threads | Time (s) | Speedup vs 1 thread |
|---|---|---|
| OpenMP off (serial build) | 3.44 | – |
| 1 thread | 3.41 | 1.0× |
| 2 threads | 1.94 | 1.8× |
| 4 threads | 1.19 | 2.9× |
| 6 threads | 0.97 | 3.5× |
| 12 threads (SMT) | 0.83 | 4.1× |

Scaling flattens beyond the 6 physical cores. The quadtree is still rebuilt serially every step, so that part of each step doesn't scale with threads. Building the tree in parallel is the next step.

## Barnes-Hut approximation vs exact mode

5,000 bodies, 10 steps, 1 thread:

| θ | Time (s) |
|---|---|
| 0 (exact: every node opened) | 2.33 |
| 0.5 | 0.118 |

θ = 0.5 runs 19.8× faster. `tests/test_simulation.cpp` checks θ = 0 forces against a direct pairwise sum to 1e-12.

## Reproduce

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel
ctest --test-dir build --output-on-failure
for t in 1 2 4 6 12; do ./build/barnes-hut --particles 20000 --steps 50 --theta 0.5 --seed 42 --threads $t; done
./build/barnes-hut --particles 5000 --steps 10 --theta 0 --seed 42 --threads 1
./build/barnes-hut --particles 5000 --steps 10 --theta 0.5 --seed 42 --threads 1
```
