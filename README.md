# Barnes-Hut 2D Simulation

A C++17 implementation of a two-dimensional gravitational N-body simulation.
The Barnes-Hut quadtree approximates distant groups of particles to reduce the
cost of force evaluation; OpenMP is used for parallel force and integration
loops when available.

## Build

Requirements: CMake 3.16 or newer and a C++17 compiler. OpenMP is optional.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

To explicitly build without OpenMP:

```sh
cmake -S . -B build-serial -DBARNES_HUT_ENABLE_OPENMP=OFF
cmake --build build-serial --parallel
```

## Run

```sh
./build/barnes-hut --particles 10000 --steps 100 --dt 0.01 \
  --theta 0.5 --seed 42 --threads 4
```

All options have defaults. Run `./build/barnes-hut --help` for the full list.
The default seed makes particle initialization repeatable. A theta of zero
disables the Barnes-Hut approximation and computes direct interactions within
quadtree leaves, which is useful for small correctness checks; larger theta
values generally increase approximation and reduce work.

## Design and validation

The public simulation and quadtree API is in `include/barnes_hut/simulation.hpp`.
The quadtree uses bucket leaves and retains particles even when many occupy the
same coordinates. Nodes containing the target particle are always opened, so a
particle is not included in its own approximated force. Forces at theta zero
are tested against a direct pairwise sum.

The tests are intentionally dependency-free. The executable reports the chosen
configuration and elapsed simulation time. For meaningful performance claims,
record compiler/version, build type, CPU, particle count, steps, theta, and
thread count, and compare against a direct-force baseline on the same machine.
