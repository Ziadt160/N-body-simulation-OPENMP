#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace barnes_hut {

inline constexpr double gravitational_constant = 1.0;
inline constexpr double softening = 1e-5;

struct Particle {
    double x = 0.0;
    double y = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double fx = 0.0;
    double fy = 0.0;
    double mass = 1.0;
};

struct SimulationConfig {
    std::size_t particle_count = 10'000;
    std::size_t steps = 100;
    double time_step = 0.01;
    double theta = 0.5;
    std::uint64_t seed = 42;
    int threads = 0;
};

class QuadTree {
public:
    explicit QuadTree(const std::vector<Particle>& particles, double theta);
    ~QuadTree();

    QuadTree(QuadTree&&) noexcept;
    QuadTree& operator=(QuadTree&&) noexcept;
    QuadTree(const QuadTree&) = delete;
    QuadTree& operator=(const QuadTree&) = delete;

    std::pair<double, double> calculate_force(const Particle& particle) const;
    std::size_t particle_count() const noexcept;
    double total_mass() const noexcept;

private:
    struct Node;

    std::unique_ptr<Node> root_;
    double theta_;
    std::size_t particle_count_ = 0;
};

void validate_config(const SimulationConfig& config);
std::vector<Particle> make_particles(const SimulationConfig& config);
void run_simulation(
    std::vector<Particle>& particles,
    const SimulationConfig& config,
    const std::function<void(std::size_t)>& on_step = {});

}  // namespace barnes_hut
