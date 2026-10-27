#include "barnes_hut/simulation.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace barnes_hut {
namespace {

constexpr std::size_t leaf_capacity = 4;
constexpr int maximum_depth = 50;

struct BoundingBox {
    double x_min;
    double x_max;
    double y_min;
    double y_max;

    bool contains(const Particle& particle) const noexcept {
        return particle.x >= x_min && particle.x <= x_max &&
               particle.y >= y_min && particle.y <= y_max;
    }
};

BoundingBox make_root_box(const std::vector<Particle>& particles) {
    if (particles.empty()) {
        throw std::invalid_argument("the particle collection must not be empty");
    }

    BoundingBox box{
        particles.front().x, particles.front().x,
        particles.front().y, particles.front().y
    };
    for (const Particle& particle : particles) {
        if (!std::isfinite(particle.x) || !std::isfinite(particle.y) ||
            !std::isfinite(particle.vx) || !std::isfinite(particle.vy) ||
            !std::isfinite(particle.mass) || particle.mass <= 0.0) {
            throw std::invalid_argument(
                "particles must have finite position/velocity and positive finite mass");
        }
        box.x_min = std::min(box.x_min, particle.x);
        box.x_max = std::max(box.x_max, particle.x);
        box.y_min = std::min(box.y_min, particle.y);
        box.y_max = std::max(box.y_max, particle.y);
    }

    const double center_x = (box.x_min + box.x_max) / 2.0;
    const double center_y = (box.y_min + box.y_max) / 2.0;
    const double extent = std::max(box.x_max - box.x_min, box.y_max - box.y_min);
    const double half_size = std::max(extent / 2.0, 1.0) + 1.0;
    return {center_x - half_size, center_x + half_size,
            center_y - half_size, center_y + half_size};
}

std::pair<double, double> pair_force(const Particle& target, const Particle& source) {
    const double dx = source.x - target.x;
    const double dy = source.y - target.y;
    const double distance_squared = dx * dx + dy * dy + softening;
    const double distance = std::sqrt(distance_squared);
    const double magnitude = gravitational_constant * target.mass * source.mass /
                             distance_squared;
    return {magnitude * dx / distance, magnitude * dy / distance};
}

}  // namespace

struct QuadTree::Node {
    explicit Node(BoundingBox node_box) : box(node_box) {}

    void subdivide() {
        const double mid_x = (box.x_min + box.x_max) / 2.0;
        const double mid_y = (box.y_min + box.y_max) / 2.0;
        children[0] = std::make_unique<Node>(BoundingBox{box.x_min, mid_x, box.y_min, mid_y});
        children[1] = std::make_unique<Node>(BoundingBox{mid_x, box.x_max, box.y_min, mid_y});
        children[2] = std::make_unique<Node>(BoundingBox{box.x_min, mid_x, mid_y, box.y_max});
        children[3] = std::make_unique<Node>(BoundingBox{mid_x, box.x_max, mid_y, box.y_max});
        is_leaf = false;
    }

    int quadrant(const Particle& particle) const noexcept {
        const double mid_x = (box.x_min + box.x_max) / 2.0;
        const double mid_y = (box.y_min + box.y_max) / 2.0;
        return (particle.x >= mid_x ? 1 : 0) | (particle.y >= mid_y ? 2 : 0);
    }

    void insert(const Particle* particle, int depth) {
        if (is_leaf) {
            if (particles.size() < leaf_capacity || depth >= maximum_depth) {
                particles.push_back(particle);
                return;
            }

            std::vector<const Particle*> existing_particles = std::move(particles);
            subdivide();
            for (const Particle* existing : existing_particles) {
                children[quadrant(*existing)]->insert(existing, depth + 1);
            }
        }
        children[quadrant(*particle)]->insert(particle, depth + 1);
    }

    void compute_mass_distribution() {
        mass = 0.0;
        com_x = 0.0;
        com_y = 0.0;

        if (is_leaf) {
            for (const Particle* particle : particles) {
                mass += particle->mass;
                com_x += particle->x * particle->mass;
                com_y += particle->y * particle->mass;
            }
        } else {
            for (const auto& child : children) {
                if (!child) {
                    continue;
                }
                child->compute_mass_distribution();
                mass += child->mass;
                com_x += child->com_x * child->mass;
                com_y += child->com_y * child->mass;
            }
        }

        if (mass > 0.0) {
            com_x /= mass;
            com_y /= mass;
        }
    }

    void accumulate_force(const Particle& target, double theta,
                          double& force_x, double& force_y) const {
        if (mass == 0.0) {
            return;
        }

        if (is_leaf) {
            for (const Particle* source : particles) {
                if (source == &target) {
                    continue;
                }
                const auto force = pair_force(target, *source);
                force_x += force.first;
                force_y += force.second;
            }
            return;
        }

        const double dx = com_x - target.x;
        const double dy = com_y - target.y;
        const double distance = std::sqrt(dx * dx + dy * dy);
        const double width = box.x_max - box.x_min;
        if (!box.contains(target) && distance > 0.0 && width / distance < theta) {
            const double distance_squared = dx * dx + dy * dy + softening;
            const double softened_distance = std::sqrt(distance_squared);
            const double magnitude = gravitational_constant * target.mass * mass /
                                     distance_squared;
            force_x += magnitude * dx / softened_distance;
            force_y += magnitude * dy / softened_distance;
            return;
        }

        for (const auto& child : children) {
            if (child) {
                child->accumulate_force(target, theta, force_x, force_y);
            }
        }
    }

    BoundingBox box;
    bool is_leaf = true;
    std::vector<const Particle*> particles;
    std::array<std::unique_ptr<Node>, 4> children{};
    double mass = 0.0;
    double com_x = 0.0;
    double com_y = 0.0;
};

QuadTree::QuadTree(const std::vector<Particle>& particles, double theta)
    : theta_(theta), particle_count_(particles.size()) {
    if (!std::isfinite(theta) || theta < 0.0) {
        throw std::invalid_argument("theta must be finite and non-negative");
    }
    root_ = std::make_unique<Node>(make_root_box(particles));
    for (const Particle& particle : particles) {
        root_->insert(&particle, 0);
    }
    root_->compute_mass_distribution();
}

QuadTree::~QuadTree() = default;
QuadTree::QuadTree(QuadTree&&) noexcept = default;
QuadTree& QuadTree::operator=(QuadTree&&) noexcept = default;

std::pair<double, double> QuadTree::calculate_force(const Particle& particle) const {
    double force_x = 0.0;
    double force_y = 0.0;
    root_->accumulate_force(particle, theta_, force_x, force_y);
    return {force_x, force_y};
}

std::size_t QuadTree::particle_count() const noexcept {
    return particle_count_;
}

double QuadTree::total_mass() const noexcept {
    return root_->mass;
}

void validate_config(const SimulationConfig& config) {
    if (config.particle_count == 0) {
        throw std::invalid_argument("particle count must be greater than zero");
    }
    if (!std::isfinite(config.time_step) || config.time_step <= 0.0) {
        throw std::invalid_argument("time step must be finite and greater than zero");
    }
    if (!std::isfinite(config.theta) || config.theta < 0.0) {
        throw std::invalid_argument("theta must be finite and non-negative");
    }
    if (config.threads < 0) {
        throw std::invalid_argument("thread count must not be negative");
    }
}

std::vector<Particle> make_particles(const SimulationConfig& config) {
    validate_config(config);
    std::mt19937_64 generator(config.seed);
    std::uniform_real_distribution<double> position_distribution(0.0, 1000.0);
    std::uniform_real_distribution<double> velocity_distribution(-1.0, 1.0);
    std::uniform_real_distribution<double> mass_distribution(1.0, 11.0);

    std::vector<Particle> particles(config.particle_count);
    for (Particle& particle : particles) {
        particle.x = position_distribution(generator);
        particle.y = position_distribution(generator);
        particle.vx = velocity_distribution(generator);
        particle.vy = velocity_distribution(generator);
        particle.mass = mass_distribution(generator);
    }
    return particles;
}

void run_simulation(std::vector<Particle>& particles, const SimulationConfig& config,
                    const std::function<void(std::size_t)>& on_step) {
    validate_config(config);
    if (particles.empty()) {
        throw std::invalid_argument("the particle collection must not be empty");
    }

#ifdef _OPENMP
    if (config.threads > 0) {
        omp_set_num_threads(config.threads);
    }
#endif

    for (std::size_t step = 0; step < config.steps; ++step) {
        const QuadTree tree(particles, config.theta);

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 32)
#endif
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(particles.size()); ++i) {
            const auto force = tree.calculate_force(particles[static_cast<std::size_t>(i)]);
            particles[static_cast<std::size_t>(i)].fx = force.first;
            particles[static_cast<std::size_t>(i)].fy = force.second;
        }

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (std::ptrdiff_t i = 0; i < static_cast<std::ptrdiff_t>(particles.size()); ++i) {
            Particle& particle = particles[static_cast<std::size_t>(i)];
            particle.vx += (particle.fx / particle.mass) * config.time_step;
            particle.vy += (particle.fy / particle.mass) * config.time_step;
            particle.x += particle.vx * config.time_step;
            particle.y += particle.vy * config.time_step;
        }

        if (on_step) {
            on_step(step + 1);
        }
    }
}

}  // namespace barnes_hut
