#include "barnes_hut/simulation.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

void require_near(double actual, double expected, double tolerance,
                  const std::string& message) {
    if (std::abs(actual - expected) > tolerance) {
        throw std::runtime_error(message);
    }
}

std::pair<double, double> direct_force(const barnes_hut::Particle& target,
                                      const std::vector<barnes_hut::Particle>& particles) {
    double force_x = 0.0;
    double force_y = 0.0;
    for (const auto& source : particles) {
        if (&source == &target) {
            continue;
        }
        const double dx = source.x - target.x;
        const double dy = source.y - target.y;
        const double distance_squared = dx * dx + dy * dy + barnes_hut::softening;
        const double distance = std::sqrt(distance_squared);
        const double magnitude = barnes_hut::gravitational_constant *
                                 target.mass * source.mass / distance_squared;
        force_x += magnitude * dx / distance;
        force_y += magnitude * dy / distance;
    }
    return {force_x, force_y};
}

void test_force_matches_direct_sum_when_theta_is_zero() {
    std::vector<barnes_hut::Particle> particles{
        {-3.0, 2.0, 0.0, 0.0, 0.0, 0.0, 2.0},
        {4.0, -1.0, 0.0, 0.0, 0.0, 0.0, 3.0},
        {1.5, 6.0, 0.0, 0.0, 0.0, 0.0, 1.5},
        {-7.0, -5.0, 0.0, 0.0, 0.0, 0.0, 4.0},
        {9.0, 3.0, 0.0, 0.0, 0.0, 0.0, 2.5},
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.0}
    };
    const barnes_hut::QuadTree tree(particles, 0.0);

    for (const auto& particle : particles) {
        const auto actual = tree.calculate_force(particle);
        const auto expected = direct_force(particle, particles);
        require_near(actual.first, expected.first, 1e-12, "x force differs from direct sum");
        require_near(actual.second, expected.second, 1e-12, "y force differs from direct sum");
    }
}

void test_all_particles_are_retained_at_identical_positions() {
    std::vector<barnes_hut::Particle> particles(12);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        particles[i].x = 5.0;
        particles[i].y = -2.0;
        particles[i].mass = static_cast<double>(i + 1);
    }

    const barnes_hut::QuadTree tree(particles, 0.5);
    require(tree.particle_count() == particles.size(), "quadtree lost particles");
    require_near(tree.total_mass(), 78.0, 1e-12, "quadtree total mass is incorrect");
    for (const auto& particle : particles) {
        const auto force = tree.calculate_force(particle);
        require(std::isfinite(force.first) && std::isfinite(force.second),
                "force for coincident particles is not finite");
        require_near(force.first, 0.0, 1e-12, "coincident particles have nonzero x force");
        require_near(force.second, 0.0, 1e-12, "coincident particles have nonzero y force");
    }
}

void test_single_particle_has_no_self_force() {
    const std::vector<barnes_hut::Particle> particles{{2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 7.0}};
    const barnes_hut::QuadTree tree(particles, 0.5);
    const auto force = tree.calculate_force(particles.front());
    require_near(force.first, 0.0, 0.0, "single particle has x self-force");
    require_near(force.second, 0.0, 0.0, "single particle has y self-force");
}

void test_target_containing_nodes_are_opened() {
    std::vector<barnes_hut::Particle> particles;
    particles.push_back({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0});
    for (int i = 0; i < 4; ++i) {
        particles.push_back({-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0});
        particles.push_back({1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0});
    }

    const barnes_hut::QuadTree tree(particles, 1e6);
    const auto force = tree.calculate_force(particles.front());
    require_near(force.first, 0.0, 1e-12,
                 "symmetric sources produced a force from a target-containing node");
    require_near(force.second, 0.0, 1e-12,
                 "symmetric sources produced a y force from a target-containing node");
}

void test_simulation_configuration_is_validated() {
    barnes_hut::SimulationConfig config;
    config.time_step = 0.0;
    bool rejected = false;
    try {
        barnes_hut::validate_config(config);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "zero time step was accepted");
}

}  // namespace

int main() {
    try {
        test_force_matches_direct_sum_when_theta_is_zero();
        test_all_particles_are_retained_at_identical_positions();
        test_single_particle_has_no_self_force();
        test_target_containing_nodes_are_opened();
        test_simulation_configuration_is_validated();
    } catch (const std::exception& error) {
        std::cerr << "Test failure: " << error.what() << '\n';
        return 1;
    }
    std::cout << "All simulation tests passed.\n";
    return 0;
}
