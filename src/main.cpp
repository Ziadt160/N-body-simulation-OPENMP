#include "barnes_hut/simulation.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

void print_help(const char* executable) {
    std::cout
        << "Usage: " << executable << " [options]\n"
        << "Barnes-Hut 2D gravitational N-body simulation.\n\n"
        << "Options:\n"
        << "  --particles N  Number of particles (default: 10000)\n"
        << "  --steps N      Number of time steps (default: 100)\n"
        << "  --dt VALUE     Time step (default: 0.01)\n"
        << "  --theta VALUE  Barnes-Hut accuracy threshold (default: 0.5)\n"
        << "  --seed N       Random seed (default: 42)\n"
        << "  --threads N    OpenMP threads; 0 uses the runtime default\n"
        << "  --help         Show this message\n";
}

std::size_t parse_size(const std::string& value, const std::string& option) {
    std::size_t consumed = 0;
    const unsigned long long parsed = std::stoull(value, &consumed);
    if (consumed != value.size() || (!value.empty() && value.front() == '-') ||
        parsed > std::numeric_limits<std::size_t>::max()) {
        throw std::invalid_argument("invalid value for " + option + ": " + value);
    }
    return static_cast<std::size_t>(parsed);
}

std::uint64_t parse_seed(const std::string& value) {
    std::size_t consumed = 0;
    const unsigned long long parsed = std::stoull(value, &consumed);
    if (consumed != value.size() || (!value.empty() && value.front() == '-')) {
        throw std::invalid_argument("invalid value for --seed: " + value);
    }
    return static_cast<std::uint64_t>(parsed);
}

double parse_double(const std::string& value, const std::string& option) {
    std::size_t consumed = 0;
    const double parsed = std::stod(value, &consumed);
    if (consumed != value.size() || !std::isfinite(parsed)) {
        throw std::invalid_argument("invalid value for " + option + ": " + value);
    }
    return parsed;
}

barnes_hut::SimulationConfig parse_arguments(int argc, char* argv[]) {
    barnes_hut::SimulationConfig config;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help") {
            print_help(argv[0]);
            std::exit(0);
        }
        if (i + 1 >= argc) {
            throw std::invalid_argument("missing value for " + option);
        }
        const std::string value = argv[++i];
        if (option == "--particles") {
            config.particle_count = parse_size(value, option);
        } else if (option == "--steps") {
            config.steps = parse_size(value, option);
        } else if (option == "--dt") {
            config.time_step = parse_double(value, option);
        } else if (option == "--theta") {
            config.theta = parse_double(value, option);
        } else if (option == "--seed") {
            config.seed = parse_seed(value);
        } else if (option == "--threads") {
            const std::size_t threads = parse_size(value, option);
            if (threads > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
                throw std::invalid_argument("thread count is too large");
            }
            config.threads = static_cast<int>(threads);
        } else {
            throw std::invalid_argument("unknown option: " + option);
        }
    }
    return config;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const barnes_hut::SimulationConfig config = parse_arguments(argc, argv);
        barnes_hut::validate_config(config);
        std::vector<barnes_hut::Particle> particles = barnes_hut::make_particles(config);

        std::cout << "Barnes-Hut 2D simulation\n"
                  << "Particles: " << config.particle_count
                  << ", steps: " << config.steps
                  << ", dt: " << config.time_step
                  << ", theta: " << config.theta
                  << ", seed: " << config.seed
                  << ", threads: "
                  << (config.threads == 0 ? std::string("runtime default")
                                          : std::to_string(config.threads))
                  << '\n';

        const auto start = std::chrono::steady_clock::now();
        barnes_hut::run_simulation(particles, config, [](std::size_t completed_steps) {
            if (completed_steps % 10 == 0) {
                std::cout << "Completed step " << completed_steps << '\n';
            }
        });
        const std::chrono::duration<double> elapsed =
            std::chrono::steady_clock::now() - start;
        std::cout << "Simulation completed in " << elapsed.count() << " seconds.\n";
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n'
                  << "Run with --help for usage.\n";
        return 1;
    }
    return 0;
}
