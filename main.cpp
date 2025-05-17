#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <chrono>
#include <cmath>
#include <sstream>
#include <algorithm>

#define M_PI 3.14159265358979323846

#define BENCHMARK_ARG "--benchmark"
#define OUTPUT_ARG "--output"

namespace
{
    bool benchmark_mode = false;
    std::string output_filename = "results";
}

/// @brief Constants of the system.
struct Params
{
    double m1;
    double m2;
    double L;
    double g;
};

// Defining a ode system as a function which returns a vector of doubles as the results and input one double one vector of doubles and the parameters from the struct above.
using ode_system = std::function<std::vector<double>(double, const std::vector<double> &, const Params &)>;

/**
 * @brief
 * Defines the equation to be solved.
 * @param t Time
 * @param y Initial state of the system
 * @param p Params struct
 * @return std::vector<double> with the results states of the system
 */
auto ode(double t, const std::vector<double> &y, const Params &p) -> std::vector<double>
{
    double theta = y[0];
    double theta_p = y[1];
    double x2 = y[2];
    double x2_p = y[3];

    double under = (p.m1 + p.m2) * p.L - p.m1 * p.L * pow(cos(theta), 2);

    double over_theta = (p.m1 + p.m2) * p.g * sin(theta) - p.m1 * p.L * pow(theta_p, 2) * sin(theta) * cos(theta);
    double theta_pp = over_theta / under;

    double over_x2 = p.m1 * p.L * pow(theta_p, 2) * sin(theta) - p.m1 * p.g * sin(theta) * cos(theta);
    double x2_pp = over_x2 / under;

    return {theta_p, theta_pp, x2_p, x2_pp};
}

/**
 * @brief
 * Solve the ode equation using the RK4 method.
 * @param ode_system The function to solve the ode system
 * @param y0 Initial state of the system
 * @param t Current time
 * @param t0 Initial time
 * @param tf Final time
 * @param h Step size
 * @param p Params struct
 * @return std::vector<std::vector<double>> with the results states of the system
 */
auto rk4_solve(
    ode_system ode_system,
    const std::vector<double> &y0,
    const double &t0,
    const double &tf,
    const double &h,
    const Params &p)
{
    using solution_t = std::vector<std::vector<double>>;
    solution_t solution = {};
    std::vector<double> y = y0;
    double t = t0;
    int n_steps = static_cast<int>((tf - t0) / h);

    solution.reserve(n_steps + 1);
    solution.push_back({t0, y[0], y[1], y[2], y[3]});

    for (size_t i = 0; i < n_steps; i++)
    {
        // Setting up the RK4 method
        std::vector<double> k1 = ode_system(t, y, p);
        std::vector<double> k2 = ode_system(t + h / 2, {y[0] + h / 2 * k1[0], y[1] + h / 2 * k1[1], y[2] + h / 2 * k1[2], y[3] + h / 2 * k1[3]}, p);
        std::vector<double> k3 = ode_system(t + h / 2, {y[0] + h / 2 * k2[0], y[1] + h / 2 * k2[1], y[2] + h / 2 * k2[2], y[3] + h / 2 * k2[3]}, p);
        std::vector<double> k4 = ode_system(t + h, {y[0] + h * k3[0], y[1] + h * k3[1], y[2] + h * k3[2], y[3] + h * k3[3]}, p);

        for (size_t j = 0; j < y.size(); j++)
        {
            y[j] += h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
        }

        t += h;
        solution.push_back({t, y[0], y[1], y[2], y[3]});
    }
    return solution;
}
// Very simple argument parser, i was going to use some std::arrays and some iterators to make it more elegant but i dont think it will matter for only 2 arguments lol
auto static handle_args(int argc, char **argv, std::ostringstream &lasterror) -> bool
{
    if (argc < 2)
    {
        return 1;
    }

    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == BENCHMARK_ARG)
        {
            benchmark_mode = true;
        }
        else if (arg == OUTPUT_ARG)
        {
            if (i + 1 < argc)
            {
                output_filename = argv[i + 1];
                i++;
            }
            else
            {
                lasterror << "Error: " << OUTPUT_ARG << " requires a filename argument.";
                return 0;
            }
        }
        else
        {
            lasterror << "Error: Unknown argument " << arg;
            return 0;
        }
    }
    return 1;
}

auto main(int argc, char **argv) -> int
{
    std::ostringstream last_error{};
    if (!handle_args(argc, argv, last_error))
    {
        std::cerr << last_error.str() << std::endl
                  << "Usage: " << argv[0] << " [" << BENCHMARK_ARG << "] "
                  << "[" << OUTPUT_ARG << " <filename> ] " << std::endl;
        return 1;
    }

    Params p = {1.0, 2.0, 1.0, 9.81}; // m1, m2, L, g

    std::vector<double> y0 = {M_PI / 6, 0.0, 0.0, 0.0}; // Initial conditions: theta, theta', x2, x2'
    double t0 = 0.0;                                    // Initial time
    double tf = 10.0;                                   // Final time
    double h = 0.01;                                    // Step size

    // We measure anyways bc the time needed to do it is negligible, but we only print if args were passed kek
    auto start_time = std::chrono::high_resolution_clock::now();
    auto results = rk4_solve(ode, y0, t0, tf, h, p);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = (end_time - start_time) * 1000;

    std::ofstream file(output_filename + ".csv");

    if (file.is_open())
    {
        file << "t,theta,theta_p',x2,x2_p'\n";
        for (const auto &result : results)
        {
            file << result[0] << "," << result[1] << "," << result[2] << "," << result[3] << "," << result[4] << "\n";
        }
        file.close();
    }
    else
    {
        std::cerr << "Unable to open file";
    }

    if (benchmark_mode)
        std::cout << "Elapsed time: " << elapsed_time.count() << " ms" << std::endl;

    std::cout << "Results saved to " << output_filename << ".csv" << std::endl;

    return 0;
}
