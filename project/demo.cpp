#include <iostream>
#include <cmath>
#include <vector>
#include <sstream> // for std::ostringstream
#include "cpgplot.h"

using namespace std;

// Define parameters for the system
const double a = 1.5;   // Adjusted price adjustment speed
const double b = 0.8;   // Adjusted demand elasticity

// Function for the derivative of P(t)
double dP_dt(double P, double P_delayed) {
    return a * (1.0 - P - b * P_delayed);
}

// Runge-Kutta 4th order method
void runge_kutta(double h, int steps, double tau, vector<double>& time, vector<double>& P) {
    double P_current = 1.0; // Initial condition for P
    double P_delayed = 1.0; // Initial delayed value
    vector<double> delay_buffer(static_cast<int>(tau / h) + 1, P_current); // Buffer for delayed values

    for (int i = 0; i < steps; ++i) {
        // Store time and state
        time.push_back(i * h);
        P.push_back(P_current);

        // Compute delayed state
        P_delayed = delay_buffer.front(); // Oldest value in the buffer

        // Runge-Kutta calculations
        double k1 = h * dP_dt(P_current, P_delayed);
        double k2 = h * dP_dt(P_current + k1 / 2.0, P_delayed);
        double k3 = h * dP_dt(P_current + k2 / 2.0, P_delayed);
        double k4 = h * dP_dt(P_current + k3, P_delayed);

        double P_next = P_current + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

        // Update buffer and state
        delay_buffer.push_back(P_next);
        delay_buffer.erase(delay_buffer.begin());

        P_current = P_next;
    }
}

// Plot results using CPGPLOT
void plot_results(const vector<vector<double> >& time_list, const vector<vector<double> >& P_list, const vector<double>& taus) {
    if (!cpgopen("/XWINDOW")) {
        cerr << "Failed to open PGPLOT window." << endl;
        return;
    }

    // Set up the plot environment
    cpgenv(0, time_list[0].back(), 0.5, 1.5, 0, 0);
    cpglab("Time", "P(t)", "Commodity Market Model - Price Dynamics");

    // Convert double vectors to float arrays and plot them with different colors
    int colors[] = {2, 3, 4, 5}; // PGPLOT color indices
    for (size_t i = 0; i < time_list.size(); ++i) {
        vector<float> time_float(time_list[i].begin(), time_list[i].end());
        vector<float> P_float(P_list[i].begin(), P_list[i].end());

        cpgsci(colors[i % 4]); // Change color
        cpgline(time_float.size(), time_float.data(), P_float.data());
    }

    // Add legend at the bottom left
    float legend_x = -25.0; // Position of the legend (more to the left)
    float legend_y = 0.55; // Initial y position of the legend
    for (size_t i = 0; i < taus.size(); ++i) {
        cpgsci(colors[i % 4]);
        stringstream ss;
        ss << "tau=" << taus[i];
        cpgptxt(legend_x, legend_y, 0, 0.0, ss.str().c_str());
        legend_y -= 0.05; // Move the next legend entry down
    }

    cpgclos();
}

int main() {
    // Time parameters
    double h = 0.1;  // Step size
    int steps = 3000; // Number of steps to observe long-term behavior

    // Different tau values to test
    double taus_array[] = {2, 4, 8, 16};
    vector<double> taus(taus_array, taus_array + sizeof(taus_array) / sizeof(taus_array[0]));

    // Vectors to store time and state for each tau
    vector<vector<double> > time_list;
    vector<vector<double> > P_list;

    // Run simulation for each tau
    for (size_t i = 0; i < taus.size(); ++i) {
        vector<double> time;
        vector<double> P;
        runge_kutta(h, steps, taus[i], time, P);
        time_list.push_back(time);
        P_list.push_back(P);
    }

    // Plot results
    plot_results(time_list, P_list, taus);

    return 0;
}

