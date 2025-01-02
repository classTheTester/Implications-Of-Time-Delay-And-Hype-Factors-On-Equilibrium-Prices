#include <iostream>
#include <cmath>
#include <vector>
#include <sstream> // for std::ostringstream
#include "cpgplot.h"

using namespace std;

// Define parameters for the system
const double a = 1.5; // Adjusted price adjustment speed
const double b = 0.8; // Adjusted demand elasticity
const double d = 0.5; // Baseline demand
const double e = 0.3; // Price elasticity of demand
const double f = 0.7; // Supply response to price
const double g = 0.6; // Decay of supply

// Function for the derivative of P(t), incorporating hype
double dP_dt(double P, double P_delayed, double hype) {
    return a * (1.0 - P - b * P_delayed) + 1.0 + hype * sin(P); // "hype" adds nonlinearity
}

// Demand dynamics function
double dQ_d_dt(double P) {
    return d - e * P;
}

// Supply response function
double dQ_s_dt(double P_delayed, double Q_s) {
    return f * P_delayed - g * Q_s;
}

// Runge-Kutta 4th order method
void runge_kutta(double h, int steps, double tau, double hype, vector<double>& time, vector<double>& P) {
    double P_current = 1.0; // Initial condition for P (shifted up)
    double P_delayed = 1.0; // Initial delayed value (shifted up)
    double Q_s_current = 1.0; // Initial condition for Q_s
    vector<double> delay_buffer(static_cast<int>(tau / h) + 1, P_current); // Buffer for delayed values

    for (int i = 0; i < steps; ++i) {
        // Store time and state
        time.push_back(i * h);
        P.push_back(P_current);

        // Compute delayed state
        P_delayed = delay_buffer.front(); // Oldest value in the buffer

        // Runge-Kutta calculations for P
        double k1 = h * dP_dt(P_current, P_delayed, hype);
        double k2 = h * dP_dt(P_current + k1 / 2.0, P_delayed, hype);
        double k3 = h * dP_dt(P_current + k2 / 2.0, P_delayed, hype);
        double k4 = h * dP_dt(P_current + k3, P_delayed, hype);

        double P_next = P_current + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

        // Runge-Kutta calculations for Q_s
        double k1_Qs = h * dQ_s_dt(P_delayed, Q_s_current);
        double k2_Qs = h * dQ_s_dt(P_delayed, Q_s_current + k1_Qs / 2.0);
        double k3_Qs = h * dQ_s_dt(P_delayed, Q_s_current + k2_Qs / 2.0);
        double k4_Qs = h * dQ_s_dt(P_delayed, Q_s_current + k3_Qs);

        double Q_s_next = Q_s_current + (k1_Qs + 2 * k2_Qs + 2 * k3_Qs + k4_Qs) / 6.0;

        // Update buffer and state
        delay_buffer.push_back(P_next);
        delay_buffer.erase(delay_buffer.begin());

        P_current = P_next;
        Q_s_current = Q_s_next;
    }
}

// Plot results using CPGPLOT
void plot_results(const vector<double>& time, const vector<double>& P, double tau, double hype) {
    if (!cpgopen("/XWINDOW")) {
        cerr << "Failed to open PGPLOT window." << endl;
        return;
    }

    // Convert double vectors to float arrays
    vector<float> time_float(time.begin(), time.end());
    vector<float> P_float(P.begin(), P.end());

    // Set up the plot environment
    cpgenv(0, time_float.back(), 0.5, 1.5, 0, 0);
    cpglab("Time", "P(t)", "Commodity Market Model - Price Dynamics with Hype");

    // Plot the data
    cpgsci(2); // Change color
    cpgline(time_float.size(), time_float.data(), P_float.data());

    // Add legend at the bottom left
    cpgsci(1); // Reset to default color for text
    stringstream ss;
    ss << "tau = " << tau << ", hype = " << hype;
    cpgptxt(0.1, 0.55, 0, 0.0, ss.str().c_str());

    cpgclos();
}

int main() {
    // Time parameters
    double h = 0.1; // Step size
    int steps = 1000; // Reduced number of steps

    // User input for tau and hype
    double tau, hype;
    cout << "Enter the value for tau: ";
    cin >> tau;
    cout << "Enter the value for hype: ";
    cin >> hype;

    // Vectors to store time and state
    vector<double> time;
    vector<double> P;

    // Run simulation
    runge_kutta(h, steps, tau, hype, time, P);

    // Plot results
    plot_results(time, P, tau, hype);

    return 0;
}
