#include <vector>
#include <iostream>
#include <cmath>
#include "writer.hpp"

double f(double t, double u) {
    return std::exp(-2*t) - 2*u;
}


/// Uses the SSP RK3 method to compute u from time 0 to time T 
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

//----------------SSPRK3Begin----------------
void SSPRK3(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {
    const unsigned int nsteps = std::round(T/dt);
    u.resize(nsteps+1);
    time.resize(nsteps+1);

    // setup
    double h = dt;

    // save inital values
    u[0] = u0;
    time[0] = 0;

    // Do the iteration for uk+1 and use uk
    for (int i = 0; i < nsteps; i++) {
        time[i] = h*i;
        double k1 = f(time[i],u[i]);
        double k2 = f(time[i]+1,u[i]+h*k1);
        double k3 = f(time[i]+0.5*h,u[i]+h*(0.25*k1 + 0.25*k2));

        u[i+1] = u[i]+h*(1.0/6.0*k1 + 1.0/6.0*k2 + 2.0/3.0*k3);
    }
}
//----------------SSPRK3End----------------

int main(int argc, char** argv) {

    double T = 10.0;

    double dt = 0.5;

    // To make some plotting easier, we take the dt parameter in as an optional
    // parameter.
    if (argc == 2) {
        dt = atof(argv[1]);
    }

    const double u0 = 0.;
    std::vector<double> time;
    std::vector<double> u;
    SSPRK3(u,time,u0,dt,T);

    writeToFile("solution.txt", u);
    writeToFile("time.txt",time);

    return 0;
}
