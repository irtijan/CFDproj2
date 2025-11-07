#include <cmath>      // header for math functions
#include <vector>     // header for vector functions
#include <print>      // header for terminal output functions
#include <algorithm>  // header for max function

int main() {
    double xi = 0.0;                                    // intial x
    double xf = 1.0;                                    // final x

    int M = 32;                                         // number of physical cells
    double dx = (xf - xi) / M;                          // grid spacing

    double a = 1.0;                                     // wave speed constant
    double T = 0.5 / a;                                 // simulation duration
    std::vector<double> v = {1.0, 0.75, 0.5, 0.25};     // cfl numbers

    std::vector<double> u(M + 2);                           // empty u at t = n
    std::vector<double> unew(u.size());                     // u at t = n + 1
    std::vector<double> uexact(unew.size());                // exact u

    double t;                                               // time
    double dt;                                              // time step
    double L1norm;                                          // L1norm
    std::vector<double> du(3);                              // du at (j + 1/2, j - 1/2, j - 3/2) 
    std::vector<double> B(2);                               // flux terms

    for (int i = 0; i < v.size(); i++) {                    // loop to compute for all v and q

        std::print("v = {:.2f}\n", v[i]);                   // print v to terminal

        dt = v[i] * dx / a;                                 // time step

        for (int x = 1; x < u.size() - 1; x++) {            // u with ics at t = 0
            if (x - 1 <= M / 4) {u[x] = 1.0;}
            if (x - 1 > M / 4) {u[x] = 0.0;}
        }

        u[0] = 1.0;                                         // dirichlet bc
        u[u.size() - 1] = u[u.size() - 2];                  // neuman bc

        t = 0.0;

        while (t + dt <= T) {                           // loop to update u with bcs until t = T
            t += dt;
            for (int j = 1; j < u.size() - 1; j++) {
                du[0] = u[j + 1] - u[j];
                du[1] = u[j] - u[j - 1];
                if (j == 1) {du[2] = du[1];}
                else {du[2] = u[j - 1] - u[j - 2];}

                for (int b = 0; b < B.size(); b++) {    // loop to compute B terms
                    if ((du[b] > 0.0 ? 1.0 : -1.0) != (du[b + 1] > 0.0 ? 1.0 : -1.0)) {B[b] = 0.0;}
                    else if (du[b] != 0.0 && du[b + 1] != 0.0 &&
                        du[b] / du[b + 1] >= 0.5 && du[b] / du[b + 1] <= 2.0) 
                    {B[b] = (du[b] > 0.0 ? 1.0 : -1.0) * std::max(std::abs(du[b]), std::abs(du[b + 1]));}
                    else {B[b] = 2.0 * (du[b] > 0.0 ? 1.0 : -1.0)
                                     * std::min(std::abs(du[b]), std::abs(du[b + 1]));}
                }

                unew[j] = u[j] - v[i] * ((u[j] - u[j - 1]) + (1.0 - v[i]) / 2.0 * (B[0] - B[1]));
            }
            unew[0] = 1.0;
            unew[unew.size() - 1] = unew[unew.size() - 2];
            u = unew;
            }

            for (int d = 1; d < uexact.size() - 1; d++) {           // loop to compute exact u
                if (dx * (d - 1) - a * T <= 0.25) {uexact[d] = 1.0;}
                else {uexact[d] = 0.0;}
            }

            L1norm = 0.0;
            for (int l = 1; l < uexact.size() - 1; l++) {           // loop to compute L1norm
                L1norm += std::abs(unew[l] - uexact[l]);
            }
            L1norm = L1norm / M;

            std::print("L1norm = {:.6f}\n", L1norm);
            // print metrics to terminal
    }
    return 0;
}