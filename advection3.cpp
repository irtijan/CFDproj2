#include <cmath>    // header for math functions
#include <vector>   // header for vector functions
#include <print>    // header for terminal output functions

int main() {
    double xi = 0.0;                                    // intial x
    double xf = 1.0;                                    // final x

    int M = 32;                                         // number of physical cells
    double dx = (xf - xi) / M;                          // grid spacing

    double a = 1.0;                                     // wave speed constant
    std::vector<double> T {0.5 / a, 1.0 / a};           // simulation duration
    std::vector<double> v = {1.0, 0.75, 0.5, 0.25};     // cfl numbers

    std::vector<double> u(M + 2);                       // empty u at t = n
    std::vector<double> unew(u.size());                 // u at t = n + 1
    std::vector<double> uexact(unew.size());            // exact u

    std::vector<double> q(4);                           // diffusion coeffs
    double t;                                           // time
    double dt;                                          // time step
    double L1norm;                                      // L1norm

    for (int i = 0; i < v.size(); i++) {                // loop to compute for all v, q, aT

        std::print("v = {:.2f}\n", v[i]);               // print v to terminal

        dt = v[i] * dx / a;                             // time step

        q = {1,
            std::abs(v[i]),
            v[i] * v[i],
            1.0 / 3.0 + 2.0 / 3.0 * v[i] * v[i]};
                      
        for (int k = 0; k < q.size(); k++) {

            for (int s = 0; s < T.size(); s++) {

                for (int x = 1; x < u.size() - 1; x++) {    // u with ics at t = 0
                    if (x - 1 <= M / 4) {u[x] = 1.0;}
                    if (x - 1 > M / 4) {u[x] = 0.0;}
                }

                u[0] = 1.0;                                 // dirichlet bc
                u[u.size() - 1] = u[u.size() - 2];          // neuman bc

                t = 0.0;

                while (t + dt <= T[s]) {                    // loop to update u with bcs until t = T
                    t += dt;
                    for (int j = 1; j < u.size() - 1; j++) {
                        unew[j] = 0.5 * (q[k] + v[i]) * u[j - 1] + (1.0 - q[k]) * u[j] + 0.5 * (q[k] - v[i]) * u[j + 1];
                    }
                    unew[0] = 1.0;
                    unew[unew.size() - 1] = unew[unew.size() - 2];
                    u = unew;
                }

                for (int d = 1; d < uexact.size() - 1; d++) {           // loop to compute exact u
                    if (dx * (d - 1) - a * T[s] <= 0.25) {uexact[d] = 1.0;}
                    else {uexact[d] = 0.0;}
                }

                L1norm = 0.0;
                for (int l = 1; l < uexact.size() - 1; l++) {           // loop to compute L1norm
                    L1norm += std::abs(unew[l] - uexact[l]);
                }
                L1norm = L1norm / M;

                std::print("q = {:.3f}, aT = {:.3f}, L1norm = {:.6f}\n", q[k],  a * T[s], L1norm);
                // print metrics to terminal
            }
        }
    }
    return 0;
}