#include <cmath>      // header for math functions
#include <vector>     // header for vector functions
#include <print>      // header for terminal output functions
#include <algorithm>  // header for max function

int main() {
    double xi = 0.0;                                    // intial x
    double xf = 1.0;                                    // final x
    double x0 = 0.25;                                   // ic
    double d = 0.1;                                     // ic

    int M = 32;                                         // number of physical cells
    double dx = (xf - xi) / M;                          // grid spacing

    double a = 1.0;                                     // wave speed constant
    double T = 0.5 / a;                                 // simulation duration
    std::vector<double> v = {1.0, 0.75, 0.5, 0.25};     // cfl numbers

    std::vector<double> u(M + 2);                       // empty u at t = n
    std::vector<double> unew(u.size());                 // u at t = n + 1
    std::vector<double> uexact(unew.size());            // exact u

    for (int i = 1; i < uexact.size() - 1; i++) {           // loop to compute exact u
        double x = dx * (i - 1.0);
        double xp = x - a * T;                              // shifted position
        if (xp < 0.0) {xp += 1.0;}                          // periodic wrap
        if (xp > 1.0) {xp -= 1.0;}

        if (xp >= x0 - d && xp <= x0 + d) {
            uexact[i] = std::cos(M_PI * (xp - x0) / (2.0 * d))
                        * std::cos(M_PI * (xp - x0) / (2.0 * d));}
        else {uexact[i] = 0.0;}
    }

    double t;                                               // time
    double dt;                                              // time step
    double L1norm;                                          // L1norm
    std::vector<double> du(3);                              // du at (j + 1/2, j - 1/2, j - 3/2) 
    std::vector<double> B(2);                               // flux terms

    for (int i = 0; i < v.size(); i++) {                    // loop to compute for all v and q

        std::print("v = {:.2f}\n", v[i]);                   // print v to terminal

        dt = v[i] * dx / a;                                 // time step

        for (int x = 1; x < u.size() - 1; x++) {    // u with ics at t = 0
            if (dx * (x - 1.0) >= x0 - d && dx * (x - 1.0) <= x0 + d) 
                {u[x] = std::cos(M_PI * (dx * (x - 1.0) - x0) / (2 * d)) 
                * std::cos(M_PI * (dx * (x - 1.0) - x0) / (2 * d));}
            else {u[x] = 0.0;}
        }


        u[0] = 1.0;                                         // dirichlet bc
        u[u.size() - 1] = u[u.size() - 2];                  // neuman bc

        t = 0.0;

        while (t + dt <= T) {                               // loop to update u with bcs until t = T
            t += dt;
            for (int j = 1; j < u.size() - 1; j++) {
                du[0] = u[j + 1] - u[j];
                du[1] = u[j] - u[j - 1];
                if (j == 1) {du[2] = du[1];}
                else {du[2] = u[j - 1] - u[j - 2];}

                for (int b = 0; b < B.size(); b++) {        // loop to compute B terms
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

            L1norm = 0.0;
            for (int l = 1; l < uexact.size() - 1; l++) {    // loop to compute L1norm
                L1norm += std::abs(unew[l] - uexact[l]);
            }
            L1norm = L1norm / M;

            std::print("L1norm = {:.6f}\n", L1norm);
            // print metrics to terminal
    }
    return 0;
}