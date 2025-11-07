#include <cmath>    // header for math functions
#include <vector>   // header for vector functions
#include <print>    // header for terminal output functions

int main() {
    double xi = 0.0;                                    // intial x
    double xf = 1.0;                                    // final x

    int M = 32;                                         // number of physical cells
    double dx = (xf - xi) / M;                          // grid spacing

    double a = 1.0;                                     // wave speed constant
    double T = 3.0 / (2.0 * a);                         // simulation duration
    std::vector<double> v = {1.0, 0.75, 0.5, 0.25};     // cfl numbers

    std::vector<double> u(M + 2);                           // empty u at t = n
    std::vector<double> unew(u.size());                     // u at t = n + 1
    std::vector<double> uexact(unew.size());                // exact u

    double Aexact = 1.0;                                    // exact amp
    double phiexact = -2.0 * M_PI * a * T;                  // exact phase

    std::vector<double> q(4);                               // diffusion coeffs
    double t;                                               // time
    double dt;                                              // time step
    double Re;                                              // real part amp
    double Im;                                              // imag part amp
    double Anum;                                            // numerical amp
    double phinum;                                          // numerical phase
    double amp_error;                                       // amp error
    double phase_error;                                     // phase error
    double L1norm;                                          // L1norm

    for (int i = 1; i < uexact.size() - 1; i++) {           // loop to compute exact u
        uexact[i] = std::sin(2.0 * M_PI * (dx * (i - 1.0) - a * T));
    }

    for (int i = 0; i < v.size(); i++) {                    // loop to compute for all v and q

        std::print("v = {:.2f}\n", v[i]);                   // print v to terminal

        dt = v[i] * dx / a;                                 // time step

        q = {1,
            std::abs(v[i]),
            v[i] * v[i],
            1.0 / 3.0 + 2.0 / 3.0 * v[i] * v[i]};
                      
        for (int k = 0; k < q.size(); k++) {

            for (int x = 1; x < u.size() - 1; x++) {            // u with ics at t = 0
                u[x] = std::sin(2.0 * M_PI * dx * (x - 1.0));
            }

            u[0] = u[u.size() - 2];                             // periodic bcs
            u[u.size() - 1] = u[1];

            t = 0.0;

            while (t + dt <= T) {                               // loop to update u with bcs until t = T
                t += dt;
                for (int j = 1; j < u.size() - 1; j++) {
                    unew[j] = 0.5 * (q[k] + v[i]) * u[j - 1] + (1.0 - q[k]) * u[j] + 0.5 * (q[k] - v[i]) * u[j + 1];
                }
                unew[0] = unew[unew.size() - 2];
                unew[unew.size() - 1] = unew[1];
                u = unew;
            }
            
            Re = 0.0;
            Im = 0.0;

            for (int r = 1; r < unew.size() - 1; r++) {          // loop to compute amp parts
                Re += unew[r] * std::cos(2.0 * M_PI * dx * (r - 1));
                Im -= unew[r] * std::sin(2.0 * M_PI * dx * (r - 1));
            }

            Anum = 2.0 / M * std::sqrt(Re * Re + Im * Im);
            phinum = std::atan2(Im, Re);

            amp_error = std::abs(Anum - Aexact);
            phase_error = std::abs(phinum - phiexact);
            
            while (phase_error > M_PI) {phase_error -= M_PI;}     // phase error wrapped to [-pi, pi]
            while (phase_error < -M_PI) {phase_error += M_PI;}
            phase_error = std::abs(phase_error - M_PI / 2.0);     // initial phase removed

            L1norm = 0.0;
            for (int l = 1; l < uexact.size() - 1; l++) {         // loop to compute L1norm
                L1norm += std::abs(unew[l] - uexact[l]);
            }
            L1norm = L1norm / M;

            std::print("q = {:.3f}, Amplitude error = {:.6f}, Phase error = {:.6f}, L1norm = {:.6f}\n", 
                q[k], amp_error, phase_error, L1norm); // print all metrics to terminal
        }
    }
    return 0;
}