#include "navier-stokes/centered.h"

#define PROF (sqrt(y + 0.1) - sqrt(0.1))

u.t[bottom] = dirichlet (0.);

int main() {
  periodic (left);
  const face vector muc[] = {1./10., 1./10.};
  mu = muc;
  run();
}

event init (t = 0) {
  foreach()
    u.x [] = PROF;
}

event forcing (i++) {
  double frac_force = 0.8, tau = 0.5;
  foreach()
    if (fabs (x - (X0 + L0/2))/L0 > frac_force/2.)
      u.x[] += dt*(PROF - u.x[])/tau;
}

event movies (t += 0.1; t < 10) {
  output_ppm (u.y, file = "uy.mp4", n = 256);
  output_ppm (u.x, file = "ux.mp4", n = 256,
	      min = 0, max = 1);
  scalar omega[];
  vorticity (u, omega);
  output_ppm (omega, file = "w.mp4", n = 256);
}
