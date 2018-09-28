#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"
#include "view.h"


int main()
{
  init_grid(128);
  double L0 = 10.
  origin(0, 0, 0);
  
  double R = 0.5;
  double w = 0.1;
  
  double xf = L0/2.;
  double yf = L0/2.;
  double zf = L0/2.;
  scalar fan[], sph[], plane[];
  fan.prolongation = fraction_refine;
  view(theta = 0.1, phi = 0.3); 

  for (double rTheta = 0.,  double rPhi = 0.; rTheta < M_PI; rTheta += M_PI/12, rPhi += M_PI/6)
  {
    double[3] rn = {sin(rTheta)*cos(-rPhi + M_PI/2.), sin(rTheta)*sin(-rPhi + M_PI/2.), np.cos(rTheta)}

    fraction(sph, sq(R) - sq(rn[0]*(x - xf)) sq(rn[1]*(y - yf)) - sq(rn[2]*(z - zf)));
    fraction(plane, rn[0]*x + rn[1]*y + rn[2]*z - w/2.);
    foreach ()
      fan[] = sph[] * plane[];
    fraction(plane, -rn[0]*x + -rn[1]*y + -rn[2]*z - w/2.);
    foreach ()
      fan[] *= plane[];
    clear();
    boundary({fan});
    cells(alpha = 0);
    draw_vof("fan");
    save("puck.mp4");
    adapt_wavelet({fan}, (double[]){0.01}, 3, 7);
  }
}
