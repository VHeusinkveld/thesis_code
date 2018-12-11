{

#line 517 "/home/vinlinux/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (val(c,0,0,0), n);
      _area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);
    }