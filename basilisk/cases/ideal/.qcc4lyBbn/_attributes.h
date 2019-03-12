typedef struct {

#line 932 "/home/vinlinux/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar);
  double (** boundary_homogeneous) (Point, Point, scalar);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;


    int z;

  } d;
  vector v;
  bool face, nodump;

#line 17 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

#line 8 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  void (* refine) (Point, scalar);

#line 94 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  void (* coarsen) (Point, scalar);

#line 76 "/home/vinlinux/basilisk/src/fractions.h"

  vector n;
