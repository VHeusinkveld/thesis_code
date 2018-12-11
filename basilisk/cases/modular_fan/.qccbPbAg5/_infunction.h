
  coord n;
  if (s.x.i < 0)
    n = mycs (point, c);
  else {
    double nn = 0.;
    {
#line 379
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 379
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
#line 379
 {
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }}
    assert (nn > 0.);
    {
#line 384

      n.x /= nn;
#line 384

      n.y /= nn;
#line 384

      n.z /= nn;}
  }
  return n;
