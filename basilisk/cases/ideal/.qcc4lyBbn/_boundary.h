#line 79 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static double _boundary0 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann(val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 79
return  neumann_homogeneous(); } return 0.; }
#line 80 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static double _boundary1 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 79 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann(-val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0)); } return 0.; } static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 79 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 80
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 80
return  neumann_homogeneous(); } return 0.; }
#line 88 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static double _boundary2 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 87 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann(val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 87 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 88
return  neumann_homogeneous(); } return 0.; }
#line 89 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static double _boundary3 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann(-val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0)); } return 0.; } static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 88 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 89
return  neumann_homogeneous(); } return 0.; }
#line 92 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static double _boundary4 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 91 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann(val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))); } return 0.; } static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 91 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); } return 0.; }
#line 93 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static double _boundary5 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 92 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann(-val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0)); } return 0.; } static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 92 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 93
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 93
return  neumann_homogeneous(); } return 0.; }
#line 37 "./physics.h"
static double _boundary6 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 36 "./physics.h"
return  dirichlet(0.); return 0.; } static double _boundary6_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 36 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 38 "./physics.h"
static double _boundary7 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 37 "./physics.h"
return  dirichlet((-1.)); return 0.; } static double _boundary7_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 37 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 39 "./physics.h"
static double _boundary8 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 38 "./physics.h"
return  dirichlet(0); return 0.; } static double _boundary8_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 38 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 40 "./physics.h"
static double _boundary9 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 39 "./physics.h"
return  dirichlet((-1.)); return 0.; } static double _boundary9_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 39 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 44 "./physics.h"
static double _boundary10 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 43 "./physics.h"
return  dirichlet(9.81*(.5)*y/273.); return 0.; } static double _boundary10_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 43 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 45 "./physics.h"
static double _boundary11 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 44 "./physics.h"
return  dirichlet(9.81*(.5)*y/273.); return 0.; } static double _boundary11_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 44 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 48 "./physics.h"
static double _boundary12 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 47 "./physics.h"
return  dirichlet((-1.)); return 0.; } static double _boundary12_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 47 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 49 "./physics.h"
static double _boundary13 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 48 "./physics.h"
return  dirichlet((-1.)); return 0.; } static double _boundary13_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 48 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 50 "./physics.h"
static double _boundary14 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 49 "./physics.h"
return  dirichlet(0.); return 0.; } static double _boundary14_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 49 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
#line 51 "./physics.h"
static double _boundary15 (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 50 "./physics.h"
return  dirichlet(0.); return 0.; } static double _boundary15_homogeneous (Point point, Point neighbor, scalar _s) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 50 "./physics.h"
return  dirichlet_homogeneous(); return 0.; }
