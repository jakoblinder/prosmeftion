#include declarations.h # tensors
#include declarations.h # declarations

Local expression = 
#include formTerm.h
;
.sort

CFunction coeff;
Bracket yu, [yu+], yd, [yd+], ye, [ye+], su2eps, su3eps, sl2Ceps, su2dK, su3dK, sl2CdK, T, gamma, H, [H+], G, W, B, GL, [GL+], WL, [WL+], BL, [BL+], D, l, lbar, L, [L+], e, ebar, eC, [eC+], u, ubar, uC, [uC+], b, bbar, dC, [dC+], q, qbar, Q, [Q+];
.sort
collect coeff;
.sort

CFunction term;
putinside term;
.sort
Format nospaces;

Print +s;
.end