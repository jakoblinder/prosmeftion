*--#[ tensors :
CFunction yu, [yu+], yd, [yd+], ye, [ye+], [su2eps], [su3eps], [sl2Ceps], [su2dK], [su3dK], [sl2CdK], [flavdK], T, gamma, TSU2, TSU3;
* Indices and functions for derivatives in SL2C notation.
CFunction sigma, sigmabar;
CFunction sigma2, sigmabar2;
* Auxiliary antisymmetric epsilons, used in combination with replace_.
CFunction [su2epsA](antisymmetric), [su3epsA](antisymmetric), [sl2CepsA](antisymmetric);

* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.
* Since two indices are also symmetric when they are cyclic and vice versa and pattern matching is not allowed for symmetric function but for cyclic it is, [sl2CdK] is declared as cyclic.
CFunction [su2dK](cyclic), [su3dK](cyclic), [sl2CdK](cyclic), [flavdK](cyclic);

*--#] tensors :

*--#[ coefficient :
Symbols d, eps, lambdah, At, g1, g2, g3, mu, lambdaphi, kappa, Ms, Mu, muM, [2L[Ms,muM]], n;
Symbols [At/Ms], [mu/Ms], [Mu/Ms], [muM/Ms];
*--#] coefficient :

*--#[ operators :
Off Statistics;
Function H, [H+], G, W, B, GL, GR, WL, WR, BL, BR;
CFunction Hc, [H+c], Gc, Wc, Bc, GLc, GRc, WLc, WRc, BLc, BRc;

Function D, l, lbar, L, [L+], e, ebar, [e_C], [e_C+], u, ubar, [u_C], [u_C+], b, bbar, [d_C], [d_C+], q, qbar, Q, [Q+];
CFunction Dc, lc, lbarc, Lc, [L+c], ec, ebarc, [e_Cc], [e_C+c], uc, ubarc, [u_Cc], [u_C+c], bc, bbarc, [d_Cc], [d_C+c], qc, qbarc, Qc, [Q+c];

Set Fieldc: WLc, WRc, BLc, BRc;
Set Field: WL, WR, BL, BR;

CFunction xi, [xi+], chi, [chi+];

Set spinors: l, e, [e_C], u, [u_C], b, [d_C], q;
Set spinorsAdj: lbar, ebar, [e_C+], ubar, [u_C+], bbar, [d_C+], qbar;
Set spinorsAll: l, e, [e_C], u, [u_C], b, [d_C], q, lbar, ebar, [e_C+], ubar, [u_C+], bbar, [d_C+], qbar;

Set spinorsc: lc, ec, [e_Cc], uc, [u_Cc], bc, [d_Cc], qc;
Set spinorsAdjc: lbarc, ebarc, [e_C+c], ubarc, [u_C+c], bbarc, [d_C+c], qbarc;
Set spinorsAllc: lc, ec, [e_Cc], uc, [u_Cc], bc, [d_Cc], qc, lbarc, ebarc, [e_C+c], ubarc, [u_C+c], bbarc, [d_C+c], qbarc;

* D2 = D_mu * D^mu:
Function D2;
* Intern abbreviation for equation of motion:
Function EOM;
* Total field strength tensor, necessary for EOM substitutions:
Function FL, FR;

*--#] operators :

*--#[ indices :
AutoDeclare Indices lor      = 4; * 4d Lorentz index.
AutoDeclare Indices Lsl      = 2; * Subscript and undotted SL2C-Index.
AutoDeclare Indices Usl      = 2; * Superscripted and undotted SL2C-Index.
AutoDeclare Indices Lsldot   = 2; * Subscript and dotted SL2C-Index.
AutoDeclare Indices Usldot   = 2; * Superscripted and dotted SL2C-Index.
AutoDeclare Indices spin     = 2; * 4 dimensional Spin index, which is summed over 2 Weyl spinors.
AutoDeclare Indices gauge    = 2; * Fundamental index of SU2.
AutoDeclare Indices gaugeadj = 3; * Adjoint index of SU2.
AutoDeclare Indices colf     = 3; * Fundamental index of SU3.
AutoDeclare Indices cola     = 8; * Adjoint index of SU3.
AutoDeclare Indices flav     = n; * Flavor index.
AutoDeclare Indices sbasis   = N; * Index for basis tensors of SU(N).

AutoDeclare Indices op; * auxiliary index for converting between commuting and noncommuting operators.

* Declare some Symbols for pattern matching
Symbols k, m;
Autodeclare Symbols i;
*--#] indices :
