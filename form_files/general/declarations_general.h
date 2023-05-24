*--#[ tensors :
CFunction yu, [yu+], yd, [yd+], ye, [ye+], [su2eps], [su3eps], [sl2Ceps], [lorentzeps], T, sigma, sigmabar, sigma2, sigmabar2, gamma, sigma2lor, TSU2, TSU3;
* Auxiliary antisymmetric epsilons, used in combination with replace_.
CFunction [su2epsA](antisymmetric), [su3epsA](antisymmetric), [sl2CepsA](antisymmetric), [lorentzepsA](antisymmetric);

* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.
* Since two indices are also symmetric when they are cyclic and vice versa, and pattern matching is not allowed for symmetric function but for cyclic, [sl2CdK] is declared as cyclic.
CFunction [su2dK](cyclic), [su3dK](cyclic), [sl2CdK](cyclic), [flavdK](cyclic);
*--#] tensors :

*--#[ NCtensors :
Function yu, [yu+], yd, [yd+], ye, [ye+], [su2eps], [su3eps], [sl2Ceps], [lorentzeps], [su2dK], [su3dK], [sl2CdK], [flavdK], T, sigma, sigmabar, sigma2, sigmabar2, gamma, sigma2lor, TSU2, TSU3;
* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.
* Since two indices are also symmetric when they are cyclic and vice versa, and pattern matching is not allowed for symmetric function but for cyclic, [sl2CdK] is declared as cyclic.
Function [su2dK](cyclic), [su3dK](cyclic), [sl2CdK](cyclic), [flavdK](cyclic);
*--#] NCtensors :

*--#[ coefficient :
Symbols d, eps, lambdah, At, g1, g2, g3, mu, lambdaphi, kappa, Ms, Mu, muM, [2L[Ms,muM]], n, N;
Symbols [At/Ms], [mu/Ms], [Mu/Ms], [muM/Ms];
*--#] coefficient :

*--#[ operators :
* Define all fields and auxiliary fields:
Function lbar, e, bbar, bC, u, [H+], qbarC, l, lbarC, qC, eC, ubar, H, qbar, B, ebarC, uC, bbarC, ebar, b, W, G, q, ubarC, lC;
CFunction Wc, [H+c], qbarc, eCc, uCc, qCc, ubarc, qc, lc, Bc, lCc, ebarCc, Gc, ec, uc, bbarc, bCc, qbarCc, lbarCc, lbarc, ubarCc, bbarCc, bc, Hc, ebarc;
CFunction Lc, [d_Cc], [u_C+c], [d_C+c], WLc, [L+c], [e_Cc], Qc, [H+c], [Q+c], BLc, GLc, [e_C+c], [u_Cc], Hc, GRc, BRc, WRc;
Function GL, [u_C+], [e_C+], Q, H, [Q+], [d_C], GR, WL, BL, BR, WR, L, [H+], [L+], [e_C], [d_C+], [u_C];

* Derivative:
Function D;
Set fieldstrengthsc: GLc, GRc, WLc, WRc, BLc, BRc;
Set fieldstrengths: GL, GR, WL, WR, BL, BR;

CFunction xi1, [xi1+], chi1, [chi1+], xi2, [xi2+], chi2, [chi2+];

Set spinors: l, e, u, b, q;
Set spinorsAdj: lbar, ebar, ubar, bbar, qbar;
Set spinorsAll: l, e, u, b, q, lbar, ebar, ubar, bbar, qbar;

Set spinorsc: lc, ec, uc, bc, qc;
Set spinorsAdjc: lbarc, ebarc, ubarc, bbarc, qbarc;
Set spinorsAllc: lc, ec, uc, bc, qc, lbarc, ebarc, ubarc, bbarc, qbarc;

* List of all possibly occurring fields bevor they are converted into Lorentz irreps
* Note: If a field (like e.g. a fieldstrength tensor B) is converted into two different Lorentz irreps (BL, BR) is has to occur twice in this list
Set AllFields: H, [H+], G, G, W, W, B, B, l, lbar, lC, lbarC, e, ebar, eC, ebarC, u, ubar, uC, ubarC, b, bbar, bC, bbarC, q, qbar, qC, qbarC;

* List of the same fields as in AllFields (in the same order!) but defined as a commuting Function
Set AllcFields: Hc, [H+c], Gc, Gc, Wc, Wc, Bc, Bc, lc, lbarc, lCc, lbarCc, ec, ebarc, eCc, ebarCc, uc, ubarc, uCc, ubarCc, bc, bbarc, bCc, bbarCc, qc, qbarc, qCc, qbarCc;

* List of the same fields (first as commutative fields) as in AllFields (in the same order!) but with possible replacements like Dirac to Weyl spinors and so on
Set AllcconFields: Hc, [H+c], GLc, GRc, WLc, WRc, BLc, BRc, Lc, [L+c], [L+c], Lc, [e_C+c], [e_Cc], [e_Cc], [e_C+c], [u_C+c], [u_Cc], [u_Cc], [u_C+c], [d_C+c], [d_Cc], [d_Cc], [d_C+c], Qc, [Q+c], [Q+c], Qc;

* Now as noncommutative fields
Set AllconFields: H, [H+], GL, GR, WL, WR, BL, BR, L, [L+], [L+], L, [e_C+], [e_C], [e_C], [e_C+], [u_C+], [u_C], [u_C], [u_C+], [d_C+], [d_C], [d_C], [d_C+], Q, [Q+], [Q+], Q;

* Set for convenient insertion of op1, op2, ... indices
Set op: op1,...,op100;

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
Autodeclare Symbols i, j;
*--#] indices :
