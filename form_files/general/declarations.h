*--#[ tensors :
CFunction yu, [yu+], yd, [yd+], ye, [ye+], su2eps, su3eps, sl2Ceps, su2dK, su3dK, sl2CdK, T, gamma;
*--#] tensors :

*--#[ declarations :
CFunction d, eps, lambdah, At, g1, g2, g3, mu, lambdaphi, kappa, Ms, Mu, muM, [2L[Ms,muM]];
Function H, [H+], G, W, B, GL, [GL+], WL, [WL+], BL, [BL+];
CFunction Hc, [H+c], Gc, Wc, Bc, GLc, [GL+c], WLc, [WL+c], BLc, [BL+c];

Function D, l, lbar, L, [L+], e, ebar, eC, [eC+], u, ubar, uC, [uC+], b, bbar, dC, [dC+], q, qbar, Q, [Q+];
CFunction Dc, lc, lbarc, Lc, [L+c], ec, ebarc, eCc, [eC+c], uc, ubarc, uCc, [uC+c], bc, bbarc, dCc, [dC+c], qc, qbarc, Qc, [Q+c];

Set Fieldc: WLc, [WL+c], BLc, [BL+c];

CFunction xi, [xi+], chi, [chi+];

Set spinors: l, e, u, b, q;
Set spinorsAdj: lbar, ebar, ubar, bbar, qbar;
Set spinorsAll: l, e, u, b, q, lbar, ebar, ubar, bbar, qbar;

Set spinorsc: lc, ec, uc, bc, qc;
Set spinorsAdjc: lbarc, ebarc, ubarc, bbarc, qbarc;
Set spinorsAllc: lc, ec, uc, bc, qc, lbarc, ebarc, ubarc, bbarc, qbarc;

AutoDeclare Indices lor      = 4; * 4d Lorentz index
AutoDeclare Indices lorA     = 4; * Auxiliary 4d Lorentz index
AutoDeclare Indices spin     = 4; * Index for Gamma matrices/ spinor index
AutoDeclare Indices spinA    = 2; * Auxiliary index for Gamma matrices/ spinor index
AutoDeclare Indices gauge    = 2; * SU(2)-index in fundamental
AutoDeclare Indices gaugeA   = 2; * Auxiliary SU(2)-index in fundamental
AutoDeclare Indices gaugeadj = 3; * SU(2)-index in adjoint
AutoDeclare Indices colf     = 3; * SU(3)-index in fundamental
AutoDeclare Indices colfA    = 3; * Auxiliary SU(3)-index in fundamental
AutoDeclare Indices cola     = 8; * SU(3)-index in adjoint
AutoDeclare Indices flav     = n; * flavor index

AutoDeclare Indices op; * auxiliary index for converting between commuting and noncommuting operators.

* Declare some Symbols for pattern matching
Symbols k,m;
*
* Indices and functions for derivatives in SL2C notation.
*
CFunction sigma, sigmabar;
CFunction sigma2, sigmabar2;
* Auxiliary antisymmtric epsilons, used in combination with replace_.
CFunction su2epsA(antisymmetric), su3epsA(antisymmetric), sl2CepsA(antisymmetric);

* Declare Kronecker Delta symbol for Sl2C Indices, because built in can not handle upper and lower (un-)dottet indices.
* Since two indices are also symmetric when they are cyclic and vice versa and pattern matching is not allowed for symmetric function but for cyclic it is, [sl2CdK] is declared as cyclic.
CFunction su2dK(cyclic), su3dK(cyclic), sl2CdK(cyclic);

Off Statistics;
*--#] declarations :
