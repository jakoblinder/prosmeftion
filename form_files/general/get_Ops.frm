Function H, [H+], B, W, G;
Function D, e, u, b, l, q, [ebar], [ubar], [bbar], [lbar], [qbar];

Function BL, BR, WL, WR, GL, GR;
Function yu, yd, ye, [yu+], [yd+], [ye+], [su2eps], [su3eps], T, gamma;

Function t;

AutoDeclare Indices lor      = 4; * 4d Lorentz index
AutoDeclare Indices spin     = 4; * Index for Gamma matrices/ spinor index
AutoDeclare Indices gauge    = 2; * SU(2)-index in fundamental
AutoDeclare Indices gaugeadj = 3; * SU(2)-index in adjoint
AutoDeclare Indices colf     = 3; * SU(3)-index in fundamental
AutoDeclare Indices cola     = 8; * SU(3)-index in adjoint
AutoDeclare Indices flav     = n; * flavor index

Set fields: H, [H+], B, W, G, D, e, u, b, l, q, [ebar], [ubar], [bbar], [lbar], [qbar];
Set tens: yu, yd, ye, [yu+], [yd+], [ye+], [su2eps], [su3eps], T, gamma;

Local expr = gamma(lor1,spin7380,spin7381)*H(gauge7150)*e(spin7381,flav7379)*D(lor1,[ebar](spin7380,flav7900))*ye(flav7658,flav7379)*[H+](gauge7150)*[ye+](flav7900,flav7658);

* Maximum dimension of Lagrangian operators
#define dim "6"
* Minimum dimension of a field
#define minFielddim "1"
* Maximum numbers of terms in one operators
#define nterms "6"

#call getOps(fields, `nterms')
#call getOps(tens, {3*`nterms'})

*Bracket H, [H+], B, W, G, D, e, u, b, l, q, [ebar], [ubar], [bbar], [lbar], [qbar], BL, BR, WL, WR, GL, GR, [d_C], [e_C], L, Q, [u_C], [d_C+], [e_C+], [L+], [Q+], [u_C+];

Format 255;
Print +ss;
.end
