#procedure fractorizeCoeff
.sort
	Symbols u, v, x, y;
* Get sign
	PolyFun sign;
	.sort
* Write positive numerical coefficient in coeff and sign in sign
	PolyFun;
	id sign(x?) = coeff(sig_(x)*x)*sign(sig_(x));
	.sort
* write everything else in frac function
	repeat;		id d?!{Ms}^n?pos_ = frac(d^n, 1);
		id d?!{Ms}^n?neg_ = frac(1, d^-n);
	endrepeat;
* combine the separate fractions
	repeat;
		id frac(u?, v?)*frac(x?, y?) = frac(u * x, v * y);
	endrepeat;
	repeat;
		id frac(u?, 1) = u;
		id sign(x?) = x;
		id coeff(1) = 1;
	endrepeat;
	.sort
	Bracket Ms;
	.sort
	Collect bbracket;
#endprocedure