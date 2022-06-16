#procedure positiveTerm
.sort
* Ensure that first term in the bracket is positive
	id bbracket(x?$inb) = bbracket(x);
* Save original expression in order to work only on the first term
	$coeffic = coefficient;
	.sort
	CFunction bsign;
	Drop coefficient;
* Get only the first term and its sign
	Local FirstTerm = firstterm_($inb);
	.sort
	PolyFun bsign;
	.sort
	PolyFun;
	id bsign(x?$s) = bsign(x);
	$sign = sig_($s);
*	Print "Sign: %$", $sign;
	.sort
	Drop FirstTerm;
* restore original expression
	Local coefficient = $coeffic;
* ensure that the first term in the bracket is positive
	id bbracket(x?) = $sign*bbracket($sign*x);
#endprocedure