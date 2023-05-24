#procedure getOps(ops)
* ops: set of non commutable functions
*   specifies which fields are printed in the table
	.sort
	Function f;
* Max is a preprocessor variable, which tells the maximum possible number of fields and tensors in one operator.
* This number is allowed to be much larger then the number of ops which really occurs.
	#define Max "1000"
	#define matched "0"
	#define matchedsomething "0"
	Table `ops'Tab(1:`Max');

	multiply left f;
	.sort
	#write "Find ops from set `ops'."
	#do i= 1, `Max'
		redefine matchedsomething "0";
		if (match(f*H?`ops'(?b)));
			id once f*H?`ops'$fac(?b$arg) = H(?b)*f;
			redefine matched "1";
		else;
*			Print "STOP";
			redefine matched "0";
		endif;
		.sort
		#if `matched' == 1
			Fill `ops'Tab(`i') = `$fac'(`$arg');
			#write "Matched `i': %$(%$)", $fac,$arg
			#toexternal "`ops'-`i': %$(%$)
", $fac,$arg
			redefine matchedsomething "1";
			goto 1;
		#endif
		if (match(f*H?!`ops'(?b)));
			id once f*H?!`ops'$fac(?b$arg) = H(?b)*f;
			redefine matchedsomething "1";
		endif;
		label 1;
		.sort
		#if ((`matched' == 0) && (`matchedsomething' == 1))
			#write "NoMatched `i': %$(%$)", $fac,$arg
		#endif

		#if `matchedsomething' == 0
			#breakdo
		#endif	
	#enddo
	id f = 1;
#endprocedure
