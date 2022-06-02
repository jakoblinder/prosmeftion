#procedure getOps(ops, Max)
*name: string
*	name of the table file which will be generated
*ops: set of non commutable functions
* 	specifies which fields are printed in the table
	.sort
	Function f;
*	#define Max "100"
	#define matched "0"
	Table `ops'Tab(1:`Max');

	multiply left f;
	.sort
	#do i= 1, `Max'
		if (match(f*H?`ops'(?b)));	
			id once f*H?`ops'$fac(?b$arg) = H(?b)*f;
			Print "Matched `i': %$(%$)", $fac,$arg;
			redefine matched "1";
		else;
*			Print "STOP";
			redefine matched "0";
		endif;
		.sort
		#if `matched' == 1
			Fill `ops'Tab(`i') = `$fac'(`$arg');
			goto 1;
		#endif
		if (match(f*H?!`ops'(?b)));
			id once f*H?!`ops'$fac(?b$arg) = H(?b)*f;
			Print "NoMatched `i': %$(%$)", $fac,$arg;
		endif;
		label 1;
	#enddo
	id f = 1;
	.sort
	printtable `ops'Tab > `ops'.t;
#endprocedure
