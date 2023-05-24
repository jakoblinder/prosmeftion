#message `PIPES_'

#setexternal `PIPE1_'
#prompt READY

#fromexternal

Local expr1=
#fromexternal
;
id a?something = 4;
.sort
#toexternal "%E\n", expr1

.end
