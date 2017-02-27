#!/bin/tcsh -f 

set file='FixingRun.txt'
@ i = 0

foreach line ("`cat $file`")
   set argv = ( $line )
   set valuea = $1
   set valueb = $2
   @ RunNumber = $valuea * 100 + $valueb
   #echo $valuea $valueb $RunNumber
   #echo $i
   @ i += 1
   ./Fitter/min2d_oct.csh ${RunNumber} > log${i} &
end
