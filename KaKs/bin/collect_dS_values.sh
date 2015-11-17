#!/bin/bash

tail -n +2 Fam.*.codmlrun/2NG.dS | grep -v "<==" | cut -f2- -d" " | tr ' ' '\n' | grep  -v "\-1.0000" | grep  . > all.dS.values.txt 
# tail -n +2 Fam.*.codmlrun/Fam.*.2NG.dS | grep -v "<==" | cut -f2- -d" " | tr ' ' '\n' | grep  -v "\-1.0000" | grep  . > all.dS.values.txt 