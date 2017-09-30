# A python Toolkit for HXMT data analysis


## the usage of HXMT software


**NOTE: Only the software for High Energy detector is introduced below.**


*Before utilizing this toolkit, a HXMT software package should be installed([download page](http://www.hxmt.org/index.php/dataan/2013-03-22-08-13-10)).*

A full description of the usage of HXMT software can be found in fhelp document while it is well installed. The complete prosedures including three steps: PI calculationg, selection of good time intervals, producing data production(light curve, background and spectrum files).

### hepical

The task will simply do the PI calculation from the  channel column of input event file.

#### examples
Process an HE event file with PI column, "infile," using the defaults to calculate the PI values. 
    
`> hepical evtfile=he_evt.fits outfile=he_pi.fits`
