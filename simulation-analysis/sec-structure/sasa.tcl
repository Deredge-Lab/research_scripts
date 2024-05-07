set outfile [open system-name_sasa.dat w]
set nf [molinfo top get numframes]
set protein [atomselect top "protein"]
set sel [atomselect top "resid 1 to 100"]
for {set i 0} {$i < $nf} {incr i} {
	molinfo top set frame $i
	### default probe radius = 1.4 angstroms
	set sasa [measure sasa 1.4 $protein -restrict $sel]  
	puts $outfile "Frame $i, SASA $sasa" 
    }
close $outfile
