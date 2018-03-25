
	capture erase lstaft.mlib
	mata: mata set matastrict off
	qui {
		do "./staft/staft_setup.mata"
		do "./staft/staft_lf2.mata"
		do "./staft/rcsgen_mata.mata"

		mata: mata mlib create lstaft, dir(.)
		mata: mata mlib add    lstaft *(), dir(.)
		mata: mata d *()  
		mata mata clear
		mata mata mlib index
	}
