
capture erase lstaft.mlib
mata: mata set matastrict off
qui {
        do "./mata/staft_setup.mata"
        do "./mata/staft_lf2.mata"
        do "./mata/rcsgen_mata.mata"

        mata: mata mlib create lstaft, dir(.)
        mata: mata mlib add    lstaft *(), dir(.)
        mata: mata d *()  
        mata mata clear
        mata mata mlib index
}
