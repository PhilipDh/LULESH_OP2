inline void CalcCourantConstraint(double *ss, double *vdov, double *arealg, double *qqc2, double *dtcourant){
    Real_t dtf = ss[0] * ss[0];

    if ( vdov[0] < Real_t(0.) ) {

        dtf = dtf
            + (*qqc2) * arealg[0] * arealg[0]
            * vdov[0] * vdov[0] ;
    }

    dtf = SQRT(dtf) ;

    dtf = arealg[0] / dtf ;

/* determine minimum timestep with its corresponding elem */
    if (vdov[0] != Real_t(0.)) {
        if ( dtf < (*dtcourant) ) {
        (*dtcourant) = dtf ;
        // courant_elem = indx ;
        }
    }
}