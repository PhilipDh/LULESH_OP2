inline void CalcHydroConstraint(double *vdov, double *dthydro){
    Real_t dtdvov = m_dvovmax / (FABS(vdov[0])+Real_t(1.e-20)) ;

    if ( (*dthydro) > dtdvov ) {
            (*dthydro) = dtdvov ;
            // hydro_elem = indx ;
    }
}