inline void CalcNewEStep4(double *delvc, double *pbvc, double *e_new, double *vnewc, double *bvc, double *p_new,
                            double *ql_old, double *qq_old, double *p_old, double *q_old, double *q_new,
                            double *pHalfStep, double *rho0, double *e_cut, double *emin){
    const Real_t sixth = Real_t(1.0) / Real_t(6.0) ;
    // Index_t ielem = regElemList[0];
    Real_t q_tilde ;

    if (delvc[0] > Real_t(0.)) {
        q_tilde = Real_t(0.) ;
    }
    else {
        Real_t ssc = ( pbvc[0] * e_new[0]
                + vnewc[0] * vnewc[0] * bvc[0] * p_new[0] ) / (*rho0) ;

        if ( ssc <= Real_t(.1111111e-36) ) {
        ssc = Real_t(.3333333e-18) ;
        } else {
        ssc = SQRT(ssc) ;
        }

        q_tilde = (ssc*ql_old[0] + qq_old[0]) ;
    }

    e_new[0] = e_new[0] - (  Real_t(7.0)*(p_old[0]     + q_old[0])
                            - Real_t(8.0)*(pHalfStep[0] + q_new[0])
                            + (p_new[0] + q_tilde)) * delvc[0]*sixth ;

    if (FABS(e_new[0]) < (*e_cut)) {
        e_new[0] = Real_t(0.)  ;
    }
    if (     e_new[0]  < (*emin) ) {
        e_new[0] = (*emin) ;
    }   
}