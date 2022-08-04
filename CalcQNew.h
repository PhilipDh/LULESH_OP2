inline void CalcQNew(double *delvc, double *pbvc, double *e_new, double *vnewc, double *bvc, double *p_new,
                    double *q_new, double *ql_old, double *qq_old, double *rho0, double *q_cut
){
    // Index_t ielem = regElemList[0];

    if ( delvc[0] <= Real_t(0.) ) {
        Real_t ssc = ( pbvc[0] * e_new[0]
                + vnewc[0] * vnewc[0] * bvc[0] * p_new[0] ) / (*rho0) ;

        if ( ssc <= Real_t(.1111111e-36) ) {
        ssc = Real_t(.3333333e-18) ;
        } else {
        ssc = SQRT(ssc) ;
        }

        q_new[0] = (ssc*ql_old[0] + qq_old[0]) ;

        if (FABS(q_new[0]) < (*q_cut)) q_new[0] = Real_t(0.) ;
    }
}