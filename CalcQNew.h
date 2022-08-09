inline void CalcQNew(const double *delvc, const double *pbvc, const double *e_new, const double *vnewc, const double *bvc, const double *p_new,
                    double *q_new, const double *ql_old, const double *qq_old, const double *rho0, const double *q_cut
){
    // Index_t ielem = regElemList[0];

    if ( delvc[0] <= double(0.) ) {
        double ssc = ( pbvc[0] * e_new[0]
                + vnewc[0] * vnewc[0] * bvc[0] * p_new[0] ) / (*rho0) ;

        if ( ssc <= double(.1111111e-36) ) {
        ssc = double(.3333333e-18) ;
        } else {
        ssc = sqrt(ssc) ;
        }

        q_new[0] = (ssc*ql_old[0] + qq_old[0]) ;

        if (fabs(q_new[0]) < (*q_cut)) q_new[0] = double(0.) ;
    }
}