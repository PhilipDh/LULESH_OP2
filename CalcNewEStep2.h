inline void CalcNewEStep2(double *compHalfStep, double *delvc, double *q_new, double *pbvc, double *e_new, double *bvc,
                            double *pHalfStep, double *ql_old, double *qq_old, double *p_old, double *q_old, double *rho0
){
    Real_t vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep[0]) ;

    if ( delvc[0] > Real_t(0.) ) {
        q_new[0] /* = qq_old[0] = ql_old[0] */ = Real_t(0.) ;
    }
    else {
        Real_t ssc = ( pbvc[0] * e_new[0]
                + vhalf * vhalf * bvc[0] * pHalfStep[0] ) / (*rho0) ;

        if ( ssc <= Real_t(.1111111e-36) ) {
        ssc = Real_t(.3333333e-18) ;
        } else {
        ssc = SQRT(ssc) ;
        }

        q_new[0] = (ssc*ql_old[0] + qq_old[0]) ;
    }

    e_new[0] = e_new[0] + Real_t(0.5) * delvc[0]
        * (  Real_t(3.0)*(p_old[0]     + q_old[0])
            - Real_t(4.0)*(pHalfStep[0] + q_new[0])) ;
}