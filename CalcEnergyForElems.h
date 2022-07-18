inline void CalcNewE(double *e_new, double *e_old, double *emin, double *delvc, double *p_old, double *q_old, double *work){
    e_new[0] = e_old[0] - Real_t(0.5) * delvc[0] * (p_old[0] + q_old[0])
        + Real_t(0.5) * work[0];

    if (e_new[0]  < (*emin) ) {
        e_new[0] = (*emin) ;
    }
}

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

inline void CalcNewEStep3(double *e_new, double *work, double *e_cut, double *emin){
    e_new[0] += Real_t(0.5) * work[0];

    if (FABS(e_new[0]) < (*e_cut)) {
        e_new[0] = Real_t(0.)  ;
    }
    if (     e_new[0]  < (*emin) ) {
        e_new[0] = (*emin) ;
    }
}

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