inline void CopyVelocityToTempArray(double *vnewc, double *vnew){
    vnewc[0] = vnew[0] ;
}

inline void ApplyLowerBoundToVelocity(double *vnewc, double *eosvmin){
    if (vnewc[0] < (*eosvmin))
        vnewc[0] = (*eosvmin) ;
}

inline void ApplyUpperBoundToVelocity(double *vnewc, double *eosvmax){
    if (vnewc[0] > (*eosvmax))
        vnewc[0] = (*eosvmax) ;
}

inline void ALE3DRelevantCheck(double *v, double *eosvmin, double *eosvmax){
    Real_t vc = v[0] ;
    if ((*eosvmin) != Real_t(0.)) {
        if (vc < (*eosvmin))
        vc = (*eosvmin) ;
    }
    if ((*eosvmax) != Real_t(0.)) {
        if (vc > (*eosvmax))
        vc = (*eosvmax) ;
    }
    if (vc <= 0.) {
        exit(VolumeError);
    }
}

inline void CopyEOSValsIntoArray(
                                double *e_old, double *e,
                                double *delvc, double *delv,
                                double *p_old, double *p,
                                double *q_old, double *q,
                                double *qq_old, double *qq,
                                double *ql_old, double *ql){
    e_old[0] = e[0] ;
    delvc[0] = delv[0] ;
    
    p_old[0] = p[0] ;
    q_old[0] = q[0] ;

    qq_old[0] = qq[0] ;
    ql_old[0] = ql[0] ;
}

inline void CalcHalfSteps(double *compression, double *vnewc, double* delvc, double *compHalfStep){
    Real_t vchalf ;
    compression[0] = Real_t(1.) / vnewc[0] - Real_t(1.);
    vchalf = vnewc[0] - delvc[0] * Real_t(.5);
    compHalfStep[0] = Real_t(1.) / vchalf - Real_t(1.);
}

inline void CheckEOSLowerBound(double *vnewc, double *compHalfStep, double *compression, double *eosvmin){
    if (vnewc[0] <= (*eosvmin)) { /* impossible due to calling func? */
        compHalfStep[0] = compression[0] ;
    }
}

inline void CheckEOSUpperBound(double *vnewc, double *compHalfStep, double *compression, double* p_old, double *eosvmax){
    if (vnewc[0] >= (*eosvmax)) { /* impossible due to calling func? */
        p_old[0]        = Real_t(0.) ;
        compression[0]  = Real_t(0.) ;
        compHalfStep[0] = Real_t(0.) ;
    }
}

inline void CalcEOSWork(double *work){
    work[0] = Real_t(0.);
}

inline void CopyTempEOSVarsBack(
                                double *p, double *p_new,
                                double *e, double *e_new,
                                double *q, double *q_new){
    p[0] = p_new[0];
    e[0] = e_new[0];
    q[0] = q_new[0];  
}