inline void CalcHalfStepBVC(double *bvc, double *compression, double *pbvc){
    Real_t c1s = Real_t(2.0)/Real_t(3.0) ;
    bvc[0] = c1s * (compression[0] + Real_t(1.));
    pbvc[0] = c1s;
}

inline void CalcPHalfstep(double *p_new, double *bvc, double *e_old, double *vnewc,double *p_cut, double *eosvmax, double *pmin){
    p_new[0] = bvc[0] * e_old[0] ;

    if    (FABS(p_new[0]) <  (*p_cut)   )
        p_new[0] = Real_t(0.0) ;

    if    ( vnewc[0] >= (*eosvmax) ) /* impossible condition here? */
        p_new[0] = Real_t(0.0) ;

    if    (p_new[0]       <  (*pmin))
        p_new[0]   = (*pmin) ;
}

inline void CalcBVC(double *bvc, double *compression, double *pbvc){
    Real_t c1s = Real_t(2.0)/Real_t(3.0) ;
    bvc[0] = c1s * (compression[0] + Real_t(1.));
    pbvc[0] = c1s;
}

inline void CalcPNew(double *p_new, double *bvc, double *e_old, double *vnewc,double *p_cut, double *eosvmax, double *pmin){
    p_new[0] = bvc[0] * e_old[0] ;

    if    (FABS(p_new[0]) <  (*p_cut)   )
        p_new[0] = Real_t(0.0) ;

    if    ( vnewc[0] >= (*eosvmax) ) /* impossible condition here? */
        p_new[0] = Real_t(0.0) ;

    if    (p_new[0]       <  (*pmin))
        p_new[0]   = (*pmin) ;
}