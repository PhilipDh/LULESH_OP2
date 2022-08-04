inline void CalcPNew(double *p_new, double *bvc, double *e_old, double *vnewc,double *p_cut, double *eosvmax, double *pmin){
    p_new[0] = bvc[0] * e_old[0] ;

    if    (FABS(p_new[0]) <  (*p_cut)   )
        p_new[0] = Real_t(0.0) ;

    if    ( vnewc[0] >= (*eosvmax) ) /* impossible condition here? */
        p_new[0] = Real_t(0.0) ;

    if    (p_new[0]       <  (*pmin))
        p_new[0]   = (*pmin) ;
}