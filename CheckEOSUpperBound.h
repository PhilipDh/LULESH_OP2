inline void CheckEOSUpperBound(double *vnewc, double *compHalfStep, double *compression, double* p_old, double *eosvmax){
    if (vnewc[0] >= (*eosvmax)) { /* impossible due to calling func? */
        p_old[0]        = Real_t(0.) ;
        compression[0]  = Real_t(0.) ;
        compHalfStep[0] = Real_t(0.) ;
    }
}