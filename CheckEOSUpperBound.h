inline void CheckEOSUpperBound(const double *vnewc, double *compHalfStep, double *compression, double* p_old, const double *eosvmax){
    if (vnewc[0] >= (*eosvmax)) { /* impossible due to calling func? */
        p_old[0]        = double(0.) ;
        compression[0]  = double(0.) ;
        compHalfStep[0] = double(0.) ;
    }
}