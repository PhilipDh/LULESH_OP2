inline void CheckEOSLowerBound(const double *vnewc, double *compHalfStep, const double *compression, const double *eosvmin){
    if (vnewc[0] <= (*eosvmin)) { /* impossible due to calling func? */
        compHalfStep[0] = compression[0] ;
    }
}