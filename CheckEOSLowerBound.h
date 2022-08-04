inline void CheckEOSLowerBound(double *vnewc, double *compHalfStep, double *compression, double *eosvmin){
    if (vnewc[0] <= (*eosvmin)) { /* impossible due to calling func? */
        compHalfStep[0] = compression[0] ;
    }
}