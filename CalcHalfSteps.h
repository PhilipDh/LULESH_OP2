inline void CalcHalfSteps(double *compression, double *vnewc, double* delvc, double *compHalfStep){
    Real_t vchalf ;
    compression[0] = Real_t(1.) / vnewc[0] - Real_t(1.);
    vchalf = vnewc[0] - delvc[0] * Real_t(.5);
    compHalfStep[0] = Real_t(1.) / vchalf - Real_t(1.);
}