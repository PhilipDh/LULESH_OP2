inline void CalcSoundSpeedForElem(double *pbvc, double *enewc, double *vnewc, double *bvc, double *pnewc, double *ss, double *rho0){
    Real_t ssTmp = (pbvc[0] * enewc[0] + vnewc[0] * vnewc[0] *
                bvc[0] * pnewc[0]) / (*rho0);
    if (ssTmp <= Real_t(.1111111e-36)) {
        ssTmp = Real_t(.3333333e-18);
    }
    else {
        ssTmp = SQRT(ssTmp);
    }
    // domain.ss(ielem) = ssTmp ;
    ss[0] = ssTmp ;
}