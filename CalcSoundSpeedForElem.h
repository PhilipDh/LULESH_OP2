inline void CalcSoundSpeedForElem(const double *pbvc, const double *enewc, const double *vnewc, const double *bvc, const double *pnewc, double *ss, const double *rho0){
    double ssTmp = (pbvc[0] * enewc[0] + vnewc[0] * vnewc[0] *
                bvc[0] * pnewc[0]) / (*rho0);
    if (ssTmp <= double(.1111111e-36)) {
        ssTmp = double(.3333333e-18);
    }
    else {
        ssTmp = sqrt(ssTmp);
    }
    // domain.ss(ielem) = ssTmp ;
    ss[0] = ssTmp ;
}