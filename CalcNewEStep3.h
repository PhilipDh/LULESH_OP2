inline void CalcNewEStep3(double *e_new, double *work, double *e_cut, double *emin){
    e_new[0] += Real_t(0.5) * work[0];

    if (FABS(e_new[0]) < (*e_cut)) {
        e_new[0] = Real_t(0.)  ;
    }
    if (     e_new[0]  < (*emin) ) {
        e_new[0] = (*emin) ;
    }
}
