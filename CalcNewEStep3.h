inline void CalcNewEStep3(double *e_new, const double *work, const double *e_cut, const double *emin){
    e_new[0] += double(0.5) * work[0];

    if (fabs(e_new[0]) < (*e_cut)) {
        e_new[0] = double(0.)  ;
    }
    if (     e_new[0]  < (*emin) ) {
        e_new[0] = (*emin) ;
    }
}
