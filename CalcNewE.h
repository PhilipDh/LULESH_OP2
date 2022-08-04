inline void CalcNewE(double *e_new, double *e_old, double *emin, double *delvc, double *p_old, double *q_old, double *work){
    e_new[0] = e_old[0] - Real_t(0.5) * delvc[0] * (p_old[0] + q_old[0])
        + Real_t(0.5) * work[0];

    if (e_new[0]  < (*emin) ) {
        e_new[0] = (*emin) ;
    }
}