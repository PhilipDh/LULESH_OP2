inline void CopyTempEOSVarsBack(
                                double *p, double *p_new,
                                double *e, double *e_new,
                                double *q, double *q_new){
    p[0] = p_new[0];
    e[0] = e_new[0];
    q[0] = q_new[0];  
}