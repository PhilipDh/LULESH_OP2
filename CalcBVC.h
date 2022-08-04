inline void CalcBVC(double *bvc, double *compression, double *pbvc){
    Real_t c1s = Real_t(2.0)/Real_t(3.0) ;
    bvc[0] = c1s * (compression[0] + Real_t(1.));
    pbvc[0] = c1s;
}