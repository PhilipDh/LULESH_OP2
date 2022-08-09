inline void CalcHalfStepBVC(double *bvc, const double *compression, double *pbvc){
    double c1s = double(2.0)/double(3.0) ;
    bvc[0] = c1s * (compression[0] + double(1.));
    pbvc[0] = c1s;
}
