inline void CalcPosForNodes(
                            double *x, double *y, double *z,
                            const double *xd, const double *yd, const double *zd,
                            const double *dt
){
//   domain.x(i) += domain.xd(i) * dt ;
    x[0] += xd[0] * (*dt) ;
//   domain.y(i) += domain.yd(i) * dt ;
    y[0] += yd[0] * (*dt) ;
//   domain.z(i) += domain.zd(i) * dt ;
    z[0] += zd[0] * (*dt) ;
}