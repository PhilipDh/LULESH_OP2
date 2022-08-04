inline void CalcVeloForNodes(
                            double *xd, double *yd, double *zd,
                            double *xdd, double *ydd, double *zdd,
                            double *dt
){
    Real_t xdtmp, ydtmp, zdtmp ;

//   xdtmp = domain.xd(i) + domain.xdd(i) * dt ;
    xdtmp = xd[0] + xdd[0] * (*dt) ;
    if( FABS(xdtmp) < m_u_cut ) xdtmp = Real_t(0.0);
//   domain.xd(i) = xdtmp ;
    xd[0] = xdtmp ;

//   ydtmp = domain.yd(i) + domain.ydd(i) * dt ;
    ydtmp = yd[0] + ydd[0] * (*dt) ;
    if( FABS(ydtmp) < m_u_cut ) ydtmp = Real_t(0.0);
//   domain.yd(i) = ydtmp ;
    yd[0] = ydtmp ;

//   zdtmp = domain.zd(i) + domain.zdd(i) * dt ;
    zdtmp = zd[0] + zdd[0] * (*dt) ;
    if( FABS(zdtmp) < m_u_cut ) zdtmp = Real_t(0.0);
//   domain.zd(i) = zdtmp ;
    zd[0] = zdtmp ;
}