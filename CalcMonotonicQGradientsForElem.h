inline void CalcMonotonicQGradientsForElem(
    const double *p_x1, const double *p_x2, const double *p_x3, const double *p_x4, const double *p_x5, const double *p_x6, const double *p_x7, const double *p_x8,
    const double *p_y1, const double *p_y2, const double *p_y3, const double *p_y4, const double *p_y5, const double *p_y6, const double *p_y7, const double *p_y8,
    const double *p_z1, const double *p_z2, const double *p_z3, const double *p_z4, const double *p_z5, const double *p_z6, const double *p_z7, const double *p_z8,
    const double *p_xd1, const double *p_xd2, const double *p_xd3, const double *p_xd4, const double *p_xd5, const double *p_xd6, const double *p_xd7, const double *p_xd8,
    const double *p_yd1, const double *p_yd2, const double *p_yd3, const double *p_yd4, const double *p_yd5, const double *p_yd6, const double *p_yd7, const double *p_yd8,
    const double *p_zd1, const double *p_zd2, const double *p_zd3, const double *p_zd4, const double *p_zd5, const double *p_zd6, const double *p_zd7, const double *p_zd8,
    const double *volo,
    const double *vnew,
    double *delx_zeta,
    double *delv_zeta,
    double *delv_xi,
    double *delx_xi,
    double *delx_eta,
    double *delv_eta
){
    const double ptiny = double(1.e-36);
    double ax,ay,az;
    double dxv,dyv,dzv;

    double x0 = p_x1[0];
    double x1 = p_x2[0];
    double x2 = p_x3[0];
    double x3 = p_x4[0];
    double x4 = p_x5[0];
    double x5 = p_x6[0];
    double x6 = p_x7[0];
    double x7 = p_x8[0]; 

    double y0 = p_y1[0];
    double y1 = p_y2[0];
    double y2 = p_y3[0];
    double y3 = p_y4[0];
    double y4 = p_y5[0];
    double y5 = p_y6[0];
    double y6 = p_y7[0];
    double y7 = p_y8[0];

    double z0 = p_z1[0];
    double z1 = p_z2[0];
    double z2 = p_z3[0];
    double z3 = p_z4[0];
    double z4 = p_z5[0];
    double z5 = p_z6[0];
    double z6 = p_z7[0];
    double z7 = p_z8[0];

    double xv0 = p_xd1[0];
    double xv1 = p_xd2[0];
    double xv2 = p_xd3[0];
    double xv3 = p_xd4[0];
    double xv4 = p_xd5[0];
    double xv5 = p_xd6[0];
    double xv6 = p_xd7[0];
    double xv7 = p_xd8[0]; 

    double yv0 = p_yd1[0];
    double yv1 = p_yd2[0];
    double yv2 = p_yd3[0];
    double yv3 = p_yd4[0];
    double yv4 = p_yd5[0];
    double yv5 = p_yd6[0];
    double yv6 = p_yd7[0];
    double yv7 = p_yd8[0];

    double zv0 = p_zd1[0];
    double zv1 = p_zd2[0];
    double zv2 = p_zd3[0];
    double zv3 = p_zd4[0];
    double zv4 = p_zd5[0];
    double zv5 = p_zd6[0];
    double zv6 = p_zd7[0];
    double zv7 = p_zd8[0];

    double vol = volo[0]*vnew[0];
    double norm = double(1.0) / ( vol + ptiny );

    double dxj = double(-0.25)*((x0+x1+x5+x4) - (x3+x2+x6+x7)) ;
    double dyj = double(-0.25)*((y0+y1+y5+y4) - (y3+y2+y6+y7)) ;
    double dzj = double(-0.25)*((z0+z1+z5+z4) - (z3+z2+z6+z7)) ;

    double dxi = double( 0.25)*((x1+x2+x6+x5) - (x0+x3+x7+x4)) ;
    double dyi = double( 0.25)*((y1+y2+y6+y5) - (y0+y3+y7+y4)) ;
    double dzi = double( 0.25)*((z1+z2+z6+z5) - (z0+z3+z7+z4)) ;

    double dxk = double( 0.25)*((x4+x5+x6+x7) - (x0+x1+x2+x3)) ;
    double dyk = double( 0.25)*((y4+y5+y6+y7) - (y0+y1+y2+y3)) ;
    double dzk = double( 0.25)*((z4+z5+z6+z7) - (z0+z1+z2+z3)) ;

    /* find delvk and delxk ( i cross j ) */
    ax = dyi*dzj - dzi*dyj ;
    ay = dzi*dxj - dxi*dzj ;
    az = dxi*dyj - dyi*dxj ;

    delx_zeta[0] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

    dxv = double(0.25)*((xv4+xv5+xv6+xv7) - (xv0+xv1+xv2+xv3)) ;
    dyv = double(0.25)*((yv4+yv5+yv6+yv7) - (yv0+yv1+yv2+yv3)) ;
    dzv = double(0.25)*((zv4+zv5+zv6+zv7) - (zv0+zv1+zv2+zv3)) ;

    // domain.delv_zeta(i) = ax*dxv + ay*dyv + az*dzv ;
    delv_zeta[0] = ax*dxv + ay*dyv + az*dzv ;

    /* find delxi and delvi ( j cross k ) */
    ax = dyj*dzk - dzj*dyk ;
    ay = dzj*dxk - dxj*dzk ;
    az = dxj*dyk - dyj*dxk ;

    // domain.delx_xi(i) = vol / sqrt(ax*ax + ay*ay + az*az + ptiny) ;
    delx_xi[0] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

    dxv = double(0.25)*((xv1+xv2+xv6+xv5) - (xv0+xv3+xv7+xv4)) ;
    dyv = double(0.25)*((yv1+yv2+yv6+yv5) - (yv0+yv3+yv7+yv4)) ;
    dzv = double(0.25)*((zv1+zv2+zv6+zv5) - (zv0+zv3+zv7+zv4)) ;

    // domain.delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;
    delv_xi[0] = ax*dxv + ay*dyv + az*dzv ;

    /* find delxj and delvj ( k cross i ) */

    ax = dyk*dzi - dzk*dyi ;
    ay = dzk*dxi - dxk*dzi ;
    az = dxk*dyi - dyk*dxi ;

    // domain.delx_eta(i) = vol / sqrt(ax*ax + ay*ay + az*az + ptiny) ;
    delx_eta[0] = vol / sqrt(ax*ax + ay*ay + az*az + ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

    dxv = double(-0.25)*((xv0+xv1+xv5+xv4) - (xv3+xv2+xv6+xv7)) ;
    dyv = double(-0.25)*((yv0+yv1+yv5+yv4) - (yv3+yv2+yv6+yv7)) ;
    dzv = double(-0.25)*((zv0+zv1+zv5+zv4) - (zv3+zv2+zv6+zv7)) ;

    // domain.delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
    delv_eta[0] = ax*dxv + ay*dyv + az*dzv ;

}