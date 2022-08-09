static inline
void SumElemFaceNormalt(double *normalX0, double *normalY0, double *normalZ0,
                       double *normalX1, double *normalY1, double *normalZ1,
                       double *normalX2, double *normalY2, double *normalZ2,
                       double *normalX3, double *normalY3, double *normalZ3,
                       const double x0, const double y0, const double z0,
                       const double x1, const double y1, const double z1,
                       const double x2, const double y2, const double z2,
                       const double x3, const double y3, const double z3)
{
   double bisectX0 = double(0.5) * (x3 + x2 - x1 - x0);
   double bisectY0 = double(0.5) * (y3 + y2 - y1 - y0);
   double bisectZ0 = double(0.5) * (z3 + z2 - z1 - z0);
   double bisectX1 = double(0.5) * (x2 + x1 - x3 - x0);
   double bisectY1 = double(0.5) * (y2 + y1 - y3 - y0);
   double bisectZ1 = double(0.5) * (z2 + z1 - z3 - z0);
   double areaX = double(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
   double areaY = double(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
   double areaZ = double(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

   *normalX0 += areaX;
   *normalX1 += areaX;
   *normalX2 += areaX;
   *normalX3 += areaX;

   *normalY0 += areaY;
   *normalY1 += areaY;
   *normalY2 += areaY;
   *normalY3 += areaY;

   *normalZ0 += areaZ;
   *normalZ1 += areaZ;
   *normalZ2 += areaZ;
   *normalZ3 += areaZ;
}

inline void IntegrateStressForElemsLoop(
                                        const double *p_x1, const double *p_x2, const double *p_x3, const double *p_x4, const double *p_x5, const double *p_x6, const double *p_x7, const double *p_x8,
                                        const double *p_y1, const double *p_y2, const double *p_y3, const double *p_y4, const double *p_y5, const double *p_y6, const double *p_y7, const double *p_y8,
                                        const double *p_z1, const double *p_z2, const double *p_z3, const double *p_z4, const double *p_z5, const double *p_z6, const double *p_z7, const double *p_z8,
                                        double *p_fx1, double *p_fx2, double *p_fx3, double *p_fx4, double *p_fx5, double *p_fx6, double *p_fx7, double *p_fx8,
                                        double *p_fy1, double *p_fy2, double *p_fy3, double *p_fy4, double *p_fy5, double *p_fy6, double *p_fy7, double *p_fy8,
                                        double *p_fz1, double *p_fz2, double *p_fz3, double *p_fz4, double *p_fz5, double *p_fz6, double *p_fz7, double *p_fz8,
                                        double *volume,
                                        const double *sigxx, const double *sigyy, const double *sigzz){
    double b[3][8] ;// shape function derivatives
    double x_local[8] ;
    double y_local[8] ;
    double z_local[8] ;

    double fx_local[8] ;
    double fy_local[8] ;
    double fz_local[8] ;

    // int nd0i = nodelist[iteration * 8 + 0] ;
    // int nd1i = elemToNode[1] ;
    // int nd2i = elemToNode[2] ;
    // int nd3i = elemToNode[3] ;
    // int nd4i = elemToNode[4] ;
    // int nd5i = elemToNode[5] ;
    // int nd6i = elemToNode[6] ;
    // int nd7i = elemToNode[7] ;

    x_local[0] = p_x1[0];
    x_local[1] = p_x2[0];
    x_local[2] = p_x3[0];
    x_local[3] = p_x4[0];
    x_local[4] = p_x5[0];
    x_local[5] = p_x6[0];
    x_local[6] = p_x7[0];
    x_local[7] = p_x8[0];
    // std::cout << p_x1[0];

    y_local[0] = p_y1[0];
    y_local[1] = p_y2[0];
    y_local[2] = p_y3[0];
    y_local[3] = p_y4[0];
    y_local[4] = p_y5[0];
    y_local[5] = p_y6[0];
    y_local[6] = p_y7[0];
    y_local[7] = p_y8[0];

    z_local[0] = p_z1[0];
    z_local[1] = p_z2[0];
    z_local[2] = p_z3[0];
    z_local[3] = p_z4[0];
    z_local[4] = p_z5[0];
    z_local[5] = p_z6[0];
    z_local[6] = p_z7[0];
    z_local[7] = p_z8[0];

    const double x0 = x_local[0] ;   const double x1 = x_local[1] ;
    const double x2 = x_local[2] ;   const double x3 = x_local[3] ;
    const double x4 = x_local[4] ;   const double x5 = x_local[5] ;
    const double x6 = x_local[6] ;   const double x7 = x_local[7] ;

    const double y0 = y_local[0] ;   const double y1 = y_local[1] ;
    const double y2 = y_local[2] ;   const double y3 = y_local[3] ;
    const double y4 = y_local[4] ;   const double y5 = y_local[5] ;
    const double y6 = y_local[6] ;   const double y7 = y_local[7] ;

    const double z0 = z_local[0] ;   const double z1 = z_local[1] ;
    const double z2 = z_local[2] ;   const double z3 = z_local[3] ;
    const double z4 = z_local[4] ;   const double z5 = z_local[5] ;
    const double z6 = z_local[6] ;   const double z7 = z_local[7] ;

    double fjxxi, fjxet, fjxze;
    double fjyxi, fjyet, fjyze;
    double fjzxi, fjzet, fjzze;
    double cjxxi, cjxet, cjxze;
    double cjyxi, cjyet, cjyze;
    double cjzxi, cjzet, cjzze;

    fjxxi = double(.125) * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
    fjxet = double(.125) * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
    fjxze = double(.125) * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

    fjyxi = double(.125) * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
    fjyet = double(.125) * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
    fjyze = double(.125) * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

    fjzxi = double(.125) * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
    fjzet = double(.125) * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
    fjzze = double(.125) * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

    /* compute cofactors */
    cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
    cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
    cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

    cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
    cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
    cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

    cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
    cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
    cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

    /* calculate partials :
        this need only be done for l = 0,1,2,3   since , by symmetry ,
        (6,7,4,5) = - (0,1,2,3) .
    */
    b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
    b[0][1] =      cjxxi  -  cjxet  -  cjxze;
    b[0][2] =      cjxxi  +  cjxet  -  cjxze;
    b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
    b[0][4] = -b[0][2];
    b[0][5] = -b[0][3];
    b[0][6] = -b[0][0];
    b[0][7] = -b[0][1];

    b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
    b[1][1] =      cjyxi  -  cjyet  -  cjyze;
    b[1][2] =      cjyxi  +  cjyet  -  cjyze;
    b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
    b[1][4] = -b[1][2];
    b[1][5] = -b[1][3];
    b[1][6] = -b[1][0];
    b[1][7] = -b[1][1];

    b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
    b[2][1] =      cjzxi  -  cjzet  -  cjzze;
    b[2][2] =      cjzxi  +  cjzet  -  cjzze;
    b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
    b[2][4] = -b[2][2];
    b[2][5] = -b[2][3];
    b[2][6] = -b[2][0];
    b[2][7] = -b[2][1];

    /* calculate jacobian determinant (volume) */
    volume[0] = double(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
    // std::cout << "Volume: " << *volume << "\n";

    for (int i = 0 ; i < 8 ; ++i) {
      b[0][i] = double(0.0);
      b[1][i] = double(0.0);
      b[2][i] = double(0.0);
    }

    //Start Sum ElemFace Normals
    double bisectX0, bisectY0, bisectZ0, bisectX1, bisectY1, bisectZ1;
    double areaX, areaY, areaZ;

    /* evaluate face one: nodes 0, 1, 2, 3 */
    // SumElemFaceNormalt(&b[0][0], &b[1][0], &b[2][0],
    //                 &b[0][1], &b[1][1], &b[2][1],
    //                 &b[0][2], &b[1][2], &b[2][2],
    //                 &b[0][3], &b[1][3], &b[2][3],
    //                 x_local[0], y_local[0], z_local[0], x_local[1], y_local[1], z_local[1],
    //                 x_local[2], y_local[2], z_local[2], x_local[3], y_local[3], z_local[3]);
    bisectX0 = double(0.5) * (x_local[3] + x_local[2] - x_local[1] - x_local[0]);
    bisectY0 = double(0.5) * (y_local[3] + y_local[2] - y_local[1] - y_local[0]);
    bisectZ0 = double(0.5) * (z_local[3] + z_local[2] - z_local[1] - z_local[0]);
    bisectX1 = double(0.5) * (x_local[2] + x_local[1] - x_local[3] - x_local[0]);
    bisectY1 = double(0.5) * (y_local[2] + y_local[1] - y_local[3] - y_local[0]);
    bisectZ1 = double(0.5) * (z_local[2] + z_local[1] - z_local[3] - z_local[0]);
    areaX = double(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    areaY = double(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    areaZ = double(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    b[0][0] += areaX;
    b[0][1] += areaX;
    b[0][2] += areaX;
    b[0][3] += areaX;

    b[1][0] += areaY;
    b[1][1] += areaY;
    b[1][2] += areaY;
    b[1][3] += areaY;

    b[2][0] += areaZ;
    b[2][1] += areaZ;
    b[2][2] += areaZ;
    b[2][3] += areaZ;
    /* evaluate face two: nodes 0, 4, 5, 1 */
    // SumElemFaceNormalt(&b[0][0], &b[1][0], &b[2][0],
    //                 &b[0][4], &b[1][4], &b[2][4],
    //                 &b[0][5], &b[1][5], &b[2][5],
    //                 &b[0][1], &b[1][1], &b[2][1],
    //                 x_local[0], y_local[0], z_local[0], x_local[4], y_local[4], z_local[4],
    //                 x_local[5], y_local[5], z_local[5], x_local[1], y_local[1], z_local[1]);
    bisectX0 = double(0.5) * (x_local[1] + x_local[5] - x_local[4] - x_local[0]);
    bisectY0 = double(0.5) * (y_local[1] + y_local[5] - y_local[4] - y_local[0]);
    bisectZ0 = double(0.5) * (z_local[1] + z_local[5] - z_local[4] - z_local[0]);
    bisectX1 = double(0.5) * (x_local[5] + x_local[4] - x_local[1] - x_local[0]);
    bisectY1 = double(0.5) * (y_local[5] + y_local[4] - y_local[1] - y_local[0]);
    bisectZ1 = double(0.5) * (z_local[5] + z_local[4] - z_local[1] - z_local[0]);
    areaX = double(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    areaY = double(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    areaZ = double(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    b[0][0] += areaX;
    b[0][4] += areaX;
    b[0][5] += areaX;
    b[0][1] += areaX;

    b[1][0] += areaY;
    b[1][4] += areaY;
    b[1][5] += areaY;
    b[1][1] += areaY;

    b[2][0] += areaZ;
    b[2][4] += areaZ;
    b[2][5] += areaZ;
    b[2][1] += areaZ;
    /* evaluate face three: nodes 1, 5, 6, 2 */
    // SumElemFaceNormalt(&b[0][1], &b[1][1], &b[2][1],
    //                 &b[0][5], &b[1][5], &b[2][5],
    //                 &b[0][6], &b[1][6], &b[2][6],
    //                 &b[0][2], &b[1][2], &b[2][2],
    //                 x_local[1], y_local[1], z_local[1], x_local[5], y_local[5], z_local[5],
    //                 x_local[6], y_local[6], z_local[6], x_local[2], y_local[2], z_local[2]);
    bisectX0 = double(0.5) * (x_local[2] + x_local[6] - x_local[5] - x_local[1]);
    bisectY0 = double(0.5) * (y_local[2] + y_local[6] - y_local[5] - y_local[1]);
    bisectZ0 = double(0.5) * (z_local[2] + z_local[6] - z_local[5] - z_local[1]);
    bisectX1 = double(0.5) * (x_local[6] + x_local[5] - x_local[2] - x_local[1]);
    bisectY1 = double(0.5) * (y_local[6] + y_local[5] - y_local[2] - y_local[1]);
    bisectZ1 = double(0.5) * (z_local[6] + z_local[5] - z_local[2] - z_local[1]);
    areaX = double(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    areaY = double(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    areaZ = double(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    b[0][1] += areaX;
    b[0][5] += areaX;
    b[0][6] += areaX;
    b[0][2] += areaX;

    b[1][1] += areaY;
    b[1][5] += areaY;
    b[1][6] += areaY;
    b[1][2] += areaY;

    b[2][1] += areaZ;
    b[2][5] += areaZ;
    b[2][6] += areaZ;
    b[2][2] += areaZ;
    /* evaluate face four: nodes 2, 6, 7, 3 */
    // SumElemFaceNormalt(&b[0][2], &b[1][2], &b[2][2],
    //                 &b[0][6], &b[1][6], &b[2][6],
    //                 &b[0][7], &b[1][7], &b[2][7],
    //                 &b[0][3], &b[1][3], &b[2][3],
    //                 x_local[2], y_local[2], z_local[2], x_local[6], y_local[6], z_local[6],
    //                 x_local[7], y_local[7], z_local[7], x_local[3], y_local[3], z_local[3]);
    bisectX0 = double(0.5) * (x_local[3] + x_local[7] - x_local[6] - x_local[2]);
    bisectY0 = double(0.5) * (y_local[3] + y_local[7] - y_local[6] - y_local[2]);
    bisectZ0 = double(0.5) * (z_local[3] + z_local[7] - z_local[6] - z_local[2]);
    bisectX1 = double(0.5) * (x_local[7] + x_local[6] - x_local[3] - x_local[2]);
    bisectY1 = double(0.5) * (y_local[7] + y_local[6] - y_local[3] - y_local[2]);
    bisectZ1 = double(0.5) * (z_local[7] + z_local[6] - z_local[3] - z_local[2]);
    areaX = double(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    areaY = double(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    areaZ = double(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    b[0][2] += areaX;
    b[0][6] += areaX;
    b[0][7] += areaX;
    b[0][3] += areaX;

    b[1][2] += areaY;
    b[1][6] += areaY;
    b[1][7] += areaY;
    b[1][3] += areaY;

    b[2][2] += areaZ;
    b[2][6] += areaZ;
    b[2][7] += areaZ;
    b[2][3] += areaZ;
    /* evaluate face five: nodes 3, 7, 4, 0 */
    // SumElemFaceNormalt(&b[0][3], &b[1][3], &b[2][3],
    //                 &b[0][7], &b[1][7], &b[2][7],
    //                 &b[0][4], &b[1][4], &b[2][4],
    //                 &b[0][0], &b[1][0], &b[2][0],
    //                 x_local[3], y_local[3], z_local[3], x_local[7], y_local[7], z_local[7],
    //                 x_local[4], y_local[4], z_local[4], x_local[0], y_local[0], z_local[0]);
    bisectX0 = double(0.5) * (x_local[0] + x_local[4] - x_local[7] - x_local[3]);
    bisectY0 = double(0.5) * (y_local[0] + y_local[4] - y_local[7] - y_local[3]);
    bisectZ0 = double(0.5) * (z_local[0] + z_local[4] - z_local[7] - z_local[3]);
    bisectX1 = double(0.5) * (x_local[4] + x_local[7] - x_local[0] - x_local[3]);
    bisectY1 = double(0.5) * (y_local[4] + y_local[7] - y_local[0] - y_local[3]);
    bisectZ1 = double(0.5) * (z_local[4] + z_local[7] - z_local[0] - z_local[3]);
    areaX = double(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    areaY = double(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    areaZ = double(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    b[0][3] += areaX;
    b[0][7] += areaX;
    b[0][4] += areaX;
    b[0][0] += areaX;

    b[1][3] += areaY;
    b[1][7] += areaY;
    b[1][4] += areaY;
    b[1][0] += areaY;

    b[2][3] += areaZ;
    b[2][7] += areaZ;
    b[2][4] += areaZ;
    b[2][0] += areaZ;
    /* evaluate face six: nodes 4, 7, 6, 5 */
    // SumElemFaceNormalt(&b[0][4], &b[1][4], &b[2][4],
    //                 &b[0][7], &b[1][7], &b[2][7],
    //                 &b[0][6], &b[1][6], &b[2][6],
    //                 &b[0][5], &b[1][5], &b[2][5],
    //                 x_local[4], y_local[4], z_local[4], x_local[7], y_local[7], z_local[7],
    //                 x_local[6], y_local[6], z_local[6], x_local[5], y_local[5], z_local[5]);
    bisectX0 = double(0.5) * (x_local[5] + x_local[6] - x_local[7] - x_local[4]);
    bisectY0 = double(0.5) * (y_local[5] + y_local[6] - y_local[7] - y_local[4]);
    bisectZ0 = double(0.5) * (z_local[5] + z_local[6] - z_local[7] - z_local[4]);
    bisectX1 = double(0.5) * (x_local[6] + x_local[7] - x_local[5] - x_local[4]);
    bisectY1 = double(0.5) * (y_local[6] + y_local[7] - y_local[5] - y_local[4]);
    bisectZ1 = double(0.5) * (z_local[6] + z_local[7] - z_local[5] - z_local[4]);
    areaX = double(0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
    areaY = double(0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
    areaZ = double(0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

    b[0][4] += areaX;
    b[0][7] += areaX;
    b[0][6] += areaX;
    b[0][5] += areaX;

    b[1][4] += areaY;
    b[1][7] += areaY;
    b[1][6] += areaY;
    b[1][5] += areaY;

    b[2][4] += areaZ;
    b[2][7] += areaZ;
    b[2][6] += areaZ;
    b[2][5] += areaZ;

    // for(int i = 0; i < 8; i++) {
    //     fx_local[i] = -( t_sigxx[0] * b[0][i] );
    //     fy_local[i] = -( t_sigxx[1] * b[1][i]  );
    //     fz_local[i] = -( t_sigxx[2] * b[2][i] );
    // }

    for(int i = 0; i < 8; i++) {
        fx_local[i] = -( sigxx[0] * b[0][i] );
        fy_local[i] = -( sigyy[0] * b[1][i]  );
        fz_local[i] = -( sigzz[0] * b[2][i] );
    }
    // for( int lnode=0 ; lnode<8 ; ++lnode ) {
    //     int gnode = elemToNode[lnode];
    //     //  domain.fx(gnode) += fx_local[lnode];
    //     //  domain.fy(gnode) += fy_local[lnode];
    //     //  domain.fz(gnode) += fz_local[lnode];
    //     m_fx[gnode] += x_local[lnode];
    //     m_fy[gnode] += y_local[lnode];
    //     m_fz[gnode] += z_local[lnode];
    //    }
    p_fx1[0] += fx_local[0];
    p_fx2[0] += fx_local[1];
    p_fx3[0] += fx_local[2];
    p_fx4[0] += fx_local[3];
    p_fx5[0] += fx_local[4];
    p_fx6[0] += fx_local[5];
    p_fx7[0] += fx_local[6];
    p_fx8[0] += fx_local[7];

    p_fy1[0] += fy_local[0];
    p_fy2[0] += fy_local[1];
    p_fy3[0] += fy_local[2];
    p_fy4[0] += fy_local[3];
    p_fy5[0] += fy_local[4];
    p_fy6[0] += fy_local[5];
    p_fy7[0] += fy_local[6];
    p_fy8[0] += fy_local[7];

    p_fz1[0] += fz_local[0];
    p_fz2[0] += fz_local[1];
    p_fz3[0] += fz_local[2];
    p_fz4[0] += fz_local[3];
    p_fz5[0] += fz_local[4];
    p_fz6[0] += fz_local[5];
    p_fz7[0] += fz_local[6];
    p_fz8[0] += fz_local[7];
}