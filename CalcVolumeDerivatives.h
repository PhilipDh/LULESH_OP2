static inline
void VoluDert(const Real_t x0, const Real_t x1, const Real_t x2,
             const Real_t x3, const Real_t x4, const Real_t x5,
             const Real_t y0, const Real_t y1, const Real_t y2,
             const Real_t y3, const Real_t y4, const Real_t y5,
             const Real_t z0, const Real_t z1, const Real_t z2,
             const Real_t z3, const Real_t z4, const Real_t z5,
             Real_t* dvdx, Real_t* dvdy, Real_t* dvdz)
{
   const Real_t twelfth = Real_t(1.0) / Real_t(12.0) ;

   *dvdx =
      (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
      (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
      (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
   *dvdy =
      - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
      (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
      (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);

   *dvdz =
      - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
      (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
      (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);

   *dvdx *= twelfth;
   *dvdy *= twelfth;
   *dvdz *= twelfth;
}

inline void CalcVolumeDerivatives(
    double *p_x1, double *p_x2, double *p_x3, double *p_x4, double *p_x5, double *p_x6, double *p_x7, double *p_x8,
    double *p_y1, double *p_y2, double *p_y3, double *p_y4, double *p_y5, double *p_y6, double *p_y7, double *p_y8,
    double *p_z1, double *p_z2, double *p_z3, double *p_z4, double *p_z5, double *p_z6, double *p_z7, double *p_z8,
    double *p_dvdx,
    double *p_dvdy,
    double *p_dvdz,
    double *p_x8n,
    double *p_y8n,
    double *p_z8n,
    double *p_v, double *p_determ, double *p_volo
){
    Real_t  x1[8],  y1[8],  z1[8] ;
    Real_t pfx[8], pfy[8], pfz[8] ;

    x1[0] = p_x1[0];
    x1[1] = p_x2[0];
    x1[2] = p_x3[0];
    x1[3] = p_x4[0];
    x1[4] = p_x5[0];
    x1[5] = p_x6[0];
    x1[6] = p_x7[0];
    x1[7] = p_x8[0];

    y1[0] = p_y1[0];
    y1[1] = p_y2[0];
    y1[2] = p_y3[0];
    y1[3] = p_y4[0];
    y1[4] = p_y5[0];
    y1[5] = p_y6[0];
    y1[6] = p_y7[0];
    y1[7] = p_y8[0];

    z1[0] = p_z1[0];
    z1[1] = p_z2[0];
    z1[2] = p_z3[0];
    z1[3] = p_z4[0];
    z1[4] = p_z5[0];
    z1[5] = p_z6[0];
    z1[6] = p_z7[0];
    z1[7] = p_z8[0];

    VoluDert(x1[1], x1[2], x1[3], x1[4], x1[5], x1[7],
            y1[1], y1[2], y1[3], y1[4], y1[5], y1[7],
            z1[1], z1[2], z1[3], z1[4], z1[5], z1[7],
            &p_dvdx[0], &p_dvdy[0], &p_dvdz[0]);
    VoluDert(x1[0], x1[1], x1[2], x1[7], x1[4], x1[6],
            y1[0], y1[1], y1[2], y1[7], y1[4], y1[6],
            z1[0], z1[1], z1[2], z1[7], z1[4], z1[6],
            &p_dvdx[3], &p_dvdy[3], &p_dvdz[3]);
    VoluDert(x1[3], x1[0], x1[1], x1[6], x1[7], x1[5],
            y1[3], y1[0], y1[1], y1[6], y1[7], y1[5],
            z1[3], z1[0], z1[1], z1[6], z1[7], z1[5],
            &p_dvdx[2], &p_dvdy[2], &p_dvdz[2]);
    VoluDert(x1[2], x1[3], x1[0], x1[5], x1[6], x1[4],
            y1[2], y1[3], y1[0], y1[5], y1[6], y1[4],
            z1[2], z1[3], z1[0], z1[5], z1[6], z1[4],
            &p_dvdx[1], &p_dvdy[1], &p_dvdz[1]);
    VoluDert(x1[7], x1[6], x1[5], x1[0], x1[3], x1[1],
            y1[7], y1[6], y1[5], y1[0], y1[3], y1[1],
            z1[7], z1[6], z1[5], z1[0], z1[3], z1[1],
            &p_dvdx[4], &p_dvdy[4], &p_dvdz[4]);
    VoluDert(x1[4], x1[7], x1[6], x1[1], x1[0], x1[2],
            y1[4], y1[7], y1[6], y1[1], y1[0], y1[2],
            z1[4], z1[7], z1[6], z1[1], z1[0], z1[2],
            &p_dvdx[5], &p_dvdy[5], &p_dvdz[5]);
    VoluDert(x1[5], x1[4], x1[7], x1[2], x1[1], x1[3],
            y1[5], y1[4], y1[7], y1[2], y1[1], y1[3],
            z1[5], z1[4], z1[7], z1[2], z1[1], z1[3],
            &p_dvdx[6], &p_dvdy[6], &p_dvdz[6]);
    VoluDert(x1[6], x1[5], x1[4], x1[3], x1[2], x1[0],
            y1[6], y1[5], y1[4], y1[3], y1[2], y1[0],
            z1[6], z1[5], z1[4], z1[3], z1[2], z1[0],
            &p_dvdx[7], &p_dvdy[7], &p_dvdz[7]);
          /* load into temporary storage for FB Hour Glass control */
    for(Index_t ii=0;ii<8;++ii){
        // Index_t jj=8*i+ii;
        // p_dvdx[ii] = pfx[ii];
        // p_dvdy[ii] = pfy[ii];
        // p_dvdz[ii] = pfz[ii];
        p_x8n[ii]  = x1[ii];
        p_y8n[ii]  = y1[ii];
        p_z8n[ii]  = z1[ii];
        // std::cout<< p_dvdx[ii] << "\n";
    }

    p_determ[0] = p_volo[0] * p_v[0];

    if(p_v[0] <= 0.0){
        exit(VolumeError);
    }
}