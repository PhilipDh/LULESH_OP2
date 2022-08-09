static inline
void VoluDert(const double x0, const double x1, const double x2,
             const double x3, const double x4, const double x5,
             const double y0, const double y1, const double y2,
             const double y3, const double y4, const double y5,
             const double z0, const double z1, const double z2,
             const double z3, const double z4, const double z5,
             double* dvdx, double* dvdy, double* dvdz)
{
   const double twelfth = double(1.0) / double(12.0) ;

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
    const double *p_x1, const double *p_x2, const double *p_x3, const double *p_x4, const double *p_x5, const double *p_x6, const double *p_x7, const double *p_x8,
    const double *p_y1, const double *p_y2, const double *p_y3, const double *p_y4, const double *p_y5, const double *p_y6, const double *p_y7, const double *p_y8,
    const double *p_z1, const double *p_z2, const double *p_z3, const double *p_z4, const double *p_z5, const double *p_z6, const double *p_z7, const double *p_z8,
    double *p_dvdx,
    double *p_dvdy,
    double *p_dvdz,
    double *p_x8n,
    double *p_y8n,
    double *p_z8n,
    const double *p_v, double *p_determ, const double *p_volo
){

        double  x1[8],  y1[8],  z1[8] ;
        double pfx[8], pfy[8], pfz[8] ;

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

        // Formerly its own function "VolumeDer", moved to be inside this function for CUDA to work
        // VoluDert(x1[1], x1[2], x1[3], x1[4], x1[5], x1[7],
        //         y1[1], y1[2], y1[3], y1[4], y1[5], y1[7],
        //         z1[1], z1[2], z1[3], z1[4], z1[5], z1[7],
        //         &p_dvdx[0], &p_dvdy[0], &p_dvdz[0]);
        const double twelfth = double(1.0) / double(12.0) ;

        p_dvdx[0] =
        (y1[2] + y1[3]) * (z1[1] + z1[2]) - (y1[1] + y1[2]) * (z1[2] + z1[3]) +
        (y1[1] + y1[5]) * (z1[4] + z1[5]) - (y1[4] + y1[5]) * (z1[1] + z1[5]) -
        (y1[3] + y1[7]) * (z1[4] + z1[7]) + (y1[4] + y1[7]) * (z1[3] + z1[7]);
        p_dvdy[0] =
        - (x1[2] + x1[3]) * (z1[1] + z1[2]) + (x1[1] + x1[2]) * (z1[2] + z1[3]) -
        (x1[1] + x1[5]) * (z1[4] + z1[5]) + (x1[4] + x1[5]) * (z1[1] + z1[5]) +
        (x1[3] + x1[7]) * (z1[4] + z1[7]) - (x1[4] + x1[7]) * (z1[3] + z1[7]);

        p_dvdz[0] =
        - (y1[2] + y1[3]) * (x1[1] + x1[2]) + (y1[1] + y1[2]) * (x1[2] + x1[3]) -
        (y1[1] + y1[5]) * (x1[4] + x1[5]) + (y1[4] + y1[5]) * (x1[1] + x1[5]) +
        (y1[3] + y1[7]) * (x1[4] + x1[7]) - (y1[4] + y1[7]) * (x1[3] + x1[7]);

        p_dvdx[0] *= twelfth;
        p_dvdy[0] *= twelfth;
        p_dvdz[0] *= twelfth;
                
        // VoluDert(x1[0], x1[1], x1[2], x1[7], x1[4], x1[6],
        //         y1[0], y1[1], y1[2], y1[7], y1[4], y1[6],
        //         z1[0], z1[1], z1[2], z1[7], z1[4], z1[6],
        //         &p_dvdx[3], &p_dvdy[3], &p_dvdz[3]);
        p_dvdx[3] =
        (y1[1] + y1[2]) * (z1[0] + z1[1]) - (y1[0] + y1[1]) * (z1[1] + z1[2]) +
        (y1[0] + y1[4]) * (z1[7] + z1[4]) - (y1[7] + y1[4]) * (z1[0] + z1[4]) -
        (y1[2] + y1[6]) * (z1[7] + z1[6]) + (y1[7] + y1[6]) * (z1[2] + z1[6]);
        p_dvdy[3] =
        - (x1[1] + x1[2]) * (z1[0] + z1[1]) + (x1[0] + x1[1]) * (z1[1] + z1[2]) -
        (x1[0] + x1[4]) * (z1[7] + z1[4]) + (x1[7] + x1[4]) * (z1[0] + z1[4]) +
        (x1[2] + x1[6]) * (z1[7] + z1[6]) - (x1[7] + x1[6]) * (z1[2] + z1[6]);
        p_dvdz[3] =
        - (y1[1] + y1[2]) * (x1[0] + x1[1]) + (y1[0] + y1[1]) * (x1[1] + x1[2]) -
        (y1[0] + y1[4]) * (x1[7] + x1[4]) + (y1[7] + y1[4]) * (x1[0] + x1[4]) +
        (y1[2] + y1[6]) * (x1[7] + x1[6]) - (y1[7] + y1[6]) * (x1[2] + x1[6]);

        p_dvdx[3] *= twelfth;
        p_dvdy[3] *= twelfth;
        p_dvdz[3] *= twelfth;
                
        // VoluDert(x1[3], x1[0], x1[1], x1[6], x1[7], x1[5],
        //         y1[3], y1[0], y1[1], y1[6], y1[7], y1[5],
        //         z1[3], z1[0], z1[1], z1[6], z1[7], z1[5],
        //         &p_dvdx[2], &p_dvdy[2], &p_dvdz[2]);
        p_dvdx[2] =
        (y1[0] + y1[1]) * (z1[3] + z1[0]) - (y1[3] + y1[0]) * (z1[0] + z1[1]) +
        (y1[3] + y1[7]) * (z1[6] + z1[7]) - (y1[6] + y1[7]) * (z1[3] + z1[7]) -
        (y1[1] + y1[5]) * (z1[6] + z1[5]) + (y1[6] + y1[5]) * (z1[1] + z1[5]);
        p_dvdy[2] =
        - (x1[0] + x1[1]) * (z1[3] + z1[0]) + (x1[3] + x1[0]) * (z1[0] + z1[1]) -
        (x1[3] + x1[7]) * (z1[6] + z1[7]) + (x1[6] + x1[7]) * (z1[3] + z1[7]) +
        (x1[1] + x1[5]) * (z1[6] + z1[5]) - (x1[6] + x1[5]) * (z1[1] + z1[5]);

        p_dvdz[2] =
        - (y1[0] + y1[1]) * (x1[3] + x1[0]) + (y1[3] + y1[0]) * (x1[0] + x1[1]) -
        (y1[3] + y1[7]) * (x1[6] + x1[7]) + (y1[6] + y1[7]) * (x1[3] + x1[7]) +
        (y1[1] + y1[5]) * (x1[6] + x1[5]) - (y1[6] + y1[5]) * (x1[1] + x1[5]);

        p_dvdx[2] *= twelfth;
        p_dvdy[2] *= twelfth;
        p_dvdz[2] *= twelfth;

        // VoluDert(x1[2], x1[3], x1[0], x1[5], x1[6], x1[4],
        //         y1[2], y1[3], y1[0], y1[5], y1[6], y1[4],
        //         z1[2], z1[3], z1[0], z1[5], z1[6], z1[4],
        //         &p_dvdx[1], &p_dvdy[1], &p_dvdz[1]);
        p_dvdx[1] =
        (y1[3] + y1[0]) * (z1[2] + z1[3]) - (y1[2] + y1[3]) * (z1[3] + z1[0]) +
        (y1[2] + y1[6]) * (z1[5] + z1[6]) - (y1[5] + y1[6]) * (z1[2] + z1[6]) -
        (y1[0] + y1[4]) * (z1[5] + z1[4]) + (y1[5] + y1[4]) * (z1[0] + z1[4]);
        p_dvdy[1] =
        - (x1[3] + x1[0]) * (z1[2] + z1[3]) + (x1[2] + x1[3]) * (z1[3] + z1[0]) -
        (x1[2] + x1[6]) * (z1[5] + z1[6]) + (x1[5] + x1[6]) * (z1[2] + z1[6]) +
        (x1[0] + x1[4]) * (z1[5] + z1[4]) - (x1[5] + x1[4]) * (z1[0] + z1[4]);
        p_dvdz[1] =
        - (y1[3] + y1[0]) * (x1[2] + x1[3]) + (y1[2] + y1[3]) * (x1[3] + x1[0]) -
        (y1[2] + y1[6]) * (x1[5] + x1[6]) + (y1[5] + y1[6]) * (x1[2] + x1[6]) +
        (y1[0] + y1[4]) * (x1[5] + x1[4]) - (y1[5] + y1[4]) * (x1[0] + x1[4]);

        p_dvdx[1] *= twelfth;
        p_dvdy[1] *= twelfth;
        p_dvdz[1] *= twelfth;

        // VoluDert(x1[7], x1[6], x1[5], x1[0], x1[3], x1[1],
        //         y1[7], y1[6], y1[5], y1[0], y1[3], y1[1],
        //         z1[7], z1[6], z1[5], z1[0], z1[3], z1[1],
        //         &p_dvdx[4], &p_dvdy[4], &p_dvdz[4]);
        p_dvdx[4] =
        (y1[6] + y1[5]) * (z1[7] + z1[6]) - (y1[7] + y1[6]) * (z1[6] + z1[5]) +
        (y1[7] + y1[3]) * (z1[0] + z1[3]) - (y1[0] + y1[3]) * (z1[7] + z1[3]) -
        (y1[5] + y1[1]) * (z1[0] + z1[1]) + (y1[0] + y1[1]) * (z1[5] + z1[1]);
        p_dvdy[4] =
        - (x1[6] + x1[5]) * (z1[7] + z1[6]) + (x1[7] + x1[6]) * (z1[6] + z1[5]) -
        (x1[7] + x1[3]) * (z1[0] + z1[3]) + (x1[0] + x1[3]) * (z1[7] + z1[3]) +
        (x1[5] + x1[1]) * (z1[0] + z1[1]) - (x1[0] + x1[1]) * (z1[5] + z1[1]);

        p_dvdz[4] =
        - (y1[6] + y1[5]) * (x1[7] + x1[6]) + (y1[7] + y1[6]) * (x1[6] + x1[5]) -
        (y1[7] + y1[3]) * (x1[0] + x1[3]) + (y1[0] + y1[3]) * (x1[7] + x1[3]) +
        (y1[5] + y1[1]) * (x1[0] + x1[1]) - (y1[0] + y1[1]) * (x1[5] + x1[1]);

        p_dvdx[4] *= twelfth;
        p_dvdy[4] *= twelfth;
        p_dvdz[4] *= twelfth;

        // VoluDert(x1[4], x1[7], x1[6], x1[1], x1[0], x1[2],
        //         y1[4], y1[7], y1[6], y1[1], y1[0], y1[2],
        //         z1[4], z1[7], z1[6], z1[1], z1[0], z1[2],
        //         &p_dvdx[5], &p_dvdy[5], &p_dvdz[5]);
        p_dvdx[5] =
        (y1[7] + y1[6]) * (z1[4] + z1[7]) - (y1[4] + y1[7]) * (z1[7] + z1[6]) +
        (y1[4] + y1[0]) * (z1[1] + z1[0]) - (y1[1] + y1[0]) * (z1[4] + z1[0]) -
        (y1[6] + y1[2]) * (z1[1] + z1[2]) + (y1[1] + y1[2]) * (z1[6] + z1[2]);
        p_dvdy[5] =
        - (x1[7] + x1[6]) * (z1[4] + z1[7]) + (x1[4] + x1[7]) * (z1[7] + z1[6]) -
        (x1[4] + x1[0]) * (z1[1] + z1[0]) + (x1[1] + x1[0]) * (z1[4] + z1[0]) +
        (x1[6] + x1[2]) * (z1[1] + z1[2]) - (x1[1] + x1[2]) * (z1[6] + z1[2]);

        p_dvdz[5] =
        - (y1[7] + y1[6]) * (x1[4] + x1[7]) + (y1[4] + y1[7]) * (x1[7] + x1[6]) -
        (y1[4] + y1[0]) * (x1[1] + x1[0]) + (y1[1] + y1[0]) * (x1[4] + x1[0]) +
        (y1[6] + y1[2]) * (x1[1] + x1[2]) - (y1[1] + y1[2]) * (x1[6] + x1[2]);

        p_dvdx[5] *= twelfth;
        p_dvdy[5] *= twelfth;
        p_dvdz[5] *= twelfth;

        // VoluDert(x1[5], x1[4], x1[7], x1[2], x1[1], x1[3],
        //         y1[5], y1[4], y1[7], y1[2], y1[1], y1[3],
        //         z1[5], z1[4], z1[7], z1[2], z1[1], z1[3],
        //         &p_dvdx[6], &p_dvdy[6], &p_dvdz[6]);

        p_dvdx[6] =
        (y1[4] + y1[7]) * (z1[5] + z1[4]) - (y1[5] + y1[4]) * (z1[4] + z1[7]) +
        (y1[5] + y1[1]) * (z1[2] + z1[1]) - (y1[2] + y1[1]) * (z1[5] + z1[1]) -
        (y1[7] + y1[3]) * (z1[2] + z1[3]) + (y1[2] + y1[3]) * (z1[7] + z1[3]);
        p_dvdy[6] =
        - (x1[4] + x1[7]) * (z1[5] + z1[4]) + (x1[5] + x1[4]) * (z1[4] + z1[7]) -
        (x1[5] + x1[1]) * (z1[2] + z1[1]) + (x1[2] + x1[1]) * (z1[5] + z1[1]) +
        (x1[7] + x1[3]) * (z1[2] + z1[3]) - (x1[2] + x1[3]) * (z1[7] + z1[3]);

        p_dvdz[6] =
        - (y1[4] + y1[7]) * (x1[5] + x1[4]) + (y1[5] + y1[4]) * (x1[4] + x1[7]) -
        (y1[5] + y1[1]) * (x1[2] + x1[1]) + (y1[2] + y1[1]) * (x1[5] + x1[1]) +
        (y1[7] + y1[3]) * (x1[2] + x1[3]) - (y1[2] + y1[3]) * (x1[7] + x1[3]);
        p_dvdx[6] *= twelfth;
        p_dvdy[6] *= twelfth;
        p_dvdz[6] *= twelfth;

        // VoluDert(x1[6], x1[5], x1[4], x1[3], x1[2], x1[0],
        //         y1[6], y1[5], y1[4], y1[3], y1[2], y1[0],
        //         z1[6], z1[5], z1[4], z1[3], z1[2], z1[0],
        //         &p_dvdx[7], &p_dvdy[7], &p_dvdz[7]);
        p_dvdx[7] =
        (y1[5] + y1[4]) * (z1[6] + z1[5]) - (y1[6] + y1[5]) * (z1[5] + z1[4]) +
        (y1[6] + y1[2]) * (z1[3] + z1[2]) - (y1[3] + y1[2]) * (z1[6] + z1[2]) -
        (y1[4] + y1[0]) * (z1[3] + z1[0]) + (y1[3] + y1[0]) * (z1[4] + z1[0]);
        p_dvdy[7] =
        - (x1[5] + x1[4]) * (z1[6] + z1[5]) + (x1[6] + x1[5]) * (z1[5] + z1[4]) -
        (x1[6] + x1[2]) * (z1[3] + z1[2]) + (x1[3] + x1[2]) * (z1[6] + z1[2]) +
        (x1[4] + x1[0]) * (z1[3] + z1[0]) - (x1[3] + x1[0]) * (z1[4] + z1[0]);

        p_dvdz[7] =
        - (y1[5] + y1[4]) * (x1[6] + x1[5]) + (y1[6] + y1[5]) * (x1[5] + x1[4]) -
        (y1[6] + y1[2]) * (x1[3] + x1[2]) + (y1[3] + y1[2]) * (x1[6] + x1[2]) +
        (y1[4] + y1[0]) * (x1[3] + x1[0]) - (y1[3] + y1[0]) * (x1[4] + x1[0]);
        p_dvdx[7] *= twelfth;
        p_dvdy[7] *= twelfth;
        p_dvdz[7] *= twelfth;
                /* load into temporary storage for FB Hour Glass control */
        for(int ii=0;ii<8;++ii){
                // int jj=8*i+ii;
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
                // exit(-1); // Volume Error
        }
}