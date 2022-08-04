inline void FBHourglassForceForElems(
    double *p_xd1, double *p_xd2, double *p_xd3, double *p_xd4, double *p_xd5, double *p_xd6, double *p_xd7, double *p_xd8,
    double *p_yd1, double *p_yd2, double *p_yd3, double *p_yd4, double *p_yd5, double *p_yd6, double *p_yd7, double *p_yd8,
    double *p_zd1, double *p_zd2, double *p_zd3, double *p_zd4, double *p_zd5, double *p_zd6, double *p_zd7, double *p_zd8,
    double *p_fx1, double *p_fx2, double *p_fx3, double *p_fx4, double *p_fx5, double *p_fx6, double *p_fx7, double *p_fx8,
    double *p_fy1, double *p_fy2, double *p_fy3, double *p_fy4, double *p_fy5, double *p_fy6, double *p_fy7, double *p_fy8,
    double *p_fz1, double *p_fz2, double *p_fz3, double *p_fz4, double *p_fz5, double *p_fz6, double *p_fz7, double *p_fz8,
    double *p_dvdx,
    double *p_dvdy,
    double *p_dvdz,
    double *p_x8n,
    double *p_y8n,
    double *p_z8n,
    double *p_determ,
    double *p_ss,
    double *p_elemMass,
    double *gamma_t
){

    Real_t hgfx[8], hgfy[8], hgfz[8] ;

    Real_t coefficient;

    Real_t hourgam[8][4];
    Real_t xd1[8], yd1[8], zd1[8] ;

    Real_t volinv=Real_t(1.0)/p_determ[0];
    Real_t ss1, mass1, volume13 ;


        for(Index_t i1=0;i1<4;++i1){
        Real_t hourmodx =
            p_x8n[0] * gamma_t[(i1*8)+0] + p_x8n[0+1] * gamma_t[(i1*8)+1] +
            p_x8n[0+2] * gamma_t[(i1*8)+2] + p_x8n[0+3] * gamma_t[(i1*8)+3] +
            p_x8n[0+4] * gamma_t[(i1*8)+4] + p_x8n[0+5] * gamma_t[(i1*8)+5] +
            p_x8n[0+6] * gamma_t[(i1*8)+6] + p_x8n[0+7] * gamma_t[(i1*8)+7];

        Real_t hourmody =
            p_y8n[0] * gamma_t[(i1*8)+0] + p_y8n[0+1] * gamma_t[(i1*8)+1] +
            p_y8n[0+2] * gamma_t[(i1*8)+2] + p_y8n[0+3] * gamma_t[(i1*8)+3] +
            p_y8n[0+4] * gamma_t[(i1*8)+4] + p_y8n[0+5] * gamma_t[(i1*8)+5] +
            p_y8n[0+6] * gamma_t[(i1*8)+6] + p_y8n[0+7] * gamma_t[(i1*8)+7];

        Real_t hourmodz =
            p_z8n[0] * gamma_t[(i1*8)+0] + p_z8n[0+1] * gamma_t[(i1*8)+1] +
            p_z8n[0+2] * gamma_t[(i1*8)+2] + p_z8n[0+3] * gamma_t[(i1*8)+3] +
            p_z8n[0+4] * gamma_t[(i1*8)+4] + p_z8n[0+5] * gamma_t[(i1*8)+5] +
            p_z8n[0+6] * gamma_t[(i1*8)+6] + p_z8n[0+7] * gamma_t[(i1*8)+7];

        hourgam[0][i1] = gamma_t[(i1*8)+0] -  volinv*(p_dvdx[0  ] * hourmodx +
                                                p_dvdy[0  ] * hourmody +
                                                p_dvdz[0  ] * hourmodz );

        hourgam[1][i1] = gamma_t[(i1*8)+1] -  volinv*(p_dvdx[0+1] * hourmodx +
                                                p_dvdy[0+1] * hourmody +
                                                p_dvdz[0+1] * hourmodz );

        hourgam[2][i1] = gamma_t[(i1*8)+2] -  volinv*(p_dvdx[0+2] * hourmodx +
                                                p_dvdy[0+2] * hourmody +
                                                p_dvdz[0+2] * hourmodz );

        hourgam[3][i1] = gamma_t[(i1*8)+3] -  volinv*(p_dvdx[0+3] * hourmodx +
                                                p_dvdy[0+3] * hourmody +
                                                p_dvdz[0+3] * hourmodz );

        hourgam[4][i1] = gamma_t[(i1*8)+4] -  volinv*(p_dvdx[0+4] * hourmodx +
                                                p_dvdy[0+4] * hourmody +
                                                p_dvdz[0+4] * hourmodz );

        hourgam[5][i1] = gamma_t[(i1*8)+5] -  volinv*(p_dvdx[0+5] * hourmodx +
                                                p_dvdy[0+5] * hourmody +
                                                p_dvdz[0+5] * hourmodz );

        hourgam[6][i1] = gamma_t[(i1*8)+6] -  volinv*(p_dvdx[0+6] * hourmodx +
                                                p_dvdy[0+6] * hourmody +
                                                p_dvdz[0+6] * hourmodz );

        hourgam[7][i1] = gamma_t[(i1*8)+7] -  volinv*(p_dvdx[0+7] * hourmodx +
                                                p_dvdy[0+7] * hourmody +
                                                p_dvdz[0+7] * hourmodz );
    }

    /* compute forces */
    /* store forces into h arrays (force arrays) */

    ss1=p_ss[0];

    mass1=p_elemMass[0];
    volume13=CBRT(p_determ[0]);

    xd1[0] = p_xd1[0];
    xd1[1] = p_xd2[0];
    xd1[2] = p_xd3[0];
    xd1[3] = p_xd4[0];
    xd1[4] = p_xd5[0];
    xd1[5] = p_xd6[0];
    xd1[6] = p_xd7[0];
    xd1[7] = p_xd8[0];

    yd1[0] = p_yd1[0];
    yd1[1] = p_yd2[0];
    yd1[2] = p_yd3[0];
    yd1[3] = p_yd4[0];
    yd1[4] = p_yd5[0];
    yd1[5] = p_yd6[0];
    yd1[6] = p_yd7[0];
    yd1[7] = p_yd8[0];

    zd1[0] = p_zd1[0];
    zd1[1] = p_zd2[0];
    zd1[2] = p_zd3[0];
    zd1[3] = p_zd4[0];
    zd1[4] = p_zd5[0];
    zd1[5] = p_zd6[0];
    zd1[6] = p_zd7[0];
    zd1[7] = p_zd8[0];

    coefficient = - m_hgcoef * Real_t(0.01) * ss1 * mass1 / volume13;

    Real_t hxx[4];
    for(Index_t i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * xd1[0] + hourgam[1][i] * xd1[1] +
                hourgam[2][i] * xd1[2] + hourgam[3][i] * xd1[3] +
                hourgam[4][i] * xd1[4] + hourgam[5][i] * xd1[5] +
                hourgam[6][i] * xd1[6] + hourgam[7][i] * xd1[7];
    }
    for(Index_t i = 0; i < 8; i++) {
        hgfx[i] = coefficient *
                    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                    hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }
    for(Index_t i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * yd1[0] + hourgam[1][i] * yd1[1] +
                hourgam[2][i] * yd1[2] + hourgam[3][i] * yd1[3] +
                hourgam[4][i] * yd1[4] + hourgam[5][i] * yd1[5] +
                hourgam[6][i] * yd1[6] + hourgam[7][i] * yd1[7];
    }
    for(Index_t i = 0; i < 8; i++) {
        hgfy[i] = coefficient *
                    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                    hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }
    for(Index_t i = 0; i < 4; i++) {
        hxx[i] = hourgam[0][i] * zd1[0] + hourgam[1][i] * zd1[1] +
                hourgam[2][i] * zd1[2] + hourgam[3][i] * zd1[3] +
                hourgam[4][i] * zd1[4] + hourgam[5][i] * zd1[5] +
                hourgam[6][i] * zd1[6] + hourgam[7][i] * zd1[7];
    }
    for(Index_t i = 0; i < 8; i++) {
        hgfz[i] = coefficient *
                    (hourgam[i][0] * hxx[0] + hourgam[i][1] * hxx[1] +
                    hourgam[i][2] * hxx[2] + hourgam[i][3] * hxx[3]);
    }

    p_fx1[0] += hgfx[0];
    p_fy1[0] += hgfy[0];
    p_fz1[0] += hgfz[0];

    p_fx2[0] += hgfx[1];
    p_fy2[0] += hgfy[1];
    p_fz2[0] += hgfz[1];

    p_fx3[0] += hgfx[2];
    p_fy3[0] += hgfy[2];
    p_fz3[0] += hgfz[2];

    p_fx4[0] += hgfx[3];
    p_fy4[0] += hgfy[3];
    p_fz4[0] += hgfz[3];

    p_fx5[0] += hgfx[4];
    p_fy5[0] += hgfy[4];
    p_fz5[0] += hgfz[4];

    p_fx6[0] += hgfx[5];
    p_fy6[0] += hgfy[5];
    p_fz6[0] += hgfz[5];

    p_fx7[0] += hgfx[6];
    p_fy7[0] += hgfy[6];
    p_fz7[0] += hgfz[6];

    p_fx8[0] += hgfx[7];
    p_fy8[0] += hgfy[7];
    p_fz8[0] += hgfz[7];
}