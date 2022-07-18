inline void CalcMonotonicQRegionForElem(
    double *delv_xi, double *delv_xi_lxim, double *delv_xi_lxip,
    double *delv_eta, double *delv_eta_letam, double *delv_eta_letap,
    double *delv_zeta, double *delv_zeta_lzetam, double *delv_zeta_lzetap,
    double *delx_xi, double *delx_eta, double *delx_zeta,
    double *elemBC,
    double *m_vdov,
    double *qq, double *ql,
    double *elemMass, double *volo, double *vnew,
    double *ptiny,
    double *monoq_limiter_mult,
    double *monoq_max_slope
){
    Real_t qlin, qquad ;
    Real_t phixi, phieta, phizeta ;
    Real_t delvm = 0.0, delvp =0.0;
    Int_t bcMask = elemBC[0] ;

    Real_t norm = Real_t(1.) / (delv_xi[0]+ (*ptiny) ) ;

    switch (bcMask & XI_M) {
    case XI_M_COMM: /* needs comm data */
    case 0:         delvm = delv_xi_lxim[0]; break ;
    case XI_M_SYMM: delvm = delv_xi[0] ;       break ;
    case XI_M_FREE: delvm = Real_t(0.0) ;      break ;
    default:          fprintf(stderr, "Error in switch at %s line %d\n",
                            __FILE__, __LINE__);
        delvm = 0; /* ERROR - but quiets the compiler */
        break;
    }
    switch (bcMask & XI_P) {
    case XI_P_COMM: /* needs comm data */
    case 0:         delvp = delv_xi_lxip[0] ; break ;
    case XI_P_SYMM: delvp = delv_xi[0] ;       break ;
    case XI_P_FREE: delvp = Real_t(0.0) ;      break ;
    default:          fprintf(stderr, "Error in switch at %s line %d\n",
                            __FILE__, __LINE__);
        delvp = 0; /* ERROR - but quiets the compiler */
        break;
    }

    delvm = delvm * norm ;
    delvp = delvp * norm ;

    phixi = Real_t(.5) * ( delvm + delvp ) ;

    delvm *= (*monoq_limiter_mult) ;
    delvp *= (*monoq_limiter_mult) ;

    if ( delvm < phixi ) phixi = delvm ;
    if ( delvp < phixi ) phixi = delvp ;
    if ( phixi < Real_t(0.)) phixi = Real_t(0.) ;
    if ( phixi > (*monoq_max_slope)) phixi = (*monoq_max_slope);

    /*  phieta     */
    norm = Real_t(1.) / ( delv_eta[0] + (*ptiny) ) ;

    switch (bcMask & ETA_M) {
        case ETA_M_COMM: /* needs comm data */
        case 0:          delvm = delv_eta_letam[0] ; break ;
        case ETA_M_SYMM: delvm = delv_eta[0] ;        break ;
        case ETA_M_FREE: delvm = Real_t(0.0) ;        break ;
        default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                __FILE__, __LINE__);
        delvm = 0; /* ERROR - but quiets the compiler */
        break;
    }
    switch (bcMask & ETA_P) {
        case ETA_P_COMM: /* needs comm data */
        case 0:          delvp = delv_eta_letap[0] ; break ;
        case ETA_P_SYMM: delvp = delv_eta[0] ;        break ;
        case ETA_P_FREE: delvp = Real_t(0.0) ;        break ;
        default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                __FILE__, __LINE__);
        delvp = 0; /* ERROR - but quiets the compiler */
        break;
    }

    delvm = delvm * norm ;
    delvp = delvp * norm ;

    phieta = Real_t(.5) * ( delvm + delvp ) ;

    delvm *= (*monoq_limiter_mult) ;
    delvp *= (*monoq_limiter_mult) ;

    if ( delvm  < phieta ) phieta = delvm ;
    if ( delvp  < phieta ) phieta = delvp ;
    if ( phieta < Real_t(0.)) phieta = Real_t(0.) ;
    if ( phieta > (*monoq_max_slope))  phieta = (*monoq_max_slope);

    /*  phizeta     */
    // norm = Real_t(1.) / ( domain.delv_zeta(ielem) + ptiny ) ;
    norm = Real_t(1.) / ( delv_zeta[0] + (*ptiny) ) ;

    switch (bcMask & ZETA_M) {
        case ZETA_M_COMM: /* needs comm data */
        // case 0:           delvm = domain.delv_zeta(domain.lzetam(ielem)) ; break ;
        case 0:           delvm = delv_zeta_lzetam[0] ; break ;
        case ZETA_M_SYMM: delvm = delv_zeta[0] ;         break ;
        case ZETA_M_FREE: delvm = Real_t(0.0) ;          break ;
        default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                __FILE__, __LINE__);
        delvm = 0; /* ERROR - but quiets the compiler */
        break;
    }
    switch (bcMask & ZETA_P) {
        case ZETA_P_COMM: /* needs comm data */
        // case 0:           delvp = domain.delv_zeta(domain.lzetap(ielem)) ; break ;
        case 0:           delvp = delv_zeta_lzetap[0] ; break ;
        case ZETA_P_SYMM: delvp = delv_zeta[0] ;         break ;
        case ZETA_P_FREE: delvp = Real_t(0.0) ;          break ;
        default:          fprintf(stderr, "Error in switch at %s line %d\n",
                                __FILE__, __LINE__);
        delvp = 0; /* ERROR - but quiets the compiler */
        break;
    }

    delvm = delvm * norm ;
    delvp = delvp * norm ;

    phizeta = Real_t(.5) * ( delvm + delvp ) ;

    delvm *= (*monoq_limiter_mult) ;
    delvp *= (*monoq_limiter_mult) ;

    if ( delvm   < phizeta ) phizeta = delvm ;
    if ( delvp   < phizeta ) phizeta = delvp ;
    if ( phizeta < Real_t(0.)) phizeta = Real_t(0.);
    if ( phizeta > (*monoq_max_slope)  ) phizeta = (*monoq_max_slope);

    /* Remove length scale */
    if ( m_vdov[0] > Real_t(0.) )  {
        qlin  = Real_t(0.) ;
        qquad = Real_t(0.) ;
    }
    else {
        Real_t delvxxi   = delv_xi[0]   * delx_xi[0]   ;
        Real_t delvxeta  = delv_eta[0]  * delx_eta[0]  ;
        Real_t delvxzeta = delv_zeta[0] * delx_zeta[0] ;

        if ( delvxxi   > Real_t(0.) ) delvxxi   = Real_t(0.) ;
        if ( delvxeta  > Real_t(0.) ) delvxeta  = Real_t(0.) ;
        if ( delvxzeta > Real_t(0.) ) delvxzeta = Real_t(0.) ;

        Real_t rho = elemMass[0] / (volo[0] * vnew[0]) ;

        qlin = -m_qlc_monoq * rho *
        (  delvxxi   * (Real_t(1.) - phixi) +
            delvxeta  * (Real_t(1.) - phieta) +
            delvxzeta * (Real_t(1.) - phizeta)  ) ;

        qquad = m_qqc_monoq * rho *
        (  delvxxi*delvxxi     * (Real_t(1.) - phixi*phixi) +
            delvxeta*delvxeta   * (Real_t(1.) - phieta*phieta) +
            delvxzeta*delvxzeta * (Real_t(1.) - phizeta*phizeta)  ) ;
    }

    qq[0] = qquad ;
    ql[0] = qlin  ;

}