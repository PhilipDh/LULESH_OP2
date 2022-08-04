inline void CalcLagrangeElemRemaining(
                                    double *dxx, double *dyy, double *dzz,
                                    double *m_vdov,
                                    double *vnew

){
    // calc strain rate and apply as constraint (only done in FB element)
    // double vdov = domain.dxx(k) + domain.dyy(k) + domain.dzz(k) ;
    double vdov = dxx[0] + dyy[0] + dzz[0] ;
    double vdovthird = vdov/double(3.0) ;

    // make the rate of deformation tensor deviatoric
    // domain.vdov(k) = vdov ;
    m_vdov[0] = vdov ;
    // domain.dxx(k) -= vdovthird ;
    // domain.dyy(k) -= vdovthird ;
    // domain.dzz(k) -= vdovthird ;

    dxx[0] -= vdovthird ;
    dyy[0] -= vdovthird ;
    dzz[0] -= vdovthird ;

// See if any volumes are negative, and take appropriate action.
    // if (domain.vnew(k) <= double(0.0))
    if (vnew[0] <= double(0.0))
    {
        exit(VolumeError);
    }
}