inline void updateVolumesForElem(
    double *vnew,
    double *v
){
    Real_t tmpV = vnew[0] ;

    if ( FABS(tmpV - Real_t(1.0)) < m_v_cut )
    tmpV = Real_t(1.0) ;

    // domain.v(i) = tmpV ;
    v[0] = tmpV ;
}