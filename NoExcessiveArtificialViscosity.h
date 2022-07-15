inline void NoExcessiveArtificialViscosity(
    double *q
){
    if ( q[0] > m_qstop ) {
        exit(QStopError);
    }
}