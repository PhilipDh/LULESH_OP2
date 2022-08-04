inline void ALE3DRelevantCheck(double *v, double *eosvmin, double *eosvmax){
    Real_t vc = v[0] ;
    if ((*eosvmin) != Real_t(0.)) {
        if (vc < (*eosvmin))
        vc = (*eosvmin) ;
    }
    if ((*eosvmax) != Real_t(0.)) {
        if (vc > (*eosvmax))
        vc = (*eosvmax) ;
    }
    if (vc <= 0.) {
        exit(VolumeError);
    }
}
