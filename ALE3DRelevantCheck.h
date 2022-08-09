inline void ALE3DRelevantCheck(const double *v, const double *eosvmin,const double *eosvmax){
    double vc = v[0] ;
    if ((*eosvmin) != double(0.)) {
        if (vc < (*eosvmin))
        vc = (*eosvmin) ;
    }
    if ((*eosvmax) != double(0.)) {
        if (vc > (*eosvmax))
        vc = (*eosvmax) ;
    }
    if (vc <= 0.) {
        // exit(-1); // Volume Error
    }
}
