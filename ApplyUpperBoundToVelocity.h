inline void ApplyUpperBoundToVelocity(double *vnewc, const double *eosvmax){
    if (vnewc[0] > (*eosvmax))
        vnewc[0] = (*eosvmax) ;
}