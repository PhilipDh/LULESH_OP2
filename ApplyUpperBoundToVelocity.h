inline void ApplyUpperBoundToVelocity(double *vnewc, double *eosvmax){
    if (vnewc[0] > (*eosvmax))
        vnewc[0] = (*eosvmax) ;
}