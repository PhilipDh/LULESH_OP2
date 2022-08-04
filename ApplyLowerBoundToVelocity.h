inline void ApplyLowerBoundToVelocity(double *vnewc, double *eosvmin){
    if (vnewc[0] < (*eosvmin))
        vnewc[0] = (*eosvmin) ;
}