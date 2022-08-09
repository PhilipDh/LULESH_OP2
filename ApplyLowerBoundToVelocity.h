inline void ApplyLowerBoundToVelocity(double *vnewc, const double *eosvmin){
    if (vnewc[0] < (*eosvmin))
        vnewc[0] = (*eosvmin) ;
}