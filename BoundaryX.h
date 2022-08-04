inline void BoundaryX(double *xdd, int *symmX){
    if(symmX[0] & BOUNDARY_NODE){
        xdd[0] = Real_t(0.0);
    }
}