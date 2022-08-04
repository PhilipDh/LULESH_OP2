inline void BoundaryZ(double *zdd, int *symmZ){
    if(symmZ[0] & BOUNDARY_NODE){
        zdd[0] = Real_t(0.0);
    }
}