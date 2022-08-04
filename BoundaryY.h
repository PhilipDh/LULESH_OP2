inline void BoundaryY(double *ydd, int *symmY){
    if(symmY[0] & BOUNDARY_NODE){
        ydd[0] = Real_t(0.0);
    }
}
