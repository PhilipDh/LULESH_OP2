inline void BoundaryX(double *xdd, int *symmX){
    if(symmX[0] & BOUNDARY_NODE){
        xdd[0] = Real_t(0.0);
    }
}

inline void BoundaryY(double *ydd, int *symmY){
    if(symmY[0] & BOUNDARY_NODE){
        ydd[0] = Real_t(0.0);
    }
}

inline void BoundaryZ(double *zdd, int *symmZ){
    if(symmZ[0] & BOUNDARY_NODE){
        zdd[0] = Real_t(0.0);
    }
}