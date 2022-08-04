inline void CalcAccelForNodes(
                            double *xdd, double *ydd, double *zdd,
                            double *fx, double *fy, double *fz,
                            double *nodalMass
                            ){
    xdd[0] = fx[0] / nodalMass[0];
    ydd[0] = fy[0] / nodalMass[0];
    zdd[0] = fz[0] / nodalMass[0];
}