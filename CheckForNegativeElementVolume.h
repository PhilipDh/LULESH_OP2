inline void CheckForNegativeElementVolume(double *determ){
    if(determ[0] <= 0.0){
        exit(VolumeError);
    }
}