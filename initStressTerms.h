inline void initStressTerms(double *sigxx, double *p, double *q)
{
    sigxx[0] = -p[0] - q[0];
    sigxx[1] = -p[0] - q[0];
    sigxx[2] = -p[0] - q[0];
}