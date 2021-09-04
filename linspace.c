double* linspace(double a, double b, int n, double u[])
{
    double c;
    int i;
    
    /* make sure number of points and array are valid */
    if(n < 2 || u == 0)
        return (void*)0;
    
    /* step size */
    c = (b - a)/(n - 1);
    
    /* fill vector */
    for(i = 0; i < n - 1; ++i)
        u[i] = a + i*c;
    
    /* fix last entry to b */
    u[n - 1] = b;
    
    /* done */
    return u;
}
