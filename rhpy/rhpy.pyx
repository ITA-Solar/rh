import numpy as np
import os
cimport numpy as np
from libc.stdlib cimport free, malloc
cdef extern from "stdlib.h":
     void *memcpy(void *dst, void *src, long n)

# we are now using Double (Float64)
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

# external routines that we are wrapping
cdef extern double **get_nstar(char *hatom_file, char *pf_file, char *abund_file,
                               int nspace, double *temp, double *nne, double *toth)


# converts from C Float **a to numpy 2D Double array (and frees a)
cdef inline np.ndarray c2npy_double(double **a, int m, int n):
     cdef np.ndarray[DTYPE_t,ndim=2]result = np.zeros((m,n),dtype=DTYPE)
     cdef double *dest
     cdef int i
     dest = <double *> malloc(m*n*sizeof(double*))	
     for i in range(m):
         memcpy(dest + i*n,a[i],n*sizeof(double*))
         #free(a[i]) # caused an error -- why?
     memcpy(result.data,dest,m*n*sizeof(double*))
     free(dest)
     free(a)
     return result     


# define our function
def nh_lte(np.ndarray[DTYPE_t, ndim=1] temp, np.ndarray[DTYPE_t,ndim=1] nelec,
           np.ndarray[DTYPE_t, ndim=1] htot):
    ''' nl_lte(temp,nelec,htot) 
        Ccalculates the 6-level Hydrogen LTE populations using Saha-Boltzmann from RH.
        
        IN:
           temp  : Temperature [K] 
           nelec : electron density [m^-3] 
           htot  : number of hydrogen atoms per volume [m^-3].
        All these arrays have to be Float64 and one-dimensional. 

        OUT:
           nh    : hydrogen populations (2d array: first index level (out of 6), 
                   second index depth. 
    '''

    cdef int N = temp.shape[0]
    cdef int NHydr = 6
    cdef int i
    cdef double *temp_c
    cdef double *nelec_c
    cdef double *htot_c
    cdef double **res

    # check if contiguous, if not force C contiguous arrays
    if not (<object>temp).flags["C_CONTIGUOUS"]:
        temp = temp.copy('C')
    if not (<object>nelec).flags["C_CONTIGUOUS"]:
        nelec = nelec.copy('C')
    if not (<object>htot).flags["C_CONTIGUOUS"]:
        htot = htot.copy('C')

    # point C variables to numpy data
    temp_c  = (<double *>temp.data)
    htot_c  = (<double *>htot.data)    
    nelec_c = (<double *>nelec.data)

    # get paths and files
    base_dir = os.getenv('RH_BASE_PATH')
    if not base_dir:
        raise IOError('Could not get environment variable RH_BASE_PATH. Is it set?')

    base_dir  += 'Atoms/'
    hatom_file = base_dir+'H_6.atom'
    pf_file    = base_dir+'pf_Kurucz.input'
    abund_file = base_dir+'abundance.input'

    res = get_nstar(hatom_file,pf_file,abund_file,N,temp_c,nelec_c,htot_c)
    result = c2npy_double(res,NHydr,N)

    return result
