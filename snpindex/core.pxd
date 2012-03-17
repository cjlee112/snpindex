


cdef extern from "string.h":
    ctypedef int size_t
    void *memset(void *,int,size_t)

cdef extern from "stdlib.h":
    void free(void *)
    void *calloc(size_t,size_t)
    void *realloc(void *,size_t)
    void qsort(void *base, size_t nmemb, size_t size,
             int (*compar)(void *,void *))



cdef extern from "math.h":
    double c_log "log" (double)
    double c_exp "exp" (double)
    double c_sqrt "sqrt" (double)


cdef extern from "cgraph.h":
    ctypedef struct CDictEntry:
        int k
        int v

    ctypedef struct CDict:
        int n
        CDictEntry *dict

    CDict *cdict_alloc(int n)
    int cdict_free(CDict *d)
    int cdict_qsort_cmp(void *void_a,void *void_b)
    CDictEntry *cdict_getitem(CDict *d,int k)

    ctypedef struct CGraphEntry:
        int k
        CDict *v

    ctypedef struct CGraph:
        int n
        CGraphEntry *dict

    CGraph *cgraph_alloc(int n)
    int cgraph_free(CGraph *d)
    int cgraph_qsort_cmp(void *void_a,void *void_b)
    CGraphEntry *cgraph_getitem(CGraph *d,int k)
    int *calloc_int(int n)


cdef class SNPDataset:
    cdef CGraph *snp
    cdef CDict *sample
    cdef readonly int nsnp,nsample,min_sample,max_sample,npos,nmut
    cdef float *p0 # PROB OF NO MUTATION AT EACH CODON
    cdef int *codon_wt # ENCODING OF WT CODONS
    cdef float *codon_log_p # LOG PROB OF ALL POSSIBLE CODONS
    cdef void save_snp_dict(self,int isnp,int snp_code,int nsample,int tmp_sample[],
                            int maxsnp)
    cdef void save_sample_dict(self,int isample,int nsnp,CDictEntry sample2snp[],
                               int maxsample)
    cdef int c_is_synonymous(self,int isnp)
    cdef int c_is_wt(self,int isnp)
    cdef int c_is_ns(self,int isnp)
    cdef int c_is_transition(self,int isnp)


cdef class SNPSamples:
    cdef SNPDataset dataset


cdef class IntSet:
    cdef CDict *d,*d_free
    cdef SNPDataset dataset
    cdef void resize(self,int newsize)
    #cdef int intersect(self,CDict *d,int n,CDictEntry *dd,CDictEntry *dd2)
    #cdef int intersect_python(self,other,int n,CDictEntry *dd,CDictEntry *dd2)
    #cdef int union(self,CDictEntry *dd2,CDictEntry *dd3)
    #cdef int difference(self,other,CDictEntry *dd3)


