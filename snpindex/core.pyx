        

def encode_codon(ipos,codon,aa,
                 nt='agtc',acode='*ACDEFGHIKLMNPQRSTVWY'):
    'encode SNP as an integer value'
    codon=codon.lower()
    codon0=nt.index(codon[0])+4*nt.index(codon[1]) \
            +16*nt.index(codon[2])
    aa0=acode.index(aa.upper()[0])
    return codon0+64*aa0+1344*ipos


def decode_codon(int snp_code,nt='agtc',acode='*ACDEFGHIKLMNPQRSTVWY'):
    'decode SNP code (not SNP index value) to obtain ipos,aa,codon'
    cdef int ipos,aa,codon
    aa=(snp_code/64)%21
    codon=snp_code%64
    ipos=snp_code/1344
    return ipos,acode[aa],nt[codon%4]+nt[(codon/4)%4]+nt[(codon/16)%4]


cdef int intersect(CDict *d,int n,CDictEntry *dd,CDictEntry *dd2):
    'save items that are in both d an dd'
    cdef int i,j
    j=0
    for i from 0 <= i < n:
        if cdict_getitem(d,dd[i].k):
            dd2[j].k=dd[i].k
            dd2[j].v=dd[i].v
            j=j+1
    return j

cdef int intersect_python(other,int n,CDictEntry *dd,CDictEntry *dd2):
    'save items that are in both other and dd'
    cdef int i,j
    j=0
    for i from 0 <= i < n:
        if dd[i].k in other:
            dd2[j].k=dd[i].k
            dd2[j].v=dd[i].v
            j=j+1
    return j




cdef int union(int n1,CDictEntry *dd1,int n2,CDictEntry *dd2,CDictEntry *dd3):
    'find union of dd1 and dd2 by linear scan'
    cdef int i,j,k
    i=0
    j=0
    k=0
    while i<n1 or j<n2:
        if i<n1 and (j>=n2 or dd1[i].k<dd2[j].k):
            dd3[k].k=dd1[i].k
            dd3[k].v=dd1[i].v
            i=i+1
        elif i<n1 and j<n2 and dd1[i].k==dd2[j].k:
            dd3[k].k=dd1[i].k
            dd3[k].v=dd1[i].v
            i=i+1
            j=j+1
        else:
            dd3[k].k=dd2[j].k
            dd3[k].v=dd2[j].v
            j=j+1
        k=k+1
    return k


cdef int difference(int n,CDictEntry *dd1,object other,CDictEntry *dd3):
    cdef int i,j,is_absent
    cdef IntSet o
    cdef CDict *d
    try:
        o=other
        d=o.d
    except TypeError:
        d=<CDict *>NULL
    j=0
    for i from 0 <= i < n:
        if d: # DIRECT ACCESS TO A CDictEntry ARRAY
            if cdict_getitem(d,dd1[i].k):
                is_absent=0
            else:
                is_absent=1
        else: # ACCESS AS PYTHON DICT
            if dd1[i].k in other:
                is_absent=0
            else:
                is_absent=1
        if is_absent:
            dd3[j].k=dd1[i].k
            dd3[j].v=dd1[i].v
            j=j+1
    return j







cdef CDict *cdict_from_dict(object d):
    'build cdict from any python dict or list'
    cdef int i,k,v
    cdef CDict *d2
    d2=cdict_alloc(len(d))
    if d2==NULL:
        raise MemoryError
    try:
        l=d.items()
    except AttributeError:
        l=[]
        for j in d:
            l.append((j,-1))
    i=0
    for k,v in l:
        try:
            d2[0].dict[i].k=k
            d2[0].dict[i].v=v
        except TypeError:
            raise TypeError('IntSet keys/values must be int!')
        i=i+1
    qsort(d2[0].dict,i,sizeof(CDictEntry),cdict_qsort_cmp)
    return d2





        




cdef class PairAnalysis:
    'performs pairwise analysis of SNPs to find associations'
    cdef SNPDataset dataset
    cdef int nsnp

    def __new__(self,SNPDataset dataset):
        self.dataset=dataset
        self.nsnp=dataset.nsnp
        #self.analyze(icap,ratio_min)

    def analyze(self,float f=2.0,int Nij_min=2):
        '''perform pairwise analysis over all mutations.
        Return tuples (Zscore,i,j,Nij,Nj-Nij,Nij0,Nj0-Nij0) that meet
        the criterion p(i|j)=Nij/Nj more than f std devs above
        p(i|j0)=Nij0/Nj0, using sample std dev given by 1/sqrt(N)
        for sample size N (i.e. Nj for p(i|j), and Nj0 for p(i|j0)).
        The tuples are sorted in descending order of Zscore.
        Arguments:
          f: detection threshold, in #std deviations (default value 2).
          Nij_min: minimum count required for any pair ij to be reported
                   (default value 2).'''
        cdef int i,j,k,l,s,isnp,n
        cdef CDictEntry *dd
        cdef CGraphEntry *gd
        cdef float g,*sigma,sigmaJ0
        cdef float Nij0,Nj0,Pij0
        cdef int *count,*N0
        f=f/2. # TWO STD DEV --> f=1
        gd=self.dataset.snp[0].dict
        count=<int *>calloc(self.nsnp*self.nsnp,sizeof(int))
        if count==NULL:
            raise MemoryError
        sigma=<float *>calloc(self.nsnp,sizeof(float))
        if sigma==NULL:
            raise MemoryError
        N0=<int *>calloc(self.dataset.npos,sizeof(float))
        if N0==NULL:
            raise MemoryError
        i=self.dataset.nsnp - 1
        n=0
        while i>=0: # COUNT #WT OBS FOR EACH CODON
            if gd[i].v: #  A REGULAR MUTATION: SUM COUNTS FOR THIS POS
                n=n+gd[i].v[0].n
            else: # END OF MUTATIONS FOR THIS POS: SAVE TOTAL COUNT
                j=gd[i].k/1344
                N0[j]=self.dataset.nsample - n # #UNMUTATED OBS
                n=0
            i=i-1
        for i from 0 <= i < self.nsnp: # COMPUTE SAMPLE STD DEV
            if gd[i].v:
                sigma[i]=f/c_sqrt(gd[i].v[0].n)
        for s from 0 <= s < self.dataset.max_sample:
            dd=self.dataset.sample[s].dict
            for i from 0 <= i < self.dataset.sample[s].n:
                isnp=dd[i].k
                l=isnp*self.nsnp
                for j from 0 <= j < self.dataset.sample[s].n:
                    if i != j:
                        k=dd[j].k
                        count[l+k] = count[l+k] +1
        # NOW COMPUTE SIGNIFICANCE OF (Nij/Nj) vs. (Nij0/Nj0)
        # Nij: count[i*self.nsnp+j]
        # Nj: gd[j].v[0].n
        # Nij0: gd[i].v[0].n - SUM(count[i*self.nsnp+k] for k in same pos as j)
        # Nj0: N0[gd[j].k/1344]
        pairList=[]
        for i from 0 <= i < self.nsnp: # FIND LINKED i,j PAIRS
            if gd[i].v==NULL: # WT CODON, NOT MUTATION, SO SKIP
                continue
            l=i*self.nsnp
            for j from 0 <= j < self.nsnp:
                if gd[j].v==NULL: # WT CODON: START OF MUTS AT THIS POS
                    Nj0=N0[gd[j].k/1344] # COUNT OF WT CODON OBS AT THIS POSITION
                    k=j+1
                    Nij0=gd[i].v[0].n # #TOTAL OBS OF MUTATION i
                    while gd[k].v: # SUBTRACT OBS OF i WITH MUTATIONS AT j POS
                        Nij0=Nij0-count[l+k]
                        k=k+1
                    sigmaJ0=f/c_sqrt(Nj0)
                    Pij0=Nij0/Nj0 + sigmaJ0 # ADJUST UPWARDS BY SAMPLE DEV
                elif i!=j and count[l+j]>=Nij_min:
                    g=count[l+j]/<float>gd[j].v[0].n - sigma[j] # ADJUST BY SAMPLE DEV
                    if g>Pij0: # OVER THRESHOLD, SAVE AS GOOD CANDIDATE
                        pairList.append(((g-Pij0+sigma[j]+sigmaJ0)
                                         *2.*f/(sigma[j]+sigmaJ0),
                                         i,j,count[l+j],gd[j].v[0].n-count[l+j],
                                         Nij0,Nj0-Nij0))

        free(sigma) # DUMP TEMPORARY ARRAYS
        free(count)
        free(N0)
        pairList.sort() # SORT BY Z SCORE, HIGHEST VALUES FIRST
        pairList.reverse()
        return pairList
                                                                

    def __dealloc__(self):
        pass






cdef class IntSet:
    '''dict interface to a set of SNPs or samples coded as integer IDs.
    In addition to the standard dict methods, it also implements set
    operations:
      A*B (intersection of A and B; *= also implemented)
      A+B (union of A and B; += also implemented)
      A-B (subtraction: subset of A not in B; -= also implemented)'''
    def __new__(self,SNPDataset dataset=None,int isnp= -1,int isample= -1,
                int nalloc=0,iterable=None):
        self.dataset=dataset
        if dataset is not None and isnp>=0:
            if isnp>=dataset.nsnp:
                raise IndexError('overflow')
            self.d=dataset.snp[0].dict[isnp].v
            if self.d==NULL:
                raise KeyError('wt codon has no sample list')
        elif dataset is not None and isample>=0:
            if isample>=dataset.max_sample:
                raise IndexError('overflow')
            self.d=dataset.sample+isample
        elif nalloc>0:
            self.d=cdict_alloc(nalloc)
            if self.d==NULL:
                raise MemoryError
            self.d_free=self.d
        elif iterable is not None:
            self.d=cdict_from_dict(iterable)
            self.d_free=self.d
        else:
            raise ValueError('must pass a valid isnp,isample,nalloc or iterable!')                

    def __len__(self):
        return self.d[0].n
    def keys(self):
        'get list of SNPs in this dataset'
        cdef int i
        l=[]
        for i from 0 <= i < self.d[0].n:
            l.append(self.d[0].dict[i].k)
        return l
    def __iter__(self):
        return iter(self.keys())
    def items(self):
        'get list of all snp:sample-list tuples'
        cdef int i
        l=[]
        for i from 0 <= i < self.nsnp:
            l.append((self.d[0].dict[i].k,self.d[0].dict[i].v))
        return l
    def iteritems(self):
        'generator for all snp:sample-list tuples'
        return iter(self.items())
    def __getitem__(self,int k):
        cdef CDictEntry *e
        e=cdict_getitem(self.d,k)
        if e==NULL:
            raise KeyError
        return e[0].v
    def __contains__(self,k):
        try:
            v=self[k]
            return True
        except KeyError:
            return False

    cdef void resize(self,int newsize):
        cdef int oldsize
        cdef CDictEntry *dd
        oldsize=self.d[0].n
        if newsize==oldsize:
            return # NOTHING TO DO
        self.d[0].n=newsize
        if newsize==0:
            newsize=1 # MINIMUM SIZE FOR REALLOC
        dd=<CDictEntry *>realloc(self.d[0].dict,newsize*sizeof(CDictEntry))
        if dd!=NULL: # REALLOC SUCCEEDED
            self.d[0].dict=dd
        elif oldsize>newsize: # UNABLE TO EXPAND AS NEEDED!
            raise MemoryError
        else: # NOT ABLE TO COMPACT, NOT REALLY A CRITICAL ERROR...
            import logging
            logging.warning('realloc failed to shorten array in snpsample.IntSet.resize?!?!')

    def __mul__(self,object other):
        cdef int j,n
        cdef IntSet o,result,q
        cdef CDict *d
        cdef CDictEntry *dd,*e
        q=self
        try:
            o=other
        except TypeError: # TREAT OTHER AS GENERIC PYTHON DICT
            result=IntSet(nalloc=q.d[0].n)
            d=q.d
            #d=self.d
            #dd=d[0].dict
            #j=testme(dd)
            j=intersect_python(other,q.d[0].n,d[0].dict,result.d[0].dict)
            result.resize(j) # SHORTEN THE LIST IF NEEDED
            return result
        if q.d[0].n<o.d[0].n:
            n=q.d[0].n
            dd=q.d[0].dict
            d=o.d
        else:
            n=o.d[0].n
            dd=o.d[0].dict
            d=q.d
        result=IntSet(nalloc=n)
        j=intersect(d,n,dd,result.d[0].dict)
        result.resize(j) # SHORTEN THE LIST IF NEEDED
        return result
    def __imul__(self,other):
        cdef int j
        cdef IntSet o
        try:
            o=other
            j=intersect(o.d,self.d[0].n,self.d[0].dict,self.d[0].dict)
        except TypeError:
            j=intersect_python(other,self.d[0].n,self.d[0].dict,
                               self.d[0].dict)
        self.resize(j) # SHORTEN THE LIST IF NEEDED
        return self # imul MUST ALWAYS RETURN self!!
    def __add__(self,other):
        cdef int k,n2
        cdef IntSet o,result,q
        cdef CDictEntry *dd2
        cdef CDict *d_free
        q=self
        try:
            o=other
            n2=o.d[0].n
            dd2=o.d[0].dict
            d_free=NULL
        except TypeError:
            d_free=cdict_from_dict(other)
            n2=d_free[0].n
            dd2=d_free[0].dict
        result=IntSet(nalloc=q.d[0].n+o.d[0].n)
        k=union(q.d[0].n,q.d[0].dict,n2,dd2,result.d[0].dict)
        result.resize(k) # SHORTEN THE LIST IF NEEDED
        if d_free: # DUMP TEMP STORAGE IF PRESENT
            cdict_free(d_free)
        return result
    def __iadd__(self,other):
        cdef int k,n2
        cdef CDict *d,*d_free
        cdef IntSet o
        cdef CDictEntry *dd2
        try:
            o=other
            n2=o.d[0].n
            dd2=o.d[0].dict
            d_free=NULL
        except TypeError:
            d_free=cdict_from_dict(other)
            n2=d_free[0].n
            dd2=d_free[0].dict
        d=cdict_alloc(self.d[0].n+o.d[0].n)
        k=union(self.d[0].n,self.d[0].dict,n2,dd2,d[0].dict) # GET UNION
        if self.d_free: # DUMP OLD MEMORY IF NEEDED
            cdict_free(self.d_free)
        self.d=d # ATTACH NEW CDICT TO THIS OBJECT
        self.d_free=d # DYNAMICALLY ALLOCATED, SO FREE LATER ON...
        self.resize(k) # SHORTEN THE LIST IF NEEDED
        if d_free: # DUMP TEMP STORAGE IF PRESENT
            cdict_free(d_free)
        return self # iadd MUST RETURN self!
    def __sub__(self,other):
        cdef int n
        cdef IntSet result,q
        q=self
        result=IntSet(nalloc=q.d[0].n)
        n=difference(q.d[0].n,q.d[0].dict,other,result.d[0].dict)
        result.resize(n)
        return result
    def __isub__(self,other):
        cdef int n
        n=difference(self.d[0].n,self.d[0].dict,other,self.d[0].dict)
        self.resize(n)
        return self # isub MUST RETURN self!
    def ndiff(self,IntSet other):
        '#distinct differences between self and other'
        cdef int i,n
        n=0
        for i from 0 <= i < self.d[0].n:
            if cdict_getitem(other.d,self.d[0].dict[i].k):
                n=n+1
        return self.d[0].n+other.d[0].n - 2*n
        
    def __dealloc__(self):
        if self.d_free:
            cdict_free(self.d_free)
    




cdef class SNPDataset:
    '''dict interface to a set of SNPs observed in a set of samples'''
    property samples:
        'dict interface to individual samples in this dataset'
        def __get__(self):
            return SNPSamples(self)

    def __new__(self,nsnp,nsample,infile):
        'create dataset from file.  See method readSortedSNPs() docstring'
        self.snp=cgraph_alloc(nsnp)
        if self.snp==NULL:
            raise MemoryError
        self.sample=<CDict *>calloc(nsample,sizeof(CDict))
        if self.sample==NULL:
            raise MemoryError
        self.readSortedSNPs(infile,nsample,nsnp)

    cdef void save_snp_dict(self,int isnp,int snp_code,int nsample,int tmp_sample[],
                            int maxsnp):
        cdef int i
        cdef CDictEntry *dd
        if isnp>=maxsnp:
            raise IndexError('overflow!  please increase nsnp')
        self.snp[0].dict[isnp].k=snp_code
        if nsample<=0: # NOTHING TO SAVE
            self.snp[0].dict[isnp].v=NULL # MARK EMPTY LIST WITH NULL POINTER
            return
        self.snp[0].dict[isnp].v=cdict_alloc(nsample)
        self.snp[0].dict[isnp].v[0].n=nsample
        dd=self.snp[0].dict[isnp].v[0].dict
        for i from 0 <= i < nsample:
            dd[i].k=tmp_sample[i]
        qsort(dd,nsample,sizeof(CDictEntry),cdict_qsort_cmp)
        self.nmut=self.nmut+1

    cdef void save_sample_dict(self,int isample,int nsnp,CDictEntry sample2snp[],
                               int maxsample):
        cdef int i
        cdef CDictEntry *dd
        if isample>=maxsample:
            raise IndexError('overflow! please increase nsample')
        self.sample[isample].n=nsnp
        if nsnp<0:
            raise IndexError('negative nsnp??')
        elif nsnp==0: # NOTHING TO DO...
            return
        dd=<CDictEntry *>calloc(nsnp,sizeof(CDictEntry))
        if dd==NULL:
            raise MemoryError
        self.sample[isample].dict=dd
        for i from 0 <= i < nsnp: # SAVE ALL THE snp_code VALUES
            dd[i].k=sample2snp[i].v
        qsort(dd,nsnp,sizeof(CDictEntry),cdict_qsort_cmp)

    cdef int c_is_synonymous(self,int isnp):
        cdef int ipos,aa,aa_wt,snp_code
        if isnp>=self.nsnp or isnp<0:
            raise IndexError('isnp out of range')
        snp_code=self.snp[0].dict[isnp].k
        aa=(snp_code/64)%21
        ipos=snp_code/1344
        aa_wt=(self.codon_wt[ipos]/64)%21
        if aa_wt==aa:
            return 1
        else:
            return 0

    cdef int c_is_wt(self,int isnp):
        cdef int snp_code
        if isnp>=self.nsnp or isnp<0:
            raise IndexError('isnp out of range')
        snp_code=self.snp[0].dict[isnp].k
        if snp_code%1344 == 0:
            return 1
        else:
            return 0

    cdef int c_is_ns(self,int isnp):
        if self.c_is_synonymous(isnp)==0 and self.c_is_wt(isnp)==0:
            return 1
        else:
            return 0

    cdef int c_is_transition(self,int isnp):
        cdef int ipos,codon,codon_wt,snp_code,xor
        if isnp>=self.nsnp or isnp<0:
            raise IndexError('isnp out of range')
        snp_code=self.snp[0].dict[isnp].k
        codon=snp_code%64
        ipos=snp_code/1344
        codon_wt=self.codon_wt[ipos]%64
        xor=codon^codon_wt
        if xor==1 or xor==4 or xor==16:
            return 1
        else:
            return 0

    def is_synonymous(self,int isnp):
        'return 1 if synonymous SNP, else 0'
        return self.c_is_synonymous(isnp)
    def is_wt(self,int isnp):
        'return 1 if WT codon (no mutation), else 0'
        return self.c_is_wt(isnp)
    def is_ns(self,int isnp):
        'return 1 if non-synonymous SNP, else 0'
        return self.c_is_ns(isnp)
    def is_transition(self,int isnp):
        'return 1 if transition SNP, else 0'
        return self.c_is_transition(isnp)
    def snp_info(self,int isnp):
        'get ipos,aa,codon,aa_wt,codon_wt for this SNP'
        cdef int ipos
        if isnp>=self.nsnp or isnp<0:
            raise IndexError('isnp out of range')
        ipos,aa,codon=decode_codon(self.snp[0].dict[isnp].k)
        ipos2,aa_wt,codon_wt=decode_codon(self.codon_wt[ipos])
        if self.c_is_wt(isnp):
            return ipos,aa_wt,codon_wt,aa_wt,codon_wt
        else:
            return ipos,aa,codon,aa_wt,codon_wt

    def similarSamples(self,int maxdiff,int mindiff=0):
        'return a list of sample pairs with mindiff <= d <= maxdiff differences'
        cdef int s1,s2,i,j,k,n1,n2
        cdef CDictEntry *dd1,*dd2
        l=[]
        for s1 from 0 <= s1 < self.max_sample:
            n1=self.sample[s1].n
            if n1>0:
##                 if s1%100==0:
##                     print s1
                dd1=self.sample[s1].dict
                for s2 from s1+1 <= s2 < self.max_sample:
                    n2=self.sample[s2].n
                    if n2>0:
                        dd2=self.sample[s2].dict
                        i=0
                        j=0
                        k=0
                        while i<n1 and j<n2:
                            if dd1[i].k<dd2[j].k:
                                i=i+1
                            elif dd1[i].k==dd2[j].k:
                                i=i+1
                                j=j+1
                                k=k+1
                            else:
                                j=j+1
                        k=n1+n2-2*k
                        if mindiff<=k and k<=maxdiff:
                            l.append((s1,s2))
        return l

        

    def readSortedSNPs(self,infile,int nsample,int maxsnp):
        '''read data, one SNP per line in tab-separated format
        (ipos, codon, aa, wt_codon, wt_aa, isample)
        The input file should be sorted so all observations of a given
        SNP (i.e. a specific codon at a specific ipos) are grouped
        together.  This can be achieved via the UNIX command
        sort <infile >outfile

        Arguments:
        - nsample: should be larger than the largest isample value
        - maxsnp: should be larger than the total # of distinct
          codon variants at all positions (including the wildtype
          codon).  E.g. allowing for 9 possible SNPs at each codon,
          plus the wildtype codon, for 100 codons use a maxsnp=1000'''
        cdef int i,j,snp_code,*tmp_sample,n,ntotal,nsnp,npos,*nmut,nsnp_pos,lastsample
        cdef int *codon_wt,isample
        cdef CDict *d
        cdef CDictEntry *dd
        tmp_sample=<int *>calloc(nsample,sizeof(int))
        nmut=<int *>calloc(maxsnp,sizeof(int))
        codon_wt=<int *>calloc(maxsnp,sizeof(int))
        if tmp_sample==NULL or nmut==NULL or codon_wt==NULL:
            raise MemoryError
        nsnp= 0
        imax=0
        npos=0
        lastpos= -1
        lastcodon='NOBODY'
        n=0
        ntotal=0
        nsnp_pos=0
        lastsample = -1
        for line in infile:
            l=line.split()
            ipos,codon,aa,codon0,aa0,isample= (int(l[0]),l[1],l[2],l[3],l[4],int(l[5]))
            if codon == codon0: # ignore wildtype codons
                continue
            if ipos!=lastpos or codon!=lastcodon:
                lastsample = -1
                if n>0: # SAVE LAST SNP LIST
                    self.save_snp_dict(nsnp,snp_code,n,tmp_sample,maxsnp)
                    nsnp=nsnp+1 # ADVANCE TO NEXT SNP SLOT
                if ipos!=lastpos: # NEW POSITION
                    if ipos>=maxsnp: # CHECK FOR OVERFLOW
                        raise IndexError('overflow! please increase nsnp')
                    if lastpos>=0: # SAVE COUNT OF "MUTATIONS" AT LAST POSITION
                        nmut[lastpos]=nsnp_pos
                    codon_wt[ipos]=encode_codon(0,codon0,aa0) # SAVE NEW WT CODON
                    self.save_snp_dict(nsnp,1344*ipos,0,NULL,maxsnp) # SAVE WT CODON
                    nsnp=nsnp+1 # ADVANCE TO NEXT SNP SLOT
                    nsnp_pos=0 # RESET COUNTERS FOR THIS POSITION
                    lastpos=ipos
                snp_code=encode_codon(ipos,codon,aa)
                if ipos>=npos:
                    npos=ipos+1
                lastcodon=codon
                n=0
            if isample == lastsample:
                continue # ignore duplicate samples
            if n>=nsample:
                raise IndexError('overflow! please increase nsample')
            tmp_sample[n]=isample
            lastsample = isample
            n= n+1
            ntotal = ntotal + 1
            nsnp_pos= nsnp_pos + 1 # COUNT MUTATIONS AT THIS CODON
        if n>0: # SAVE LAST SNP LIST
            self.save_snp_dict(nsnp,snp_code,n,tmp_sample,maxsnp)
            nsnp = nsnp + 1
            if ipos>=maxsnp:
                raise IndexError('overflow! please increase nsnp')
            nmut[ipos]=nsnp_pos
        free(tmp_sample)
        qsort(self.snp[0].dict,nsnp,sizeof(CGraphEntry),cgraph_qsort_cmp)
        self.snp[0].n=nsnp # SAVE TOTAL COUNT OF ALL MUTATION TYPES
        self.nsnp=nsnp
        self.npos=npos
        # WT CODON IS 1ST snp ENTRY FOR EACH POSITION

        d=cdict_alloc(ntotal) # CREATE TEMP STORAGE FOR SORTING BY SAMPLE ID
        if d==NULL:
            raise MemoryError
        n=0
        for i from 0 <= i < nsnp: # SAVE (SAMPLE ID,isnp)
            if self.snp[0].dict[i].v: # REAL MUTATION WITH SAMPLES TO SAVE
                dd=self.snp[0].dict[i].v[0].dict
                for j from 0 <= j < self.snp[0].dict[i].v[0].n:
                    d[0].dict[n].k=dd[j].k
                    d[0].dict[n].v=i
                    n= n+1
        dd=d[0].dict
        if ntotal < n:
            raise IndexError('n > ntotal  Debug!')
        elif ntotal > n:
            print 'Some duplicate snp-sample entries were filtered.'
        qsort(dd,n,sizeof(CDictEntry),cdict_qsort_cmp) # SORT BY SAMPLE ID
        isample = -1
        n=0
        j=0
        for i from 0 <= i < ntotal:
            if dd[i].k!=isample:
                if n>0: # SAVE LAST SAMPLE LIST
                    self.save_sample_dict(isample,n,dd+i-n,nsample)
                isample=dd[i].k # NEW SAMPLE ID
                n=0
                j=j+1 # COUNT #SAMPLES
            n=n+1 # COUNT #SNP IN THIS SAMPLE
        if n>0: # SAVE LAST SAMPLE LIST
            self.save_sample_dict(isample,n,dd+ntotal-n,nsample)
            self.min_sample=dd[0].k # 1ST SAMPLE ID
            self.max_sample=isample+1 # END OF VALID SAMPLES
        cdict_free(d)
        self.nsample=j
        self.p0=<float *>calloc(npos,sizeof(float))
        for i from 0 <= i < npos:
            self.p0[i]=(j-nmut[i])/<float>j
        free(nmut)
        self.codon_wt=<int *>realloc(codon_wt,npos*sizeof(int))
        if self.codon_wt==NULL:
            free(codon_wt)
            raise MemoryError


    def __len__(self):
        '#snps in this dataset'
        return self.nmut
    def keys(self):
        'get list of SNPs in this dataset'
        cdef int i
        l=[]
        for i from 0 <= i < self.nsnp:
            if self.snp[0].dict[i].v:
                l.append(i)
        return l
    def __iter__(self):
        return iter(self.keys())
    def items(self):
        'get list of all snp:sample-list tuples'
        cdef int i,j
        cdef CDictEntry *dd
        l=[]
        for i from 0 <= i < self.nsnp:
            if self.snp[0].dict[i].v:
                l.append((i,IntSet(self,i)))
        return l
    def iteritems(self):
        'generator for all snp:sample-list tuples'
        return iter(self.items())
    def __getitem__(self,int k):
        if k<0 or k>self.nsnp:
            raise KeyError('out of bounds')
        return IntSet(self,k)

    def __dealloc__(self):
        cdef int i
        if self.snp:
            cgraph_free(self.snp)
        if self.sample:
            for i from 0 <= i < self.nsample:
                if self.sample[i].dict!=NULL:
                    free(self.sample[i].dict)
            free(self.sample)
        if self.p0:
            free(self.p0)
        if self.codon_wt:
            free(self.codon_wt)


cdef class SNPSamples:
    def __new__(self,SNPDataset dataset):
        self.dataset=dataset

    def __len__(self):
        return self.dataset.nsample
    def keys(self):
        cdef int i
        l=[]
        for i from 0 <= i < self.dataset.max_sample:
            if self.dataset.sample[i].n>0:
                l.append(i)
        return l
    def __iter__(self):
        return iter(self.keys())
    def items(self):
        cdef int i
        l=[]
        for i from 0 <= i < self.dataset.max_sample:
            if self.dataset.sample[i].n>0:
                l.append((i,IntSet(self.dataset,isample=i)))
        return l
    def iteritems(self):
        return iter(self.items())
    def __getitem__(self,k):
        if k>=0 and k<self.dataset.max_sample and self.dataset.sample[k].n>0:
            return IntSet(self.dataset,isample=k)
        else:
            raise KeyError('invalid isample:%d' %isample)
        
