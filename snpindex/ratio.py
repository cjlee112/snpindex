from rpy2 import robjects

# get some standard R classes
r = robjects.r
fet = r['fisher.test']
iv = robjects.IntVector
rmat = r['matrix']

def ratio_interval(l, maxP=0.0001):
    '''Expects list of tuples of the form
    (Zscore,i,j,Nij,Nj-Nij,Nij0,Nj0-Nij0) e.g. from PairAnalysis.analyze()
    Returns sorted list of tuples of the form
    (pval, isnp1, isnp2, ratioMLE, ratioLowerBound)'''
    output = []
    for t in l:
        result = fet(rmat(iv(t[3:]), 2, 2), alternative='greater')
        p = result.r['p.value'][0][0]
        if p <= maxP:
            output.append((p, t[1], t[2], result.r['estimate'][0][0],
                           result.r['conf.int'][0][0]))
    output.sort()
    return output

