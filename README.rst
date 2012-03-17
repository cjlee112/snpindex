#######################################
Fast SNP Indexing and Pairwise Analysis
#######################################

A python extension module for analyzing statistical 
associations between SNP observations from a set of samples.

Functionality Overview
----------------------

* store and index SNP observations from a large dataset of samples,
  using compact integer representation for sample and SNP IDs,
  and fast C code for indexing. It loads and indexes the full
  50000 sample protease + RT dataset in a few seconds,
  using about 30 MB for data storage (on my laptop).
* convenience functions for getting information about individual
  SNPs.
* dict-like interface to individual SNPs or samples.
* sets of samples (i.e. observations of a specific SNP)
  and sets of SNPs (observed in an individual sample) are
  both stored using the IntSet class, providing a dict-like
  interface plus set operations (intersection A*B, union A+B,
  difference A-B).
* the PairAnalysis class performs fast identification of
  pairs with significant statistical association. It uses
  a fast approximation for finding significant i,j pairs, and
  returns the counts needed for futher statistical tests like
  the Fisher Exact Test. It performs complete pairwise analysis
  of the full 50000 sample protease + RT dataset in about 1
  second, using about 40 MB for temporary data storage.

Installation
------------

You need Pyrex to compile this module.  Build it as follows::

  python setup.py build

Install it as follows::

  python setup.py install

Tutorial
--------

Currently, the module can load data from a tab-separated text
file of the following format::

  101	acc	T	ccc	P	1116
  101	cac	H	ccc	P	24931
  101	cac	H	ccc	P	25304
  101	cac	H	ccc	P	33851
  101	cca	P	ccc	P	10673
  101	cca	P	ccc	P	1073
  101	cca	P	ccc	P	10814
  101	cca	P	ccc	P	119
  101	cca	P	ccc	P	12268
  101	cca	P	ccc	P	1237

The tab-separated fields must be, in order:
ipos, codon, aa, codon_wt, aa_wt, sampleID,  

This file should be sorted in order of (ipos,codon) 
so that observations of a given SNP are grouped together. 
This can be achieved using the UNIX command::

  sort <infile >outfile

Once you have the data file formatted, loading it and 
using snpindex is incredibly easy::

  import snpindex
  ifile = open('hiv_snp_sample_obs.txt')
  s = snpindex.SNPDataset(4000, 51000, ifile)
  ifile.close()

The SNPDataset constructor takes three arguments:

1. the maximum bound on number of SNPs in the dataset.
   An easy rule of thumb is to multiply the number of codons
   in your dataset by 10 (i.e. 9 possible SNPs per codon,
   plus the wildtype codon, which is also indexed by SNPDataset).
   If you give too low a bound, it will raise an IndexError
   exception on overflow.
2. the maximum bound on sampleID. If you give too low a bound,
   it will raise an IndexError exception on overflow.
3. file object to read the dataset from.

A SNPDataset object acts like a dictionary whose keys are SNP 
IDs and whose values are IntSet objects containing the list of 
samples for each SNP. These IntSet objects in turn act like 
dictionaries whose keys are sample IDs, and which additionally 
provide set intersection, union and difference operations A*B, 
A+B, A-B. The SNPDataset object also has a samples attribute, 
which acts like a dictionary whose keys are sampleIDs, and 
whose values are IntSet objects containing a list of SNP IDs 
for each sample.

Examples::

  >>> samples=s.samples
  >>> seq1=samples[2000] # GET SAMPLE 2000
  >>> seq1.keys() # LIST ITS SNPs
  [62, 70, 73, 91, 211, 217, 474, 525, 692, 881, 908, 941, 1002, 1161, 1250, 1295, 1346, 1442, 1526, 1577, 1614
  >>> mut62=s[62] # GET SET OF SAMPLES CONTAINING SNP 62
  >>> len(mut62) # NUMBER OF SAMPLES CONTAINING THIS SNP
  2334
  >>> mut62.keys()[:10] # SHOW THE FIRST TEN
  [3, 18, 27, 34, 62, 93, 94, 122, 144, 162]
  >>> seq2=samples[mut62.keys()[0] ] # GET SAMPLE 3
  >>> len(seq2) # NUMBER OF SNPs IN THIS SAMPLE
  54
  >>> join=seq1*seq2 # INTERSECTION OF THE TWO SETS
  >>> len(join) # NUMBER OF SNPs SHARED IN THE TWO SAMPLES
  19
  >>> join.keys() # LIST THE SHARED SNPs
  [62, 70, 91, 211, 881, 908, 1002, 1295, 1526, 1921, 1990, 2068, 2411, 2461, 2497, 2546, 2615, 2648, 2727]
  >>> seq2.keys()
  [28, 54, 62, 70, 91, 211, 270, 332, 465, 483, 587, 601, 674, 700, 866, 877, 881, 890, 901, 908, 1002, 1035,]
  >>> sum=seq1+seq2 # UNION OF THE TWO SETS
  >>> len(sum) 78 >>> len(join)
  19
  >>> len(seq1)+len(seq2)-len(join)
  78
  >>> diff=seq1-seq2 # SET OF SNPS IN seq1 BUT NOT seq2
  >>> len(diff)
  24
  >>> len(seq1)
  43

SNP information methods
-----------------------

* SNPDataset.is_ns(isnp): return 1 if non-synonymous SNP, else 0
* SNPDataset.is_synonymous(isnp): return 1 if synonymous SNP, else 0
* SNPDataset.is_transition(isnp): return 1 if transition SNP, else 0
* SNPDataset.is_wt(isnp): return 1 if WT codon (no mutation), else 0
* SNPDataset.snp_info(isnp): get ipos,aa,codon,aa_wt,codon_wt for this SNP

PairAnalysis
------------

This class provides fast identification of SNP pairs that show 
statistically significant association.
To avoid the very slow R fisher_test() computation 
(which takes about 15 ms per pair on my Powerbook; this calculation 
alone would take over 38 hours for a single pass through the 
HIV protease + RT dataset), I am using a fast approximation 
that is about 100,000-fold faster:
the Fisher test computes P(Theta1/Theta2>1 | n1, N1, n2,N2), 
where Theta1 and Theta2 are the (hidden) binomial probabilities 
for dataset1 and dataset2 (whose maximum-likelihood estimates 
are n1/N1 and n2/N2 respectively). For a given n1,N1,n2,N2 
dataset, the posterior distributions for Theta1 and Theta2 are 
themselves binomial distributions (which become tighter for 
increasing N1 or N2).
For N>10 these binomial distributions become well-approximated 
by a normal distribution. Thus, computing a confidence interval 
on Theta1 (or Theta2) is equivalent to determining the error 
distribution of the sample mean (of a sample of size N) from a 
normal distribution.
For a normal distribution, the distribution of the means of 
random samples of size N is itself normally distributed around 
the true (hidden) mean, with a standard deviation proportional 
to 1/sqrt(N). Thus, we can approximate the confidence interval 
for Theta1 as [sample_mean +/- 1/sqrt(N1)]. To test this 
approximation, you can use the R epitools package that 
provides exact computation of the confidence interval for binomials::

  >>> from rpy import r
  >>> r.library('epitools')
  ['epitools', 'methods', 'stats', 'graphics', 'grDevices', 'utils', 'datasets', 'base']
  >>> r.binom_exact(50,100,0.95) # COMPUTE 95% CONF INTERVAL FOR n=50, N=100
  {'conf.level': 0.94999999999999996, 'upper': 0.60167887049669888,
  'lower': 0.39832112950330106, 'proportion': 0.5, 'n': 100, 'x': 50}
  >>> import math
  >>> r.binom_exact(5,10,0.95)['upper']-.5
  0.31291397155260148
  >>> 1/math.sqrt(10)
  0.31622776601683794
  >>> r.binom_exact(50,100,0.95)['upper']-.5
  0.10167887049669888
  >>> 1/math.sqrt(100)
  0.10000000000000001
  >>> r.binom_exact(500,1000,0.95)['upper']-.5
  0.031450827028208561
  >>> 1/math.sqrt(1000)
  0.031622776601683791
  >>> r.binom_exact(5000,10000,0.95)['upper']-.5
  0.0098486194100191327
  >>> 1/math.sqrt(10000)
  0.01


* To evaluate the significance of association
  between a pair of SNPs i,j, 
  we can compute a Z score: 
  Z=(P(i | j)-P(i | j0))*2/(sigma(j)+sigma(j0)), 
  where sigma(j)=1/sqrt(Nj) and sigma(j0)=1/sqrt(Nj0). 
  This is just standard Z score with the standard deviation 
  in the denominator computed as the average of the two 
  separate sample standard deviations.

PairAnalysis implements this as a fast C computaion. 
Specifically, it returns a list of SNPs i,j such that 
P(i | j) > P(i | j0) + f*sigma(j)+sigma(j0), 
where P(i | j)=Nij/Nj, P(i | j0)=Nij0/Nj0, 
and sigma(j) is the estimated standard deviation of 
the sample mean for a sample of size Nj drawn from a 
normal distribution. It returns a list of tuples of 
the form (Zscore, i, j, Nij, Nj - Nij, Nij0, Nj0 - Nij0), 
sorted in descending order of Zscore::

  >>> pa = snpindex.PairAnalysis(s)
  >>> l = pa.analyze()
  >>> len(l)
  33353
  >>> l[:10]
  [(103.13453484160702, 688, 670, 30580, 7040, 987.0, 11243.0), (94.755925779817673, 670, 688, 30580, 988, 5694.0, 11184.0), (76.63764376031807, 2378, 2380, 14515, 4206, 7787.0, 21081.0), (75.205579500615912, 211, 670, 24696, 12924, 1499.0, 10731.0), (74.616678483068739, 2380, 2378, 14515, 8498, 4203.0, 22590.0), (65.296536822600387, 670, 211, 24696, 1500, 11638.0, 10644.0), (65.141497670227167, 1063, 2042, 2598, 709, 2395.0, 38619.0), (64.92690574498549, 402, 670, 21308, 16312, 1290.0, 10940.0), (63.952537074124649, 1037, 1063, 4148, 2283, 3754.0, 38033.0), (63.469163331578784, 163, 579, 1769, 272, 735.0, 45088.0)]


If you have the rpy2 python interface to R installed, you can run the
convenience module snpindex.ratio to compute Fisher Exact Test
p-values and ratio estimators::

  >>> from snpindex import ratio
  >>> ratios = ratio.ratio_interval(l[600:650])
  >>> ratios[:10]
  [(0.0, 23, 336, 11.076943849005225, 10.383011675171922),
  (0.0, 23, 540, 9.4580710019886869, 8.8927194575171882),
  (0.0, 23, 597, 8.2182451817711399, 7.8163332825344316),
  (0.0, 163, 579, 397.01451579416636, 348.12879286664923),
  (0.0, 211, 200, 6.9639008504851594, 6.528372829088962),
  (0.0, 211, 402, 2.8830836139892249, 2.7931755682356463),
  (0.0, 211, 670, 13.678218225229886, 13.020703551029637),
  (0.0, 211, 688, 5.1962107599476663, 5.0188359353635068),
  (0.0, 336, 540, 54.538614686114933, 50.536948684820928),
  (0.0, 394, 200, 4.5278731838454567, 4.3139943490009092)]

It returns a list of tuples of the form
(pvalue, isnp1, isnp2, ratioMLE, ratioLowerBound)
i.e. the maximum likelihood estimator for the ratio and
the 95% confidence lower bound.

