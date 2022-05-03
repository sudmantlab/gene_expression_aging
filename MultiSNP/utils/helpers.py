from bisect import bisect_right, bisect_left
from hashlib import sha256
import numpy as np
 
def BinarySearchLeft(a, x):
    i = bisect_left(a, x)
    if i != len(a):
        return i
    else:
        return -1

def BinarySearchRight(a, x):
    i = bisect_right(a, x)
    if i != 0:
        return (i-1)
    else:
        return -1

class Seed:
    # seed used for each model/data shuffle 
    # for each new seeded process i, use seed str(metaseed)+str(i)
    # Generate replicatable random seed for each iteration of data
    def __init__(self, tissue, metaseed=423671):
        # initialize seed with tissue type and metaseed
        self.metaseed = metaseed
        self.tissue = tissue

    def get(self, gene=0, boot_iter=0, t="model", old=True):
        data = [self.metaseed, self.tissue, str(gene), "boot"+str(boot_iter),
                t, old]
        h = sha256(np.array(data))
        seed = np.frombuffer(h.digest(), dtype='uint32')
        return seed

class GTArray:
    # Genotype data structure
    # Takes numpy arrays of genotype dosage, chromosome, positions, 
    # reference and alt alleles
    def __init__(self, gt, indiv, chrom, pos):
        self.gt = gt 
        self.indiv = indiv
        self.pos = pos 
        self.chrom = chrom
        #self.ref = ref
        #self.alt = alt
        self.chrom_start = {c:BinarySearchLeft(self.chrom, c) for c in set(self.chrom)}
        self.chroms = list(set(self.chrom))

    # Index range of genetic positions
    # returns genotype of positions as pandas dataframe
    #def range(self, chrom, start, end):
        
