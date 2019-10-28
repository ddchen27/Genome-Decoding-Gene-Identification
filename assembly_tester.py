from compsci260lib import *
import re
import sys

def assembly_tester(reads, contigs, supercontigs):
#     contig0 = get_fasta_dict('/Users/David/CS260/Problem Set 7/contig0.fasta').get('contig0')
#     contig1 = get_fasta_dict('/Users/David/CS260/Problem Set 7/contig1.fasta').get('contig1')
#     supercontig0 = get_fasta_dict('/Users/David/CS260/Problem Set 7/supercontig.fasta').get('super_contig0')
#     print "Lengths in nucleotides of contig 0, contig 1, and supercontig:", len(contig0), len(contig1), len(supercontig0) - len(" <---bridged gap---> ")
    
    contig0 = contigs.get('contig0')
    contig1 = contigs.get('contig1')
    read = sorted(reads.keys())
    for i in range(0, len(read) - 1, 2):
        paired_read1 = reads.get(read[i])
        paired_read2 = reads.get(read[i + 1])
        
        # 1. Every read appears somewhere in the set of contigs in the correct orientation
        if(re.search(paired_read1, contig0) == None and re.search(paired_read1, contig1) == None):
            print "First condition failed for paired read:", read[i]
        if(re.search(paired_read2, contig0) == None and re.search(paired_read2, contig1) == None):
            print "First condition failed for paired read:", read[i + 1]
        
        # 2. If mated pair of reads are in the same contig, the distance between is >=1980 and <=2020
        if(re.search(paired_read1, contig0) and re.search(paired_read2, contig0)):
            distance = abs(re.search(paired_read1, contig0).start() - (re.search(paired_read2, contig0).start() + len(paired_read2) - 1))
            if(distance < 1980 or distance > 2020):
                print "Second condition failed for paired reads in contig 0", read[i], read[i + 1]
        if(re.search(paired_read1, contig1) and re.search(paired_read2, contig1)):
            distance = abs(re.search(paired_read1, contig1).start() - (re.search(paired_read2, contig1).start() + len(paired_read2) - 1))
            if(distance < 1980 or distance > 2020):
                print "Second condition failed for paired reads in contig 1", read[i], read[i + 1]
                
        # 3. When a mated pair of reads in different contigs, the first read in the pair appears before the second read in the pair
        pair1 = paired_read1
        pair2 = paired_read2
        if(re.search(paired_read1, supercontigs) and re.search(paired_read2, supercontigs)): #checks which pair is first
            if(re.search(paired_read1, supercontigs).start() > re.search(paired_read2, supercontigs).start()):
                pair1 = paired_read2
                pair2 = paired_read1
        if((re.search(pair1, contig0) and re.search(pair2, contig1)) or (re.search(pair1, contig1) and re.search(pair2, contig0))):
            if(re.search(pair1, supercontigs).start() > re.search(pair2, supercontigs).start()):
                print "Third condition failed for paired reads:", read[i], read[i + 1]
        
if __name__ == '__main__':
    reads = get_fasta_dict('/Users/David/CS260/Problem Set 7/paired.reads.fasta')
    contigs = {}
    contigs['contig0'] = get_fasta_dict('/Users/David/CS260/Problem Set 7/contig0.fasta').get('contig0')
    contigs['contig1'] = get_fasta_dict('/Users/David/CS260/Problem Set 7/contig1.fasta').get('contig1')
    supercontigs = get_fasta_dict('/Users/David/CS260/Problem Set 7/supercontig.fasta').get('super_contig0')
    assembly_tester(reads, contigs, supercontigs)
