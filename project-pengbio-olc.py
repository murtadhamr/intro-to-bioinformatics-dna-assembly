# -*- coding: utf-8 -*-
"""
Created on Sun Dec 17 16:52:55 2017

@author: Lenovo
"""

# from Bio import SeqIO
# record_dict = SeqIO.parse("pengbio.fasta_2", "fasta")
# print(record_dict.id)

dict_new = {}
dict_n = []
index = []

sample = {
    '1': 'GGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTCGTCCAGACCCCTAGC',
     '3': 'GTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGTCGTGAAC',
     '2': 'CTTTACCCGGAAGAGCGGGACGCTGCCCTGCGCGATTCCAGGCTCCCCACGGG',
     '5': 'CGATTCCAGGCTCCCCACGGGGTACCCATAACTTGACAGTAGATCTC',
     '4': 'TGCGAGGGAAGTGAAGTATTTGACCCTTTACCCGGAAGAGCG',
     '6': 'TGACAGTAGATCTCGTCCAGACCCCTAGCTGGTACGTCTTCAGTAGAAAATTGTTTTTTTCTTCCAAGAGGTCGGAGT'
}

#david fasta 2 memuat 5 sequence sample
from Bio import SeqIO
for record in SeqIO.parse("pengbio.fasta_2", "fasta"):
    index = record.id.split('.')
    dict_new.update({index[1]: str(record.seq)})
#     dict_n.append({index[1]:str(record.seq)})
    

def meanLength(data): 
    s = 0
    tot = 0
    for x in data:
        s += len(data[x])
        tot += 1
    return s/float(tot)

def getOverlap(left, right):
    for i in range(len(left)):
        if left[i:] == right[:len(left)-i]:
            return left[i:]
    return ''

def getAllOverlaps(reads):
    d = dict()
    for name1, seq1 in reads.items():
        for name2, seq2 in reads.items():
            if name1 == name2:
                continue
            if name1 not in d:
                d[name1] = dict()
            d[name1][name2] = len(getOverlap(seq1, seq2))

    return d

def findFirstRead(overlaps):
    for i in overlaps:
        print(i)
        signifOverlaps = False
        for j in d[i]:
            if d[j][i] > 3:
                signifOverlaps = True
        if not signifOverlaps:
            return i
        
def findKeyForLargestValue(d):
#     return sorted(d.items(), key=lambda x: x[1])[-1][0]
    m = max(d.values())
    for k in d:
        if d[k] == m:
            return k
        
def findOrder(first, d):
    if max(d[first].values()) < 3:
        return [first]
    else:
        nextRead = findKeyForLargestValue(d[first])
        return [first] + findOrder(nextRead, d)
    
def assembleGenome(readOrder, reads, overlaps):
    genome = ''
    for readName in readOrder[:-1]:
        rightOverlap = max(x for x in overlaps[readName].values() if x >= 3)
        genome += reads[readName][:-rightOverlap]
    genome += reads[readOrder[-1]]
    return genome

###############################################################################

d=getAllOverlaps(dict_new)

findFirstRead(d)

order=findOrder(findFirstRead(d), d)

reference = assembleGenome(order,dict_new,d)
print (reference)
print ("Reference gen terpanjang",len(reference))
