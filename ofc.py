import sys
import math
import argparse
from Bio.Seq import reverse_complement, transcribe, back_transcribe, translate
from Bio.Alphabet import generic_dna

parser=argparse.ArgumentParser()
parser.add_argument("sequence", help="enter your sequence", type=argparse.FileType('r'))#file buffer
parser.add_argument("-f", "--frame", help="enter from which frame you want to read", type=int, action="store", default=1) #integer num of frame
parser.add_argument("-o", "--orf", help="use it if you want to read it in Open Reading Frame", action="store_true") #boolean for ORF
parser.add_argument("-r", "--reverse", help="if you want to read it as reverse complement seq from 3' to 5'", action="store_true") #boolean for reverse complement
args=parser.parse_args()

seq_buf=args.sequence.readlines() #sequence string here
num=0
length=0
fr=args.frame-1 #number of frame

seq_buf.pop(0)
for i in range(0, len(seq_buf)):
	seq_buf[i]=seq_buf[i].strip('\n') #omiting possible spaces 
seq_buf=''.join(seq_buf)

if args.reverse==False: #printing from 5' to 3'
    coding_dna1=seq_buf #main seq string
    coding_dna2=coding_dna1[fr:] #in case we want other frame
    amin_seq=translate(coding_dna2, to_stop=args.orf) #2dim vector of translated aminoacids divided with '*' when STOP codon
    print "Sequence from requested frame(",(fr+1) ,"): ", coding_dna2
    print "Complementary sequence:", reverse_complement(coding_dna2)[::-1]
    print "Aminoacids sequence:", amin_seq
    amino=amin_seq.split('*')
        
        
else: #printing from 3' to 5'
    coding_dna1=reverse_complement(seq_buf) #main seq string, reverse complement of course
    coding_dna2=coding_dna1[fr:] #in case we want other frame
    amin_seq=translate(coding_dna2, to_stop=args.orf) #2dim vector of translated aminoacids divided with '*' when STOP codon
    print "Reverse complement sequence from requested frame(",(fr+1) ,"): ", coding_dna2
    print "Complementary sequence:", reverse_complement(coding_dna2)
    print "Aminoacids sequence:", amin_seq
    amino=amin_seq.split('*')

if args.orf==False: #when we want the whole seq
    for j in range(len(amino)):
            if amino[j]=="":
                num=num+3
                continue
            length=length+len(amino[j])*3
            print amino[j], "from position", length+num-len(amino[j])*3 + 1,"to position", length+num
            num=num+3
