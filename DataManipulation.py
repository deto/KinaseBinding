# -*- coding: utf-8 -*-
"""
Created on Sat Jul 06 16:51:20 2013

@author: daved_000
"""
import numpy as np;
import re;

aa2i_dict = dict([['A',0], ['C',1], ['D',2], ['E',3], ['F',4], 
               ['G',5], ['H',6], ['I',7], ['K',8], ['L',9], 
                ['M',10], ['N',11], ['P',12], ['Q',13], ['R',14], 
                ['S',15], ['T',16], ['V',17], ['W',18], ['Y',19]]);

i2aa_dict = dict((v,k) for k,v in aa2i_dict.iteritems());                

def aa2i(aa_letter):
    return aa2i_dict[aa_letter];
    

"""
Form A:

For a sequence of N residues, create a sequence that is 20xN long composed
of all zeros and ones.  Location M is equal to one if position A in the
input sequence has amino acid X.  M = 15*A + X   
   
   
"""


#Take a list of sequences and process each sequence into form A
#For K sequences, produce a matrix that is 300xK
def list_to_formA(indata):
    outdata = np.zeros(shape = (len(indata),300), dtype=np.uint32);
    for i,row in enumerate(indata):
        outdata[i,:] = seq_to_formA(row);
    
    return outdata;

def seq_to_formA(inseq):
    outseq = np.zeros(shape = (1,300), dtype=np.uint32);    
    for i,res in enumerate(inseq.upper()):
        if(res != '_' and res != 'X'):
            try:
                pos = aa2i_dict[res] + i*20;
                outseq[0,pos] = 1;
            except KeyError:
                pass;
    
    return outseq;

#Output the 20x15 count matrix
def formA_to_count(formA):
    m1 = np.sum(formA, axis=0);
    m2 = np.reshape(m1, newshape=(20,15), order='F');
    return m2;

#uses a list of the 15 char string representation
#Ouputs the 20x15 count matrix
#This naive implementation doesn't take advantage of FormA.  Much slower!
def chunks_to_count(chunks):
    out = np.zeros(shape=(20,15), dtype=np.uint32);
    for row in chunks:
        for i,letter in enumerate(row):
            if(aa2i_dict.has_key(letter)):
                out[aa2i_dict[letter], i] += 1;
            
        
    return out;

    

def print_count(count_matrix):
    positions = range(-7,8);
    col_headers = ','.join(['{:6d}'.format(p) for p in positions]);
    print '      ', col_headers;
    for i in range(0,count_matrix.shape[0]):
        res = i2aa_dict[i];
        print '   ', res, ':', ','.join(['{:6d}'.format(n) for n in count_matrix[i,:]])

def print_binom(binomial_matrix):
    positions = range(-7,8);
    col_headers = ','.join(['{:6d}'.format(p) for p in positions]);
    print '     ', col_headers;
    for i in range(0,binomial_matrix.shape[0]):
        res = i2aa_dict[i];
        print '   ', res, ':', ','.join(['{:-6.1f}'.format(n) for n in binomial_matrix[i,:]])



"""
subset_formA returns a modified formA that fixes a residue at a position

pos represents the -7 to +7 location of the residue to fix
res is a character representing the residue
"""
def subset_formA(formA, pos, res):
    #Using the condition, create a mask by which to multiply each row.
    #then just use matrix multiplication.
    index = aa2i_dict[res] + (pos+7) * 20;
    new_rows = np.where(formA[:,index] == 1)[0]; #returns a tuple, need first element
    
    if(len(new_rows) == 0):
        raise Exception("Error in subset_formA: No matches. Position: " + str(pos) + ", Residue: " + res);
        
    return formA[new_rows,:];


#Takes two (20x15) count matrices as input. Returns a 20x15 probability matrix
def Binomial_Probability_Matrix(data, background):
    #Check if the data/background are in FormA, and if so, convert to count matrix
    if(data.shape[1] == 300):
        data = formA_to_count(data);
    
    if(background.shape[1] == 300):
        background = formA_to_count(background);
    
    from scipy.stats import binom;
    M = np.sum(data, axis=0);  #sum for each position
    M = M[0]; #Just take the first sum for now, examine this later
    
    num_bg = np.sum(background, axis=0)[0];
    prob_bg = background / float(num_bg);
    
    out = binom.sf(data, M, prob_bg, loc=1);
    out = -1*np.log10(out);
    return out;


#sequence is a string.  Each is 15 characters long
#Motifs are a list of lists of chars.  Each is 15 characters long
def motif_matches(sequence, motifs):
    sequence = sequence.upper();
    total_matches = 0;
    for motif in motifs:
        reg = re.compile(''.join(motif));
        if(reg.match(sequence)):
            total_matches += 1;
            print "-------------------";
            print ''.join(motif);
            print sequence;
    
    return total_matches;


