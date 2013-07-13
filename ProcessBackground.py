# -*- coding: utf-8 -*-
"""
Created on Sat Jul 06 15:41:04 2013

@author: daved_000

This will contain vital steps for processing the background information

"""

"""
Takes background data and groups it into as many N-length sequences as possible

'fasta_in' is a list of 2-element lists.  For each 2-element list, the first
element is the FASTA ID and the second is the sequence

'N' determines the lenght of sequences to create
"""
import FileIO;
import DataManipulation;


def fasta_to_chunks(fasta_in, N=15):
    out = list();
    for row in fasta_in:
        new_list = [row[1][i:i+N] for i in range(0,len(row[1])-N+1)];
        out.extend(new_list);
    
    return out;

#Selects only sequences whose center residue is 'res'
def filter_chunks(chunks_in, res):
    N = len(chunks_in[0]);
    center = N/2;
    out = [x.upper() for x in chunks_in if x[center] == res];
    
    return out;

def load_fasta_background(filename, center='.'):
    fasta_list = FileIO.read_fasta(filename);
    result = fasta_to_chunks(fasta_list);
    
    if(center != '.'):
        result = filter_chunks(result, center);
        
    result = result[0:300000]; #Concatenate because I don't have enough ram    
    formA = DataManipulation.list_to_formA(result);
    
    return dict([('seq',result),
                 ('formA',formA)]);