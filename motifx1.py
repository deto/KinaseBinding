# -*- coding: utf-8 -*-
"""
Created on Sat Jul 06 23:40:12 2013

@author: daved_000

First attempt at MOTIF-X
"""
import FileIO;
import ProcessBackground;
import copy;
from DataManipulation import *;

def motif_x_iter(data_A, bg_A, prev_result):
    print ''.join(prev_result);
    data_count = formA_to_count(data_A);
    bg_count = formA_to_count(bg_A);
    
    bpm = Binomial_Probability_Matrix(data_count, bg_count);
    
    
    new_results = list();
    for pos in xrange(0,14):
        for res_i in xrange(0,20):
            score = bpm[res_i, pos];
            if score > 5 and (not np.isinf(score)):
                new_fix = copy.deepcopy(prev_result);
                char_to_fix = i2aa_dict[res_i];
                new_fix[pos] = char_to_fix;
                
                data_A_sub = subset_formA(data_A, pos-7, char_to_fix);
                bg_A_sub = subset_formA(bg_A, pos-7, char_to_fix);
                
                temp_results = motif_x_iter(data_A_sub, bg_A_sub, new_fix);
                new_results.extend(temp_results);
                
            
    
    
    if len(new_results) == 0:
        return [prev_result];
    else:
        return new_results;
    

if(__name__ == '__main__'):
    kinase = 'PKCD';
    start_center = 'S';
    starting_result = list('.......' + start_center + '.......');

    kinase_data = FileIO.load_specific_kinase(kinase, start_center);
        
    background = ProcessBackground.load_fasta_background('uniprot-%28organism%3A9606+keyword%3A1185%29+AND+reviewed%3Ayes.fasta',
                                                      center = start_center);

    out = motif_x_iter(kinase_data['formA'], background['formA'], starting_result);
    print '-----------------';
    for a in out:
        print ''.join(a);

