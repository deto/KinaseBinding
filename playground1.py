# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 22:03:43 2013

@author: daved_000
"""
import FileIO;
import ProcessBackground;
from DataManipulation import *;

from motifx1 import motif_x_iter;

match_kinase = "PKCD";

mismatch_kinases = ['CDK5','GSK3B','p38-alpha'];

start_center = 'S';         
starting_result = list('.......' + start_center + '.......');

background = ProcessBackground.load_fasta_background('uniprot-%28organism%3A9606+keyword%3A1185%29+AND+reviewed%3Ayes.fasta',
                                                      center = start_center);


mismatch_motifs = list();
for kinase in mismatch_kinases:
    kinase_data = FileIO.load_specific_kinase(kinase, start_center);
    motifs = motif_x_iter(kinase_data['formA'], background['formA'], starting_result);
    mismatch_motifs.extend(motifs);
    


match_kinase_data = FileIO.load_specific_kinase(match_kinase, start_center);

match_counts = np.zeros(len(match_kinase_data['seq']));
for i,seq in enumerate(match_kinase_data['seq']):
    match_counts[i] = motif_matches(seq, mismatch_motifs);

