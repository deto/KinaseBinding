# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 14:20:36 2013

@author: daved_000


Attempt to solve for the count matrix using Weave

"""
from scipy.weave import converters
from scipy import weave
import numpy as np;

#Create a giant list of strings
ss = np.random.randint(0,20, (1000000,15));

def slowcount(list_strings):
    cm = np.zeros((20,15));
    for s in list_strings:
        for i,c in enumerate(s):
            cm[c,i] += 1;
    
    return cm


def fastercount(list_strings):
    cm = np.zeros((20,15));
    for i in range(0,15):
        for c in range(0,20):
            cm[c,i] += np.count_nonzero(ss[:,i] == c);
        
    return cm;


def fastestcount(list_strings):
    cm = np.zeros((20,15), 'i'); 
    
    code = """
            #line 29 "CountMatrixAttempt.py"
            int test;
            test = 0;
            for(int s=0; s<list_strings.length()(0); s++)
            {
                for(int i=0; i<15; i++)
                {
                    test = list_strings(s,i);
                    cm(test,i)++;
                }
            }
            """;
    
    out = weave.inline(code, ['cm','list_strings'], 
                       type_converters=converters.blitz,
                       compiler = 'gcc');
                       
    return cm;

