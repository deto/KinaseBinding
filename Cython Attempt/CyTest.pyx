# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 23:38:42 2013

@author: daved_000
"""
import numpy as np
cimport numpy as np
DTYPE = np.int
ctypedef np.int_t DTYPE_t

def primes(int kmax):
    cdef int n, k, i
    cdef int p[1000]
    result = []
    if kmax > 1000:
        kmax = 1000
    k = 0
    n = 2
    while k < kmax:
        i = 0
        while i < k and n % p[i] != 0:
            i = i + 1
        if i == k:
            p[k] = n
            k = k + 1
            result.append(n)
        n = n + 1
    return result

def test_loop(n):
    cdef int total
    total = 0
    while total < n:
        total = total + 1
        
def list_to_count(seqs):
    #cdef np.ndarray cm = np.zeros((26,15), dtype=DTYPE)
    
    cdef int i, N, j, index
    i = 0
    N = len(seqs)
    while i < N:
        j = 0
        while j < 15:
            index = seqs[i][j];     
            #cm[index][j] += 1
            j = j + 1
        
        i = i + 1
    
    #return cm


        