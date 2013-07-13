# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 23:38:42 2013

@author: daved_000
"""

def primes(kmax):
    p = [0]*1000;
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
    """Fibonacci series test"""
    total = 0;
    while total < n:
        total = total + 1;
    

