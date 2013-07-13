# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 23:11:18 2013

@author: daved_000
"""
import random;
import numpy as np;
import pyximport;
import time;
pyximport.install(setup_args={"script_args":["--compiler=mingw32"]}, reload_support=True)

import CyTest;
import CyTest2;

#Create a giant list of strings
ss = np.random.randint(0,20, (300000,15));

#Create the same structure in FormA
ss_big = np.random.randint(0,2, (300000,300));

cm = np.zeros((20,15));
def test1():
    start = time.time();
    for i in range(0,20):
        cm[i,:] = np.sum(ss==i, axis=0);
    print time.time() - start;

def test2():
    start = time.time();
    out = sum(ss_big,axis=1);
    print time.time() - start;