#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 09:14:54 2018

@author: colosu
"""

import matplotlib.pyplot as plt

for J in range(50):
    FEPfile = "./FEP/fep" + str(J+1) + ".txt"
    Sqfile = "./Sq/Sq" + str(J+1) + ".txt"
    PSqfile = "./PSq/PSq" + str(J+1) + ".txt"
    FEPSqfile = "./FEPSq/fepSq" + str(J+1) + ".png"
    FEPPSqfile = "./FEPPSq/fepPSq" + str(J+1) + ".png"
    
    FEPFile = open(FEPfile, "r")
    SqFile = open(Sqfile, "r")
    PSqFile = open(PSqfile, "r")
    
    
    FEP = []
    for line in FEPFile:
        FEP.append(float(line))
    Sq = []
    for line in SqFile:
        Sq.append(float(line))
    PSq = []
    for line in PSqFile:
        PSq.append(float(line))
    
    FEPFile.close()
    SqFile.close()
    PSqFile.close()
    
    # plot patterns
    fig1 = 'FEP vs Sq'
    plt.figure(fig1)
    plt.title(fig1)
    plt.xlabel("FEP")
    plt.ylabel("Sq")
    #plt.xlim(-0.1,1.1)
    #plt.ylim(-0.1,2)
    plt.scatter(FEP, Sq, marker='o', color='blue', linewidths=0.1, s=20)
    plt.savefig(FEPSqfile)
    plt.clf()
    
    fig2 = 'FEP vs PSq'
    plt.figure(fig2)
    plt.title(fig2)
    plt.xlabel("FEP")
    plt.ylabel("PSq")
    #plt.xlim(-0.1,1.1)
    #plt.ylim(-0.1,1.1)
    plt.scatter(FEP, PSq, marker='o', color='blue', linewidths=0.1, s=20)
    plt.savefig(FEPPSqfile)
    plt.clf()