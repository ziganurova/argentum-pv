#!/usr/bin/python

hapNum = 100
hapLen = 10000*2500 #Haplotype length: 
#Effective population sizes
N_0 = 10000.0

mu = 10**(-8) #Mutation rate per site per generation
mu1 = 10**(-8) #Recombination rate per site per generation

print 'SCRM command generator for two splitted populations'
command = './scrm ' + `int(hapNum)` + ' 1 '

theta = 4*mu*N_0*hapLen
rho = 4*mu1*N_0*hapLen

command += '-t ' + `theta` + ' -r ' + `rho` + ' ' + `hapLen` + ' -l 10 '

command += '> test.txt'
print command
