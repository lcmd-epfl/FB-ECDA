#!/usr/bin/python
import sys
import re
import math
import numpy as np
import argparse
from scipy.linalg import fractional_matrix_power

###### Input Parser ######
parser = argparse.ArgumentParser(description='FB-ECDA is a decomposition analysis tool of electronic coupling for charge transfer. Currently, it only supports charge transfer between the same type of molecules')
parser.add_argument('lmon1', help='Gaussian log file of molecule 1')
parser.add_argument('lmon2', help='Gaussian log file of molecule 2')
parser.add_argument('ldim', help='Gaussian log file of dimer (molecule 1+2)')
parser.add_argument('pmon1', help='Gaussian pun file of molecule 1')
parser.add_argument('pmon2', help='Gaussian pun file of molecule 2')
parser.add_argument('pdim', help='Gaussian pun file of dimer (molecule 1+2)')
parser.add_argument('Na', help='The index of MO of molecule')
#parser.add_argument('Nb', help='The index of MO of molecule 2')
parser.add_argument('ao', help='# of AO of each element for the chosen basis set')
parser.add_argument('fragid', help='File specifying the way of fragmentization')

args = parser.parse_args()
file_dict = vars(args)
lmon1 = file_dict['lmon1']
lmon2 = file_dict['lmon2']
ldim = file_dict['ldim']
pmon1 = file_dict['pmon1']
pmon2 = file_dict['pmon2']
pdim = file_dict['pdim']
Na = int(file_dict['Na'])
Nb = Na
#Nb = int(file_dict('Nb'))
ao = file_dict['ao']
fragid = file_dict['fragid']


def dict_AO(filename):
    lines = filename.readlines()
    NB_ele = dict()
    for l in lines:
        tmp = l.split()
        NB_ele.update({tmp[0]:int(tmp[1])})
    return NB_ele

def fragment_parser(fragid):
    with open(fragid, 'r') as f:
        line = f.readlines()
    fragidlist = dict()
    for fragment in line:
        fragtmp = []
        tmp = re.split(r'[,\n]', fragment)
        fragname = tmp[0]
        for i in range(1, len(tmp)-1):
            fragtmp.extend(range(int(tmp[i].split('-')[0]),int(tmp[i].split('-')[-1])+1))
        fragidlist.update({fragname:fragtmp})
    return fragidlist

def chunkstring(string, length):
    string = string.strip('\n')
    return list(string[0+i:length+i] for i in range(0, len(string), length))

def Read_Num_Basis(f):
    for l in f:
        if 'primitive' in l:
        ######## Store number of basis in a monomer
            N_basis = int(l.split()[0])
            break
    return N_basis

def Read_Overlap(f,N_basis):
    l = f.readlines()
    tmp_alpha = []
    for i in range(len(l)):
        if '*** Overlap ***' in l[i]:
           tmp_alpha.append(i)
    S = np.zeros((N_basis,N_basis))
    row = int(math.ceil(float(N_basis)/float(5)))
    count = 0
    count2 = 0
    for i in range(row):
        for j in range(N_basis-i*5):
            tmp = [x.replace('D','E') for x in l[tmp_alpha[-1]+count+2+j].split()[1:6]]
            tmp = np.array(tmp,dtype='float')
            #print tmp[-1]
            for k in range(len(tmp)):
                #print count2+j,count2+k
                S[count2+j][count2+k] = tmp[k]
        count += j+2
        count2 += 5
    S = S + S.T - np.diag(S.diagonal())
    return S

def Read_Coeff_Eigen(f,N_basis,mon,N_basis2):
    Coeff = []
    Eigen = []
    l = f.readlines()
    row = int(math.ceil(float(N_basis)/float(5)))
    #print row
    cpbasis = [0]*N_basis2
    for i in range(N_basis):
        Eigen.append(l[(row+1)*i+1][18:33].replace('D','E'))
        tmp2 = []
        for j in range(1,row+1):
            tmp = [x.replace('D','E') for x in chunkstring(l[(row+1)*i+j+1],15)]
            tmp2.extend(tmp)
        if mon == 1:
           Coeff.extend(tmp2)
           Coeff.extend(cpbasis)
        elif mon == 2:
           Coeff.extend(cpbasis)
           Coeff.extend(tmp2)
        else:
           Coeff.extend(tmp2) 

    Coeff = np.array(Coeff,dtype=float)
    Eigen = np.diag(np.array(Eigen,dtype=float))
    return Coeff, Eigen

def Construct_AObasis(f,NB_ele):
    l = f.readlines()
    B_A = []
    for i in range(len(l)):
        if 'Charge' in l[i]:
            for j in range(len(l)):
                if 'Symmetry' in l[i+j+2]:
                    break
                else:
                    B_A.append(NB_ele[l[i+j+1].split()[0]])
            break
    Natom = len(B_A)/2
    return B_A,Natom

def AO_index(B_A,Atom_index):
    start = 1
    for i in range(Atom_index-1):
        start += B_A[i]
    end = start + B_A[Atom_index-1] 
    AO = range(start,end)
    return AO

######## Read number of AO basis of monomer1
with open (ao,'r') as f:
     NB_ele = dict_AO(f)
f.close()
######## Read number of AO basis of monomer1
with open (lmon1,'r') as f:
     N_basis_mon1 = Read_Num_Basis(f)
f.close()
######## Read number of AO basis of monomer2
with open (lmon2,'r') as f:
     N_basis_mon2 = Read_Num_Basis(f)
f.close()
######## Read number of AO basis and overlap matrix in dimer
with open (ldim,'r') as f:
     N_basis_dim = Read_Num_Basis(f)
f.close()
with open (ldim,'r') as f:
     B_A,Natom = Construct_AObasis(f,NB_ele)
f.close()
with open (ldim,'r') as f:
     S = Read_Overlap(f,N_basis_dim)
f.close()
######## Read number of MO coefficient and eigen value of monomer1
with open (pmon1,'r') as f:
     Coeff_mon1,Eigen_mon1 = Read_Coeff_Eigen(f,N_basis_mon1,1,N_basis_mon2)
     Coeff_mon1 = np.reshape(Coeff_mon1,(N_basis_mon1,N_basis_dim))
f.close()
######## Read number of MO coefficient and eigen value of monomer2
with open (pmon2,'r') as f:
     Coeff_mon2,Eigen_mon2 = Read_Coeff_Eigen(f,N_basis_mon2,2,N_basis_mon1)
     Coeff_mon2 = np.reshape(Coeff_mon2,(N_basis_mon2,N_basis_dim))
f.close()
######## Read number of MO coefficient and eigen value of dimer
with open (pdim,'r') as f:
     Coeff_dim,Eigen_dim = Read_Coeff_Eigen(f,N_basis_dim,0,N_basis_dim)
     Coeff_dim = np.reshape(Coeff_dim,(N_basis_dim,N_basis_dim))
f.close()
######## gamma1 with (N_basis_mon1xN_basis_dim) dimesion
Orb1 = Coeff_mon1[Na-1]
Orb2 = Coeff_mon2[Nb-1]
S12 = 0
for i in range(N_basis_dim):
    for j in range(N_basis_dim):
        S12 += Orb1[i]*Orb2[j]*S[i][j]
Smat = np.array([[1.0,S12],[S12,1.0]])
U = fractional_matrix_power(Smat,-0.5)
Orb1_Ld = U[0][0]*Orb1+ U[0][1]*Orb2
Orb2_Ld = U[1][0]*Orb1+ U[1][1]*Orb2

###### Fragmentization ######
fm1_A = []
fm2_A = []
fid1 = []
fid2 = []
frag = fragment_parser(fragid)
for i in frag:
    fid1.append(i)
    fid2.append(i)
    fm1_A.append(frag[i])
    tmp = [x+Natom for x in frag[i]]
    fm2_A.append(tmp)
fm1=[]
fm2=[]
for frag in fm1_A:
    tmp = []
    for index in frag:
        tmp.extend(AO_index(B_A,index))
    fm1.append(tmp)
for frag in fm2_A:
    tmp = []
    for index in frag:
        tmp.extend(AO_index(B_A,index))
    fm2.append(tmp)
Orbm1 = []
Orbm2 = []
for frag in fm1:
    tmp = np.zeros(N_basis_dim)
    for i in range(N_basis_dim):
        if i+1 in frag:
           tmp[i] = Orb1_Ld[i]
        if i+1 > N_basis_mon1:
           tmp[i] = Orb1_Ld[i]/len(fm1)
    Orbm1.append(tmp) 
for frag in fm2:
    tmp = np.zeros(N_basis_dim)
    for i in range(N_basis_dim):
        if i+1 in frag:
           tmp[i] = Orb2_Ld[i]
        if i+1 <= N_basis_mon1:
           tmp[i] = Orb2_Ld[i]/len(fm2)
    Orbm2.append(tmp)
J = []
for i,frag1 in enumerate(Orbm1):
    for j,frag2 in enumerate(Orbm2):
        gamma1 = np.matmul(frag1.T, np.matmul(S, Coeff_dim.T))
        gamma2 = np.matmul(frag2.T, np.matmul(S, Coeff_dim.T))
        tmp = np.matmul(gamma1,np.matmul(Eigen_dim,gamma2.T))*27.2114
        J.append(tmp)
        print "V_"+fid1[i]+'-'+fid2[j]+":", tmp 
print "V_tot:", np.sum(J)
