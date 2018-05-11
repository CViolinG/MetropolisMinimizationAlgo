from sys import exit
from prody import *
import numpy as np
import scipy as sp
from math import exp,log,sqrt
from random import *

dt = 1
k = 0
total=0
N=parseArray('N.txt', delimiter=", ")
n=N.shape[0]
m=N.shape[1]
print n
if n!=m:
    exit("Square matrix is required")

def writePtxt(): #writes a naive probability for each state
    file = open("P.txt", "w+")
    total = 0
    for i in range(N.shape[0]):
        total+=N[i,i]
    Pwrite = np.zeros(N.shape[0])
    for i in range(N.shape[0]):
        Pwrite[i] = N[i,i] / float(total)
        write = "%s\n"%Pwrite[i]
        file.write(write)
    file.close()

def writePtxtEvenProbs():
    file = open("P.txt", "w+")

    for i in range(N.shape[0]):
        file.write("0.04167\n")
    file.close()

writePtxt()
#writePtxtEvenProbs()

#use Naieve P.txt matrix by histogram estimation for P values
#Also try with P.txt from Detailed balance


####Is the Offf Diagonal ratio the same for any polynomial or just the exponential
### Literature search this. Properties of Tridiagonal Matrices


P=parseArray('P.txt')
m=P.shape[0]
if n!=m:
    exit("P and N are inconsistent")

R=np.zeros((n,n))
for i in range(n):
    for j in range(n):
        if i!=j:
            R[i,j] = P[i] * (N[i,j] + N[j,i]) / (N[i,i]*P[j] + N[j,j]*P[i])
        else:
            R[i,i]=0
            for k in range(n):
                if(i!=k):
                    R[i,i]-=P[k]*(N[k,i]+N[i,k]) / (N[i,i]*P[k] + N[k,k]*P[i])

writeArray("R0.txt",R,format="%f")
lmda=np.linalg.eigvals(R)
writeArray("lmda0.txt",lmda,format="%f")
r=sp.linalg.expm(R)
writeArray('testest.txt',r,format="%f")
D=np.zeros((n-1))
F=np.zeros((n-1))
for i in range(n-1):
    D[i]=R[i,i+1]*sqrt(P[i+1]/P[i])
writeArray("D0.txt",D,format="%f")


logl=0
loglt=0
for i in range(n):
    for j in range(n):
        logl+=log(r[i,j])*N[i,j]

print 0,logl,1

seed()
loglt=0
burnInReject=0
burnInAccept=0
ratio=0
for t in range(1,1000000):
    for i in range(1, n-1):
        dr=(random()-0.5)*0.05
        R[i+1,i]+=dr
        #R[i-1,i]-=dr # This is an extra criterion
        R[i,i] = -1 * (R[i+1,i] + R[i-1,i])
        r=sp.linalg.expm(R)
        loglt=0
    for ii in range(n):
        for jj in range(n):
            if(r[ii,jj]<0.000000000000000000000001):
                loglt+=0
            else:
                loglt+=log(r[ii,jj])*N[ii,jj]
    if(t==1000):
        ratio = min(float(burnInAccept)/burnInReject, float(burnInReject),burnInAccept)
        print "Burn In Completed. Acceptance Ratio set to : %s" % ratio
    if loglt<logl and random()> ratio:#exp((loglt-logl)):#reject Change
        R[i+1,i]+=dr
        R[i-1,i]-=dr*P[j]/P[i]
        if(t>1000):
            if t%1000 == 0:
                print t,loglt,0
            #k+=1
        else:
            burnInReject+=1
        if(k==80):
            break
    else:#Accept Change
        if(t>1000):
            if t%1000 == 0:
                print t,loglt,k, logl
            k=0
        else:
            burnInAccept+=1
        logl=loglt
writeArray("R.txt",R,format="%f")
lmda=np.linalg.eigvals(R)
writeArray("lmda.txt",lmda,format="%f")
writeArray("r.txt",r,format="%f")
r=sp.linalg.expm(R)
for i in range(n-1):
    D[i]=R[i,i+1]*sqrt(P[i+1]/P[i])
writeArray("D.txt",D,format="%f")
for i in range(1, n-1):
    F[i]  +=    F[i-1] - -0.6 * np.log(abs(R[i,i+1]/R[i+1,i]))
writeArray("Free.txt",F,format="%f")
exit("Done!")

