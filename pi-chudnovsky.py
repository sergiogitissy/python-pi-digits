#! /usr/bin/python

import math
from gmpy import mpz
from time import time
#from sympy.mpmath.lib import pi_chudnovsky, pi_agm, fpi
#from sympy.mpmath import pi, mp
from  mpmath import  mp
import sys
import os

NN=0
NN_max=0
max_bs_time=0
def writing_file (is_writ, file_name, buffer):
    if is_writ:
      t=time()
      dec = open(file_name,"w")
      dec.write(buffer)
      dec.close()
      print "Writing ",digits,"digits in ",(time()-t),"secs"
      statinfo = os.stat(file_name)
      print statinfo.st_size,"bytes"


def pi_chudnovsky_gmpy_mpz(digits):
    """
    Compute int(pi * 10**digits)
    
    This is done using Chudnovsky's series with binary splitting
    """
    C = 640320
    C3_OVER_24 = C**3 // 24
    cpt=1
    def bs(a, b):
        """
        Computes the terms for binary splitting the Chudnovsky infinite series

        a(a) = +/- (13591409 + 545140134*a)
        p(a) = (6*a-5)*(2*a-1)*(6*a-1)
        b(a) = 1
        q(a) = a*a*a*C3_OVER_24

        returns P(a,b), Q(a,b) and T(a,b)
        """
        #print "a=",a,"b=",b 
               global NN,NN_max
        global max_bs_time
        NN+=1
        t_start=time()
        if b - a == 1:
            #print "END",
            # Directly compute P(a,a+1), Q(a,a+1) and T(a,a+1)
            if a == 0:
                Pab = Qab = mpz(1)
            else:
                Pab = mpz((6*a-5)*(2*a-1)*(6*a-1))
                Qab = mpz(a*a*a*C3_OVER_24)

            Tab = Pab * (13591409 + 545140134*a) # a(a) * p(a)
            if a & 1:
                Tab = -Tab
            else:
              # Recursively compute P(a,b), Q(a,b) and T(a,b)
              # m is the midpoint of a and b
              m = (a + b) // 2
              # Recursively calculate P(a,m), Q(a,m) and T(a,m)
              t_time=time()-t_start
              if t_time>max_bs_time:
                max_bs_time=t_time
                NN_max=NN

              print "A",
              Pam, Qam, Tam = bs(a, m)

              # Recursively calculate P(m,b), Q(m,b) and T(m,b)
              t_time=time()-t_start
              if t_time>max_bs_time:
                max_bs_time=t_time
                NN_max=NN
              print "B",
              Pmb, Qmb, Tmb = bs(m, b)

              print "Last step"
              # Now combine
              Pab = Pam * Pmb
              Qab = Qam * Qmb
              Tab = Qmb * Tam + Pam * Tmb

            return Pab, Qab, Tab
            
    # how many terms to compute
    DIGITS_PER_TERM = math.log10(C3_OVER_24/6/2/6)
    N = int(digits/DIGITS_PER_TERM + 1)
    print "DIGITS_PER_TERM=",DIGITS_PER_TERM ,"N=",N
    # Calclate P(0,N) and Q(0,N)

    #call bs fucntion
    P, Q, T = bs(0, N)

    one_squared = mpz(10)**(2*digits)
    sqrtC = (10005*one_squared).sqrt()
    print NN,"call to bs, max time in bs",max_bs_time,"for NN=",NN_max
    return (Q*426880*sqrtC)
    
    
 
####################
#   MAIN
####################
if (len(sys.argv)>1):
  prec= int (sys.argv[1])
else:
  prec=10**3


if __name__ == "__main__":

    digits=prec
    t=time()
    pi_chudnovsky_gmpy_str= str(pi_chudnovsky_gmpy_mpz(digits))
    print "pi_chudnovsky: ",digits,"digits (5 last digits)",pi_chudnovsky_gmpy_str[-5:],"time",(time()-t)
    file_name= "pi_"+str(prec)+".txt"
    writing_file (1, file_name, pi_chudnovsky_gmpy_str)
    print pi_chudnovsky_gmpy_str

    #for i in range (1,100):
     #for j in range (10*i-9, 10*i+1):
      #print pi_chudnovsky_gmpy_str[j],
     #print 


raw_input("Hit a key")
