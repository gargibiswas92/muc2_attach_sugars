#!/usr/bin/python

import sys
import types
import math
import numpy as np
from math import sqrt
from sympy.solvers import nonlinsolve
from sympy.solvers import solve
from sympy import symbols

infile1 = open("/home_b/gargi/project_muc/attach_sugar_mucin/muc2_smog20_tail.gro", "r")
infile2 = open("/home_b/gargi/project_muc/attach_sugar_mucin/glycan_pr.txt", "r")
outfile1 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_coord_muc.txt", "w")
outfile2 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_gro_muc.txt", "w")
outfile3 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_def_muc.txt", "w")
outfile4 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_bond_muc.txt", "w")

res_num = []
x_cord_prot = []
y_cord_prot = []
z_cord_prot = []
count_line = 0

for lines in infile1:
    a = lines[0:5]
    b = lines[5:8]
    c = lines[20:28]
    d = lines[28:36]
    e = lines[36:44]
    res_num.append(int(a))
    x_cord_prot.append(float(c))
    y_cord_prot.append(float(d))
    z_cord_prot.append(float(e))
    count_line = count_line + 1
    
res_prot_s = []
count_prot_s = 0
count_no_sug = 0

for lines2 in infile2:
    res_prot_s.append(int(lines2))
    count_prot_s = count_prot_s + 1
    
#print(count_line, count_prot_s)
x,y,z = symbols('x, y, z', real=True)

for i in range(count_prot_s):
    for j in range(count_line):
        if res_prot_s[i] == res_num[j] and res_prot_s[i] <= 40391 and res_prot_s[i] >= 10001:
            x2 = x_cord_prot[j]
            y2 = y_cord_prot[j]
            z2 = z_cord_prot[j]
            x3 = x_cord_prot[j+1]
            y3 = y_cord_prot[j+1]
            z3 = z_cord_prot[j+1]
            x33 = x_cord_prot[j-1]
            y33 = y_cord_prot[j-1]
            z33 = z_cord_prot[j-1]
            vec1 = math.sqrt(((x2 - x3)*(x2 - x3)) + ((y2 - y3)*(y2 - y3)) + ((z2 - z3)*(z2 - z3)))
            vec2 = math.sqrt(((x2 - x33)*(x2 - x33)) + ((y2 - y33)*(y2 - y33)) + ((z2 - z33)*(z2 - z33)))
            mult1 = vec1*math.cos((62.8155*math.pi)/180)*0.37783
            mult2 = vec2*math.cos((88.0857*math.pi)/180)*0.37783
            eq1 = ((x - x2)**2 + (y - y2)**2 + (z - z2)**2) - 0.143
            eq2 = ((x - x2)*(x2 - x3)) + ((y - y2)*(y2 - y3)) + ((z - z2)*(z2 - z3)) - mult1
            eq3 = ((x - x2)*(x2 - x33)) + ((y - y2)*(y2 - y33)) + ((z - z2)*(z2 - z33)) - mult2
            sol1, sol2 = nonlinsolve([eq1, eq2, eq3], [x,y,z])
            
            if type(sol1[0]) == float and type(sol1[1]) == float and type(sol1[2]) == float:
               s1 = float(sol1[0])
               s2 = float(sol1[1])
               s3 = float(sol1[2])
               
               k = int(i) + 40392
               ss = 'HETATM {:>5s}  C1 0VA {:>5s} {:>10.2f} {:>7.2f} {:>7.2f}{}'.format(str(k), str(k), s1, s2, s3, "\n")
               ss2 = '{:>5d}0VA     C1{:>5d} {:>7.3f} {:>7.3f} {:>7.3f}{}'.format(k, k, s1, s2, s3, "\n")
               ss3 = ' {:>5s}       NB_2  {:>5s}    0VA     C1  {:>5s}{}'.format(str(k), str(k), str(k), "\n")

               ss4 = ' {:>5s}   {:>5s}   1        3.778300432e-01 2.000000000e+04{}'.format(str(k), str(res_prot_s[i]), "\n")
               print(s1, s2, s3)

            elif type(sol2[0]) == float and type(sol2[1]) == float and type(sol2[2]) == float:
                 s1 = float(sol2[0])
                 s2 = float(sol2[1])
                 s3 = float(sol2[2])
                 
                 k = int(i) + 40392
                 ss = 'HETATM {:>5s}  C1 0VA {:>5s} {:>10.2f} {:>7.2f} {:>7.2f}{}'.format(str(k), str(k), s1, s2, s3, "\n")
                 ss2 = '{:>5d}0VA     C1{:>5d} {:>7.3f} {:>7.3f} {:>7.3f}{}'.format(k, k, s1, s2, s3, "\n")
                 ss3 = ' {:>5s}       NB_2  {:>5s}    0VA     C1  {:>5s}{}'.format(str(k), str(k), str(k), "\n")

                 ss4 = ' {:>5s}   {:>5s}   1        3.778300432e-01 2.000000000e+04{}'.format(str(k), str(res_prot_s[i]), "\n")
                 print(s1, s2, s3)
        
            else:
                count_no_sug = count_no_sug + 1
                print("No real root", "\n")
#            if np.iscomplex(sol1[0]) == False and np.iscomplex(sol1[1]) == False and np.iscomplex(sol1[2]) == False:
#               s1 = float(sol1[0])
#               s2 = float(sol1[1])
#               s3 = float(sol1[2])

#            elif np.iscomplex(sol2[0]) == False and np.iscomplex(sol2[1]) == False and np.iscomplex(sol2[2]) == False:
#                 s1 = float(sol2[0])
#                 s2 = float(sol2[1])
#                 s3 = float(sol2[2])
        
#            else:
#                count_no_sug = count_no_sug + 1
#                print("No real root", "\n")
            k = int(i) + 40392
            ss = 'HETATM {:>5s}  C1 0VA {:>5s} {:>10.2f} {:>7.2f} {:>7.2f}{}'.format(str(k), str(k), s1, s2, s3, "\n")
            ss2 = '{:>5d}0VA     C1{:>5d} {:>7.3f} {:>7.3f} {:>7.3f}{}'.format(k, k, s1, s2, s3, "\n")
            ss3 = ' {:>5s}       NB_2  {:>5s}    0VA     C1  {:>5s}{}'.format(str(k), str(k), str(k), "\n")

            ss4 = ' {:>5s}   {:>5s}   1        3.778300432e-01 2.000000000e+04{}'.format(str(k), str(res_prot_s[i]), "\n")
            print(s1, s2, s3)
            outfile1.writelines(ss)
            outfile2.writelines(ss2)
            outfile3.writelines(ss3)
            outfile4.writelines(ss4)
            
outfile1.close()
infile1.close()
infile2.close()
outfile2.close()
outfile3.close()
    
