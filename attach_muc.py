#!/usr/bin/python

import sys
import math
from math import sqrt
from sympy.solvers import nonlinsolve
from sympy.solvers import solve
from sympy import symbols

infile1 = open("/home_b/gargi/project_muc/attach_sugar_mucin/muc2_smog20_tail.gro", "r")
infile2 = open("/home_b/gargi/project_muc/attach_sugar_mucin/glycan_pr.txt", "r")
outfile1 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_coordinates.txt", "w")
outfile2 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_grofile.txt", "w")
outfile3 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_mol_def.txt", "w")
outfile4 = open("/home_b/gargi/project_muc/attach_sugar_mucin/sug_bond_def.txt", "w")
outfile5 = open("/home_b/gargi/project_muc/attach_sugar_mucin/left_out.txt", "w")

res_num = []
x_cord_prot = []
y_cord_prot = []
z_cord_prot = []
count_line = 0

for line in infile1:
    a = line[0:5]
    b = line[5:8]
    c = line[20:28]
    d = line[28:36]
    e = line[36:44]
    
    res_num.append(int(a))
    x_cord_prot.append(float(c))
    y_cord_prot.append(float(d))
    z_cord_prot.append(float(e))
    count_line = count_line + 1
    
res_prot_s = []
count_prot_s = 0

for lines2 in infile2:
    res_prot_s.append(int(lines2))
    count_prot_s = count_prot_s + 1
    
#print(count_line, count_prot_s)
x,y,z = symbols('x, y, z', real=True)
no_gly_count = 0
no_gly = []

for i in range(count_prot_s):
    for j in range(count_line):
        if res_prot_s[i] == res_num[j] and res_prot_s[i] <= 40391 and res_prot_s[i] >= 16884:
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
            m = str(sol1[0])
            n = str(sol1[1])
            q = str(sol1[2])
      
            if (('I' in m) == False) and (('I' in n) == False) and (('I' in q) == False):
                s1 = float(m)
                s2 = float(n)
                s3 = float(q)
            else:
                no_gly_count = no_gly_count + 1
                p = res_num[j]
                no_gly.append(p)
                sp = [str(res_num[j]), "\n"]
                outfile5.writelines(sp)

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
    
