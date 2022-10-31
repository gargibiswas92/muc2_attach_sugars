#!/usr/bin/python

import sys
import math
from math import sqrt
from sympy.solvers import nonlinsolve
from sympy.solvers import solve
from sympy import symbols
import random

infile1 = open("/home_b/gargi/project_muc/attach_sugar_mucin/muc2_smog20_tail.gro", "r")
outfile3 = open("/home_b/gargi/project_muc/attach_sugar_mucin/glycan_pr2.txt", "w")

res_num = []
count_line = 0
res_name = []

for lines in infile1:
    a1 = int(lines[0:5])
    a2 = str(lines[5:8])
    res_num.append(a1)
    res_name.append(a2)
    count_line = count_line + 1

#for lines in infile1:
#    a, b, c, d, e, f, g, h, i = lines.split()
#    res_num.append(int(b))
#    res_name.append(d)
#    count_line = count_line + 1
    
count_gly_site = 0
gly_site = []

for i in range(count_line):
    if res_num[i] >= 16884:
       if res_name[i] == 'SER' or res_name[i] == 'THR':
          #print(res_num[i])
          gly_site.append(res_num[i])
          count_gly_site = count_gly_site + 1
#print(count_gly_site)

perc = count_gly_site*0.5
perc_int = int(perc)
#print(perc)

li2 = random.sample(gly_site,perc_int)
#print(li2)
 
for i in range(len(li2)):
    ss = [str(i), " ",str(li2[i]), "\n"]
    outfile3.writelines(ss)      
    
