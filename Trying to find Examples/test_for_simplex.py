import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import ians_numers
import math
import itertools


bf = fieldmath.Zp(5)
f25 = fieldmath.FieldExtension(bf, [1,1,1], involution_pow=5)
d = "-1"

#print(f"(1+1a)*(1+1a) = {f25.print_elm(f25.multiply(d,d))}")

hessa_sic = fieldmath.create_matrix(
    [
        ["1","1","1",         "-1","-1","-1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",     "-1","-1a","-1a^2"],
        ["-1","-1a^2","-1a",  "0","0","0",         "1","1a^2","1a"   ]
    ], f25)*"1"

emilys_fav_frame = fieldmath.create_matrix(
    [
        ["1","1","1",         "1","1","1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",  "1","1a","1a^2"   ],
        ["1","1a^2","1a",     "0","0","0",      "1","1a^2","1a"   ]
    ], f25)

print(frame_thy.is_etf(hessa_sic))

print(frame_thy.contains_simplex(2, hessa_sic, get_all=True))
# My code is wrong! a = 2? 
print(frame_thy.contains_simplex(3, hessa_sic, get_all=True))