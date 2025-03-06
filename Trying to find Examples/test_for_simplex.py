import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import ians_numers
import math
import itertools


bf = fieldmath.Zp(5)
f25 = fieldmath.FieldExtension(bf, [1,1,1], involution_pow=5)

#print(f"(1+1a)*(1+1a) = {f25.print_elm(f25.multiply(d,d))}")

hessa_sic = fieldmath.create_matrix(
    [
        ["1","1","1",         "-1","-1","-1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",     "-1","-1a","-1a^2"],
        ["-1","-1a^2","-1a",  "0","0","0",         "1","1a^2","1a"   ]
    ], f25)

emilys_fav_frame = fieldmath.create_matrix(
    [
        ["1","1","1",         "1","1","1",      "0","0","0"       ],
        ["0","0","0",         "1","1a","1a^2",  "1","1a","1a^2"   ],
        ["1","1a^2","1a",     "0","0","0",      "1","1a^2","1a"   ]
    ], f25)

idk = fieldmath.create_matrix(
    [
        [ 1,  1,  1,  1,  1,  1,  1,  1,  1,  1],
        [ 1,  1,  1,  1, -1, -1, -1, -1, -1, -1],
        [ 1, -1, -1, -1,  1,  1,  1, -1, -1, -1],
        [-1,  1, -1, -1,  1, -1, -1,  1,  1, -1],
        [-1, -1,  1, -1, -1,  1,  -1, 1, -1,  1],
        [-1, -1, -1,  1, -1, -1,  1, -1,  1,  1]
    ], fieldmath.Zp(3))

f9 = fieldmath.FieldExtension(fieldmath.Zp(3), [2, 2, 1], 3)
d = "1a"
d2 = fieldmath.pow_over_field(d, 2, f9)
d3 = fieldmath.pow_over_field(d, 3, f9)
d5 = fieldmath.pow_over_field(d, 5, f9)
d6 = fieldmath.pow_over_field(d, 6, f9)
d7 = fieldmath.pow_over_field(d, 7, f9)

idk2 = fieldmath.create_matrix(
    [
        [ 1, 2, 0, 0, 0, 0, 0, 0, d, d, d, d,d3,d3,d3,d3],
        [ d, d,d5,d5,d5,d5,d5,d5, 1, 1, 1, 1,d6,d6,d6,d6],
        [ 0, 0, 0, 0, 0, 0, 1, 2, d, d,d5,d5,d3,d3,d7,d7],
        [ 0, 0, 0, 0, 1, 2, 0, 0, d,d5, d,d5,d3,d7,d3,d7],
        [ 0, 0,d2,d6, 0, 0 ,0 ,0,d7,d3,d3,d7, d,d5,d5, d]
    ], f9)

print(frame_thy.is_frame(idk2))
print(frame_thy.is_etf(idk2))

print(frame_thy.contains_simplex(3, idk2))

#print(frame_thy.contains_simplex(17-3, idk2))