import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import fieldmath
import frame_thy
import ians_numers
import math


f125 = fieldmath.FieldExtension(fieldmath.Zp(5),[1,1,0,1])
f25 = fieldmath.FieldExtension(fieldmath.Zp(5),[3,0,1])

def iter_numbers_2():
    keys = ["p", "d", "n", "s", "k", "sqrt(2d-1)", "n-1="]
    results = []
    for s in range(6, 100):
        for p in [x for x in range(3,100) if ians_numers.isprime(x)]:
            if s % p != 0:
                one_half = fieldmath.Zp(p).reciprocal(2)
                smallest_d = ((s**2+1)*one_half) % p
                if smallest_d % p != 0: 
                    for k in range(1, 9):
                        d = smallest_d+p*k
                        n = 2*d 
                        if s < d:
                            prime_fac = ians_numers.prime_factors(n-1)
                            if len(list(set(prime_fac))) == 1:
                                if ians_numers.is_sum_of_squares(n, prime_fac):
                                    if math.comb(n,s+1) < 2100000:
                                        yield (p, d, n, s, k, math.sqrt(n-1), f"{prime_fac[0]}^{len(prime_fac)}")


def iter_numbers():
    keys = ["p", "d", "n", "s", "k", "sqrt(2d-1)", "n-1="]
    results = []
    for s in range(2, 100):
        for k in range(2, 5):
            for p in [x for x in range(100) if ians_numers.isprime(x)]:
                if ((s*s)-(k*s)-k+1) % p == 0 and (k*(s+1)) % 2 == 0 and 2*s < k*(s+1) and s % p != 0:
                    #print(f"s: {s}, k: {k}, p:{p}, so n=2d={k*(s+1)}, d={k*(s+1)//2} and 1/mu^2 = s^2 mod p: {True if ((k*(s+1)-1) -(s*s)) % p == 0 else False}")
                    d = k*(s+1)//2
                    n = k*(s+1)
                    # Is constructable conference matrix
                    prime_fac = ians_numers.prime_factors(n-1)
                    if len(list(set(prime_fac))) == 1:
                        if ians_numers.is_sum_of_squares(n, prime_fac):
                            yield (p, d, n, s, k, math.sqrt(n-1), f"{prime_fac[0]}^{len(prime_fac)}")


for (p, d, n, s, k, sqrt_real, prime_decomp) in iter_numbers_2():
    print(f"Testing: p={p}, d={d}, n={n}, s={s}")
    if prime_decomp == "5^2":
        conf_fld = f25
    elif prime_decomp == "5^3":
        conf_fld = f125
    elif prime_decomp.split("^")[-1] == "1":
        conf_fld = fieldmath.Zp(n-1)
    else:
        print(f"Prime decomp is hard: {prime_decomp}")
        continue

    fld = fieldmath.Zp(p)
    G = frame_thy.create_gram_of_d_2d_etf_from_conference_mat(
        frame_thy.create_conference_matrix(
            conf_fld, 
            fieldmath.Zp(p)
        ), 
        fld, 
        s, 
        True
        )
    
    compute_discrim = False

    is_G_frm, (rnk_G, disc_G) = frame_thy.is_frame(None, G, compute_discrim)
    isetf, (a,b,c_G) = frame_thy.is_etf(None, G, is_G_frm)

    assert is_G_frm and isetf, "yikes"

    # We can test some possible frame vectors
    its = 0
    print("There are: ", (math.comb(G.rows, s+1)), "things to check")

    # Note I do believe $s=13 as well
    for sub_frame_G, sub_frame_cols in fieldmath.iter_sub_mats(G, s+1, randomize=False):
        sub_frame_G_cols = G.get_sub_matrix_from_cols(sub_frame_cols)
        if its % 1000 == 999:
            print(".", end='', flush=True)
        its+=1

        if sub_frame_G_cols.rank() != s:
            continue
        
        is_Gsub_frm, (rnk_Gsub, disc_Gsub) = frame_thy.is_frame(None, sub_frame_G, compute_discrim)
        is_tgt, c = frame_thy.is_tight(None, sub_frame_G)
        if is_tgt and is_Gsub_frm and sub_frame_G_cols.rank() == rnk_Gsub:
            print(f"\nYo! the frame vectors {list(sub_frame_cols)} is a ({a},{b},{c})-ETF, therefore a regular {s}-simplex!")
            break
    print("\nSad stuff, no regular simplex here, moving on to the next")