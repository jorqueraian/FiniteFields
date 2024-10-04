import re
import math
import pandas as pd

# https://stackoverflow.com/questions/567222/simple-prime-number-generator-in-python
def isprime(n):
    return re.compile(r'^1?$|^(11+)\1+$').match('1' * n) is None

def prime_factors(n):
    i = 2
    factors = []
    while i * i <= n:
        if n % i:
            i += 1
        else:
            n //= i
            factors.append(i)
    if n > 1:
        factors.append(n)
    return factors


def is_sum_of_squares(n, prime_facs, get_vals= False):
    dist_prime_fac = list(set(prime_facs))
    for prime in dist_prime_fac:
        if (prime - 3) % 4 == 0 and dist_prime_fac.count(prime) % 2 != 0:
            return False
    return True


if __name__ == "__main__":
    # give some s, dimension of simplex, 
    # k, number of simplices that make up entire ETF (column- wise)
    # gives the number of vectors n = 2d = k(s+1), and 
    keys = ["p", "d", "n", "s", "k", "sqrt(2d-1)", "n-1="]
    results = []
    for s in range(1, 100):
        for k in range(2, 100):
            for p in [x for x in range(100) if isprime(x)]:
                if ((s*s)-(k*s)-k+1) % p == 0 and (k*(s+1)) % 2 == 0 and 2*s < k*(s+1):
                    #print(f"s: {s}, k: {k}, p:{p}, so n=2d={k*(s+1)}, d={k*(s+1)//2} and 1/mu^2 = s^2 mod p: {True if ((k*(s+1)-1) -(s*s)) % p == 0 else False}")
                    d = k*(s+1)//2
                    n = k*(s+1)
                    # Is constructable conference matrix
                    prime_fac = prime_factors(n-1)
                    if len(list(set(prime_fac))) == 1:
                        if is_sum_of_squares(n, prime_fac):
                            results.append([p, d, n, s, k, math.sqrt(n-1), f"{prime_fac[0]}^{len(prime_fac)}"])


    nums = pd.DataFrame(results, columns=keys)

    nums.to_csv('numbers.csv', index=False) 