import numpy as np
import sympy
import math
from tqdm import tqdm


def wt(x):
    n = 0
    for i in x:
        n += i[0]
    return n

def mod2(x):
    numer, denom = x.as_numer_denom()
    if numer % 2 == 0:
        return 0
    else:
        return 1

def gaussian(A, s, n, k, l):
    M = np.concatenate((A,s), axis = 1)
    M_rref = sympy.Matrix(M).rref(iszerofunc=lambda x: x % 2==0)
    M_rref = M_rref[0].applyfunc(lambda x: mod2(x))
    R = np.array(M_rref.tolist())

    Test = R[0:n-k-l , 0:n-k-l]
    H1 = R[0:n-k-l, n-k-l:n]
    H2 = R[n-k-l:n-k, n-k-l:n]
    s1 = R[0:n-k-l, n:n+1]
    s2 = R[n-k-l:n-k, n:n+1]
    return Test, H1, H2, s1, s2

"""
reference: 
https://stackoverflow.com/questions/31190182/sympy-solving-matrices-in-a-finite-field by Chris Chudzicki
"""

def kbits(n, k):
    limit=1<<n
    val=(1<<k)-1
    while val<limit:
        yield "{0:0{1}b}".format(val,n)
        minbit=val&-val #rightmost 1 bit
        fillbit = (val+minbit)&~val  #rightmost 0 to the left of that bit
        val = val+minbit | (fillbit//(minbit<<1))-1

"""
reference: 
https://stackoverflow.com/questions/58069431/find-all-binary-strings-of-certain-weight-has-fast-as-possible by Matt Timmermans
"""

#project to first l component
def proj(v, l): 
    w = []
    for i in range(l):
        w.append(v[i][0])
    w = np.array(w)[:, None]
    return w


#### Initialize (n, w, k, H, s) with input.txt
f = open('input.txt', 'r')
f.readline()
n = int(f.readline())
f.readline()
f.readline()
f.readline()
W = int(f.readline())
f.readline()

Data_H = []
while True:
    x = f.readline()
    if(x[0] == '#'):
        break
    Data_H.append(x)

k = len(Data_H)
H = []

for D in Data_H:
    V = []
    for i in range(len(D)-1):
        V.append(int(D[i]))
    H.append(V)
H = np.array(H)
H = np.transpose(H)
I = np.eye(n-k, dtype = np.int64)
H = np.concatenate((I, H), axis = 1)


Data_s = f.readline()
s = []
for i in range(len(Data_s)):
    s.append(int(Data_s[i]))
s = np.array(s)
s = s[:, None]

f.close()
#### Initialize Done

# main function
ans = []

p = 8   #choose between 0 <= p <= W, (should it be multiple of 4 ? )
l = 20   #choose between 0 <= l <= n - k, should be k == l (mod 2)

p1 = int(p/2)
p2 = int(p1/2) 
l1 = int(np.log2(math.comb(p, p1)))

#iterate with P
while True:
    flag = False
    rng = np.random.default_rng()
    I = np.eye(n, dtype = np.int64)
    Permute = rng.permutation(I)
    H_bar = np.matmul(H, Permute)

    #we have to check linearly independency of H_bar and skip if not.
    #we have to do this in modulo 2

    ##if first_col_full_rank(H_bar, n-k-l) == False:
    ##    continue

    #gaussian elimination
    Test, H1, H2, s1, s2 = gaussian(H_bar, s, n, k, l)

    if np.array_equal(Test, np.eye(n-k-l)) == False:
        continue

    L1 = []
    L2 = []
    L3 = []
    L4 = []
    LL1 = []
    LL2 = []
    L = []
    y_set = kbits(int((k+l)/2), p2)
    y1_set = []
    y2_set = []

    for y in y_set:
        y1 = []
        y2 = []
        for i in y:
            y1.append(int(i))
            y2.append(0)
        for i in y:
            y1.append(0)
            y2.append(int(i))
        y1 = np.array(y1)[:, None]
        y2 = np.array(y2)[:, None]
        y1_set.append(y1)
        y2_set.append(y2)


    for y1 in y1_set: 
        L1.append((y1, np.matmul(H2, y1)%2))
        L3.append((y1, np.matmul(H2, y1)%2))

    for y2 in y2_set:
        L2.append((y2, np.matmul(H2, y2)%2))
        L4.append((y2, (np.matmul(H2, y2) + s2)%2))

    t = np.random.randint(0, 2, (l1, 1))

    L1_bar = []
    L2_bar = []
    for v in L1:
        if np.array_equal(proj(v[1], l1)[:min((int)((k+l)/2), l1)], t[:min((int)((k+l)/2), l1)]) :
            L1_bar.append(v)

    if (int)((k+l)/2) < l1:
        for v in L2:
            if np.array_equal(proj(v[1], l1)[(int)((k+l)/2):l1], t[(int)((k+l)/2):l1]) :
                L2_bar.append(v)
    else:
        L2_bar = L2

    for v in (L1_bar):
        for w in (L2_bar):
            z = ((v[0]+w[0])%2, (v[1]+w[1])%2)
            LL1.append(z)


    L3_bar = L4_bar = []
    for v in L3:
        if np.array_equal(proj(v[1], l1)[:min((int)((k+l)/2), l1)], t[:min((int)((k+l)/2), l1)]) :
            L3_bar.append(v)

    if (int)((k+l)/2) < l1:
        for v in L4:
            if np.array_equal(proj(v[1], l1)[(int)((k+l)/2):l1], t[(int)((k+l)/2):l1]) :
                L4_bar.append(v)
    else:
        L4_bar = L4


    for v in (L3_bar):
        for w in (L4_bar):
            z = ((v[0]+w[0])%2, (v[1]+w[1])%2)
            LL2.append(z)

    for x1 in LL1:
        for x2 in LL2:
            if np.array_equal(x1[1], x2[1]) :
                L.append((x1[0] + x2[0])%2)

    for e2 in L:
        e1 = (np.matmul(H1, e2) + s1)%2
        if wt(e1) <= W - p:
            ans = np.matmul(Permute, np.concatenate((e1, e2)))
            flag = True
            break

    if flag == True:
        break

for i in ans:
    print(i[0], end = "")
print()