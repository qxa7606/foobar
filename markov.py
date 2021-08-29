from copy import deepcopy
from fractions import Fraction, gcd

# getting cofactor of matrix (used to get determinant and adj)
def get_cofactor(m, ii, jj):
    return [x[:jj] + x[jj + 1:] for x in m[:ii] + m[ii + 1:]]

def get_determinant(m):
    # didn't add error catching because existent is guaranteed
    if len(m) == 1:
        return m[0][0]
    if(len(m) == 2):
        value = m[0][0] * m[1][1] - m[1][0] * m[0][1]
        return value
    total = 0
    for col in range(len(m)):
        sign = (-1) ** (col)
        tmp_det = get_determinant(get_cofactor(m, 0, col))
        total += (sign * m[0][col] * tmp_det)
    return total

# get adjoint of matrix
def get_adj(m):
    l = len(m)
    return zip(*[[get_determinant(get_cofactor(m, i, j)) * (1 if (i + j) % 2 == 0 else -1) for j in range(l)] for i in range(l)])

# get matrix invert
def matrix_inversion(m):
    d = get_determinant(m)
    adj = get_adj(m)
    l = len(adj)
    return [[adj[i][j] * Fraction(1, d) for j in range(l)] for i in range(l)]

def matrix_subtraction(a, b):
    num_rows = len(a)
    num_cols = len(a[0])
    return [[a[i_row][i_col] - b[i_row][i_col] for i_col in range(num_cols)] for i_row in range(num_rows)]

def matrix_multiplication(a, b):
    return [[sum([cols[i] * row[i] for i in range(len(a[0]))]) for cols in zip(*b)] for row in a]

# return a common denominator for each probability
def get_common_denom(l):
    lcm = 1
    for i in [x.denominator for x in l]:
        lcm = lcm*i//gcd(lcm, i)
    return [x.numerator * (lcm / x.denominator) for x in l] + [lcm]

'''
Approach:
Given m, and t meaning transient states and a meaning absorbing states

     t     a
t    Q  |  R
     -------
a    O  |  I

                    -1
The answer is (I - Q)  * R

^^ Will give us a matrix of probabilities starting from any of the transient states,
Since our problem states we always start from S0, the solution is the first row
'''
def solution(m):
    if len(m) < 2:
        return [1,1]
    m = [[Fraction(e, sum(s)) if e else 0 for e in s] for s in m] # convert to fractions
    ab_states = [sum(x) == x[i] for i, x in enumerate(m)] # get absorbing/terminal states
    qq = [[e for i, e in enumerate(s) if not ab_states[i]] for i, s in enumerate(m) if not ab_states[i]] # get Q
    r = [[e for i, e in enumerate(s) if ab_states[i]] for i, s in enumerate(m) if not ab_states[i]] # get R
    ii = [[0 if i != j else 1 for j in range(len(qq))] for i in range(len(qq))] # get I
    iq = matrix_subtraction(ii, qq) # (I - Q)
    iq = matrix_inversion(iq) # (I - Q) inverse
    res = matrix_multiplication(iq, r) # (I - Q)^-1 * R 
    return get_common_denom(res[0])


print(

#     matrix_inversion([
#     [5, -2, 2, 7],
#     [1, 0, 0, 3],
#     [-3, 1, 5, 0],
#     [3, -1, -9, 4]
# ])

solution([
            [1, 1, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 1, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 1, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [1, 0, 1, 0, 1, 0, 1, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            ])

)