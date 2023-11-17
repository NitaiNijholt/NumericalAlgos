import numpy as np

def householder_qr(A):
    # Following Algorithm 3.1 in Heath
    A = A.copy()
    m,n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))

    # Loop over columns
    for k in range(min(n,m-1)):

        # Compute Householder vector for current col
        a_k = -np.sign(A[k,k]) * np.linalg.norm(A[k,k:])
        v_k = np.zeros((1,m))
        v_k[0,k:] = A[k,k:]
        e_k = np.zeros((m, 1))      
        e_k[k,:] = 1
        v_k = v_k.T - a_k * e_k
        b_k = v_k.T @ v_k

        # Skip current column if its already zero
        if b_k != 0:

            # Apply trnasformation to remaining submatrix
            for j in range(k,n):
                gamma_j = v_k.T @ A[:,j]
                
                # BUG HERE
                A[:,j] = A[:,j] - (2*gamma_j / b_k)*v_k.T
                Q[:,k:] = Q[:,k:] - (2*gamma_j / b_k)*v_k
                R[k:,k:] = A[k:,k:]
    
    return Q,R


# https://www.tutorialexample.com/best-practice-to-numpy-create-hilbert-matrix-for-beginners-numpy-tutorial/
def hilbert(n):
    x = np.arange(1, n+1) + np.arange(0, n)[:, np.newaxis]
    return 1.0/x



A = hilbert(5)
Q , R = householder_qr(A)

print(A)
print(np.multiply(Q,R))
print('END')
