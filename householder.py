import numpy as np


# https://www.tutorialexample.com/best-practice-to-numpy-create-hilbert-matrix-for-beginners-numpy-tutorial/
def hilbert(n):
    x = np.arange(1, n+1) + np.arange(0, n)[:, np.newaxis]
    return 1.0/x


def householder_qr_old(A):
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
                print(gamma_j)

            # P_k = np.eye(m)
            # v_k_v_k_T = v_k @ v_k.T
            # v_k_v_k_T = v_k_v_k_T[k:,k:]
            # P_k[k:, k:] = (np.eye(m-k) - (2 / b_k) * v_k_v_k_T)
            
            # Update Q and R matrices
            # Q = Q @ P_k.T
            # R = P_k @ R

                # BUG HERE
                A[:,j] = A[:,j] - (2*gamma_j / b_k)*v_k.T
                Q[:,k:] = Q[:,k:] - (2*gamma_j / b_k)*v_k
                R[k:,k:] = A[k:,k:]
    
    #print('A loop', A)

    return Q,R

#A = hilbert(5)
#Q, R = householder_qr_old(A)

#print('OLD')

#print(Q)
#print('OLD')

def householder_qr2(A):
    A = A.copy()
    m, n = A.shape
    Q = np.eye(m)
    R = A.copy()
    Hs = []

    for k in range(min(n, m - 1)):
        # Compute Householder vector for current column
        a_k = -np.sign(R[k, k]) * np.linalg.norm(R[k:, k])
        v_k = np.zeros((m, 1))
        v_k[k:, 0] = R[k:, k]
        e_k = np.zeros((m, 1))
        e_k[k, 0] = 1
        v_k = v_k - a_k * e_k
        b_k = v_k.T @ v_k

        # Skip current column if it's already zero
        if b_k != 0:
            # Apply transformation to remaining submatrix
            P_k = np.eye(m)

            v_k_v_k_T = v_k @ v_k.T
            v_k_v_k_T = v_k_v_k_T[k:,k:]
            P_k[k:, k:] = (np.eye(m-k) - (2 / b_k) * v_k_v_k_T)
            
            # Update Q and R matrices
            Q = Q @ P_k.T
            R = P_k @ R

    return Q, R



def householder_qr(A):
    """"Householder QR decomposition."""
    A = A.copy()
    m, n = A.shape
    Q = np.eye(m) #np.eye(m, n)
    R = A.copy()
    minimum = min(m-1,n)

    for k in range(minimum):
        alpha_k = -np.sign(A[k,k]) * np.linalg.norm(A[k:,k])
        e_k = np.zeros(m-k)
        e_k[0] = alpha_k
        # print(e_k)
        # print(e_k.shape)
        
        v_k = A[k:,k] - e_k
        #print(v_k)
        beta_k = v_k.T @ v_k

        if beta_k == 0:
            continue

        # Loop over row instead...
        # V_k obtained from first column shuld be applied
        # to other elements from row 1
        for j in range(k,n):
            if j >= k:
                #print('Inner loop k j',k,j)
                # Apply same v_k to other columns...

                # print('v_k.shape', v_k.shape)
                # print('A[k:,j].shape', A[k:,j].shape)
                # print('v_k.T',  v_k.T.shape)
                # print(v_k.T,A[:,j])

                gamma_j = v_k.T @ A[k:,j]
                
                A[k:,j] = A[k:,j] - (2 * gamma_j/beta_k) * v_k

                #print((2 * gamma_j/beta_k) * v_k)
                #print(v_k)
                #print('IS V_K NORML',np.linalg.norm(v_k))
                
                # First column and main diagonal of Q are OK!
                #print('1',Q,'2',Q[k:,j],'3',v_k.T)
        for j in range(m):     
            gamma_j = v_k.T @ Q[k:,j]
#                gamma_j = v_k.T @ Q[k:,j]
#                Q[k:,j] = Q[k:,j] - (2 * gamma_j/beta_k) * v_k

            Q[k:,j] = Q[k:,j] - (2 * gamma_j/beta_k) * v_k
                #A[k:,j] = A[k:,j] - (2 * gamma_j/beta_k) * v_k
                #Q[k:,j] = Q[k:,j] - (2 * gamma_j/beta_k) * v_k # Fixed line
                #Q[k:,j] = Q[k:,j] - (2 * gamma_j/beta_k) * v_k
                #print(k,j,Q[k:,j])
                #print(R)#Q)
                #print(Q[k:,j])
                #gamma_j = v_k.T @ Q #[k:,j]

                #H = np.eye(m-k)
                #H = H[k:,j] - (2 * gamma_j/beta_k) * v_k
                #print(H)
                #Q = Q @ H                

                #H = np.eye(m)
                #H[k:, k:] = np.eye(m-k) - (2 * np.outer(v_k, v_k.T)) / beta_k
                #Q = Q @ H
                #H = np.eye(m-k) - (2 * np.outer(v_k, v_k.T)) / beta_k
                #H = np.eye(m-k) - (2 * np.outer(v_k, v_k.T)) / beta_k
                #H = H - (2 * np.outer(v_k, v_k)) / beta_k
                #print(H[:,0])
                #print(Q[k:,j].shape,H[:,j].shape)
                #print(Q[k:, k:],H)
                #Q[k:, k:] = Q[k:, k:] @ H

                
                #Q[k:, j] = Q[k:, j] @ H[:,j]
                #Q[j,k:] = Q[j,k:] @ H[:,0]
                #H = np.eye(m-k) - 2 * v_k * v_k.T
                #H = np.outer(v_k,v_k.T)
                #print('V',V)
                #Q[k:,j] = Q[k:,j] @ H
                # Reflection
                
                #Q[k:,j] = Q[k:,j] - (2 * gamma_j/beta_k) * v_k 

                # print('v_k original hilbert?',(2 * gamma_j/beta_k) * v_k)
                # print((2 * gamma_j/beta_k) * v_k)

                # Add extra zeros to v_k to it is always same size=column
                #v_k2 = np.zeros(m)
                #v_k2[k:] = v_k
                # print(v_k2)
                #v_k = A[k:,k] - e_k
                
                # Get H which is square, then v_k square thing...
                
                #H = np.eye(m-k)
                #H = H[k:,j] - (2 * gamma_j/beta_k) * v_k
    R = A
    return Q.T, R


A = hilbert(5)
print('A',A)

Q , R = householder_qr(A)


print('Q',Q)
#print('R',R)
#print(Q @ R)
#print('END')


def modified_gram_schmidt(A):
    """Returns QR decomposition of matrix A obtained by classic gram schmidt orgthonanlization procedure. 
    Q and R are returned as a tuple. Method used specified in Heath Scientific Computing: An introductory survey: 5th edition p.132"""
    A = A.copy()
    m, n = A.shape
    Q = np.zeros((m, n))
    R = np.zeros((n, n))
    for k in range(n):
        R[k,k] = np.linalg.norm(A[:,k])
        if R[k,k] == 0:
            break
        Q[:,k] = A[:, k]/R[k,k]
        # starting at 1, to adhere to condition in Heath
        for j in range(k+1, n):
            R[k,j] = Q[:,k].T @ A[:, j]
            A[:, j] = A[:, j] - R[k,j]*Q[:,k]
    return Q, R


Q , R = modified_gram_schmidt(A)

print('Q GM',Q)
#print('R GM',R)
#print('QR GM',Q @ R)
print('END')