import numpy as np

def est_homography(X, Y):
    """
    Calculates the homography H of two planes such that Y ~ H*X
    If you want to use this function for hw5, you need to figure out
    what X and Y should be.
    Input:
        X: 4x2 matrix of (x,y) coordinates
        Y: 4x2 matrix of (x,y) coordinates
    Returns:
        H: 3x3 homogeneours transformation matrix s.t. Y ~ H*X

    """

    ##### STUDENT CODE START #####
    x_prime = [Y[0][0], Y[1][0], Y[2][0], Y[3][0]]
    y_prime = [Y[0][1], Y[1][1], Y[2][1], Y[3][1]]
    x = [X[0][0], X[1][0], X[2][0], X[3][0]]
    y = [X[0][1], X[1][1], X[2][1], X[3][1]]

    A = []
    for i in range(4):
        A.append([-x[i], -y[i], -1, 0, 0, 0, x[i] * x_prime[i], y[i] * x_prime[i], x_prime[i]])
        A.append([0, 0, 0, -x[i], -y[i], -1, x[i] * y_prime[i], y[i] * y_prime[i], y_prime[i]])

    [U, S, Vt] = np.linalg.svd(A)

    H = np.array(Vt[-1])
    H = np.array(H, dtype=np.float32)
    H = H.reshape(3, 3)
    ##### STUDENT CODE END #####

    return H