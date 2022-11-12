import numpy as np

def est_pixel_world(pixels, R_wc, t_wc, K):
    """
    Estimate the world coordinates of a point given a set of pixel coordinates.
    The points are assumed to lie on the x-y plane in the world.
    Input:
        pixels: N x 2 coordiantes of pixels
        R_wc: (3, 3) Rotation of camera in world
        t_wc: (3, ) translation from camera to world
        K: 3 x 3 camara intrinsics
    Returns:
        Pw: N x 3 points, the world coordinates of pixels
    """

    ##### STUDENT CODE START #####
    n = len(pixels)
    Pw = np.zeros((n, 3))

    for num, pixel in enumerate(pixels):
        R = np.linalg.inv(R_wc)
        t = (-1 * R @ t_wc.reshape(3, 1)).reshape(3, 1)

        Numerator = (R_wc @ t).reshape(3, 1)
        Denominator = (R_wc @ np.linalg.inv(K) @ np.transpose(np.hstack((pixel, [1])))).reshape(3, 1)

        lmbda = Numerator[2, 0] / Denominator[2, 0]

        Wc = (lmbda*Denominator - Numerator).reshape(3,)
        Pw[num] = Wc

    ##### STUDENT CODE END #####
    return Pw
