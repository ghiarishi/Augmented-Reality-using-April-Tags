import numpy as np

def est_Pw(s):
    """
    Estimate the world coordinates of the April tag corners, assuming the world origin
    is at the center of the tag, and that the xy plane is in the plane of the April
    tag with the z axis in the tag's facing direction. See world_setup.jpg for details.
    Input:
        s: side length of the April tag

    Returns:
        Pw: 4x3 numpy array describing the world coordinates of the April tag corners
            in the order of a, b, c, d for row order. See world_setup.jpg for details.

    """

    ##### STUDENT CODE START #####
    S = s*0.5
    a = np.array([-S, -S, 0])
    b = np.array([S, -S, 0])
    c = np.array([S, S, 0])
    d = np.array([-S, S, 0])
    Pw = np.array([a, b, c, d])
    ##### STUDENT CODE END #####

    return Pw