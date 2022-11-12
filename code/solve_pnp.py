from est_homography import est_homography
import numpy as np

def PnP(Pc, Pw, K=np.eye(3)):
    """
    Solve Perspective-N-Point problem with collineation assumption, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3, ) numpy array describing camera translation in the world (t_wc)

    """

    ##### STUDENT CODE START #####
    # Homography Approach

    Pw = Pw[:,0:2]
    H = est_homography(Pw, Pc)

    H = H/H[2,2]

    H = np.linalg.inv(K) @ H

    h1 = H[:,0].reshape(3,1)
    h2 = H[:,1].reshape(3,1)
    h3 = np.cross(h1, h2, axis = 0).reshape(3,1)
    h3_alt = H[:,2]

    H = np.hstack((h1, h2, h3))

    [U, S, Vt] = np.linalg.svd(H)
    

    # Following slides: Pose from Projective Transformation

    det_uvt = np.linalg.det(U @ Vt)
    diag_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, det_uvt]])

    R_old = U @ diag_mat @ Vt
    R_new = R_old.T

    t_old = h3_alt/np.linalg.norm(h1)
    t_new = -R_old.T @ t_old

    ##### STUDENT CODE END #####

    return R_new, t_new