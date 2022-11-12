import numpy as np

def P3P(Pc, Pw, K=np.eye(3)):
    """
    Solve Perspective-3-Point problem, given correspondence and intrinsic

    Input:
        Pc: 4x2 numpy array of pixel coordinate of the April tag corners in (x,y) format
        Pw: 4x3 numpy array of world coordinate of the April tag corners in (x,y,z) format
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3,) numpy array describing camera translation in the world (t_wc)

    """
    ##### STUDENT CODE START #####
    # Invoke Procrustes function to find R, t
    # You may need to select the R and t that could transoform all 4 points correctly. 
    # R,t = Procrustes(Pc_3d, Pw[1:4])
    ##### STUDENT CODE END #####
    f = (K[0,0] + K[1,1])/2
    ux = K[0,2]
    uy = K[1,2]
    
    remove_k = np.array([ux, uy])

    Pc_new = Pc - remove_k

    calibrated_coords = (np.hstack((Pc_new, f*np.ones((4,1))))).T

    w1 = Pw[0,:].reshape(3,1)
    w2 = Pw[1,:].reshape(3,1)
    w3 = Pw[2,:].reshape(3,1)
    w4 = Pw[3,:].reshape(3,1)

    q1 = calibrated_coords[:,0]
    q2 = calibrated_coords[:,1]
    q3 = calibrated_coords[:,2]
    q4 = calibrated_coords[:,3]


    j1 = q1/np.linalg.norm(q1)
    j2 = q2/np.linalg.norm(q2)
    j3 = q3/np.linalg.norm(q3)

    cos_alpha = np.dot(j2, j3)
    cos_beta = np.dot(j1, j3)
    cos_gamma = np.dot(j1, j2)

    a = np.linalg.norm(w3-w2)
    b = np.linalg.norm(w3-w1)
    c = np.linalg.norm(w2-w1)

    brac1 = (a**2 - c**2)/b**2
    brac2 = (a**2 + c**2)/b**2
    brac3 = (b**2 - c**2)/b**2
    brac4 = (b**2 - a**2)/b**2

    A4 = ((brac1 - 1)**2) - ((4*(c**2)/(b**2))*((cos_alpha)**2))

    A3 = 4*((brac1*(1-brac1)*cos_beta)
    - ((1-brac2)*cos_alpha*cos_gamma)
    + ((2*(c**2)/(b**2))*((cos_alpha)**2)*cos_beta))
    
    A2 = 2*((brac1**2) 
    - (1) 
    + (2*(brac1**2)*((cos_beta)**2))
    + (2*brac3*((cos_alpha)**2))
    - (4*brac2*cos_alpha*cos_beta*cos_gamma)
    + (2*brac4*((cos_gamma)**2)))

    A1 = 4*(((-1)*brac1*((1+brac1)*cos_beta)) 
    + ((2*(a**2)/(b**2))*(((cos_gamma)**2)*cos_beta))
    - ((1-brac2)*cos_alpha*cos_gamma))

    A0 = (((1+brac1)**2) - (4*(a**2)/(b**2))*((cos_gamma)**2))

    coefficients = [A4, A3, A2, A1, A0]
    roots = np.roots(coefficients)
    
    bracket_stuff = ((a**2 - c**2)/b**2)
    final_roots = []
    for root in roots:
        if root.imag == 0 and root.real > 0:
            u = (((-1+bracket_stuff)*(root.real**2)) - (2*(bracket_stuff)*cos_beta*root.real) + 1 + bracket_stuff)/(2*(cos_gamma - (root.real*cos_alpha)))
            if u > 0:
                final_roots.append((u,root.real))

    error = []
    R_arr = []
    T_arr = []

    for u,v in final_roots:
        s1 = np.sqrt((c**2)/(1 + (u**2) - (2*u*cos_gamma)))
        s2 = u*s1
        s3 = v*s1
        
        p1 = s1*j1
        p2 = s2*j2
        p3 = s3*j3

        Pc_new = np.vstack((p1, p2, p3))

        R,t = Procrustes(Pc_new, Pw[0:3])

        R_arr.append(R)
        T_arr.append(t)
        R_new = (R.T)
        T_new = (-R_new@t).reshape(3,1)
        
        cmp = K@(R_new@w4 + T_new)
        print(cmp)
        error.append(np.linalg.norm(cmp[:2,:] - Pc[3,:].reshape(2,1)))

    index = np.argmin(error)

    R = R_arr[index]
    t = T_arr[index]
    
    return R, t

def Procrustes(X, Y):
    """
    Solve Procrustes: Y = RX + t

    Input:
        X: Nx3 numpy array of N points in camera coordinate (returned by your P3P)
        Y: Nx3 numpy array of N points in world coordinate
    Returns:
        R: 3x3 numpy array describing camera orientation in the world (R_wc)
        t: (3,) numpy array describing camera translation in the world (t_wc)

    """

    ##### STUDENT CODE START #####

    ##### STUDENT CODE END #####
    A = Y.T
    B = X.T

    # centroid
    A_bar = (np.mean(A, axis = 1)).reshape(3,1)
    B_bar = (np.mean(B, axis = 1)).reshape(3,1)

    A = A - A_bar
    B = B - B_bar

    R_hat = A @ B.T

    [U, S, Vt] = np.linalg.svd(R_hat)

    det_uvt = np.linalg.det(U @ Vt)
    diag_mat = np.array([[1, 0, 0], [0, 1, 0], [0, 0, det_uvt]])
    R = U @ diag_mat @ Vt

    t = (A_bar - (R @ B_bar)).reshape(3,)

    return R, t

# X = np.array([[ 0.12291957,  0.06189272,  0.7015967],
#  [ 0.05758767,  0.00441186,  0.81126754],
#  [-0.06603083,  0.02759465,  0.74977749]])
# Y = np.array([[ 0.07, -0.07,  0.  ],
#  [ 0.07,  0.07,  0.  ],
#  [-0.07,  0.07,  0.  ]])

Pc = np.array([[304.28405762 ,346.36758423],
 [449.04196167 ,308.92901611],
 [363.24179077 ,240.77729797],
 [232.29425049 ,266.60055542]])

Pw = np.array([[-0.07 ,-0.07 , 0],
 [ 0.07 ,-0.07 ,0],
 [ 0.07  ,0.07 , 0],
 [-0.07  ,0.07 , 0]])

K = np.array([[823.8,0,304.8],
             [0,822.8,236.3],
             [0,0,1]])

P3P(Pc, Pw, K)

