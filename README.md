# Augmented-Reality-using-April-Tags

In this repository, I have implemented an augmented reality application to output a video that contains several virtual object models as if they exist in the real
world. Furthermore, I have specified pixel positions to place an arbitrary object.

In order to do this, I began by recovering the camera poses through two different approaches: 

1) solving the Perspective-N-Point (PnP) problem with coplanar assumption

2) solving the Persepective-three-point (P3P) and the Procrustes problem

After retrieving the 3D relationship between the camera and world, we can place an arbitrary objects in the scene.

Note: Although the Open-CV library has been used, I have implemented the above algorithms without calling the OpenCV functions as a black box. 

Output: 

![gifoutput](https://user-images.githubusercontent.com/72302800/201510512-3d5b2ee7-5981-4bd9-afca-39401af07fae.gif)
