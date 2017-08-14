# 3D-face-feature-etraction-base-on-2D-pic
using 2D pic to learn the facial features then map them to 3D mesh model, extracting the feature from vertex nor. Also trying curvature filter &amp; cluster and spheric harmonics way  

This project is built by vs2015.

Libs:
  
  use opengl, boost, eigen3. You must download the glut , the boost and the eigen3 for drawing, math and matrix calculation.


Files:
  
  in 2d_feature dir, there are 11 3D face. 
  
  
    mesh.bin ，stores the triangle mesh model with structure :

      nv(an integer ，describes the number of vertexs) 

      nf(an integer ，describes the number of meshs) 

      vb(nv * 3 floats ，describe the coordinate of the nv vertex)

      ib(nf * 3 integers ，describe the vertex index of each triangle mesh)
   
   
    marks.txt : 

      the feature points learned from 2D pic(no 2D pic upload cause the privacy problem) 

      从二维平面图中提取的特征点在三维空间上对应的坐标


    sharp_point.txt :

      the coordinate of the vertexs which pass curvature check

      满足曲率一定大小的点的坐标
  
  
    DBScan_result.txt :

      clustering the sharp_points by DBScan, then find the center of each classes. This file stores the coordinate of the centers and something for drawing the model.

      对过滤后的点用DBscan方法进行聚类，这个文件存储了这些中心的坐标以及一些用于绘图的数据
  
  
    DBScan_feature.txt :  

      mark the center as eyes, ears, nose tip artificially. the  variance of the  vertex normals of the vertexs around the facial features are the final feature of a 3D face.

      人为地对聚类后的点标注为眼睛，耳朵，鼻头。以这些特征点的邻域的点的法向量方差作为最后的特征点。
  
  
    spheric_harmonic_feature.txt ：

      the method refers to the paper "Rotation Invariant Spherical Harmonic Representation of 3D Shape Descriptors" (not sure about the corectness)

      采用文章《Rotation Invariant Spherical Harmonic Representation of 3D Shape Descriptors》的方法（不能保证其正确性）
  
  
    feature_from_2d.txt : 

      the  variance of the  vertex normals of the vertexs around the 2D facial features are the final feature of a 3D face.

      利用2D特征点周围的点的法向量方差作为特征点


  in dir test, there are the code and the project files(Mesh.h cpp, drawface.cpp, Vertex.h cpp , the 5 files are useless)
