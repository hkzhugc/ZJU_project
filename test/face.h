#pragma once
#include <vector>
#include "new_vertex.h"
#include <Eigen\Dense>
using namespace Eigen;
typedef struct index3
{
	int idx1;
	int idx2;
	int idx3;
}index3;

typedef struct Curvature
{
	double max_curvature;
	double min_curvature;
}Curvature;

class Face
{
public:
	Face();
	Face(int vertexNum);
	~Face();
	void read_vertex(Vertex v1);
	void readVertex2Matrix(Vertex  v1);
	void read_face(index3  index);
	void read_feature(double feature);
	void read_right_point(Position point);
	void cal_curvature();
	void cal_varOfVertexNNor();//计算方差，并使得点云均值为0
	
	void cal_2d_feature_2_3d_feature(float * vb, size_t size, int face_index);//use the 2d features to cal the point nor var in those areas 
	void rotate_by_2d_feature(float * vb, size_t nv, float * fvb);//rotate the face

	void zeroOneNormalize(double max, double min);//01标准化
	void feature_normalize(std::vector<double> max, std::vector<double> min);//0-1 normalize
	double feature_L2_norm() const;//return the L2 norm of the feature vector
	void DBScan(double eps, int minPts, int face_index);
	void feature_extract(float *vb, size_t size, double r, int face_index);
	static double euler_dis(const Face &f1, const Face &f2);//return euler distance
	static double L1_dis(const Face & f1, const Face & f2);//return L1 dis
	static double inner_product(const Face &f1, const Face &f2);//return inner product
	static double cos_dis(const Face &f1, const Face &f2);//return cos distance
	void add_D2_feature(int dime_num);
	double ICP_rotate(Face & basedFace, int rotate_nor, double tolerance = 0.001);//以baseFace作为基础，将this->face旋转,如果rotate_nor为1,那么以法向量作为定义的距离
	void print_feature();
	
	
	
	Vertex left_eye;	
	Vertex right_eye;
	Vertex nose;	
	Vertex left_ear;	
	Vertex right_ear;
	
	MatrixXd face_position;


	//all from left to right
	Matrix<double, 3, 15> face_curve;
	Matrix<double, 3, 12> eyebrow;
	Matrix<double, 3, 8> eye_curve;
	Matrix<double, 3, 11> nose_curve;
	Matrix<double, 3, 18> lips;
	Matrix<double, 3, 1>  nose_tip;
	Matrix<double, 3, 8> eyeball;
	Matrix<double, 3, 2> eyes;
	Matrix<double, 3, 2> ear_tip;




	static Position default_left_eye;
	static Position default_right_eye;
	static Position default_left_ear;
	static Position default_right_ear;
	static Position default_nose;
	
	int vertex_num;
	int face_num;
	MatrixXd vertex;
	MatrixXd nor;
	int *cluster_set;//并查集，用来表明根据密度聚类后的从属关系。
	int * class_no;
	std::vector<int> cluster_parent;
	std::vector<Position> class_mean;
	std::vector <Position> feature_point;
	std::vector <Vertex> v;//store vertex
	std::vector <index3> f;// store faces
	Curvature * c;//store points curvature
	std::vector <double> acc_face_area;//acc_face_area[i] - acc_face_area[i-1] = area of f[i]
	std::vector <double> feature_vec;//as the name shows
private:
	int nor_num;
};