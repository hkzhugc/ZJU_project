#include "face.h"
#include <time.h>
#include <boost/random.hpp>
#include <boost/random/random_device.hpp>
using namespace boost;
kreutzer1986 engine;

Face::Face()
{
	vertex_num = 0;
	face_num = 0;
	nor_num = 0;
	v.clear();
	f.clear();
	acc_face_area.clear();
	feature_vec.clear();
}

Face::Face(int vertexNum)
{
	vertex_num = 0;
	face_num = 0;
	nor_num = 0;
	face_position = MatrixXd(3, 5);
	vertex = MatrixXd(3, vertexNum);
	nor = MatrixXd(3, vertexNum);
	c = (Curvature *)malloc(vertexNum * sizeof(Curvature));
	memset(c, 0, vertexNum * sizeof(Curvature));
	cluster_set = NULL;
	class_no = NULL;
	feature_point.clear();
	v.clear();
	f.clear();
	class_mean.clear();
	cluster_parent.clear();
	acc_face_area.clear();
	feature_vec.clear();
}

Face::~Face()
{
}

void Face::read_vertex( Vertex  v1)
{
	v.push_back(v1);
}

void Face::readVertex2Matrix(Vertex  v1)
{
	Position p = v1.getPos();
	Position n = v1.getNor();
	Vector3d v(p.getX(), p.getY(), p.getZ());
	Vector3d nor(n.getX(), n.getY(), n.getZ());
	vertex.block<3, 1>(0, vertex_num) = v;
	this->nor.block<3, 1>(0, vertex_num) = nor;
	vertex_num++;
}


void Face::read_face(index3  index)
{
	double area = Position::cal_triangle_area(v[index.idx1].getPos(), v[index.idx2].getPos(), v[index.idx3].getPos());
	if (face_num == 0)
	{	
		acc_face_area.push_back(area);
	}
	else
	{
		double pre_area = acc_face_area[face_num - 1];
		acc_face_area.push_back(area + pre_area);
	}
	f.push_back(index);
	face_num++;
}

void Face::read_feature(double feature)
{
	feature_vec.push_back(feature);
}

void Face::read_right_point(Position point)
{
	feature_point.push_back(point);
}

void Face::cal_curvature()
{
	printf("enter 79\n");
	bool *vis = (bool *)malloc(v.size() * sizeof(bool));
	memset(vis, false, v.size() * sizeof(bool));
	for (int i = 0; i < f.size(); i++)
	{
		int idx[3] = { f[i].idx1, f[i].idx2, f[i].idx3 };
		for (int j = 0; j < 3; j++)
		{		
			int v_index = idx[j];
			Position v_nor = v[v_index].getNor();
			int index[2] = { idx[(j + 1) % 3], idx[(j + 2) % 3] };
			Position v1, v2, v3;
			v1 = v[v_index].getPos(); v2 = v[index[0]].getPos(); v3 = v[index[1]].getPos();
			double k1 = Position::inner_product(v_nor, v2 - v1) / ((v2 - v1).norm());
			double k2 = Position::inner_product(v_nor, v3 - v1) / ((v3 - v1).norm());
			if (!vis[v_index])
			{
				vis[v_index] = true;
				c[v_index].max_curvature = (k1 > k2) ? k1 : k2;
				c[v_index].min_curvature = k1 + k2 - c[v_index].max_curvature;
			}
			else
			{
				if (k1 > c[v_index].max_curvature || k2 > c[v_index].max_curvature)
					c[v_index].max_curvature = k1 > k2 ? k1 : k2;
				if (k1 < c[v_index].min_curvature || k2 < c[v_index].min_curvature)
					c[v_index].min_curvature = k1 < k2 ? k1 : k2;
			}
		}
	}
}

void Face::cal_varOfVertexNNor()
{
	Vector3d pos_mean(vertex.row(0).mean(), vertex.row(1).mean(), vertex.row(2).mean());
	Vector3d nor_mean(nor.row(0).mean(), nor.row(1).mean(), nor.row(2).mean());
	Vector3d pos_var = Vector3d::Zero();
	Vector3d nor_var = Vector3d::Zero();
	for (int i = 0; i < vertex.cols(); i++)
	{
		Vector3d pos_diff_i = vertex.col(i) - pos_mean;
		Vector3d nor_diff_i = nor.col(i) - nor_mean;
		vertex.col(i) -= pos_mean;
		nor.col(i) -= nor_mean;
		Vector3d p = pos_diff_i.array().square();
		Vector3d n = nor_diff_i.array().square();
		pos_var += p;
		nor_var += n;
	}
	pos_var /= vertex.cols();
	nor_var /= vertex.cols();
	/*
	//位置方差
	for (int i = 0; i < 3; i++)
	{
		read_feature(pos_var[i]);
	}
	*/
	/*
	//法向量方差
	for (int i = 0; i < 3; i++)
	{
		read_feature(nor_var[i]);
	}
	*/
	
}


bool inCube(Vector3d p, Vector3d min, Vector3d max)
{
	for (int i = 0; i < 3; i++)
	{
		if (p[i] < min[i] || p[i] > max[i])
			return false;
	}
	return true;
}

void Face::cal_2d_feature_2_3d_feature(float * vb, size_t size, int face_index)
{
	//cal the var of point's nor in a cube(max_x, min_x, max_y, min_y, max_z, min_z) 
	Vector3d min, max;//record the cube
	std::vector<Vector3d> nors_in_cube;//records the points's nors
	Vector3d mean;//records the mean of the nors
	Vector3d var;//records the var
	int ii;
	
	char file_path[50];
	FILE * fp;
	sprintf_s(file_path, "..\\2d_feature\\%d\\feature_from_2d.txt", face_index);
	fopen_s(&fp, file_path, "r");
	if (fp)
	{
		printf("feature_from_2d has already calculated\n");
		for (int i = 0; i < 12; i++)
		{
			double var;
			fscanf_s(fp, "%lf", &var);
			read_feature(var);
		}
	}
	else
	{
		printf("feature_from_2d hasn't  calculated yet\n");
		fopen_s(&fp, file_path, "w");
#define ADD2DFEATURE(FACIAL_FEATURE)\
min = max = FACIAL_FEATURE.col(0); \
nors_in_cube.clear(); \
mean = Vector3d::Zero(); \
var = Vector3d::Zero(); \
ii = 0;\
for (int i = 1; i < FACIAL_FEATURE.cols(); i++)\
{\
Vector3d p = FACIAL_FEATURE.col(i); \
for (int j = 0; j < 3; j++)\
{\
if (min[j] > p[j])\
min[j] = p[j]; \
if (max[j] < p[j])\
	max[j] = p[j]; \
}\
}\
\
for (size_t i = 0; i < size; ++i)\
{\
Vector3d p(vb[3 * i], vb[3 * i + 1], vb[3 * i + 2]); \
if (inCube(p, min, max))\
{\
mean += nor.col(i); \
nors_in_cube.push_back(nor.col(i)); \
}\
}\
\
mean /= nors_in_cube.size(); \
for (auto k = nors_in_cube.begin(); ii < nors_in_cube.size(); ++ii, ++k)\
{\
Vector3d diff = (*k - mean).array().square(); \
var += diff; \
}\
var /= nors_in_cube.size(); \
for (int i = 0; i < 3; i++)\
{\
read_feature(var[i]); \
fprintf_s(fp, "%lf ", var[i]);\
}
		ADD2DFEATURE(eyebrow);//cal the eyebrow feature
		ADD2DFEATURE(nose_curve);
		ADD2DFEATURE(lips);
		Matrix<double, 3, 18> eye_set;
		for (int i = 0; i < 18; i++)
		{
			if (i < 8)
				eye_set.col(i) = eye_curve.col(i);
			else if (i < 16)
				eye_set.col(i) = eyeball.col(i - 8);
			else
				eye_set.col(i) = eyes.col(i - 16);
		}
		ADD2DFEATURE(eye_set);
	}
	fclose(fp);
}

void Face::rotate_by_2d_feature(float * vb, size_t size, float * fvb)
{
	//rotate to make the nose tips on the z axis around the y axis
	double r = sqrt((nose_tip.x() * nose_tip.x()) + (nose_tip.z() * nose_tip.z()));
	double sintheta = nose_tip.x() /  r;
	double costheta = nose_tip.z() / r;
	Matrix3d Rt;

	//around the y axis
	Rt << costheta, 0, -sintheta, 0, 1, 0, sintheta, 0, costheta;

	std::cout << Rt * nose_tip << std::endl;
	face_curve = Rt *  face_curve;
	eyebrow = Rt *  eyebrow;
	eye_curve = Rt *  eye_curve;
	nose_curve = Rt *  nose_curve;
	lips = Rt *  lips;
	nose_tip = Rt *   nose_tip;
	eyeball = Rt *  eyeball;
	eyes = Rt *  eyes;
	ear_tip = Rt *  ear_tip;

	for (int i = 0; i < size; ++i)
	{
		Vector3d p(vb[3 * i], vb[3 * i + 1], vb[3 * i + 2]);
		p = Rt * p;
		for (int j = 0; j < 3; ++j)
		{
			vb[3 * i + j] = p[j];
		}
	}

	for (int i = 0; i < 77; ++i)
	{
		Vector3d p(fvb[3 * i], fvb[3 * i + 1], fvb[3 * i + 2]);
		p = Rt * p;
		for (int j = 0; j < 3; ++j)
		{
			fvb[3 * i + j] = p[j];
		}
	}


	//rotate the make the ear top aligned to eyebrow on y axis

	double y = ear_tip.row(1).mean();
	double z = ear_tip.row(2).mean();

	//in XoY, the ordinary is (x, y)
	Vector3d ordinary(0, y, z);
	
	//shift to ordinary
	{
		for (int i = 0; i < face_curve.cols(); ++i)
		{
			face_curve.col(i) -= ordinary;
		}

		for (int i = 0; i < eyebrow.cols(); ++i)
		{
			eyebrow.col(i) -= ordinary;
		}

		for (int i = 0; i < eye_curve.cols(); ++i)
		{
			eye_curve.col(i) -= ordinary;
		}

		for (int i = 0; i < nose_curve.cols(); ++i)
		{
			nose_curve.col(i) -= ordinary;
		}

		for (int i = 0; i < lips.cols(); ++i)
		{
			lips.col(i) -= ordinary;
		}

		for (int i = 0; i < nose_tip.cols(); ++i)
		{
			nose_tip.col(i) -= ordinary;
		}

		for (int i = 0; i < eyeball.cols(); ++i)
		{
			eyeball.col(i) -= ordinary;
		}

		for (int i = 0; i < eyes.cols(); ++i)
		{
			eyes.col(i) -= ordinary;
		}

		for (int i = 0; i < ear_tip.cols(); ++i)
		{
			ear_tip.col(i) -= ordinary;
		}

		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				vb[3 * i + j] -= ordinary[j];
			}
		}

		for (int i = 0; i < 77; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				fvb[3 * i + j] -= ordinary[j];
			}
		}
	}
	
	
	//cal the mean positon of eyebrow
	Vector3d eyebrow_mean(eyes.row(0).mean(), eyes.row(1).mean(), eyes.row(2).mean());


	//rotate the face around the x axis
	r = sqrt((eyebrow_mean.y() * eyebrow_mean.y()) + (eyebrow_mean.z() * eyebrow_mean.z()));
	sintheta = eyebrow_mean.y() / r;
	costheta = eyebrow_mean.z() / r;
	Rt << 1, 0, 0, 0, costheta, -sintheta, 0, sintheta, costheta;

	std::cout << Rt * eyebrow_mean << std::endl;
	face_curve = Rt *  face_curve;
	eyebrow = Rt *  eyebrow;
	eye_curve = Rt *  eye_curve;
	nose_curve = Rt *  nose_curve;
	lips = Rt *  lips;
	nose_tip = Rt *   nose_tip;
	eyeball = Rt *  eyeball;
	eyes = Rt *  eyes;
	ear_tip = Rt *  ear_tip;
	
	for (int i = 0; i < size; ++i)
	{
		Vector3d p(vb[3 * i], vb[3 * i + 1], vb[3 * i + 2]);
		p = Rt * p;
		for (int j = 0; j < 3; ++j)
		{
			vb[3 * i + j] = p[j];
		}
	}

	for (int i = 0; i < 77; ++i)
	{
		Vector3d p(fvb[3 * i], fvb[3 * i + 1], fvb[3 * i + 2]);
		p = Rt * p;
		for (int j = 0; j < 3; ++j)
		{
			fvb[3 * i + j] = p[j];
		}
	}

	//shift back
	{
		for (int i = 0; i < face_curve.cols(); ++i)
		{
			face_curve.col(i) += ordinary;
		}

		for (int i = 0; i < eyebrow.cols(); ++i)
		{
			eyebrow.col(i) += ordinary;
		}

		for (int i = 0; i < eye_curve.cols(); ++i)
		{
			eye_curve.col(i) += ordinary;
		}

		for (int i = 0; i < nose_curve.cols(); ++i)
		{
			nose_curve.col(i) += ordinary;
		}

		for (int i = 0; i < lips.cols(); ++i)
		{
			lips.col(i) += ordinary;
		}

		for (int i = 0; i < nose_tip.cols(); ++i)
		{
			nose_tip.col(i) += ordinary;
		}

		for (int i = 0; i < eyeball.cols(); ++i)
		{
			eyeball.col(i) += ordinary;
		}

		for (int i = 0; i < eyes.cols(); ++i)
		{
			eyes.col(i) += ordinary;
		}

		for (int i = 0; i < ear_tip.cols(); ++i)
		{
			ear_tip.col(i) += ordinary;
		}

		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				vb[3 * i + j] += ordinary[j];
			}
		}

		for (int i = 0; i < 77; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				fvb[3 * i + j] += ordinary[j];
			}
		}
	}
	ordinary[1] = 30 + nose_tip[1];
	{
		for (int i = 0; i < face_curve.cols(); ++i)
		{
			face_curve.col(i) -= ordinary;
		}

		for (int i = 0; i < eyebrow.cols(); ++i)
		{
			eyebrow.col(i) -= ordinary;
		}

		for (int i = 0; i < eye_curve.cols(); ++i)
		{
			eye_curve.col(i) -= ordinary;
		}

		for (int i = 0; i < nose_curve.cols(); ++i)
		{
			nose_curve.col(i) -= ordinary;
		}

		for (int i = 0; i < lips.cols(); ++i)
		{
			lips.col(i) -= ordinary;
		}

		for (int i = 0; i < nose_tip.cols(); ++i)
		{
			nose_tip.col(i) -= ordinary;
		}

		for (int i = 0; i < eyeball.cols(); ++i)
		{
			eyeball.col(i) -= ordinary;
		}

		for (int i = 0; i < eyes.cols(); ++i)
		{
			eyes.col(i) -= ordinary;
		}

		for (int i = 0; i < ear_tip.cols(); ++i)
		{
			ear_tip.col(i) -= ordinary;
		}

		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				vb[3 * i + j] -= ordinary[j];
			}
		}

		for (int i = 0; i < 77; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				fvb[3 * i + j] -= ordinary[j];
			}
		}
	}
}

void Face::zeroOneNormalize(double max, double min)
{
	MatrixXd minVertex = MatrixXd::Ones(3, vertex.cols()) * min;
	vertex = (vertex - minVertex) / max;
}

void Face::feature_normalize(std::vector<double> max, std::vector<double> min)
{
	for (int i = 0; i < feature_vec.size(); i++)
	{
		if (i > max.size() || i > min.size() || i > feature_vec.size())
		{
			printf("wrong\n");
			system("PAUSE");
			return;
		}
		if (max[i] == min[i])
			feature_vec[i] = 0.0;
		else
			feature_vec[i] = (feature_vec[i] - min[i]) / (max[i] - min[i]);
	}
}

double Face::feature_L2_norm() const
{
	double sum = 0.0;
	for (int i = 0; i < feature_vec.size(); i++)
	{
		sum += feature_vec[i] * feature_vec[i];
	}
	return sqrt(sum);
}

void Face::DBScan(double eps, int minPts, int face_index)
{
	if (cluster_set)
		free(cluster_set);
	if (class_no)
		free(class_no);
	if (cluster_parent.size() > 0)
		cluster_parent.clear();
	if (class_mean.size() > 0)
		class_mean.clear();
	cluster_set = (int *)malloc(sizeof(int) * feature_point.size());
	class_no = (int *)malloc(sizeof(int) * feature_point.size());
	memset(class_no, 0, sizeof(int) * feature_point.size());//0表示为噪音点
	char file_path[50];
	FILE * fp;
	sprintf_s(file_path, "..\\2d_feature\\%d\\DBScan_result.txt", face_index);
	fopen_s(&fp, file_path, "r");
	if (fp)
	{
		printf("DBScan already done\n");
		double point[3];
		Position p;
#define READ_FACIAL_FEATURE(SPOT) \
		for (int i = 0; i < 3; i++)\
		{\
			fscanf_s(fp, "%lf", &point[i]);\
		}\
		p = Position(point[0], point[1], point[2]);\
		SPOT = Vertex(p, Position())

		printf("read facial feature\n");
		READ_FACIAL_FEATURE(left_eye);
		READ_FACIAL_FEATURE(right_eye);
		READ_FACIAL_FEATURE(left_ear);
		READ_FACIAL_FEATURE(right_ear);
		READ_FACIAL_FEATURE(nose);

		printf("read class_mean\n");
		int class_size;
		fscanf_s(fp, "%d", &class_size);
		for (int i = 0; i < class_size; i++)
		{
			double point[3];
			for (int j = 0; j < 3; j++)
			{
				fscanf_s(fp, "%lf", &point[j]);
			}
			class_mean.push_back(Position(point[0], point[1], point[2]));
		}


		printf("read class no and cluser set\n");
		for (int i = 0; i < feature_point.size(); i++)
		{
			fscanf_s(fp, "%d", &class_no[i]);
			fscanf_s(fp, "%d", &cluster_set[i]);
		}

	}
	else
	{
		fopen_s(&fp, file_path, "w");
		printf("file not exist\n糖糖把%d点记住\n", feature_point.size());
		for (int i = 0; i < feature_point.size(); i++)
		{
			cluster_set[i] = i;//point self
		}
		bool * vis = (bool *)malloc(sizeof(bool) * feature_point.size());
		memset(vis, false, sizeof(bool) * feature_point.size());
		int cnt = 0;
		for (int i = 0; i < feature_point.size(); ++i)
		{
			if (!vis[i])
			{
				int class_cnt = 0;
				Position sum_class = Position();
				vis[i] = true;
				Position point = feature_point[i];
				std::vector<int> neighbors_index;
				neighbors_index.clear();
				for (int j = 0; j < feature_point.size(); ++j)
				{
					if (i == j) continue;
					if ((point - feature_point[j]).norm() < eps)
					{
						neighbors_index.push_back(j);
					}
				}
				if (neighbors_index.size() < minPts)//拓展点
				{
					class_no[i] = 0;//为噪音

				}
				else
				{
					printf("我觉得%d点还OK\n", i);
					//system("PAUSE");
					++class_cnt;
					cnt++;
					class_no[i] = 1;
					sum_class += feature_point[i];
					cluster_parent.push_back(i);
					for (int ii = 0; ii < neighbors_index.size(); ++ii)
					{
						//printf("!!!!!!!!!\n");
						if (!vis[neighbors_index[ii]])
						{
							vis[neighbors_index[ii]] = true;
							Position point2 = feature_point[neighbors_index[ii]];
							std::vector<int> neighbors_index2;
							neighbors_index2.clear();
							for (int jj = 0; jj < feature_point.size(); ++jj)
							{
								if (ii == jj) continue;
								if ((point2 - feature_point[jj]).norm() < eps)
								{
									neighbors_index2.push_back(jj);
								}
							}
							if (neighbors_index2.size() > minPts)
							{
								//printf("%dpoint也还可以\n", neighbors_index[ii]);
								neighbors_index.insert(neighbors_index.end(), neighbors_index2.begin(), neighbors_index2.end());
								//printf("crashed\n");
							}
						}
						if (class_no[neighbors_index[ii]] != 1)
						{
							cnt++;
							sum_class += feature_point[neighbors_index[ii]];
							++class_cnt;
							class_no[neighbors_index[ii]] = 1;
							cluster_set[neighbors_index[ii]] = i;//指向父亲
						}
						//printf("end loop\n");
					}
					sum_class.mean(class_cnt);
					Position var = Position();
					for (int ii = 0; ii < neighbors_index.size(); ++ii)
					{
						Position point2 = feature_point[neighbors_index[ii]];
						var += (point2 - sum_class) * (point2 - sum_class);
					}
					var.mean((int)(neighbors_index.size()));
					if (sum_class.getY() > -10 && sum_class.getY() < 20 && sum_class.getZ() > 85 && fabs(sum_class.getX()) < 40)
					{
						if (sum_class.getX() < 0)
						{
							left_eye = Vertex(sum_class, var);
							printf("found left eye at %lf %lf %lf\n", sum_class.getX(), sum_class.getY(), sum_class.getZ());
							//system("PAUSE");
						}
						else
						{
							right_eye = Vertex(sum_class, var);
							printf("found right eye at %lf %lf %lf\n", sum_class.getX(), sum_class.getY(), sum_class.getZ());
							//system("PAUSE");
						}
					}
					else if (fabs(sum_class.getX()) < 10 && sum_class.getY() < -30 && sum_class.getY() > -60 && sum_class.getZ() > 70)
					{
						//更靠近下面的
						if (nose.getPos() == Position() || (nose.getPos() != Position() && (nose.getPos().getY() > sum_class.getY())))
							nose = Vertex(sum_class, var);
						printf("found nose  at %lf %lf %lf\n", sum_class.getX(), sum_class.getY(), sum_class.getZ());
						//system("PAUSE");
					}
					else if (fabs(sum_class.getX()) > 80 && sum_class.getY() < 10 && sum_class.getY() > -30 && fabs(sum_class.getZ()) < 40)
					{
						if (sum_class.getX() < 0)
						{
							if (left_ear.getPos() == Position() || (left_ear.getPos() != Position() && sum_class.getX() < left_ear.getPos().getX()))//更靠近两端的
							{
								left_ear = Vertex(sum_class, var);
								printf("found left ear at %lf %lf %lf\n", sum_class.getX(), sum_class.getY(), sum_class.getZ());
								//system("PAUSE");
							}
						}
						else//更靠近两端的
						{
							if (right_ear.getPos() == Position() || (right_ear.getPos() != Position() && sum_class.getX() > right_ear.getPos().getX()))
							{
								right_ear = Vertex(sum_class, var);
								printf("found right ear at %lf %lf %lf\n", sum_class.getX(), sum_class.getY(), sum_class.getZ());
								//system("PAUSE");
							}

						}
					}
					class_mean.push_back(sum_class);
				}
			}
		}
		printf("哇, 一共有%d点\n", cnt);

		//record the facial feature
		fprintf_s(fp, "%lf %lf %lf\n", left_eye.getPos().getX(), left_eye.getPos().getY(), left_eye.getPos().getZ());
		fprintf_s(fp, "%lf %lf %lf\n", right_eye.getPos().getX(), right_eye.getPos().getY(), right_eye.getPos().getZ());
		fprintf_s(fp, "%lf %lf %lf\n", left_ear.getPos().getX(), left_ear.getPos().getY(), left_ear.getPos().getZ());
		fprintf_s(fp, "%lf %lf %lf\n", right_ear.getPos().getX(), right_ear.getPos().getY(), right_ear.getPos().getZ());
		fprintf_s(fp, "%lf %lf %lf\n", nose.getPos().getX(), nose.getPos().getY(), nose.getPos().getZ());


		//record the class_mean
		fprintf_s(fp, "%d\n", (int)(class_mean.size()));
		for (size_t i = 0; i < class_mean.size(); ++i)
		{
			Position p = class_mean[i];
			fprintf_s(fp, "%lf %lf %lf\n", p.getX(), p.getY(), p.getZ());
		}

		int *cluster_cnt = (int *)malloc(sizeof(int) * feature_point.size());
		memset(cluster_cnt, 0, sizeof(int) * feature_point.size());
		for (int i = 0; i < feature_point.size(); i++)
		{
			//record the farther
			fprintf_s(fp, "%d %d\n", class_no[i], cluster_set[i]);

			if (class_no[i])
				cluster_cnt[cluster_set[i]]++;
		}
		printf("class_size is %d\n", cluster_parent.size());
		for (int i = 0; i < cluster_parent.size(); i++)
		{
			printf("set%d has %d children\n", cluster_parent[i], cluster_cnt[cluster_parent[i]]);
		}
	}
	fclose(fp);
}

void Face::feature_extract(float *vb, size_t size, double r, int face_index)
{
	if (left_eye.getPos() == Position())
	{
		if (right_eye.getPos() == Position());
			//left_eye = Vertex(default_left_eye, left_eye.getNor());
		else
			left_eye = Vertex(Position(-right_eye.getPos().getX(), right_eye.getPos().getY(), right_eye.getPos().getZ()), right_eye.getNor());
	}
	if (right_eye.getPos() == Position())
	{
		if (left_eye.getPos() == Position());
			//right_eye = Vertex(default_right_eye, right_eye.getNor());
		else
			right_eye = Vertex(Position(-left_eye.getPos().getX(), left_eye.getPos().getY(), left_eye.getPos().getZ()), left_eye.getNor());
	}
	if (left_ear.getPos() == Position())
	{
		if (right_ear.getPos() == Position());
			//left_ear = Vertex(default_left_ear, left_ear.getNor());
		else
			left_ear = Vertex(Position(-right_ear.getPos().getX(), right_ear.getPos().getY(), right_ear.getPos().getZ()), right_ear.getNor());
	}
	if (right_ear.getPos() == Position())
	{
		if (left_ear.getPos() == Position());
			//right_ear = Vertex(default_right_ear, right_ear.getNor());
		else
			right_ear = Vertex(Position(-left_ear.getPos().getX(), left_ear.getPos().getY(), left_ear.getPos().getZ()), left_ear.getNor());
	}
	if (nose.getPos() == Position());
		//nose = Vertex(default_nose, nose.getNor());

	char file_path[50];
	FILE * fp;
	sprintf_s(file_path, "..\\2d_feature\\%d\\DBScan_feature.txt", face_index);
	fopen_s(&fp, file_path, "r");
	if (fp)
	{
		for (int i = 0; i < 15; i++)
		{
			double var;
			fscanf_s(fp, "%lf", &var);
			read_feature(var);
		}
	}
	else
	{
		fopen_s(&fp, file_path, "w");
		Position mean_position[5] = { left_eye.getPos(), right_eye.getPos(), left_ear.getPos(), right_ear.getPos(), nose.getPos() };
		for (int i = 0; i < 5; i++)
		{
			std::vector <Position> temp;
			temp.clear();
			printf("%d feature at %lf, %lf, %lf\n", i, mean_position[i].getX(), mean_position[i].getY(), mean_position[i].getZ());
			if (mean_position[i] == Position())
			{
				for (int k = 0; k < 3; k++)
				{
					read_feature(0);// read 0 as feature
				}
				continue;
			}
			for (size_t j = 0; j < size; ++j)
			{
				Position p = Position(vb[3 * j], vb[3 * j + 1], vb[3 * j + 2]);
				double distance = (p - mean_position[i]).norm();

				if (distance <= r)
				{
					temp.push_back(p);
				}
			}
			if (temp.size() == 0)
			{
				for (int k = 0; k < 3; k++)
				{
					read_feature(0);// read var as feature
				}
				continue;
			}
			MatrixXd m(3, temp.size());
			int cnt = 0;
			for (auto k = temp.begin(); cnt < temp.size(); cnt++, ++k)
			{
				Position p = *k;
				m.block<3, 1>(0, cnt) = Vector3d(p.getX(), p.getY(), p.getZ());
			}
			Vector3d var = Vector3d::Zero();
			Vector3d mean(m.row(0).mean(), m.row(1).mean(), m.row(2).mean());
			face_position.block<3, 1>(0, i) = mean;
			for (size_t j = 0; j < temp.size(); ++j)
			{
				Vector3d diff_i = m.col(j) - mean;
				Vector3d p = diff_i.array().square();
				var += p;
			}
			var /= m.cols();
			for (int k = 0; k < 3; k++)
			{
				fprintf_s(fp, "%lf ", var[k]);
				read_feature(var[k]);// read var as feature
			}
		}
	}
}

double Face::euler_dis(const Face & f1, const Face & f2)
{
	int min_feature_num = (f1.feature_vec.size() < f2.feature_vec.size()) ? f1.feature_vec.size() : f2.feature_vec.size();
	printf("min_feature_num = %d\n", min_feature_num);
	double sum = 0;
	for (int i = 0; i < min_feature_num; i++)
	{
		double L1_dis_i = f1.feature_vec[i] - f2.feature_vec[i];
		sum += L1_dis_i * L1_dis_i;
	}
	return sqrt(sum);
}

double Face::L1_dis(const Face & f1, const Face & f2)
{
	int min_feature_num = (f1.feature_vec.size() < f2.feature_vec.size()) ? f1.feature_vec.size() : f2.feature_vec.size();
	printf("min_feature_num = %d\n", min_feature_num);
	double sum = 0;
	for (int i = 0; i < min_feature_num; i++)
	{
		double L1_dis_i = f1.feature_vec[i] - f2.feature_vec[i];
		sum += fabs(L1_dis_i);
	}
	return sqrt(sum);
}

double Face::inner_product(const Face & f1, const Face & f2)
{
	int min_feature_num = (f1.feature_vec.size() < f2.feature_vec.size()) ? f1.feature_vec.size() : f2.feature_vec.size();
	double sum = 0;
	for (int i = 0; i < min_feature_num; i++)
	{
		double inner_product_i = f1.feature_vec[i] * f2.feature_vec[i];
		sum += inner_product_i;
	}
	return sum;
}

double Face::cos_dis(const Face & f1, const Face & f2)
{
	printf("min_feature_num = %d\n", f2.feature_vec.size());
	return inner_product(f1, f2) / (f1.feature_L2_norm() * f2.feature_L2_norm());
}

void Face::add_D2_feature(int dime_num)
{
	int count = dime_num;
	while (dime_num--)
	{
		Position randomPosition[2];
		//select 2 point randomly by the area
		for (int i = 0; i < 2; i++)
		{
			uniform_real<double> real(0, acc_face_area.back());

			double area = real(engine);
			int mid = acc_face_area.size() / 2;
			int left = 0;
			int right = acc_face_area.size();
			//printf("%d faces, area = %lf, max_area = %lf\n", acc_face_area.size(), area, acc_face_area.back());
			while (1)
			{
				//printf("mid = %d, acc_face_area[mid] = %lf, acc_face_area[mid + 1] = %lf\n", mid, acc_face_area[mid], acc_face_area[mid + 1]);
				//system("PAUSE");
				if (acc_face_area[mid] >= area && mid == 0)//the first triangle
					break;
				if (acc_face_area[mid] <= area && mid == acc_face_area.size() - 1)//the last triangle
					break;
				if (acc_face_area[mid] <= area && acc_face_area[mid + 1] >= area)
					break;
				if (acc_face_area[mid] > area)
				{
					right = mid;
					mid = (left + right) / 2;
				}
				else if (acc_face_area[mid + 1] < area)
				{
					left = mid;
					mid = (left + right) / 2;
				}
			}//binsearch to find the triangle randomly by the area
			;
			Position A(vertex.col(f[mid].idx1)[0], vertex.col(f[mid].idx1)[1], vertex.col(f[mid].idx1)[2]);
			Position B(vertex.col(f[mid].idx2)[0], vertex.col(f[mid].idx2)[1], vertex.col(f[mid].idx2)[2]);
			Position C(vertex.col(f[mid].idx3)[0], vertex.col(f[mid].idx3)[1], vertex.col(f[mid].idx3)[2]);
			uniform_real<double> r1_gen(0, 1);
			uniform_real<double> r2_gen(0, 1);
			double r1 = r1_gen(engine);
			double r2 = r2_gen(engine);
			randomPosition[i] = A * (1.0 - sqrt(r1)) + B * (sqrt(r1) * (1 - r2)) + C * (r2 * sqrt(r1));
		}
		double euler_dis = Position::distance(randomPosition[0], randomPosition[1]);
		printf("(%lf, %lf, %lf) & (%lf, %lf, %lf) euler_dis%d = %lf\n", 
			randomPosition[0].getX(), randomPosition[0].getY(), randomPosition[0].getZ(),
			randomPosition[1].getX(), randomPosition[1].getY(), randomPosition[1].getZ(),
			count - dime_num, euler_dis);
		//system("PAUSE");
		feature_vec.push_back(euler_dis / (count));
	}
}


double Face::ICP_rotate(Face & basedFace, int rotate_nor, double tolerance)//不总是可以正确求解
{
	int length = 500;
	MatrixXd src(3, length);
	MatrixXd dst(3, length);
	VectorXd vis_src = VectorXd::Zero(vertex.cols());
	VectorXd vis_dst = VectorXd::Zero(basedFace.vertex.cols());
	printf("line243\n");
	for (int i = 0; i < length; i++)//generate a random vertex list
	{
		srand((unsigned)time(NULL));
		int index1 = rand() % (vertex.cols());
		while(vis_src[index1] == 1)
			index1 = rand() % (vertex.cols());
		vis_src[index1] = 1;
		int index2 = rand() % basedFace.vertex.cols();
		while (vis_dst[index2] == 1)
			index2 = rand() % (basedFace.vertex.cols());
		vis_dst[index2] = 1;
		//printf("index1 = %d, index2 = %d\n", index1, index2);
		if (!rotate_nor)
		{
			src.block<3, 1>(0, i) = vertex.col(index1);
			dst.block<3, 1>(0, i) = basedFace.vertex.col(index2);
		}
		else
		{
			src.block<3, 1>(0, i) = nor.col(index1);
			dst.block<3, 1>(0, i) = basedFace.nor.col(index2);
		}
	}
	printf("line263\n");
	double pre_error = 0;
	while (1)
	{
		VectorXd distances(length);
		int indecies[500];
		for (int i = 0; i < length; i++)
		{
			double min_dis = 999;//由于已经归一化，点的距离不超过根号3
			for (int j = 0; j < length; j++)
			{
				double dis = (src.col(i) - dst.col(i)).norm();
				if (dis < min_dis)
				{
					indecies[i] = j;
					min_dis = dis;
					distances[i] = dis;
				}
			}
		}
		double error = distances.mean();
		printf("error is %lf\n", error);
		//system("PAUSE");
		//cal centroid
		Vector3d centroid_src(src.row(0).mean(), src.row(1).mean(), src.row(2).mean());
		Vector3d centroid_dst(dst.row(0).mean(), dst.row(1).mean(), dst.row(2).mean());

		for (int i = 0; i < length; i++)
		{
			src.col(i) -= centroid_src;
			dst.col(i) -= centroid_dst;
		}//shift to centroid

		Matrix3d W = src * dst.transpose();
		JacobiSVD<MatrixXd> svd(W, ComputeThinU | ComputeThinV);
		Matrix3d R = svd.matrixU() * svd.matrixV().transpose();
		src = R * src;
		if (fabs(error - pre_error) < 0.001)
		{
			printf("pre_error = %lf, error = %lf, %lf\n ", pre_error, error, fabs(error - pre_error));
			break;
		}
			
		pre_error = error;
	}
	return pre_error;
}

void Face::print_feature()
{
	printf("the  %d features are:\n", feature_vec.size());
	for (int i = 0; i < feature_vec.size(); i++)
	{
		printf("%lf  ", feature_vec[i]);
	}
	printf("\n");
}
