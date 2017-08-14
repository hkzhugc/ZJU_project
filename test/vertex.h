#pragma once
#include "Position.h"
#include <math.h>
class Vertex
{
public:
	Vertex();
	Vertex(float x, float y, float z, float x_nor, float y_nor, float z_nor);
	Vertex(const Position &p1, const Position &p2);
	~Vertex();
	void normalizeNor();
	static Vertex  diffV(Vertex &a, Vertex &b);
	static Vertex  cal_nor(Vertex &a, Vertex &b);
	static double eurler_dis_of_feature(Vertex &a, Vertex &b);
	static Vertex  Combine(Vertex &v1, Vertex &v2);
	static double cos_dis(Vertex &a, Vertex &b);

	Position pos;
	Position nor;

	double x;
	double y;
	double z;
	double x_nor;
	double y_nor;
	double z_nor;
};

Vertex::Vertex()
{
	x = 0; y = 0; z = 0; x_nor = 0; y_nor = 0; z_nor = 0;
}

Vertex::Vertex(float x, float y, float z, float x_nor, float y_nor, float z_nor)
{
	this->x = x; this->y = y; this->z = z; this->x_nor = x_nor; this->y_nor = y_nor; this->z_nor = z_nor;
}

Vertex::Vertex(const Position & p1, const Position & p2)
{
	pos = p1;
	nor = p2;
}


Vertex::~Vertex()
{
}

inline void Vertex::normalizeNor()
{
	double r = sqrt(x_nor*x_nor + y_nor*y_nor + z_nor*z_nor);
	if (fabs(r - 0.0) < 1e-8) return;
	x_nor /= r;
	y_nor /= r;
	z_nor /= r;
}

inline Vertex  Vertex::diffV(Vertex & a, Vertex &b)
{
	Vertex v = Vertex(a.x - b.x, a.y - b.y, a.z - b.z, a.x_nor - b.x_nor, a.y_nor - b.y_nor, a.z_nor - b.z_nor);
	return v;
}

inline Vertex  Vertex::cal_nor(Vertex & a, Vertex & b)
{
	Vertex v = Vertex(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x, 0, 0, 0);
	double r = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	Vertex v1 =  Vertex(v.x / r, v.y / r, v.z / r, 0, 0, 0);
	return v1;
}

inline double Vertex::eurler_dis_of_feature(Vertex & a, Vertex & b)
{
	Vertex v = diffV(a, b);
	double r = sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.x_nor*v.x_nor + v.y_nor*v.y_nor + v.z_nor*v.z_nor);
	return r;
}

inline Vertex  Vertex::Combine(Vertex & v1, Vertex & v2)
{
	Vertex v = Vertex(v1.x_nor, v1.y_nor, v1.z_nor, v2.x_nor, v2.y_nor, v2.z_nor);
	return v;
}

inline double Vertex::cos_dis(Vertex & a, Vertex & b)
{
	double inner_product = a.x * b.x + a.y * b.y + a.z * b.z + a.x_nor * b.x_nor + a.y_nor * b.y_nor + a.z_nor * b.z_nor;
	double eurler_dis = eurler_dis_of_feature(a, Vertex()) * eurler_dis_of_feature(b, Vertex());
	return inner_product / eurler_dis;
}

