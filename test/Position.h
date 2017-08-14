#pragma once
#include <Eigen\Dense>
using namespace Eigen;
class Position
{
public:
	Position();
	Position(const Position &p);
	Position(double x, double y, double z);
	static double distance(Position &p1, Position &p2);
	static double cos(Position &p1, Position &p2);
	static double inner_product(Position &p1, Position &p2);
	static Position cal_nor(Position &p1, Position &p2);
	static double cal_triangle_area( Position &p1,  Position &p2,  Position &p3);
	Position operator + (const Position &p);
	Position operator - (const Position &p);
	Position operator * (const Position &p);
	Position operator * (const double num);
	Position operator - ();
	Position & operator = (const Position &p);
	Position & operator += (const Position &p);
	bool operator != (const Position &p) const;
	bool operator == (const Position &p) const;
	void mean(int num);
	void mean(double num);
	void normalize();
	void print_pos();
	double norm();
	double getX();
	double getY();
	double getZ();
	~Position();

private:
	double x;
	double y;
	double z;
};

