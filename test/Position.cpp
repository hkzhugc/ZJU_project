#include "Position.h"
#include <math.h>
#include <cstdio>

Position::Position()
{
	x = 0; y = 0; z = 0;
}

Position::Position(const Position & p)
{
	x = p.x; y = p.y; z = p.z;
}

Position::Position(double x, double y, double z)
{
	this->x = x; this->y = y; this->z = z;
}

Position::~Position()
{
}

double Position::distance(Position &p1, Position & p2)
{
	return (p1 - p2).norm();
}

double Position::cos(Position & p1, Position & p2)
{
	return inner_product(p1, p2) / (p1.norm() * p2.norm());
}

double Position::inner_product(Position & p1, Position & p2)
{
	return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

Position Position::cal_nor(Position & a, Position & b)
{
	Position v = Position(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
	double r = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	if (fabs(r - 0) < 1e-15)
		return Position();
	return Position(v.x / r, v.y / r, v.z / r);
}

double Position::cal_triangle_area(Position & p1, Position & p2, Position & p3)
{
	double a = (p1 - p2).norm();
	double b = (p1 - p3).norm();
	double c = (p2 - p3).norm();
	double p = (a + b + c) / 2.0;
	return sqrt(p*(p - a)*(p - b)*(p - c));
}

Position Position::operator+(const Position & p)
{
	return Position(this->x + p.x, this->y + p.y, this->z + p.z);
}

Position Position::operator-(const Position & p)
{
	return Position(this->x - p.x, this->y - p.y, this->z - p.z);
}

Position Position::operator*(const Position & p)
{
	return Position(x * p.x, y * p.y, z * p.z);
}

Position Position::operator*(const double num)
{
	return Position(x * num, y * num, z * num);
}

Position Position::operator-()
{
	return Position(-x, -y, -z);
}

Position & Position::operator=(const Position & p)
{
	x = p.x; y = p.y; z = p.z;
	return *this;
}

Position & Position::operator+=(const Position & p)
{
	*this = *this + p;
	return *this;
}

bool Position::operator!=(const Position & p) const
{
	return !(x == p.x && y == p.y && z == p.z);
}

bool Position::operator==(const Position & p) const
{
	return (x == p.x && y == p.y && z == p.z);
}

void Position::mean(int num)
{
	x /= num; y /= num; z /= num;
}

void Position::mean(double num)
{
	if (fabs(num - 0) < 1e-8) return;
	x /= num; y /= num; z /= num;
}

void Position::normalize()
{
	double r = norm();
	if (fabs(r - 0) < 1e-8) return;
	mean(r);
}

void Position::print_pos()
{
	printf("%lf, %lf, %lf", x, y, z);
}

double Position::norm()
{
	return sqrt(x * x + y * y + z * z);
}

double Position::getX()
{
	return x;
}

double Position::getY()
{
	return y;
}

double Position::getZ()
{
	return z;
}
