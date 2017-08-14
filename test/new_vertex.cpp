#include "new_vertex.h"
void Vertex::addNor(Position & nor, const double degree)
{
	this->nor += (nor * Position(degree, degree, degree));//add weighted  vertex nor
	total_degree += degree;
}

void Vertex::mean(int num)
{
	pos.mean(num);
}

void Vertex::normalize()
{
	nor.mean(total_degree);
	nor.normalize();
}

Position Vertex::getPos()
{
	return pos;
}

Position Vertex::getNor()
{
	return nor;
}

Vertex & Vertex::operator=(const Vertex & v)
{
	pos = v.pos;
	nor = v.nor;
	total_degree = v.total_degree;
	return *this;
}
