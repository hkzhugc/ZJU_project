#pragma once
#include "Position.h"
class Vertex
{
public:
	Vertex()
	{
		pos = Position();
		nor = Position();
		total_degree = 0;
	}
	Vertex(const Position &pos, const Position &nor)
	{
		this->pos = pos; this->nor = nor;
		total_degree = 0;
	}
	Vertex(const Vertex &v)
	{
		pos = v.pos;
		nor = v.nor;
		total_degree = v.total_degree;
	}
	void addNor(Position &nor, const double degree);
	void mean(int num);
	void normalize();
	Position getPos();
	Position getNor();
	Vertex & operator = (const Vertex &v);
	Vertex::~Vertex()
	{
	}

private:
	Position pos;
	Position nor;
	double total_degree;
};

