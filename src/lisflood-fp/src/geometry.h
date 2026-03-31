#pragma once
#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "lisflood.h"

typedef struct Geometry
{
	int xsz;
	int ysz;
	NUMERIC_TYPE blx;
	NUMERIC_TYPE bly;
	NUMERIC_TYPE tly;
	NUMERIC_TYPE dx;
	NUMERIC_TYPE dy;
} Geometry;

#endif // GEOMETRY_H
