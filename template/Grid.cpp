#include "precomp.h"
#include "Grid.h"

Grid::Grid( int columns, int rows )
{
	( *this ).columns = columns;
	( *this ).rows = rows;

	( *this ).cells = new GridCell[columns * rows];
}

Grid::~Grid()
{
	delete cells;
}

void Grid::ConstructGrid( const vector<Boid> &b )
{
	ComputeBoundingBox( b );
	StoreInCells( b );
	ReorderData();
}

void Grid::ConstructGrid( Vec3 min, Vec3 max, const vector<Boid> &b )
{
	( *this ).min = min;
	( *this ).max = max;

	StoreInCells( b );
	ReorderData();
}

void Grid::QueryGrid( const Boid &b, const int r, vector<Boid> &out )
{
	// todo: query it
}

void Grid::QueryGrid( const Vec3 &p, const int r, vector<Boid> &out )
{
	// todo: same as the other query.
}

void Grid::ComputeBoundingBox( const vector<Boid> &b )
{
	// todo: compute it
}

void Grid::StoreInCells( const vector<Boid> &b )
{
	// todo: store 'dem.
}

void Grid::ReorderData()
{
	// todo: reorder the data to make it more friendly.
}

void Grid::ComputeBoundingBoxOfCell( int x, int y, Vec3 *min, Vec3 *max )
{
	// todo: compute 'dem boxes.
}

GridCell::GridCell()
{
	boids = new Boid[NUMBER_OF_ELEMENTS_IN_CELL];
}
