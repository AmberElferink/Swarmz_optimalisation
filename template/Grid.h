#pragma once

#define NUMBER_OF_ELEMENTS_IN_CELL 50

class Grid
{

  public:
	// represents the number of columns / rows.
	int columns, rows;

	// represents the boundingbox.
	Vec3 min, max;

	// regular constructor / deconstructor
	Grid( int columns, int rows );
	~Grid();

	// constructs the grid. If no min / max
	// is provided the bounding box will be
	// computed dynamically.
	void ConstructGrid( const vector<Boid> &b );
	void ConstructGrid( Vec3 min, Vec3 max, const vector<Boid> &b, Boid *buffer );

	// queries the grid, stores the result in
	// the out vector. Take note: reuse the vector.
	void QueryGrid( const Boid &b, const int r, vector<Boid> &out );
	void QueryGrid( const Vec3 &p, const int r, vector<Boid> &out );

  private:
	GridCell *cells;

	// computes the bounding box, dynamically.
	void ComputeBoundingBox( const vector<Boid> &b );

	// stores all boids in their corresponding cells.
	void StoreInCells( const vector<Boid> &b );

	// reorders the data to make it more cache-friendly.
	void ReorderData();

	// given that we have the min / max values and
	// the number of columns and row: compute a
	// bounding box.
	void ComputeBoundingBoxOfCell( int x, int y, Vec3 *min, Vec3 *max );
};

struct GridCell
{
  public:
	// the default constructor :))
	GridCell();

	// the buffer for the boids in this cell.
	Boid boids[NUMBER_OF_ELEMENTS_IN_CELL];

	// the number of boids in this cell.
	int count;

	// Adds the given boid to this cell. If
	// the cell is full, a random boid is
	// evicted.
	void AddBoid( const Boid &b );

	// 'Clears'  this 
	void Clear();
};