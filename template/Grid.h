#pragma once

#define NUMBER_OF_ELEMENTS_IN_CELL 50

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

class Grid
{

  public:
	// represents the number of columns / rows.
	int xcells, ycells, zcells;

	// represents the boundingbox.
	Vec3 min, max;

	// regular constructor / deconstructor
	Grid( int xcells, int ycells, int zcells );
	~Grid();

	inline int CalculateGridCellIndex( int celX, int celY, int celZ );
	inline bool CheckInsideGrid( int celX, int celY, int celZ );


	void ComputeGridIndex( const Boid &b, int &celX, int &celY, int &celZ );

	// constructs the grid. If no min / max
	// is provided the bounding box will be
	// computed dynamically.
	void ConstructGrid( const vector<Boid> &b );
	void ConstructGrid( Vec3 min, Vec3 max, const vector<Boid> &b);

	// queries the grid, stores the result in
	// the out vector. Take note: reuse the vector.
	void QueryGrid( const Boid &b, const int r, vector<NearbyBoid> &out, float PerceptionRadius, float BlindspotAngleDeg, int celX, int celY, int celZ );

  private:
	GridCell *cells;

	// computes the bounding box, dynamically.
	void ComputeBoundingBox( const vector<Boid> &b );

	// stores all boids in their corresponding cells.
	void StoreInCells( const vector<Boid> &b );

	// reorders the data to make it more cache-friendly.
	void ReorderData();

};