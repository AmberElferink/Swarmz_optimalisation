
#include "precomp.h"

#include "Grid.h"
#include "precomp.h"

Grid::Grid( int nx, int ny, int nz )
{
	// receive 'dem dimensions.
	( *this ).nx = nx;
	( *this ).ny = ny;
	( *this ).nz = nz;

	( *this ).cells = new GridCell[nx * ny * nz];
}

Grid::~Grid()
{
	// huuurrrrr!
	delete cells;
}

inline int Grid::CalculateGridCellIndex( int ix, int iy, int iz )
{
	// compute the right index. Take note: this function
	// is inline. It will be placed in the code.
	return ix + iy * nx + iz * ( nx * ny );
}

inline bool Grid::CheckInsideGrid( int ix, int iy, int iz )
{
	// check for all dimensions whether it's
	// within the grid.
	if ( ix < 0 || ix >= nx )
		return false;
	if ( iy < 0 || iy >= ny )
		return false;
	if ( iz < 0 || iz >= nz )
		return false;

	return true;
}

void sw::Grid::ClearGrid()
{
	// compute the total number of cells.
	int l = nx * ny * nz;

	// for every cell, clear it out.
	for ( int j = 0; j < l; j++ )
		cells[j].Clear();
}

void Grid::ConstructGrid( const vector<Boid> &b )
{
	// clear out the grid of any
	// previous usage
	ClearGrid();

	// compute the bounding box
	ComputeBoundingBox( b );

	// store the data
	StoreInCells( b );
}

void Grid::QueryGrid( const Boid &b, const int r, vector<NearbyBoid> &out, float PerceptionRadius, float BlindspotAngleDeg, int celX, int celY, int celZ )
{
	// if the location is not inside
	// the grid, skip it.
	if ( !CheckInsideGrid( celX, celY, celZ ) )
		return;

	// retrieve the cell
	GridCell gridCell = cells[CalculateGridCellIndex( celX, celY, celZ )];

	// do cell computations, this is copied (for now)
	// - fix expensive operations
	// - apply if statements sooner (e.g., compute 
	// distance -> check, compute angle -> check, 
	// etc), allows for early-opt out.
	for ( int i = 0; i < gridCell.count; i++ )
	{
		//compute distance between b and test
		Boid target = gridCell.boids[i];
		const Vec3 &p1 = b.Position;
		const Vec3 &p2 = target.Position;
		Vec3 vec = p2 - p1;
		float distance = vec.Length();
		float blindAngle = b.Velocity.Negative().AngleTo( vec );
		// check if they are the same or not ( todo: this is broken at this point)
		if ( b.Position.DistanceToSqr( target.Position ) > 0.0001f )
		{
			// check if the distance is nearby enough
			if ( distance <= PerceptionRadius )
			{
				// check if we can 'see it'
				if ( BlindspotAngleDeg <= blindAngle || b.Velocity.Length() == 0 )
				{
					NearbyBoid nb;
					nb.boid = &target;
					nb.distance = distance;
					nb.direction = vec;
					out.push_back( nb );
				}
			}
		}
	}
}

void Grid::ComputeBoundingBox( const vector<Boid> &b )
{
	// default values
	float minX = FLT_MAX;
	float maxX = -FLT_MAX;
	float minY = FLT_MAX;
	float maxY = -FLT_MAX;
	float minZ = FLT_MIN;
	float maxZ = -FLT_MAX;

	// loop over the boids to find
	// the actual value's.

	// todo: use SOA and SiMD instructions
	// to prevent branching.
	for ( int i = 0; i < b.size(); i++ )
	{
		if ( b[i].Position.X < minX ) minX = b[i].Position.X;
		if ( b[i].Position.X > maxX ) maxX = b[i].Position.X;
		if ( b[i].Position.Y < minY ) minY = b[i].Position.Y;
		if ( b[i].Position.Y > maxY ) maxY = b[i].Position.Y;
		if ( b[i].Position.Z < minZ ) minZ = b[i].Position.Z;
		if ( b[i].Position.Z > maxZ ) maxZ = b[i].Position.Z;
	}

	// store the final min / max values.
	( *this ).min = Vec3( minX, minY, minZ );
	( *this ).max = Vec3( maxX, maxY, maxZ );
}

void Grid::ComputeGridIndex( const Boid &b, int &celX, int &celY, int &celZ )
{
	// place the boid into 'grid-space'
	Vec3 boidPosRelative = b.Position - min;

	// on default, we take the 0'th index
	// this happens if one of the dimensions
	// is not being used.
	celX = 0;
	celY = 0;
	celZ = 0;

	// compute the grid-space of the grid
	// hurr
	Vec3 gridExtend = max - min;

	// check if the dimension is not 'flat'.
	if ( abs( gridExtend.X ) >= 0.0001f )
	{
		float cellSizeX = gridExtend.X / ( nx - 1 );
		celX = (int)( boidPosRelative.X / cellSizeX );
	}

	if ( abs( gridExtend.Y ) >= 0.0001f )
	{
		float cellSizeY = gridExtend.Y / ( ny - 1 );
		celY = (int)( boidPosRelative.Y / cellSizeY );
	}

	if ( abs( gridExtend.Z ) >= 0.0001f )
	{
		float cellSizeZ = gridExtend.Z / ( nz - 1 );
		celZ = (int)( boidPosRelative.Z / cellSizeZ );
	}
}

void Grid::StoreInCells( const vector<Boid> &vb )
{
	for ( Boid b : vb )
	{
		// retrieve the index
		int ix, iy, iz;
		ComputeGridIndex( b, ix, iy, iz );

		// add to the correct cell
		int i = CalculateGridCellIndex( ix, iy, iz );
		//printf( "ix: %i, iy: %i, iz: %i, i: %i\r\n", ix, iy, iz, i );
		cells[i].AddBoid( b );
	}
}

GridCell::GridCell()
{
	count = 0;
}

void GridCell::AddBoid( const Boid &b )
{
	if ( count < NUMBER_OF_ELEMENTS_IN_CELL )
	{
		boids[count] = b;
		count++;
	}
	else
	{
		// multiple policies available:
		// - do nothing (lose the boid)
		// - evict with any cache eviction policy, RR preferred
	}
}

void GridCell::Clear()
{
	// best clear evor.
	count = 0;
}
