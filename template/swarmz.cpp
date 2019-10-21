
#include "precomp.h"

#include "Grid.h"
#include "precomp.h"

Grid::Grid( int xcells, int ycells, int zcells )
{
	( *this ).xcells = xcells;
	( *this ).ycells = ycells;
	( *this ).zcells = zcells;

	( *this ).cells = new GridCell[xcells * ycells * zcells];
}

Grid::~Grid()
{
	delete cells;
}

inline int Grid::CalculateGridCellIndex( int celX, int celY, int celZ )
{
	return celX + celY * xcells + celZ * ( xcells * ycells );
}

inline bool Grid::CheckInsideGrid( int celX, int celY, int celZ )
{
	if ( celX < 0 || celX >= xcells )
		return false;
	if ( celY < 0 || celY >= ycells )
		return false;
	if ( celZ < 0 || celZ >= zcells )
		return false;

	return true;
}

void Grid::ConstructGrid( const vector<Boid> &b )
{
	int l = xcells * ycells * zcells;

	for ( int j = 0; j < l; j++ )
	{
		cells[j].Clear();
	}

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
	throw new exception( "unimplemented" );
}

void Grid::QueryGrid( const Boid &b, const int r, vector<NearbyBoid> &out, float PerceptionRadius, float BlindspotAngleDeg, int celX, int celY, int celZ )
{
	if ( !CheckInsideGrid( celX, celY, celZ ) )
		return;

	GridCell gridCell = cells[CalculateGridCellIndex( celX, celY, celZ )];

	for ( int i = 0; i < gridCell.count; i++ )
	{
		//compute distance between b and test
		Boid test = gridCell.boids[i];
		const Vec3 &p1 = b.Position;
		const Vec3 &p2 = test.Position;
		Vec3 vec = p2 - p1;
		float distance = vec.Length();
		float blindAngle = b.Velocity.Negative().AngleTo( vec );
		if (
			b.Position.DistanceToSqr( test.Position ) < 0.0001f &&			  // check if they are the same or not ( todo: check on some id)
			distance <= PerceptionRadius &&									  // check if the distance is nearby enough
			( BlindspotAngleDeg <= blindAngle || b.Velocity.Length() == 0 ) ) // check if we can 'see it'
		{
			NearbyBoid nb;
			nb.boid = &test;
			nb.distance = distance;
			nb.direction = vec;
			out.push_back( nb );
		}

		//compute the angle between b and test
		// if both okay, add to vector (and length == 0 to see if the boid is standing still and can look around)
	}
}

void Grid::ComputeBoundingBox( const vector<Boid> &b )
{
	float minX = FLT_MAX;
	float maxX = -FLT_MAX;
	float minY = FLT_MAX;
	float maxY = -FLT_MAX;
	float minZ = FLT_MIN;
	float maxZ = -FLT_MAX;

	for ( int i = 0; i < b.size(); i++ )
	{
		if ( b[i].Position.X < minX ) minX = b[i].Position.X;
		if ( b[i].Position.X > maxX ) maxX = b[i].Position.X;
		if ( b[i].Position.Y < minY ) minY = b[i].Position.Y;
		if ( b[i].Position.Y > maxY ) maxY = b[i].Position.Y;
		if ( b[i].Position.Z < minZ ) minZ = b[i].Position.Z;
		if ( b[i].Position.Z > maxZ ) maxZ = b[i].Position.Z;
	}
	( *this ).min = Vec3( minX, minY, minZ );
	( *this ).max = Vec3( maxX, maxY, maxZ );
}

void Grid::ComputeGridIndex( const Boid &b, int &celX, int &celY, int &celZ )
{
	celX = 0;
	celY = 0;
	celZ = 0;

	Vec3 boidPosRelative = b.Position - min;
	float sizeX = max.X - min.X;
	if ( abs( sizeX ) >= 0.0001f )
	{
		float cellSizeX = sizeX / ( xcells - 1 );
		celX = (int)( boidPosRelative.X / cellSizeX );
	}

	float sizeY = max.Y - min.Y;
	if ( abs( sizeY ) >= 0.0001f )
	{
		float cellSizeY = sizeY / (ycells - 1);
		celY = (int)( boidPosRelative.Y / cellSizeY );
	}

	float sizeZ = max.Z - min.Z;
	if ( abs( sizeZ ) >= 0.0001f )
	{
		float cellSizeZ = sizeZ / (zcells - 1);
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
		GridCell *cell = &cells[i];
		cell->AddBoid( b );
	}
}

void Grid::ReorderData()
{
	// todo: reorder the data to make it more friendly.
}

GridCell::GridCell()
{
	count = 0;
}

void GridCell::AddBoid( const Boid &b )
{
	if ( count >= NUMBER_OF_ELEMENTS_IN_CELL )
	{
		//eventueel random evicten
	}
	else
	{
		boids[count] = b;
		count++;
	}
}

void GridCell::Clear()
{
	count = 0;
}
