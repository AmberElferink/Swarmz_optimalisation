#include "precomp.h"
#include "Grid.h"

Grid::Grid( int xcells , int ycells, int zcells )
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

inline int Grid::CalculateGridCellIndex(int celX, int celY, int celZ)
{
	return celX + celY * xcells + celZ * ( xcells * ycells );
}

inline bool Grid::CheckInsideGrid(int celX, int celY, int celZ)
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
	if(!CheckInsideGrid( celX, celY, celZ ))
		return;

	GridCell gridCell = cells[CalculateGridCellIndex( celX, celY, celZ )];

	for (int i = 0; i < gridCell.count; i++)
	{
		//compute distance between b and test
		Boid test = gridCell.boids[i];
		const Vec3 &p1 = b.Position;
		const Vec3 &p2 = test.Position;
		Vec3 vec = p2 - p1;
		float distance = vec.Length();
		float blindAngle = b.Velocity.Negative().AngleTo( vec );
		if ( ( b.id ) != test.id && distance <= PerceptionRadius && ( BlindspotAngleDeg <= blindAngle || b.Velocity.Length() == 0 ) )
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

	for (int i = 0; i < b.size(); i++)
	{
		if ( b[i].Position.X < minX ) minX = b[i].Position.X;
		if ( b[i].Position.X > maxX ) maxX = b[i].Position.X;
		if ( b[i].Position.Y < minY ) minY = b[i].Position.Y;
		if ( b[i].Position.Y > maxY ) maxY = b[i].Position.Y;
		if ( b[i].Position.Z < minZ ) minZ = b[i].Position.Z;
		if ( b[i].Position.Z > maxZ ) maxZ = b[i].Position.Z;
	}
	( *this ).min = Vec3( min.X, min.Y, min.Z );
	( *this ).max = Vec3( max.X, max.Y, max.Z );
}

void Grid::ComputeGridIndex(const Boid &b, int &celX, int &celY,  int &celZ)
{

	Vec3 boidPosRelative = b.Position - min;
	float sizeX = max.X - min.X;
	float cellSizeX = sizeX / xcells;
	celX = (int)( boidPosRelative.X / cellSizeX );

	float sizeY = max.Y - min.Y;
	float cellSizeY = sizeY / ycells;
	 celY = (int)( boidPosRelative.Y / cellSizeY );

	float sizeZ = max.Z - min.Z;
	float cellSizeZ = sizeZ / zcells;
	celZ = (int)( boidPosRelative.Z / cellSizeZ );

	//check for the custom constructor
	//if implemented
}

void Grid::StoreInCells( const vector<Boid> &b )
{
	Vec3 boidPosRelative = b[0].Position - min;
	// todo: store 'dem.
	float sizeX = max.X - min.X;
	float cellSizeX = sizeX / xcells;
	int celX = (int) (boidPosRelative.X / cellSizeX);

	float sizeY = max.Y - min.Y;
	float cellSizeY = sizeY / ycells;
	int celY = (int)( boidPosRelative.Y / cellSizeY );

	float sizeZ = max.Z - min.Z;
	float cellSizeZ = sizeZ / zcells;
	int celZ = (int)( boidPosRelative.Z / cellSizeZ );

	GridCell *cell = &cells[celX + celY * xcells + celZ * ( xcells * ycells )];
	cell->AddBoid( b[0] );
}

void Grid::ReorderData()
{
	// todo: reorder the data to make it more friendly.
}

void GridCell::AddBoid( const Boid &b )
{
	if (count >= NUMBER_OF_ELEMENTS_IN_CELL)
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
