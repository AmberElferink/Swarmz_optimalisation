
#include "precomp.h"
using namespace sw;

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

void Grid::ClearGrid()
{
	// compute the total number of cells.
	int l = nx * ny * nz;

	// for every cell, clear it out.
	for ( int j = 0; j < l; j++ )
		cells[j].Clear();
}

void Grid::ConstructGrid( const vector<Boid> &b, float perceptionRadius )
{
	// clear out the grid of any
	// previous usage
	ClearGrid();

	// compute the bounding box
	ComputeBoundingBox( b, perceptionRadius );

	// store the data
	StoreInCells( b );
}

void Grid::QueryGrid( const Boid &b, const int r, vector<NearbyBoid> &out, float PerceptionRadius, float BlindspotAngleDeg, int ix, int iy, int iz )
{
	// if the location is not inside
	// the grid, skip it.
	if ( !CheckInsideGrid( ix, iy, iz ) )
		return;

	// retrieve the cell
	const GridCell& gridCell = cells[CalculateGridCellIndex( ix, iy, iz )];

	// do cell computations, this is copied (for now)
	// - fix expensive operations
	// - apply if statements sooner (e.g., compute
	// distance -> check, compute angle -> check,
	// etc), allows for early-opt out.

	//for ( const Boid &target : gridCell.boids ) ;
	for ( int i = 0; i < gridCell.count; i++ )
	{
		//compute distance between b and test
		const Boid &target = gridCell.boids[i];
		const Vec3 &p1 = b.Position;
		const Vec3 &p2 = target.Position;


		Vec3 distanceVec = p2 - p1;
		float distance = distanceVec.Length();
		Vec3 distanceVecNorm = distanceVec / distance;
		//float dstSqr = b.Position.DistanceToSqr( p2 );

		Vec3 bNegVelocity = b.Velocity.Negative();
		float bNegVelocityLength = bNegVelocity.Length();



		// check if they are the same or not ( todo: this is broken at this point)
		if ( distance > 0.001f )
		{
			// check if the distance is nearby enough
			if ( distance <= PerceptionRadius )
			{
				float blindAngle = 0;
				if ( bNegVelocityLength > 0.000001f )
				{
					Vec3 bNegVelocityNorm = bNegVelocity / bNegVelocityLength;
					blindAngle = bNegVelocityNorm.AngleToNorm( distanceVecNorm );
				}

				// check if we can 'see it'
				if ( BlindspotAngleDeg <= blindAngle || bNegVelocityLength == 0 )
				{
					NearbyBoid nb;
					// was: nb.boid = &target
					nb.boid = target;
					nb.distance = distance;
					nb.direction = distanceVec;
					out.emplace_back( nb ); //TODO dit is vaag, moet met mov toch?
				}
			}
		}
	}
}




// todo: remove this function
void Grid::DrawGrid( Surface *surface, Pixel density )
{
	// find the maximum density over the z
	// dimension
	int max = 0;
	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			// gather over the z dimension
			int count = 0;
			for ( int z = 0; z < nz; z++ )
			{
				int index = CalculateGridCellIndex( x, y, z );
				count += cells[index].count;

			}

			// keep track of the largest
			if ( count > max )
				max = count;
		}
	}

	float epsilon = 1.0f;
	// draw all the boxes, stacking on
	// the z dimension
	for ( int x = 0; x < nx; x++ )
	{
		for ( int y = 0; y < ny; y++ )
		{
			// stack on the z dimension
			int count = 0;
			for ( int z = 0; z < nz; z++ )
			{
				int index = CalculateGridCellIndex( x, y, z );
				count += cells[index].count;
			}

			// draw the box
			float factor = (float)count / max;
			Pixel output = ScaleColor( density, (int)( 255 * factor ) );
			surface->Box(
				( SCRWIDTH >> 1 ) + minbb.X + x * step.X + epsilon,
				( SCRHEIGHT >> 1 ) + minbb.Y + y * step.Y + epsilon,
				( SCRWIDTH >> 1 ) + minbb.X + ( x + 1 ) * step.X - epsilon,
				( SCRHEIGHT >> 1 ) + minbb.Y + ( y + 1 ) * step.Y - epsilon, output );
		}
	}

	// draw the global box
	surface->Box(
		( SCRWIDTH >> 1 ) + minbb.X,
		( SCRHEIGHT >> 1 ) + minbb.Y,
		( SCRWIDTH >> 1 ) + maxbb.X,
		( SCRHEIGHT >> 1 ) + maxbb.Y, density );
}

void Grid::ComputeBoundingBox( const vector<Boid> &b, float perceptionRadius )
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

	float epsilon = 0.001f;
	Vec3 avgStep = Vec3(
		( maxX - minX ) / nx + epsilon,
		( maxY - minY ) / ny + epsilon,
		( maxZ - minZ ) / nz + epsilon );

	( *this ).step = Vec3(
		max( perceptionRadius, avgStep.X ),
		max( perceptionRadius, avgStep.Y ),
		max( perceptionRadius, avgStep.Z ) );

	Vec3 stepDiff = step - avgStep;
	Vec3 bbOffset = Vec3(
		stepDiff.X * nx * 0.5f,
		stepDiff.Y * ny * 0.5f,
		stepDiff.Z * nz * 0.5f );
	// store the final min / max values.
	( *this ).minbb = Vec3( minX, minY, minZ ) - bbOffset;
	( *this ).maxbb = Vec3( maxX, maxY, maxZ ) + bbOffset;

	// store the step size for this bounding box
	// add a small bit due to floating point inprecision
	( *this ).step = Vec3(
		max( perceptionRadius, ( maxX - minX ) / nx ) + epsilon,
		max( perceptionRadius, ( maxY - minY ) / ny ) + epsilon,
		max( perceptionRadius, ( maxZ - minZ ) / nz ) + epsilon );

	// todo: center when one (or more) dimension(s) of step size is too small
}

void Grid::ComputeGridIndex( const Boid &b, int &celX, int &celY, int &celZ )
{
	// place the boid into 'grid-space'
	Vec3 boidPosRelative = b.Position - minbb;

	// on default, we take the 0'th index
	// this happens if one of the dimensions
	// is not being used.
	celX = 0;
	celY = 0;
	celZ = 0;

	// check if the dimension is not 'flat'.
	celX = (int)( boidPosRelative.X / step.X );
	celY = (int)( boidPosRelative.Y / step.Y );
	celZ = (int)( boidPosRelative.Z / step.Z );
}

void Grid::StoreInCells( const vector<Boid> &vb )
{
	for ( const Boid &b : vb )
	{
		// retrieve the index
		int ix, iy, iz;
		ComputeGridIndex( b, ix, iy, iz );

		// add to the correct cell
		const int i = CalculateGridCellIndex( ix, iy, iz );
		GridCell &cell = cells[i];		
		cell.AddBoid( b );
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
