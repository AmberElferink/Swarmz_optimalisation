
#include "precomp.h"
using namespace sw;

Grid::Grid( int numberOfBuckets, int sizeOfBuckets, int nx, int ny, int nz )
{
	// receive 'dem dimensions.
	( *this ).nx = nx;
	( *this ).ny = ny;
	( *this ).nz = nz;

	// construct an rng
	std::random_device rd;
	eng = std::mt19937( rd() );

	// construct a bucket pool
	bp = new BucketPool( numberOfBuckets, sizeOfBuckets );

	// construct the cells
	int n = nx * ny * nz;
	cells.resize( n );
	for ( int j = 0; j < n; j++ )
	{
		cells[j] = new GridCell( numberOfBuckets );
	}
}

Grid::~Grid()
{
	// huuurrrrr!
	cells.clear();

	// harrrrrr!!
	bp->~BucketPool();
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
		cells[j]->Clear();

	bp->Clear();
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

// perhaps: add boid index number?
void Grid::QueryGrid( const Boid &b, SumVectors &s, const float PerceptionRadius, const float BlindspotAngleDeg, const int ix, const int iy, const int iz, const DistanceType SeparationType )
{
	// if the location is not inside
	// the grid, skip it.
	if ( !CheckInsideGrid( ix, iy, iz ) )
		return;

	float epsilon = 0.0001f;

	// --------------------------------------------------------------------------------
	// Loop hoisting

	// compute negative direction
	float bVelocityNegX = -b.Velocity.X;
	float bVelocityNegY = -b.Velocity.Y;
	float bVelocityNegZ = -b.Velocity.Z;

	// compute length
	float bVelocityLength = FloatVCalc::Length( b.Velocity.X, b.Velocity.Y, b.Velocity.Z );

	// compute negative normalized direction
	float bVelocityLengthRecpr = 1.0f / bVelocityLength;
	float bVelocityNegNormX = bVelocityNegX * bVelocityLengthRecpr;
	float bVelocityNegNormY = bVelocityNegY * bVelocityLengthRecpr;
	float bVelocityNegNormZ = bVelocityNegZ * bVelocityLengthRecpr;

	// --------------------------------------------------------------------------------
	// Prepare the structures

	// represents the relevant data to work with.
	__declspec( align( 64 ) ) float directionToBoidX[NUMBER_OF_ELEMENTS];
	__declspec( align( 64 ) ) float directionToBoidY[NUMBER_OF_ELEMENTS];
	__declspec( align( 64 ) ) float directionToBoidZ[NUMBER_OF_ELEMENTS];
	__declspec( align( 64 ) ) float distanceToBoid[NUMBER_OF_ELEMENTS];
	__declspec( align( 64 ) ) float anglesToBoid[NUMBER_OF_ELEMENTS];

	// represents the relevant indices within a bucket
	int relevantIndices[NUMBER_OF_ELEMENTS];

	// --------------------------------------------------------------------------------
	// Perform the computations

	// retrieve the cell
	const GridCell *gridCell = cells[CalculateGridCellIndex( ix, iy, iz )];
	for ( int bi = 0; bi < gridCell->numberOfBuckets; bi++ )
	{
		// retrieve the bucket
		const Bucket *bucket = bp->GetBucket( gridCell->bpi[bi] );

		int phaseOne = bucket->count;

		// phase 1:
		//  - compute all the distances
		//  - store the relevant indices

		for ( int i = 0; i < phaseOne; i++ )
		{
			directionToBoidX[i] = bucket->posX[i] - b.Position.X;
			directionToBoidY[i] = bucket->posY[i] - b.Position.Y;
			directionToBoidZ[i] = bucket->posZ[i] - b.Position.Z;

			distanceToBoid[i] = sqrtf(
				directionToBoidX[i] * directionToBoidX[i] +
				directionToBoidY[i] * directionToBoidY[i] +
				directionToBoidZ[i] * directionToBoidZ[i] );

			relevantIndices[i] = i;
		}

		// check for phase 2:
		//  - distance too small check
		//		(if true, change the SumVectors directly with the random direction)
		//  - distance too large check
		//  - id check

		int phaseTwo = 0;
		for ( int i = 0; i < phaseOne; i++ )
		{
			if ( s.index == bucket->indx[i] )
				continue;

			if ( distanceToBoid[i] > PerceptionRadius )
				continue;

			if ( distanceToBoid[i] < epsilon )
			{
				// put a random factor on SumVectors - two boids
				// that are not the same are on top of one
				// another!

				// (the following line is _old_)
				// separationSum += Vec3::GetRandomUniform( eng ) * 1000;

				continue;
			}

			// move the data, same cache line
			directionToBoidX[phaseTwo] = directionToBoidX[i];
			directionToBoidY[phaseTwo] = directionToBoidY[i];
			directionToBoidZ[phaseTwo] = directionToBoidZ[i];

			distanceToBoid[phaseTwo] = distanceToBoid[i];

			relevantIndices[phaseTwo] = relevantIndices[i];

			phaseTwo++;
		}

		// phase 2: compute all angles
		//  - skip this phase if velocity is smaller than epsilon

		// check for phase 3:
		//  - skip this check if velocity is smaller than epsilon
		//  - angle too large check

		int phaseThree = 0;
		if ( bVelocityLength < epsilon )
		{
			// we're standing still, we can safely look
			// around without bumping our head into
			// some window

			phaseThree = phaseTwo;
		}
		else
		{
			// we're moving, find out which boids
			// we can see without moving our
			// head too much

			for ( int i = 0; i < phaseTwo; i++ )
			{
				// compute the angles
				float recd = 1.0f / distanceToBoid[i];
				float distanceVecNormX = directionToBoidX[i] * recd;
				float distanceVecNormY = directionToBoidY[i] * recd;
				float distanceVecNormZ = directionToBoidZ[i] * recd;

				anglesToBoid[i] = FloatVCalc::AngleToNorm( bVelocityNegNormX, bVelocityNegNormY, bVelocityNegNormZ, distanceVecNormX, distanceVecNormY, distanceVecNormZ );
			}

			for ( int i = 0; i < phaseTwo; i++ )
			{
				// if it's behind us, continue!
				if ( BlindspotAngleDeg > anglesToBoid[i] )
					continue;

				// move the data, same cache line
				directionToBoidX[phaseThree] = directionToBoidX[i];
				directionToBoidY[phaseThree] = directionToBoidY[i];
				directionToBoidZ[phaseThree] = directionToBoidZ[i];

				distanceToBoid[phaseThree] = distanceToBoid[i]; 

				relevantIndices[phaseThree] = relevantIndices[i];

				phaseThree++;
			}
		}

		// phase 3: update the SumVectors values with the remaining values

		s.count += phaseThree;
		for (int i = 0; i < phaseThree; i++)
		{
			//calculate the sumVecs based on this neighbour
			float separationFactor = SWARMZ_TransformDistance( distanceToBoid[i], SeparationType );
			//separationSum += closeBoid.direction.Negative() * separationFactor;
			s.separationSumX += ( -directionToBoidX[i] ) * separationFactor;
			s.separationSumY += ( -directionToBoidY[i] ) * separationFactor;
			s.separationSumZ += ( -directionToBoidZ[i] ) * separationFactor;

			//headingSum += closeBoid.boid.Velocity;
			s.headingSumX += bucket->velX[relevantIndices[i]];
			s.headingSumY += bucket->velY[relevantIndices[i]];
			s.headingSumZ += bucket->velZ[relevantIndices[i]];

			//positionSum += closeBoid.boid.Position;
			s.positionSumX += bucket->posX[relevantIndices[i]];
			s.positionSumY += bucket->posY[relevantIndices[i]];
			s.positionSumZ += bucket->posZ[relevantIndices[i]];
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
				const GridCell *gridCell = cells[CalculateGridCellIndex( x, y, z )];
				for ( int bi = 0; bi < gridCell->numberOfBuckets; bi++ )
				{
					const Bucket *bucket = bp->GetBucket( gridCell->bpi[bi] );
					count += bucket->count;
				}
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
				const GridCell *gridCell = cells[CalculateGridCellIndex( x, y, z )];
				for ( int bi = 0; bi < gridCell->numberOfBuckets; bi++ )
				{
					const Bucket *bucket = bp->GetBucket( gridCell->bpi[bi] );
					count += bucket->count;
				}
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
	int index = 0;
	for ( const Boid &b : vb )
	{

		// retrieve the index
		int ix, iy, iz;
		ComputeGridIndex( b, ix, iy, iz );

		// add to the correct cell
		const int i = CalculateGridCellIndex( ix, iy, iz );
		GridCell *cell = cells[i];
		cell->AddBoid( bp, b, index );
		index++;
	}
}

GridCell::GridCell( int n )
{
	numberOfBuckets = 0;
	bpi = new int[n];
}

GridCell::~GridCell()
{
	delete bpi;
}

void GridCell::AddBoid( BucketPool *bp, const Boid &boid, int index )
{
	// base case: the first boid will always trigger a bucket.
	if ( numberOfBuckets == 0 )
	{
		bpi[numberOfBuckets] = bp->ReserveBucket();
		numberOfBuckets++;
	}

	Bucket *b = bp->GetBucket( bpi[numberOfBuckets - 1] );

	// consequent case: bucket may be full
	if ( b->count >= b->maximum )
	{
		// reserve a new bucket
		bpi[numberOfBuckets] = bp->ReserveBucket();

		// retrieve the bucket as the one we're using right now.
		b = bp->GetBucket( bpi[numberOfBuckets] );

		// keep track of the amount of buckets we're
		// using.
		numberOfBuckets++;
	}

	// we know the bucket _must_ have space for this
	// boid, so we can just add it in.
	int bCount = b->count;
	b->posX[bCount] = boid.Position.X;
	b->posY[bCount] = boid.Position.Y;
	b->posZ[bCount] = boid.Position.Z;
	b->velX[bCount] = boid.Velocity.X;
	b->velY[bCount] = boid.Velocity.Y;
	b->velZ[bCount] = boid.Velocity.Z;
	b->indx[bCount] = index;
	b->count++;
}

void GridCell::Clear()
{
	// best clear evor.
	numberOfBuckets = 0;
}
