
#include "precomp.h"
using namespace sw;

//https://stackoverflow.com/questions/46974513/code-for-acos-with-avx256
__m256 acos( __m256 x )
{
	__m256 xp = _mm256_and_ps( x, _mm256_castsi256_ps( _mm256_set1_epi32( 0x7FFFFFFF ) ) );
	// main shape
	__m256 one = _mm256_set1_ps( 1.0 );
	__m256 t = _mm256_sqrt_ps( _mm256_sub_ps( one, xp ) );
	// polynomial correction factor based on xp
	__m256 c3 = _mm256_set1_ps( -0.02007522 );
	__m256 c2 = _mm256_fmadd_ps( xp, c3, _mm256_set1_ps( 0.07590315 ) );
	__m256 c1 = _mm256_fmadd_ps( xp, c2, _mm256_set1_ps( -0.2126757 ) );
	__m256 c0 = _mm256_fmadd_ps( xp, c1, _mm256_set1_ps( 1.5707963267948966 ) );
	// positive result
	__m256 p = _mm256_mul_ps( t, c0 );
	// correct for negative x
	__m256 n = _mm256_sub_ps( _mm256_set1_ps( 3.14159265359 ), p );
	return _mm256_blendv_ps( p, n, x );
}

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
void Grid::QueryGrid(
	const Boid &b, const int boidIndex, NearbyBoidsData &s,
	const float PerceptionRadius, const float BlindspotAngleDeg,
	const int ix, const int iy, const int iz,
	const DistanceType SeparationType )
{
	// if the location is not inside
	// the grid, skip it.
	if ( !CheckInsideGrid( ix, iy, iz ) )
		return;

	float epsilon = 0.0001f;

	// info about the _mm256_cmp_ps function
	// https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_mm256_cmp_ps&expand=744
	// https://stackoverflow.com/questions/16988199/how-to-choose-avx-compare-predicate-variants
	// https://stackoverflow.com/questions/8627331/what-does-ordered-unordered-comparison-mean

	// --------------------------------------------------------------------------------
	// Loop hoisting

	//---------------------------------------------------------------------
	// Widening constants

	__m256 PerceptionRadius4 = _mm256_set1_ps( PerceptionRadius );
	__m256 BlindspotAngleDeg4 = _mm256_set1_ps( BlindspotAngleDeg );
	__m256 toRadian4 = _mm256_set1_ps( toRadian );

	__m256i bIndex4 = _mm256_set1_epi32( boidIndex );
	__m256i ones = _mm256_set1_epi32( 0xffffffff );

#pragma region Loop hoisting

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

#pragma endregion

	//---------------------------------------------------------------------
	// Widening boids

#pragma region Widening of the boid

	__m256 bPositionX4 = _mm256_set1_ps( b.Position.X );
	__m256 bPositionY4 = _mm256_set1_ps( b.Position.Y );
	__m256 bPositionZ4 = _mm256_set1_ps( b.Position.Z );

	__m256 bVelocityLengthRecpr4 = _mm256_set1_ps( bVelocityLengthRecpr );
	__m256 bVelocityNegNormX4 = _mm256_set1_ps( bVelocityNegNormX );
	__m256 bVelocityNegNormY4 = _mm256_set1_ps( bVelocityNegNormY );
	__m256 bVelocityNegNormZ4 = _mm256_set1_ps( bVelocityNegNormZ );

#pragma endregion

	// --------------------------------------------------------------------------------
	// Prepare the structures

#pragma region Structure preparation

	// __declspec( align( 64 ) )
	// means that the start address of these arrays are
	// a multiple of 64. SiMD requires the start address to
	// be a multiple of 16 (SSE) or 32 (AVX).

	// represents the relevant data to work with.
	__declspec( align( 64 ) ) union {
		float directionToBoidX[ELEMENTS_IN_BUCKET];
		__m256 directionToBoidX4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float directionToBoidY[ELEMENTS_IN_BUCKET];
		__m256 directionToBoidY4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};
	__declspec( align( 64 ) ) union {
		float directionToBoidZ[ELEMENTS_IN_BUCKET];
		__m256 directionToBoidZ4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float distanceToBoid[ELEMENTS_IN_BUCKET];
		__m256 distanceToBoid4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float distanceToBoidInv[ELEMENTS_IN_BUCKET];
		__m256 distanceToBoidInv4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float anglesToBoid[ELEMENTS_IN_BUCKET];
		__m256 anglesToBoid4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		int mask[ELEMENTS_IN_BUCKET];
		__m256 mask4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float tempA[ELEMENTS_IN_BUCKET];
		__m256 tempA4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float tempB[ELEMENTS_IN_BUCKET];
		__m256 tempB4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float tempC[ELEMENTS_IN_BUCKET];
		__m256 tempC4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

#pragma endregion

	// --------------------------------------------------------------------------------
	// Perform the computations

	// retrieve the cell
	int index = CalculateGridCellIndex( ix, iy, iz );
	const GridCell *gridCell = cells[index];

	for ( int bi = 0, bl = gridCell->numberOfBuckets; bi < bl; bi++ )
	{
		// retrieve the bucket
		const Bucket *bucket = bp->GetBucket( gridCell->bpi[bi] );
		const int iterations = ( bucket->count + SIMDSIZE - 1 ) / SIMDSIZE;

		// --------------------------------------------------------------------------------
		// Attempt to prefetch

		if ( bi + 1 < bl )
		{
			const Bucket *bucketPrefetch = bp->GetBucket( gridCell->bpi[bi + 1] );
			PREFETCH( bucketPrefetch->velX4 );
			PREFETCH( bucketPrefetch->velY4 );
			PREFETCH( bucketPrefetch->velZ4 );
			PREFETCH( bucketPrefetch->posX4 );
			PREFETCH( bucketPrefetch->posY4 );
			PREFETCH( bucketPrefetch->posZ4 );
			PREFETCH( &bPositionX4 );
			PREFETCH( &bPositionY4 );
			PREFETCH( &bPositionZ4 );
		}

		// --------------------------------------------------------------------------------
		// Compute direction
		// Compute (inverse) distance

		for ( int i = 0; i < iterations; i++ )
		{
			directionToBoidX4[i] = _mm256_sub_ps( bucket->posX4[i], bPositionX4 );
			directionToBoidY4[i] = _mm256_sub_ps( bucket->posY4[i], bPositionY4 );
			directionToBoidZ4[i] = _mm256_sub_ps( bucket->posZ4[i], bPositionZ4 );
		}

		for ( int i = 0; i < iterations; i++ )
		{
			// compute the distance
			distanceToBoid4[i] = _mm256_sqrt_ps(
				_mm256_add_ps(
					_mm256_add_ps(
						_mm256_mul_ps( directionToBoidX4[i], directionToBoidX4[i] ),
						_mm256_mul_ps( directionToBoidY4[i], directionToBoidY4[i] ) ),
					_mm256_mul_ps( directionToBoidZ4[i], directionToBoidZ4[i] ) ) );

			distanceToBoidInv4[i] = _mm256_rcp_ps( distanceToBoid4[i] );
		}

		// --------------------------------------------------------------------------------
		// Compute index mask
		// Compute distance mask

		for ( int i = 0; i < iterations; i++ )
		{
			// todo: inverse mask

			// ((uint*)(mask4[0].m256_f32))[0]
			mask4[i] = _mm256_castsi256_ps( 
				// inverse
				_mm256_xor_si256( 
					_mm256_cmpeq_epi32( bIndex4, bucket->indx4[i] ), 
					ones ) );
		}

		for ( int i = 0; i < iterations; i++ )
		{
			mask4[i] =
				_mm256_and_ps(
					_mm256_cmp_ps( distanceToBoid4[i], PerceptionRadius4, _CMP_LT_OQ ),
					mask4[i] );
		}

		// --------------------------------------------------------------------------------
		// Compute the normalized direction
		// Compute the dot products
		// Compute the cosines
		// Compute the angle
		// Compute the mask for the angle

		if ( bVelocityLength >= epsilon )
		{
			for ( int i = 0; i < iterations; i++ )
			{
				tempA4[i] = _mm256_mul_ps( directionToBoidX4[i], distanceToBoidInv4[i] );
				tempB4[i] = _mm256_mul_ps( directionToBoidY4[i], distanceToBoidInv4[i] );
				tempC4[i] = _mm256_mul_ps( directionToBoidZ4[i], distanceToBoidInv4[i] );
			}

			for ( int i = 0; i < iterations; i++ )
			{
				tempA4[i] = _mm256_mul_ps( tempA4[i], tempA4[i] );
				tempB4[i] = _mm256_mul_ps( tempB4[i], tempB4[i] );
				tempC4[i] = _mm256_mul_ps( tempC4[i], tempC4[i] );
			}

			for ( int i = 0; i < iterations; i++ )
			{
				tempA4[i] = _mm256_add_ps( tempA4[i], _mm256_add_ps( tempB4[i], tempC4[i] ) );
				tempA4[i] = _mm256_min_ps( tempA4[i], _mm256_set1_ps( 1.0f ) );
				tempA4[i] = _mm256_max_ps( tempA4[i], _mm256_set1_ps( -1.0f ) );
			}

			// prefetch the radians / perception radius
			PREFETCH( &toRadian4 );
			PREFETCH( &PerceptionRadius4 );

			for ( int i = 0; i < iterations; i++ )
			{
				tempC4[i] = acos( tempA4[i] );
			}

			for ( int i = 0; i < iterations; i++ )
			{
				anglesToBoid4[i] = _mm256_mul_ps( toRadian4, tempC4[i] );
			}

			for ( int i = 0; i < iterations; i++ )
			{
				mask4[i] = _mm256_and_ps(
					_mm256_cmp_ps( BlindspotAngleDeg4, anglesToBoid4[i], _CMP_GT_OQ ),
					mask4[i] );
			}
		}

		// --------------------------------------------------------------------------------
		// Compute the number of valid boids
		// Compute the seperation factors
		// Compute the resulting forces in a masked manner

		int boids = 0;
		for ( int i = 0; i < iterations; i++ )
		{
			for ( int j = 0; j < SIMDSIZE; j++ )
			{
				if ( mask[i] )
					boids++;
			}
		}
		s.count = boids;

		switch ( SeparationType )
		{
		case DistanceType::LINEAR:
			for ( int i = 0; i < iterations; i++ )
			{ tempB4[i] = distanceToBoid4[i]; }
			break;

		case DistanceType::INVERSE_LINEAR:
			for ( int i = 0; i < iterations; i++ )
			{
				tempB4[i] = _mm256_rcp_ps( distanceToBoid4[i] );
				tempB4[i] = _mm256_blendv_ps( tempB4[i], _mm256_setzero_ps(),
											  _mm256_cmp_ps( distanceToBoid4[i], _mm256_setzero_ps(), _CMP_EQ_OQ ) );
			}
			break;

		case DistanceType::QUADRATIC:
			for ( int i = 0; i < iterations; i++ )
			{ tempB4[i] = _mm256_mul_ps( distanceToBoid4[i], distanceToBoid4[i] ); }
			break;

		case DistanceType::INVERSE_QUADRATIC:
			for ( int i = 0; i < iterations; i++ )
			{
				tempB4[i] = _mm256_rcp_ps( _mm256_mul_ps( distanceToBoid4[i], distanceToBoid4[i] ) );
				tempB4[i] = _mm256_blendv_ps( tempB4[i], _mm256_setzero_ps(),
											  _mm256_cmp_ps( distanceToBoid4[i], _mm256_setzero_ps(), _CMP_EQ_OQ ) );
			}
			break;
		}

		for ( int i = 0; i < iterations; i++ )
		{
			if ( _mm256_movemask_ps( mask4[i] ) )
			{
				s.separationSumX4 = _mm256_blendv_ps(
					s.separationSumX4,
					_mm256_add_ps( s.separationSumX4, _mm256_mul_ps( _mm256_sub_ps( _mm256_setzero_ps(), directionToBoidX4[i] ), tempB4[i] ) ),
					mask4[i] );

				s.separationSumY4 = _mm256_blendv_ps(
					s.separationSumY4,
					_mm256_add_ps( s.separationSumY4, _mm256_mul_ps( _mm256_sub_ps( _mm256_setzero_ps(), directionToBoidY4[i] ), tempB4[i] ) ),
					mask4[i] );

				s.separationSumZ4 = _mm256_blendv_ps(
					s.separationSumZ4,
					_mm256_add_ps( s.separationSumZ4, _mm256_mul_ps( _mm256_sub_ps( _mm256_setzero_ps(), directionToBoidZ4[i] ), tempB4[i] ) ),
					mask4[i] );

				s.headingSumX4 = _mm256_blendv_ps(
					s.headingSumX4,
					_mm256_add_ps( s.headingSumX4, bucket->velX4[i] ),
					mask4[i] );

				s.headingSumY4 = _mm256_blendv_ps(
					s.headingSumY4,
					_mm256_add_ps( s.headingSumY4, bucket->velY4[i] ),
					mask4[i] );

				s.headingSumZ4 = _mm256_blendv_ps(
					s.headingSumZ4,
					_mm256_add_ps( s.headingSumZ4, bucket->velZ4[i] ),
					mask4[i] );

				s.positionSumX4 = _mm256_blendv_ps(
					s.positionSumX4,
					_mm256_add_ps( s.positionSumX4, bucket->posX4[i] ),
					mask4[i] );

				s.positionSumY4 = _mm256_blendv_ps(
					s.positionSumY4,
					_mm256_add_ps( s.positionSumY4, bucket->posY4[i] ),
					mask4[i] );

				s.positionSumZ4 = _mm256_blendv_ps(
					s.positionSumZ4,
					_mm256_add_ps( s.positionSumZ4, bucket->posZ4[i] ),
					mask4[i] );
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

	// todo: use SOA and __m256 instructions
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
