
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
void Grid::QueryGrid(
	const Boid &b, SumVectors &s,
	const float PerceptionRadius, const float BlindspotAngleDeg,
	const int ix, const int iy, const int iz,
	const DistanceType SeparationType )
{
	// if the location is not inside
	// the grid, skip it.
	if ( !CheckInsideGrid( ix, iy, iz ) )
		return;

	float epsilon = 0.0001f;

	// --------------------------------------------------------------------------------
	// Loop hoisting

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

	__m256 toRadian4 = _mm256_set1_ps( toRadian );
	__m256 ones = _mm256_set1_ps( 1.0f );

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

	__m256i bIndex4 = _mm256_set1_epi32( s.index );

#pragma endregion

	//---------------------------------------------------------------------
	// Widening boids

#pragma region Accumulation prepartion

	union {
		float separationSumX[8];
		__m256 separationSumX4;
	};

	union {
		float separationSumY[8];
		__m256 separationSumY4;
	};

	union {
		float separationSumZ[8];
		__m256 separationSumZ4;
	};

	union {
		float headingSumX[8];
		__m256 headingSumX4;
	};

	union {
		float headingSumY[8];
		__m256 headingSumY4;
	};

	union {
		float headingSumZ[8];
		__m256 headingSumZ4;
	};

	union {
		float positionSumX[8];
		__m256 positionSumX4;
	};

	union {
		float positionSumY[8];
		__m256 positionSumY4;
	};

	union {
		float positionSumZ[8];
		__m256 positionSumZ4;
	};

	// initialize it all to 0.
	separationSumX4 = _mm256_setzero_ps();
	separationSumY4 = _mm256_setzero_ps();
	separationSumZ4 = _mm256_setzero_ps();

	headingSumX4 = _mm256_setzero_ps();
	headingSumY4 = _mm256_setzero_ps();
	headingSumZ4 = _mm256_setzero_ps();

	positionSumX4 = _mm256_setzero_ps();
	positionSumY4 = _mm256_setzero_ps();
	positionSumZ4 = _mm256_setzero_ps();

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
		float anglesToBoid[ELEMENTS_IN_BUCKET];
		__m256 anglesToBoid4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	// represents the relevant data of the bucket to work with.
	__declspec( align( 64 ) ) union {
		float bucketPositionX[ELEMENTS_IN_BUCKET];
		__m256 bucketPositionX4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float bucketPositionY[ELEMENTS_IN_BUCKET];
		__m256 bucketPositionY4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float bucketPositionZ[ELEMENTS_IN_BUCKET];
		__m256 bucketPositionZ4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float bucketVelocityX[ELEMENTS_IN_BUCKET];
		__m256 bucketVelocityX4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float bucketVelocityY[ELEMENTS_IN_BUCKET];
		__m256 bucketVelocityY4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		float bucketVelocityZ[ELEMENTS_IN_BUCKET];
		__m256 bucketVelocityZ4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};

	__declspec( align( 64 ) ) union {
		int mask[ELEMENTS_IN_BUCKET];
		__m256i mask4[ELEMENTS_IN_BUCKET / SIMDSIZE];
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
		int phaseOne = bucket->count;

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

		// phase 1:
		//  - compute all the distances
		//  - store the relevant indices

		memcpy( bucketVelocityX4, bucket->velX4, ELEMENTS_IN_BUCKET * ( sizeof( float ) ) );
		memcpy( bucketVelocityY4, bucket->velY4, ELEMENTS_IN_BUCKET * ( sizeof( float ) ) );
		memcpy( bucketVelocityZ4, bucket->velZ4, ELEMENTS_IN_BUCKET * ( sizeof( float ) ) );
		memcpy( bucketPositionX4, bucket->posX4, ELEMENTS_IN_BUCKET * ( sizeof( float ) ) );
		memcpy( bucketPositionY4, bucket->posY4, ELEMENTS_IN_BUCKET * ( sizeof( float ) ) );
		memcpy( bucketPositionZ4, bucket->posZ4, ELEMENTS_IN_BUCKET * ( sizeof( float ) ) );

		int phaseOneSimd = ( phaseOne + SIMDSIZE - 1 ) / SIMDSIZE;
		for ( int i = 0; i < phaseOneSimd; i++ )
		{
			directionToBoidX4[i] = _mm256_sub_ps( bucketPositionX4[i], bPositionX4 );
			directionToBoidY4[i] = _mm256_sub_ps( bucketPositionY4[i], bPositionY4 );
			directionToBoidZ4[i] = _mm256_sub_ps( bucketPositionZ4[i], bPositionZ4 );

			// compute the distance
			distanceToBoid4[i] = _mm256_sqrt_ps(
				_mm256_add_ps(
					_mm256_add_ps(
						_mm256_mul_ps( directionToBoidX4[i], directionToBoidX4[i] ),
						_mm256_mul_ps( directionToBoidY4[i], directionToBoidY4[i] ) ),
					_mm256_mul_ps( directionToBoidZ4[i], directionToBoidZ4[i] ) ) );
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
				// todo
				// put a random factor on SumVectors - two boids
				// that are not the same are on top of one
				// another!

				// (the following line is _old_)
				// separationSum += Vec3::GetRandomUniform( eng ) * 1000;

				continue;
			}

			if ( phaseTwo != i )
			{
				// move the data, same cache line
				directionToBoidX[phaseTwo] = directionToBoidX[i];
				directionToBoidY[phaseTwo] = directionToBoidY[i];
				directionToBoidZ[phaseTwo] = directionToBoidZ[i];

				distanceToBoid[phaseTwo] = distanceToBoid[i];

				// transfer our local copy of the bucket contents
				bucketVelocityX[phaseTwo] = bucketVelocityX[i];
				bucketVelocityY[phaseTwo] = bucketVelocityY[i];
				bucketVelocityZ[phaseTwo] = bucketVelocityZ[i];

				bucketPositionX[phaseTwo] = bucketPositionX[i];
				bucketPositionY[phaseTwo] = bucketPositionY[i];
				bucketPositionZ[phaseTwo] = bucketPositionZ[i];
			}

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
			int phaseTwoSimd = ( phaseTwo + SIMDSIZE - 1 ) / SIMDSIZE;
			for ( int i = 0; i < phaseTwoSimd; i++ )
			{
				// compute the normalized direction
				__m256 recd4 = _mm256_div_ps( ones, distanceToBoid4[i] );
				__m256 directionToBoidNormX4 = _mm256_mul_ps( directionToBoidX4[i], recd4 );
				__m256 directionToBoidNormY4 = _mm256_mul_ps( directionToBoidY4[i], recd4 );
				__m256 directionToBoidNormZ4 = _mm256_mul_ps( directionToBoidZ4[i], recd4 );

				// compute dotproduct
				__m256 dotProduct4 = _mm256_add_ps(
					_mm256_add_ps(
						_mm256_mul_ps( bVelocityNegNormX4, directionToBoidNormX4 ),
						_mm256_mul_ps( bVelocityNegNormY4, directionToBoidNormY4 ) ),
					_mm256_mul_ps( bVelocityNegNormZ4, directionToBoidNormZ4 ) );

				// min / max
				anglesToBoid4[i] = _mm256_mul_ps( toRadian4, _mm256_acos_ps( dotProduct4 ) );
			}

			for ( int i = 0; i < phaseTwo; i++ )
			{
				// if it's behind us, continue!
				if ( BlindspotAngleDeg > anglesToBoid[i] )
					continue;

				if ( phaseThree != i )
				{
					// move the data, same cache line
					directionToBoidX[phaseThree] = directionToBoidX[i];
					directionToBoidY[phaseThree] = directionToBoidY[i];
					directionToBoidZ[phaseThree] = directionToBoidZ[i];

					distanceToBoid[phaseThree] = distanceToBoid[i];

					// transfer our local copy of the bucket contents
					bucketVelocityX[phaseThree] = bucketVelocityX[i];
					bucketVelocityY[phaseThree] = bucketVelocityY[i];
					bucketVelocityZ[phaseThree] = bucketVelocityZ[i];

					bucketPositionX[phaseThree] = bucketPositionX[i];
					bucketPositionY[phaseThree] = bucketPositionY[i];
					bucketPositionZ[phaseThree] = bucketPositionZ[i];
				}

				phaseThree++;
			}
		}

		// phase 3:
		// erase the few elements that will be loaded
		// and processed by avx instructions, but should
		// not be used. For example, say of the total of
		// 16 boids only 13 are valid. Then we want to turn

		// [d, d, d, d, d, d, d, d, d, d, d, d, d, d, d, d]
		// into									  ----------
		// [d, d, d, d, d, d, d, d, d, d, d, d, d, 0, 0, 0]
		//										  ----------

		// where d is some data. This is to ensure that the last
		// feel elements can be safely loaded in with the vector
		// instructions, but do not influence the result.

		int elementsToErase = ( ( ELEMENTS_IN_BUCKET - phaseThree ) % SIMDSIZE ) * sizeof( float );
		if ( elementsToErase )
		{
			memset( directionToBoidX + phaseThree, 0, elementsToErase );
			memset( directionToBoidY + phaseThree, 0, elementsToErase );
			memset( directionToBoidZ + phaseThree, 0, elementsToErase );
			memset( bucketVelocityX + phaseThree, 0, elementsToErase );
			memset( bucketVelocityY + phaseThree, 0, elementsToErase );
			memset( bucketVelocityZ + phaseThree, 0, elementsToErase );
			memset( distanceToBoid + phaseThree, 0, elementsToErase );
			memset( bucketPositionX + phaseThree, 0, elementsToErase );
			memset( bucketPositionY + phaseThree, 0, elementsToErase );
			memset( bucketPositionZ + phaseThree, 0, elementsToErase );
		}

		// phase 4: update the SumVectors values with the remaining values

		s.count += phaseThree;
		int phaseThreeSimd = ( phaseThree + SIMDSIZE - 1 ) / SIMDSIZE;
		for ( int i = 0; i < phaseThreeSimd; i++ )
		{
			// --------------------------------------------------------------------------------
			// Compute the seperation factor (bootiful)

			// info about the _mm256_cmp_ps function
			// https://software.intel.com/sites/landingpage/IntrinsicsGuide/#text=_mm256_cmp_ps&expand=744
			// https://stackoverflow.com/questions/16988199/how-to-choose-avx-compare-predicate-variants
			// https://stackoverflow.com/questions/8627331/what-does-ordered-unordered-comparison-mean

			__m256 seperationFactor4 = _mm256_setzero_ps();
			switch ( SeparationType )
			{
			case DistanceType::LINEAR:
				seperationFactor4 = distanceToBoid4[i];
				break;

			case DistanceType::INVERSE_LINEAR:
				seperationFactor4 = _mm256_div_ps( ones, distanceToBoid4[i] );
				//                                    false case           true case
				seperationFactor4 = _mm256_blendv_ps( seperationFactor4, _mm256_setzero_ps(),
													  _mm256_cmp_ps( distanceToBoid4[i], _mm256_setzero_ps(), _CMP_EQ_OQ ) );
				break;

			case DistanceType::QUADRATIC:
				seperationFactor4 = _mm256_mul_ps( distanceToBoid4[i], distanceToBoid4[i] );
				break;

			case DistanceType::INVERSE_QUADRATIC:
				seperationFactor4 = _mm256_div_ps( ones, _mm256_mul_ps( distanceToBoid4[i], distanceToBoid4[i] ) );
				//                                    false case		   true case
				seperationFactor4 = _mm256_blendv_ps( seperationFactor4, _mm256_setzero_ps(),
													  _mm256_cmp_ps( distanceToBoid4[i], _mm256_setzero_ps(), _CMP_EQ_OQ ) );
				break;
			}

			// Compute everything!
			// suuuupparr fast! (or not)
			separationSumX4 = _mm256_add_ps( separationSumX4, _mm256_mul_ps( _mm256_sub_ps( _mm256_setzero_ps(), directionToBoidX4[i] ), seperationFactor4 ) );
			separationSumY4 = _mm256_add_ps( separationSumY4, _mm256_mul_ps( _mm256_sub_ps( _mm256_setzero_ps(), directionToBoidY4[i] ), seperationFactor4 ) );
			separationSumZ4 = _mm256_add_ps( separationSumZ4, _mm256_mul_ps( _mm256_sub_ps( _mm256_setzero_ps(), directionToBoidZ4[i] ), seperationFactor4 ) );

			headingSumX4 = _mm256_add_ps( headingSumX4, bucketVelocityX4[i] );
			headingSumY4 = _mm256_add_ps( headingSumY4, bucketVelocityY4[i] );
			headingSumZ4 = _mm256_add_ps( headingSumZ4, bucketVelocityZ4[i] );

			positionSumX4 = _mm256_add_ps( positionSumX4, bucketPositionX4[i] );
			positionSumY4 = _mm256_add_ps( positionSumY4, bucketPositionY4[i] );
			positionSumZ4 = _mm256_add_ps( positionSumZ4, bucketPositionZ4[i] );
		}
	}

	// todo: a lot of times 0 is added?
	// sum up all the results horizontally ( :| )
	for ( int i = 0, l = min( SIMDSIZE, s.count ); i < l; i++ )
	{
		s.separationSumX += separationSumX[i];
		s.separationSumY += separationSumY[i];
		s.separationSumZ += separationSumZ[i];

		s.headingSumX += headingSumX[i];
		s.headingSumY += headingSumY[i];
		s.headingSumZ += headingSumZ[i];

		s.positionSumX += positionSumX[i];
		s.positionSumY += positionSumY[i];
		s.positionSumZ += positionSumZ[i];
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
