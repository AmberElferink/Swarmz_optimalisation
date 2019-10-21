/*
* Author: Michael Galetzka, 2016
*
* This is free and unencumbered software released into the public domain.
*
* Anyone is free to copy, modify, publish, use, compile, sell, or
* distribute this software, either in source code form or as a compiled
* binary, for any purpose, commercial or non-commercial, and by any
* means.
*
* In jurisdictions that recognize copyright laws, the author or authors
* of this software dedicate any and all copyright interest in the
* software to the public domain. We make this dedication for the benefit
* of the public at large and to the detriment of our heirs and
* successors. We intend this dedication to be an overt act of
* relinquishment in perpetuity of all present and future rights to this
* software under copyright law.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
* OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
* OTHER DEALINGS IN THE SOFTWARE.
*
* For more information, please refer to <http://unlicense.org>
*/

#pragma once

#include <cmath>
#include <random>
#include <unordered_map>
#include <vector>

#include <atomic>

#define PI2 6.28318530717958647692
#define NUMBER_OF_ELEMENTS_IN_CELL 50

namespace sw
{

enum class DistanceType
{
	LINEAR,
	INVERSE_LINEAR,
	QUADRATIC,
	INVERSE_QUADRATIC
};

struct Vec3
{
	float X;
	float Y;
	float Z;

	explicit Vec3( float x = 0, float y = 0, float z = 0 ) : X( x ), Y( y ), Z( z )
	{
	}

	float DistanceTo( const Vec3 &other ) const
	{
		return std::sqrt( std::pow( other.X - X, 2 ) + std::pow( other.Y - Y, 2 ) + std::pow( other.Z - Z, 2 ) );
	}

	float DistanceToSqr( const Vec3 &other ) const
	{
		return std::pow( other.X - X, 2 ) + std::pow( other.Y - Y, 2 ) + std::pow( other.Z - Z, 2 );
	}

	float Length() const
	{
		return std::sqrt( std::pow( X, 2 ) + std::pow( Y, 2 ) + std::pow( Z, 2 ) );
	}

	float DotProduct( const Vec3 &v ) const
	{
		return X * v.X + Y * v.Y + Z * v.Z;
	}

	float AngleTo( const Vec3 &v ) const
	{
		float l1 = Length();
		float l2 = v.Length();
		if ( l1 == 0 || l2 == 0 )
		{
			return 0;
		}
		return static_cast<float>( std::acos( DotProduct( v ) / ( l1 * l2 ) ) * 360 / PI2 );
	}

	Vec3 Normalized() const
	{
		float length = Length();
		if ( length == 0 )
		{
			return Vec3();
		}
		return Vec3( X / length, Y / length, Z / length );
	}

	Vec3 Negative() const
	{
		return Vec3( -X, -Y, -Z );
	}

	Vec3 &operator+=( const Vec3 &other )
	{
		X += other.X;
		Y += other.Y;
		Z += other.Z;
		return *this;
	}

	Vec3 operator/( float scalar ) const
	{
		return Vec3( X / scalar, Y / scalar, Z / scalar );
	}

	Vec3 operator*( float scalar ) const
	{
		return Vec3( X * scalar, Y * scalar, Z * scalar );
	}

	Vec3 operator+( const Vec3 &other ) const
	{
		return Vec3( X + other.X, Y + other.Y, Z + other.Z );
	}

	Vec3 operator-( const Vec3 &other ) const
	{
		return Vec3( X - other.X, Y - other.Y, Z - other.Z );
	}

	bool operator==( const Vec3 &rhs ) const
	{
		return X == rhs.X && Y == rhs.Y && Z == rhs.Z;
	}

	Vec3 ClampLength( float length ) const
	{
		float l = Length();
		if ( l > length )
		{
			return Normalized() * length;
		}
		return *this;
	}

	static Vec3 GetRandomUniform( std::mt19937 &engine )
	{
		std::uniform_real_distribution<float> thetaRange( 0.0f, PI2 );
		std::uniform_real_distribution<float> oneRange( 0, 1 );
		float theta = thetaRange( engine );
		float r = sqrt( oneRange( engine ) );
		float z = sqrt( 1.0f - r * r ) * ( oneRange( engine ) > 0.5f ? -1.0f : 1.0f );
		return Vec3( r * cos( theta ), r * sin( theta ), z );
	}
};

struct Vec3Hasher
{

	typedef std::size_t result_type;

	result_type operator()( sw::Vec3 const &v ) const
	{
		result_type const h1( std::hash<float>()( v.X ) );
		result_type const h2( std::hash<float>()( v.Y ) );
		result_type const h3( std::hash<float>()( v.Z ) );
		return ( h1 * 31 + h2 ) * 31 + h3;
	}
};

struct Boid
{
  public:
	Vec3 Position;
	Vec3 Velocity;
	Vec3 Acceleration;

	int id;

	Boid()
	{
		Position = Vec3( 0 );
		Velocity = Vec3( 0 );
		Acceleration = Vec3( 0 );

		id = 0;
	}

	explicit Boid( Vec3 pos, Vec3 vel ) : Position( pos ), Velocity( vel )
	{
		id = 0;
	}
};

struct NearbyBoid
{
	Boid *boid;
	Vec3 direction;
	float distance;
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
	void ConstructGrid( Vec3 min, Vec3 max, const vector<Boid> &b );

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

class Swarm
{
  public:
#ifdef DEBUG_PERFORMANCE
	Graph *tBuildVoxelCache = new Graph( "Build voxel cache", 100, 0x00FF0000, 0, 1.0f );
	Graph *tGetNearbyBoids = new Graph( "Get nearby Boids", 100, 0x00FF0000, 0, 100.0f );
#endif

	float PerceptionRadius = 30;

	float SeparationWeight = 120;
	DistanceType SeparationType = DistanceType::INVERSE_QUADRATIC;

	float AlignmentWeight = 1;
	float CohesionWeight = 1;

	float SteeringWeight = 0.1f;
	std::vector<Vec3> SteeringTargets;
	DistanceType SteeringTargetType = DistanceType::LINEAR;

	float BlindspotAngleDeg = 20.0f;
	float MaxAcceleration = 10.0f;
	float MaxVelocity = 20.0f;

	explicit Swarm( std::vector<Boid> *entities ) : boids( entities )
	{
		std::random_device rd;
		eng = std::mt19937( rd() );

		grid = new Grid( 20, 20, 20 );
	}

	void Update( float delta )
	{
		UpdateAcceleration();

		for ( auto &b : *boids )
		{
			b.Velocity = ( b.Velocity + b.Acceleration * delta ).ClampLength( MaxVelocity );
			b.Position += b.Velocity * delta;
		}
	}

	void UpdateAcceleration()
	{
		grid->ConstructGrid( (*boids) );

		if ( PerceptionRadius == 0 )
		{
			PerceptionRadius = 1;
		}

#ifdef DEBUG_PERFORMANCE
		tBuildVoxelCache->Start();
#endif
		buildVoxelCache();
#ifdef DEBUG_PERFORMANCE
		tBuildVoxelCache->StopAndStore();
#endif

		for ( auto &b : *boids )
		{
			updateBoid( b );
		}

#ifdef DEBUG_PERFORMANCE
		tGetNearbyBoids->Store();
#endif
	}

  private:
	Grid* grid;
	std::vector<Boid> *boids;
	std::unordered_map<Vec3, std::vector<Boid *>, Vec3Hasher> voxelCache;
	std::mt19937 eng;

	void updateBoid( Boid &b )
	{
		Vec3 separationSum;
		Vec3 headingSum;
		Vec3 positionSum;
		Vec3 po = b.Position;

#ifdef DEBUG_PERFORMANCE
		tGetNearbyBoids->Start();
#endif

		auto vnb = getNearbyBoids( b );


#ifdef DEBUG_PERFORMANCE
		tGetNearbyBoids->Stop();
#endif

		for ( NearbyBoid &closeBoid : vnb )
		{
			if ( closeBoid.distance == 0 )
			{
				separationSum += Vec3::GetRandomUniform( eng ) * 1000;
			}
			else
			{
				float separationFactor = TransformDistance( closeBoid.distance, SeparationType );
				separationSum += closeBoid.direction.Negative() * separationFactor;
			}
			headingSum += closeBoid.boid->Velocity;
			positionSum += closeBoid.boid->Position;
		}

		Vec3 steeringTarget = b.Position;
		float targetDistance = -1;
		for ( auto &target : SteeringTargets )
		{
			float distance = TransformDistance( target.DistanceTo( b.Position ), SteeringTargetType );
			if ( targetDistance < 0 || distance < targetDistance )
			{
				steeringTarget = target;
				targetDistance = distance;
			}
		}

		// Separation: steer to avoid crowding local flockmates
		Vec3 separation = vnb.size() > 0 ? separationSum / vnb.size() : separationSum;

		// Alignment: steer towards the average heading of local flockmates
		Vec3 alignment = vnb.size() > 0 ? headingSum / vnb.size() : headingSum;

		// Cohesion: steer to move toward the average position of local flockmates
		Vec3 avgPosition = vnb.size() > 0 ? positionSum / vnb.size() : b.Position;
		Vec3 cohesion = avgPosition - b.Position;

		// Steering: steer towards the nearest target location (like a moth to the light)
		Vec3 steering = ( steeringTarget - b.Position ).Normalized() * targetDistance;

		// calculate boid acceleration
		Vec3 acceleration;
		acceleration += separation * SeparationWeight;
		acceleration += alignment * AlignmentWeight;
		acceleration += cohesion * CohesionWeight;
		acceleration += steering * SteeringWeight;
		b.Acceleration = acceleration.ClampLength( MaxAcceleration );
	}

	std::vector<NearbyBoid> getNearbyBoids( const Boid &b ) const
	{
		vector<NearbyBoid> vnb;

		// retrieve the index
		int ix, iy, iz;
		grid->ComputeGridIndex( b, ix, iy, iz );

		// loop over 'dem shizzles
		for ( int x = ix -1; x < ix +1; x++ )
		{
			for ( int y = iy -1; y < iy + 1; y++ )
			{
				for ( int z = iz -1; z < iz + 1; z++ )
				{
					grid->QueryGrid( b, 0, vnb, PerceptionRadius, BlindspotAngleDeg, x,y,z );
				}
			}
		}

		return vnb;
	}

	void checkVoxelForBoids( const Boid &b, std::vector<NearbyBoid> &result, const Vec3 &voxelPos ) const
	{
		auto iter = voxelCache.find( voxelPos );
		if ( iter != voxelCache.end() )
		{
			for ( Boid *test : iter->second )
			{
				const Vec3 &p1 = b.Position;
				const Vec3 &p2 = test->Position;
				Vec3 vec = p2 - p1;
				float distance = vec.Length();
				float blindAngle = b.Velocity.Negative().AngleTo( vec );
				if ( ( &b ) != test && distance <= PerceptionRadius && ( BlindspotAngleDeg <= blindAngle || b.Velocity.Length() == 0 ) )
				{
					NearbyBoid nb;
					nb.boid = test;
					nb.distance = distance;
					nb.direction = vec;
					result.push_back( nb );
				}
			}
		}
	}

	void buildVoxelCache()
	{
		voxelCache.clear();
		voxelCache.reserve( boids->size() );
		for ( auto &b : *boids )
		{
			voxelCache[getVoxelForBoid( b )].push_back( &b );
		}
	}

	Vec3 getVoxelForBoid( const Boid &b ) const
	{
		float r = std::abs( PerceptionRadius );
		const Vec3 &p = b.Position;
		Vec3 voxelPos;
		voxelPos.X = static_cast<int>( p.X / r );
		voxelPos.Y = static_cast<int>( p.Y / r );
		voxelPos.Z = static_cast<int>( p.Z / r );
		return voxelPos;
	}

	float TransformDistance( float distance, DistanceType type )
	{
		if ( type == DistanceType::LINEAR )
		{
			return distance;
		}
		else if ( type == DistanceType::INVERSE_LINEAR )
		{
			return distance == 0 ? 0 : 1 / distance;
		}
		else if ( type == DistanceType::QUADRATIC )
		{
			return std::pow( distance, 2 );
		}
		else if ( type == DistanceType::INVERSE_QUADRATIC )
		{
			float quad = std::pow( distance, 2 );
			return quad == 0 ? 0 : 1 / quad;
		}
		else
		{
			return distance; // throw exception instead?
		}
	}
};
} // namespace sw
