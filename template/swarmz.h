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
namespace sw
{
#define PI2 6.28318530717958647692
#define PI 3.14159265
#define toRadian 57.29577951308f
#define NUMBER_OF_ELEMENTS_IN_CELL 500
#define GRIDSIZE 20

#define indexToAcosRange 0.0078125f //this is 2/256. The acos table was filled with acos[i] = std::acos( ( 2.0f / 256.0f ) * (float) i - 1 ); \
									//this means you can calculate the lookup index by: int i = (acosinput + 1) /  indexToAcosRange.
const float acosTable[256] = {3.141593, 3.016511, 2.964585, 2.924661, 2.890937, 2.861166, 2.834198, 2.809348, 2.786171, 2.764360, 2.743688, 2.723987, 2.705124, 2.686994, 2.669514, 2.652613, 2.636232, 2.620323, 2.604842, 2.589755, 2.575028, 2.560635, 2.546551, 2.532753, 2.519224, 2.505945, 2.492901, 2.480078, 2.467462, 2.455043, 2.442809, 2.430750, 2.418859, 2.407125, 2.395542, 2.384102, 2.372799, 2.361627, 2.350579, 2.339651, 2.328837, 2.318133, 2.307534, 2.297036, 2.286634, 2.276326, 2.266108, 2.255976, 2.245928, 2.235960, 2.226068, 2.216252, 2.206508, 2.196833, 2.187225, 2.177683, 2.168203, 2.158784, 2.149423, 2.140120, 2.130872, 2.121677, 2.112533, 2.103440, 2.094395, 2.085397, 2.076445, 2.067537, 2.058671, 2.049848, 2.041064, 2.032320, 2.023613, 2.014943, 2.006309, 1.997709, 1.989143, 1.980609, 1.972107, 1.963635, 1.955193, 1.946780, 1.938394, 1.930036, 1.921704, 1.913397, 1.905114, 1.896856, 1.888620, 1.880407, 1.872215, 1.864044, 1.855893, 1.847761, 1.839648, 1.831554, 1.823477, 1.815416, 1.807372, 1.799343, 1.791330, 1.783330, 1.775345, 1.767372, 1.759413, 1.751465, 1.743529, 1.735604, 1.727689, 1.719784, 1.711889, 1.704002, 1.696124, 1.688254, 1.680391, 1.672534, 1.664684, 1.656840, 1.649001, 1.641167, 1.633337, 1.625511, 1.617689, 1.609869, 1.602051, 1.594236, 1.586422, 1.578609, 1.570796, 1.562984, 1.555171, 1.547357, 1.539541, 1.531724, 1.523904, 1.516082, 1.508256, 1.500426, 1.492592, 1.484753, 1.476908, 1.469058, 1.461202, 1.453339, 1.445469, 1.437590, 1.429704, 1.421808, 1.413903, 1.405989, 1.398064, 1.390128, 1.382180, 1.374220, 1.366248, 1.358262, 1.350263, 1.342249, 1.334221, 1.326177, 1.318116, 1.310039, 1.301944, 1.293831, 1.285700, 1.277549, 1.269378, 1.261186, 1.252973, 1.244737, 1.236478, 1.228196, 1.219889, 1.211557, 1.203198, 1.194813, 1.186400, 1.177958, 1.169486, 1.160984, 1.152450, 1.143884, 1.135284, 1.126650, 1.117980, 1.109273, 1.100529, 1.091745, 1.082921, 1.074056, 1.065148, 1.056195, 1.047198, 1.038153, 1.029059, 1.019916, 1.010721, 1.001473, 0.992169, 0.982809, 0.973390, 0.963910, 0.954367, 0.944760, 0.935085, 0.925341, 0.915524, 0.905633, 0.895665, 0.885616, 0.875484, 0.865266, 0.854958, 0.844557, 0.834059, 0.823460, 0.812756, 0.801942, 0.791014, 0.779966, 0.768794, 0.757491, 0.746051, 0.734468, 0.722734, 0.710842, 0.698784, 0.686550, 0.674130, 0.661515, 0.648692, 0.635647, 0.622369, 0.608839, 0.595042, 0.580958, 0.566564, 0.551838, 0.536750, 0.521270, 0.505360, 0.488980, 0.472079, 0.454598, 0.436469, 0.417606, 0.397904, 0.377233, 0.355421, 0.332245, 0.307395, 0.280426, 0.250656, 0.216931, 0.177008, 0.125082};

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
	//return std::sqrt( std::pow( X, 2 ) + std::pow( Y, 2 ) + std::pow( Z, 2 ) );
	float Length() const
	{
		return std::sqrtf( X * X + Y * Y + Z * Z );
	}

	float Length2() const
	{
		return X * X + Y * Y + Z * Z;
	}

	float DotProduct( const Vec3 &v ) const
	{
		return X * v.X + Y * v.Y + Z * v.Z;
	}

	//lookup value in lookup table. x should be between -1 and 1
	inline float ACOS( float x )
	{
		assert( x >= -1.1 && x <= 1.1 ); //there can be a bit error due to floats, and it is clamped.
		x = std::min( x, 1.0f );
		x = std::max( x, -1.0f );

		int i = ( x + 1 ) / indexToAcosRange;
		return acosTable[i];
	}

	// return static_cast<float>(std::acos(DotProduct(v) / (l1 * l2)) * 360 / PI2);
	float AngleTo( const Vec3 &v ) const
	{
		float l1 = Length();
		float l2 = v.Length();
		if ( l1 == 0 || l2 == 0 )
		{
			return 0;
		}
		return static_cast<float>( std::acos( DotProduct( v ) / ( l1 * l2 ) ) * toRadian );
	}

	float AngleToNorm( const Vec3 &v ) //faster if lengths had to be calculated in outer scope anyway
	{
		float dotpr = DotProduct( v );
		return ACOS( DotProduct( v ) ) * toRadian; // number is: 180 / PI
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

	int numberOfNearbyBoids;

	Boid()
	{
	}

	//Boid( const Boid &b ) = delete;

	explicit Boid( Vec3 pos, Vec3 vel ) : Position( pos ), Velocity( vel )
	{
		numberOfNearbyBoids = 0;
	}
};

struct NearbyBoid
{
	Boid boid;
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
	// represents the number of cells in the
	// given dimension.
	int nx, ny, nz;

	// represents the boundingbox of the grid.
	Vec3 minbb, maxbb;

	// represents the step sizes within the grid.
	Vec3 step;

	// regular constructor / deconstructor
	Grid( int nx, int ny, int nz );
	~Grid();

	// computes the index to use for within the grid.
	inline int CalculateGridCellIndex( int ix, int iy, int iz );

	// computes whether all the provided indices are within the grid.
	inline bool CheckInsideGrid( int ix, int iy, int iz );

	// computes the indices for all dimensions.
	void ComputeGridIndex( const Boid &b, int &ix, int &iy, int &iz );

	// constructs the grid. If no min / max
	// is provided the bounding box will be
	// computed dynamically.
	void ConstructGrid( const vector<Boid> &b, float perceptionRadius );

	// queries the grid, stores the result in
	// the out vector. Take note: reuse the vector.
	void QueryGrid( const Boid &b, const int r, vector<NearbyBoid> &out, float PerceptionRadius, float BlindspotAngleDeg, int ix, int iy, int iz );

	void DrawGrid( Surface *surface, Pixel density );

  private:
	GridCell *cells;

	// clears out the grid.
	void ClearGrid();

	// computes the bounding box, dynamically.
	void ComputeBoundingBox( const vector<Boid> &b, float perceptionRadius );

	// stores all boids in their corresponding cells.
	void StoreInCells( const vector<Boid> &b );
};

class Swarm
{
  public:
	Grid *grid;

	std::vector<Vec3> SteeringTargets;
	DistanceType SteeringTargetType = DistanceType::LINEAR;
	DistanceType SeparationType = DistanceType::INVERSE_QUADRATIC;

	float AlignmentWeight = 1;
	float CohesionWeight = 1;
	float SteeringWeight = 0.1f;
	float SeparationWeight = 1;

	float PerceptionRadius = 30.0f;
	float BlindspotAngleDeg = 20.0f;

	float MaxAcceleration = 10.0f;
	float MaxVelocity = 20.0f;

	explicit Swarm( std::vector<Boid> *entities ) : boids( entities )
	{
		std::random_device rd;
		eng = std::mt19937( rd() );
		grid = new Grid( GRIDSIZE, GRIDSIZE, GRIDSIZE );
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
		if ( PerceptionRadius == 0 )
			PerceptionRadius = 1;

		grid->ConstructGrid( ( *boids ), PerceptionRadius );

		for ( auto &b : *boids )
			updateBoid( b );
	}

  private:
	// not thread safe
	std::vector<NearbyBoid> vnb;
	std::mt19937 eng;

	// read-only thread safe
	std::vector<Boid> *boids;

	void updateBoid( Boid &b )
	{
		Vec3 separationSum;
		Vec3 headingSum;
		Vec3 positionSum;
		Vec3 po = b.Position;

		vnb.clear();
		getNearbyBoids( b, vnb );
		b.numberOfNearbyBoids = vnb.size();

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
			headingSum += closeBoid.boid.Velocity;
			positionSum += closeBoid.boid.Position;
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

	void getNearbyBoids( const Boid &b, std::vector<NearbyBoid> &vnb ) const
	{

		// retrieve the index
		int ix, iy, iz;
		grid->ComputeGridIndex( b, ix, iy, iz );

		// the offsets
		const int sx = 1;
		const int sy = 1;
		const int sz = 1;

		// loop over 'dem shizzles
		for ( int x = ix - sx, lx = ix + sx; x <= lx; x++ )
		{
			for ( int y = iy - sy, ly = iy + sy; y <= ly; y++ )
			{
				for ( int z = iz - sz, lz = iz + sz; z <= lz; z++ )
				{
					grid->QueryGrid( b, 0, vnb, PerceptionRadius, BlindspotAngleDeg, x, y, z );
				}
			}
		}
	}
	//void checkVoxelForBoids( const Boid &b, std::vector<NearbyBoid> &result, const Vec3 &voxelPos ) const
	//{
	//	auto iter = voxelCache.find( voxelPos );
	//	if ( iter != voxelCache.end() )
	//	{
	//		for ( Boid *test : iter->second )
	//		{
	//			const Vec3 &p1 = b.Position;
	//			const Vec3 &p2 = test->Position;

	//			Vec3 distanceVec = p2 - p1;
	//			float distance = distanceVec.Length();
	//			Vec3 distanceVecNorm = distanceVec / distance;

	//			Vec3 bNegVelocity = b.Velocity.Negative();
	//			float bNegVelocityLength = bNegVelocity.Length();

	//			float blindAngle = 0;
	//			if ( bNegVelocityLength > 0.000001f && distance > 0.000001f )
	//			{
	//				Vec3 bNegVelocityNorm = bNegVelocity / bNegVelocityLength;
	//				blindAngle = bNegVelocityNorm.AngleToNorm( distanceVecNorm );
	//			}

	//			if ( ( &b ) != test &&
	//				distance <= PerceptionRadius &&
	//				( BlindspotAngleDeg <= blindAngle || bNegVelocityLength == 0 ) )
	//			{
	//				NearbyBoid nb;
	//				nb.boid = test;
	//				nb.distance = distance;
	//				nb.direction = distanceVec;
	//				result.push_back( nb );
	//			}
	//		}
	//	}
	//}
	//void buildVoxelCache()
	//{
	//	voxelCache.clear();
	//	voxelCache.reserve( boids->size() );
	//	for ( auto &b : *boids )
	//	{
	//		voxelCache[getVoxelForBoid( b )].push_back( &b );
	//	}
	//}

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
