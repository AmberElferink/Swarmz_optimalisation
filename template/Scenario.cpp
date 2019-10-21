#include "precomp.h"

void Scenario::Update( float dt )
{
	// construct the camera position vector
	Vec3 position = Vec3(
		camera_x - ( SCRWIDTH >> 1 ),
		camera_y - ( SCRHEIGHT >> 1 ),
		0.0f );

	// set the targets, update each boid
	swarm->SteeringTargets = targets;
	swarm->Update( dt );

	// manually update the position
	for ( Boid b : boids )
	{ b.Position += b.Velocity * dt; }
}

void Scenario::Draw( Surface *screen )
{
	screen->Clear( 0x000000 );

	// construct the camera position vector
	Vec3 position = Vec3(
		camera_x + ( SCRWIDTH >> 1 ),
		camera_y + ( SCRHEIGHT >> 1 ),
		0.0f );

	// find the total number of nearby boids
	int maximumNearbyBoids = 0;
	for ( Boid b : boids )
	{
		if ( b.numberOfNearbyBoids > maximumNearbyBoids )
			maximumNearbyBoids = b.numberOfNearbyBoids;
	}

	// draw each boid
	for ( Boid b : boids )
	{
		// determine the plot positions
		Vec3 plot_position = position + b.Position * camera_scale;
		Vec3 plot_velocity = plot_position + b.Velocity.Normalized() * draw_velocity_distance;
		Vec3 plot_acceleration = plot_velocity + b.Acceleration.Normalized() * draw_acceleration_distance;

		// plottin' dem
		Pixel factorColor = 0x00ff00;

		// deal with it
		float factor = ( (float)b.numberOfNearbyBoids / (maximumNearbyBoids + 1));
		Pixel output = ScaleColor( factorColor, (factor) * ( 255 ) );

		screen->PlotSafe( plot_acceleration.X, plot_acceleration.Y, boidAcceleration );
		screen->PlotSafe( plot_velocity.X, plot_velocity.Y, boidVelocity );
		screen->PlotSafe( plot_position.X, plot_position.Y, output );
	}

	// draw each target
	for ( Vec3 v : targets )
	{
		screen->Box(
			position.X + ( v.X - 4.0f ) * camera_scale,
			position.Y + ( v.Y - 4.0f ) * camera_scale,
			position.X + ( v.X + 4.0f ) * camera_scale,
			position.Y + ( v.Y + 4.0f ) * camera_scale,
			0xffffff );
	}
}

void Scenario::DrawVoxelDensity( Surface *screen )
{
	// requires some functionality to become public / available.
}

void Scenario::ChangeScale( float scale )
{

	// given that:
	//  - our zoom changes from 5 to 10
	//  - our viewport was 100 / 5 = 20 to 100 / 5 = 20, (20, 20)
	// and now is 100 / 10 = 10 by 100 / 10 = 10. (10, 10)
	//  - our camera is at (10, 10), which must move to (5, 5)

	// move the camera accordingly
	float factor = camera_scale / scale;
	camera_x *= factor;
	camera_y *= factor;

	printf( "Factor: %f \r\n", factor );

	// aaannnddd don't forget to update the scale
	camera_scale = scale;
}

void ScenarioRandom::Init( int count )
{
	for ( int j = 0; j < count; j++ )
	{
		Vec3 pos = Vec3(
			SPAWN_ORIGIN_X + ( RandomFloat() - 0.5f ) * ( SPAWN_WIDTH ),
			SPAWN_ORIGIN_Y + ( RandomFloat() - 0.5f ) * ( SPAWN_HEIGHT ),
			0.0f );

		Vec3 acc = Vec3(
			( RandomFloat() - 0.5f ) * 2,
			( RandomFloat() - 0.5f ) * 2,
			0.0f );

		Boid boid = Boid( pos, acc );
		boids.push_back( boid );
	}

	targets.push_back( Vec3( SPAWN_ORIGIN_X, SPAWN_ORIGIN_Y, 0 ) );
	swarm->SteeringWeight = 0.02f;
}

void ScenarioCircle::Init( int count )
{
	// x = r cos(t)    y = r sin(t)
	for ( int j = 0; j < count; j++ )
	{
		float t = ( (float)j / count ) * PI2;
		Vec3 pos = Vec3(
			SPAWN_ORIGIN_X + SPAWN_RADIUS * cos( t ),
			SPAWN_ORIGIN_Y + SPAWN_RADIUS * sin( t ),
			0 );
		Vec3 acc = Vec3(
			0.0f,
			0.0f,
			0.0f );

		Boid boid = Boid( pos, acc );
		boids.push_back( boid );
	}

	targets.push_back( Vec3( SPAWN_ORIGIN_X, SPAWN_ORIGIN_Y, 0 ) );
	swarm->SteeringWeight = 0.02f;
}
