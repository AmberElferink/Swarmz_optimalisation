#include "precomp.h"

void Scenario::MouseDown( int key )
{
	// 1 = left
	// 2 = middle
	// 3 = right

	switch ( key )
	{
	case 1:
		camera_scale = min( camera_scale_max, camera_scale * ( 2.0f - camera_scale_factor ) );
		break;

	case 3:
		camera_scale = max( camera_scale_min, camera_scale * camera_scale_factor );
		break;
	}
}

void Scenario::KeyDown( int key )
{
	// 4 = a
	// 22 = s
	// 7 = d
	// 26 = w

	switch ( key )
	{
	case 4:
		camera_x += camera_speed;
		break;

	case 7:
		camera_x -= camera_speed;
		break;

	case 22:
		camera_y -= camera_speed;
		break;

	case 26:
		camera_y += camera_speed;
		break;
	}
}

void Scenario::Update( float dt )
{
	// construct the camera position vector
	Vec3 position = Vec3( camera_x, camera_y, 0.0f );

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
	Vec3 position = Vec3( camera_x, camera_y, 0.0f );

	// draw each boid
	for ( Boid b : boids )
	{
		// determine the plot positions
		Vec3 plot_position = position + b.Position * camera_scale;
		Vec3 plot_velocity = plot_position + b.Velocity.Normalized() * draw_velocity_distance;
		Vec3 plot_acceleration = plot_velocity + b.Acceleration.Normalized() * draw_acceleration_distance;

		// plottin' dem
		screen->PlotSafe( plot_position.X, plot_position.Y, boidPosition );
		screen->PlotSafe( plot_velocity.X, plot_velocity.Y, boidVelocity );
		screen->PlotSafe( plot_acceleration.X, plot_acceleration.Y, boidAcceleration );
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
