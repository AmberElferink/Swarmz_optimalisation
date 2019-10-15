#include "precomp.h"

// represents the number of boids.
const int COUNT = 1000;

// represents the box in which boids can spawn.
const int SPAWN_WIDTH = 800;
const int SPAWN_HEIGHT = 500;
const int SPAWN_ORIGIN_X = 400;
const int SPAWN_ORIGIN_Y = 250;

// represents a few debugging colors
const Pixel boidPosition = 0xffffff;
const Pixel boidVelocity = 0xffff00;
const Pixel boidAcceleration = 0xff0000;

// used for debugging purposes
float const draw_velocity_distance = 2.0f;
float const draw_acceleration_distance = 2.0f;

// represents a basic camera

float camera_x = 0.0f;
float camera_y = 0.0f;
float camera_speed = 5.0f;

float camera_scale = 1.0f;
float camera_scale_min = 0.1f;
float camera_scale_max = 10.0f;
float camera_scale_factor = 0.9f;

// represents the set of boids
vector<Boid> boids;
vector<Vec3> targets;

// represents the swarm
Swarm *swarm;

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	for ( int j = 0; j < COUNT; j++ )
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

	swarm = new Swarm( &boids );
	swarm->SteeringWeight = 0.01f;
}

void Game::KeyDown( int key )
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

void Game::MouseDown( int key )
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

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick( float deltaTime )
{
	float dt = deltaTime / 100.0f;
	screen->Clear( 0x000000 );

	// add the mouse as a target
	int mx, my;
	SDL_GetMouseState( &mx, &my );
	Vec3 vm = Vec3(
			mx + camera_x,
			my + camera_y,
			0 );

	screen->Box( vm.X - 4.0f, vm.Y - 4.0f, vm.X + 4.0f, vm.Y + 4.0f, 0xffffff );

	targets.clear();
	targets.push_back( vm );

	swarm->SteeringTargets = targets;

	// update each boid
	swarm->Update( dt );

	// construct the camera position vector
	Vec3 position = Vec3( camera_x, camera_y, 0.0f );

	//printf( "x:%f, y:%f", camera_x, camera_y );

	// try to draw each boid
	for ( Boid b : boids )
	{
		// update the position
		b.Position += b.Velocity * dt;

		// determine the plot positions
		Vec3 plot_position = position + b.Position * camera_scale;
		Vec3 plot_velocity = plot_position + b.Velocity.Normalized() * draw_velocity_distance;
		Vec3 plot_acceleration = plot_velocity + b.Acceleration.Normalized() * draw_acceleration_distance;

		screen->PlotSafe( plot_position.X, plot_position.Y, boidPosition );
		screen->PlotSafe( plot_velocity.X, plot_velocity.Y, boidVelocity );
		screen->PlotSafe( plot_acceleration.X, plot_acceleration.Y, boidAcceleration );
	}
}