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

vector<Graph *> graphs;
Graph g1( "boids loop", 100, 0x00FF0000, 0, 100 );
Graph g2( "boids draw loop", 100, 0x00FF0000, 0, 1 );

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	graphs.push_back( &g1 );
	graphs.push_back( &g2 );

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

void Game::KeyDown( SDL_Scancode key )
{
	// 4 = a
	// 22 = s
	// 7 = d
	// 26 = w

	switch ( key )
	{
	case SDL_SCANCODE_A:
		camera_x += camera_speed;
		break;

	case SDL_SCANCODE_D:
		camera_x -= camera_speed;
		break;

	case SDL_SCANCODE_S:
		camera_y -= camera_speed;
		break;

	case SDL_SCANCODE_W:
		camera_y += camera_speed;
		break;
	}
}

void Game::MouseDown( int key )
{
	//	// 1 = left
	//	// 2 = middle
	//	// 3 = right
	//
	//	switch ( key )
	//	{
	//	case 1:
	//		camera_scale = min( camera_scale_max, camera_scale * ( 2.0f - camera_scale_factor ) );
	//		break;
	//
	//	case 3:
	//		camera_scale = max( camera_scale_min, camera_scale * camera_scale_factor );
	//		break;
	//	}
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

	screen->Box( int( vm.X - 4.0f ), int( vm.Y - 4.0f ), int( vm.X + 4.0f ), int( vm.Y + 4.0f ), 0xffffff );

	targets.clear();
	targets.push_back( vm );

	g1.Start();
	swarm->SteeringTargets = targets;

	// update each boid
	swarm->Update( dt );

	// construct the camera position vector
	Vec3 position = Vec3( camera_x, camera_y, 0.0f );

	//printf( "x:%f, y:%f", camera_x, camera_y );

	g2.Start();
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
	g2.StopAndStore();
	g1.StopAndStore();

	//ImGui::ShowDemoWindow();
	ImGui::Text( "Hello, world %d", 123 );

	if ( ImGui::CollapsingHeader( "Controls", ImGuiTreeNodeFlags_DefaultOpen) )
	{
		
		ImGui::SliderFloat( "zoom level", &camera_scale, camera_scale_min, camera_scale_max );
	}
	//TreeNodeEx gives indent 
	if ( ImGui::CollapsingHeader( "Graphs", ImGuiTreeNodeFlags_DefaultOpen ) )
	{
		for (int i = 0; i < graphs.size(); i++)
		{
			char gNr[3];
			snprintf( gNr, sizeof( gNr ), "g%i", i );
			//(sc = 100ms)
			char title[40];
			snprintf( title, sizeof( title ), "%s(scale=%.1fms)", graphs[i]->m_name, graphs[i]->m_scaleMax );
			if ( graphs[i]->m_showGraph )
				ImGui::PlotHistogram( gNr , graphs[i]->m_graphData, graphs[i]->m_graphWidth, 0, title, graphs[i]->m_scaleMin, graphs[i]->m_scaleMax, ImVec2( 0, 80 ) );
		}
	}
}