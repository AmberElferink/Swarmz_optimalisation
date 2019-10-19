#include "precomp.h"

// represents the number of boids.
const int COUNT = 400;

Scenario *scenario;

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	scenario = new ScenarioCircle();
	scenario->Init( COUNT );
}

void Game::KeyDown( int key )
{
	printf( "Keyboard pressed: %i\r\n", key );
	scenario->KeyDown( key );

	SwitchScenario( key );
}

void Tmpl8::Game::SwitchScenario( int key )
{
	switch ( key )
	{
	// key 0 = random
	case 30:
		scenario->~Scenario();
		scenario = new ScenarioRandom( );
		scenario->Init( COUNT );
		break;

	// key 1 = circle
	case 31:		
		scenario->~Scenario();
		scenario = new ScenarioCircle( );
		scenario->Init( COUNT );
		break;
	}
}

void Game::MouseDown( int key )
{
	printf( "Mousekey pressed: %i\r\n", key );
	scenario->MouseDown( key );
}

// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick( float deltaTime )
{
	scenario->Update( deltaTime * 0.01f );
	scenario->Draw( screen );
}