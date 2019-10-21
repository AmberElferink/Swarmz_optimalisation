#include "precomp.h"

// represents the number of boids.
const int COUNT = 1000;

Scenario *scenario;

//new graphs will automatically be added to graphs
vector<Graph *> graphs;
Graph g0( "total loop", 100, 0x00ff0000, 0, 100 );
Graph g1( "boids loop", 100, 0x00FF0000, 0, 100 );
Graph g2( "boids draw loop", 100, 0x00FF0000, 0, 1 );

float min_global, max_global;

void DrawGUI()
{
	//ImGui::ShowDemoWindow();
	ImGui::Text( "Hello, world %d", 123 );

	if ( ImGui::CollapsingHeader( "Controls", ImGuiTreeNodeFlags_DefaultOpen ) )
	{
		float camera_scale = scenario->camera_scale;
		if ( ImGui::SliderFloat( "zoom level", &camera_scale, scenario->camera_scale_min, scenario->camera_scale_max ) )
		{
			scenario->ChangeScale( camera_scale );
		}
	}

	if ( ImGui::CollapsingHeader( "Statistics", ImGuiTreeNodeFlags_DefaultOpen ) )
	{
		float min, max, avg, std;
		min = minimum( g0.m_graphData, 100 );
		max = maximum( g0.m_graphData, 100 );
		avg = average( g0.m_graphData, 100 );
		std = stdev( g0.m_graphData, 100 );

		if ( min_global > min )
			min_global = min;

		if ( max_global < max )
			max_global = max;

		char buffer[40];
		snprintf( buffer, sizeof( buffer ), "Recent / overall performance measurements (in ms)" );
		ImGui::Text( buffer );

		snprintf( buffer, sizeof( buffer ), "min: %7.4f         /         min: %7.4f", min, min_global );
		ImGui::Text( buffer );

		snprintf( buffer, sizeof( buffer ), "max: %7.4f         /         max: %7.4f", max,  max_global );
		ImGui::Text( buffer );

		snprintf( buffer, sizeof( buffer ), "avg: %7.4f", avg );
		ImGui::Text( buffer );

		snprintf( buffer, sizeof( buffer ), "std: %7.4f", std );
		ImGui::Text( buffer );
	}

	//TreeNodeEx gives indent
	if ( ImGui::CollapsingHeader( "Graphs", ImGuiTreeNodeFlags_DefaultOpen ) )
	{
		for ( int i = 0; i < graphs.size(); i++ )
		{
			char gNr[3];
			snprintf( gNr, sizeof( gNr ), "g%i", i );
			//(sc = 100ms)
			char title[40];
			snprintf( title, sizeof( title ), "%s(scale=%.1fms)", graphs[i]->m_name, graphs[i]->m_scaleMax );
			if ( graphs[i]->m_showGraph )
				ImGui::PlotHistogram( gNr, graphs[i]->m_graphData, graphs[i]->m_graphWidth, 0, title, graphs[i]->m_scaleMin, graphs[i]->m_scaleMax, ImVec2( 0, 80 ) );
		}
	}
}

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	scenario = new ScenarioCircle();
	scenario->Init( COUNT );
}

void Game::KeyDown( SDL_Scancode key )
{
	SwitchScenario( key );
}

void Tmpl8::Game::SwitchScenario( SDL_Scancode key )
{
	switch ( key )
	{
	// key 1 = random
	case SDL_SCANCODE_1:
		scenario->~Scenario();
		scenario = new ScenarioRandom();
		scenario->Init( COUNT );
		break;

	// key 2 = circle
	case SDL_SCANCODE_2:
		scenario->~Scenario();
		scenario = new ScenarioCircle();
		scenario->Init( COUNT );
		break;
	}
}

void Game::MouseDown( int key )
{
	printf( "Mousekey pressed: %i\r\n", key );
}
// -----------------------------------------------------------
// Main application tick function
// -----------------------------------------------------------
void Game::Tick( float deltaTime )
{
	g0.Start();
	g1.Start();
	scenario->Update( deltaTime * 0.01f );
	g2.Start();
	scenario->Draw( screen );

	g2.StopAndStore();
	g1.StopAndStore();
	g0.StopAndStore();

	DrawGUI();
}