#include "precomp.h"

// represents the number of boids.
const int COUNT = 500;

Scenario *scenario;

//new graphs will automatically be added to graphs
vector<Graph *> graphs;
Graph g1( "boids loop", 100, 0x00FF0000, 0, 100 );
Graph g2( "boids draw loop", 100, 0x00FF0000, 0, 1 );

void DrawGUI()
{
	//ImGui::ShowDemoWindow();
	ImGui::Text( "Hello, world %d", 123 );

	if ( ImGui::CollapsingHeader( "Controls", ImGuiTreeNodeFlags_DefaultOpen ) )
	{
		float camera_scale = scenario->camera_scale;
		if (ImGui::SliderFloat("zoom level", &camera_scale, scenario->camera_scale_min, scenario->camera_scale_max))
		{
			scenario->ChangeScale( camera_scale );
		}
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
	g1.Start();
	scenario->Update( deltaTime * 0.01f );
	g2.Start();
	scenario->Draw( screen );

	g2.StopAndStore();
	g1.StopAndStore();

	DrawGUI();
}