#include "precomp.h"
using namespace sw;

// represents the number of boids.
const int COUNT = 5000;

Scenario *scenario;

//new graphs will automatically be added to graphs
vector<Graph *> graphs;
Graph graphTotal( "total loop", 100, 0x00ff0000, 0, 600 );
Graph graphUpdate( "boids loop", 100, 0x00FF0000, 0, 600 );
Graph graphDraw( "boids draw loop", 100, 0x00FF0000, 0, 10 );

bool paused = false;
bool step = false;
float min_global, max_global;

//this was used to generate the hardcoded lookup table
void GenerateLookupTable()
{
	float acos[256];
	//input of acos is between -1 and 1.
	for ( int i = 0; i < 256; i++ ) //this is dividing the cirlce in 256 pieces instead of 360. It's faster since it's a power of two
	{
		acos[i] = std::acos( ( 2.0f / 256.0f ) * (float)i - 1 );
		printf( "%f, ", acos[i] );
	}
	int y = 0;
}

void DrawGUI()
{
	//ImGui::ShowDemoWindow();

	ImGui::NewLine();
	ImGui::Text( "Debug window - contains various controls to assist in debugging." );
	ImGui::NewLine();

	if ( ImGui::BeginTabBar( "##tabs", ImGuiTabBarFlags_None ) )
	{

		if ( ImGui::BeginTabItem( "Statistics" ) )
		{
			float min, max, avg, std;
			min = minimum( graphTotal.m_graphData, 100 );
			max = maximum( graphTotal.m_graphData, 100 );
			avg = average( graphTotal.m_graphData, 100 );
			std = stdev( graphTotal.m_graphData, 100 );

			if ( min_global > min )
				min_global = min;

			if ( max_global < max )
				max_global = max;

			char buffer[50];
			snprintf( buffer, sizeof( buffer ), "Performance measurements (in ms)" );
			ImGui::Text( buffer );

			float current = 0.50f * graphTotal.m_graphData[graphTotal.m_graphWidth - 1] +
							0.25f * graphTotal.m_graphData[graphTotal.m_graphWidth - 2] +
							0.15f * graphTotal.m_graphData[graphTotal.m_graphWidth - 3] +
							0.10f * graphTotal.m_graphData[graphTotal.m_graphWidth - 4];

			snprintf( buffer, sizeof( buffer ), "Current: %7.4f", current );
			ImGui::Text( buffer );

			snprintf( buffer, sizeof( buffer ), "Recent               /         Overall " );
			ImGui::Text( buffer );

			snprintf( buffer, sizeof( buffer ), "min: %7.4f         /         min: %7.4f", min, min_global );
			ImGui::Text( buffer );

			snprintf( buffer, sizeof( buffer ), "max: %7.4f         /         max: %7.4f", max, max_global );
			ImGui::Text( buffer );

			snprintf( buffer, sizeof( buffer ), "avg: %7.4f", avg );
			ImGui::Text( buffer );

			snprintf( buffer, sizeof( buffer ), "std: %7.4f", std );
			ImGui::Text( buffer );

			ImGui::EndTabItem();
		}

		if ( ImGui::BeginTabItem( "Controls" ) )
		{
			ImGui::NewLine();
			ImGui::Text( "Simulation controls" );

			ImGui::Checkbox( "Pause ", &paused );
			ImGui::SameLine( 120 );
			if ( ImGui::Button( "Step ", ImVec2( 200, 20 ) ) )
			{ step = true; }

			float camera_scale = scenario->camera_scale;
			if ( ImGui::SliderFloat( "zoom level", &camera_scale, scenario->camera_scale_min, scenario->camera_scale_max ) )
			{ scenario->ChangeScale( camera_scale ); }

			ImGui::NewLine();
			ImGui::Text( "Swarm controls" );

			ImGui::SliderFloat( "Weight alignment", &scenario->swarm->AlignmentWeight, 0, 10 );
			ImGui::SliderFloat( "Weight separation", &scenario->swarm->SeparationWeight, 0, 10 );
			ImGui::SliderFloat( "Weight cohesion", &scenario->swarm->CohesionWeight, 0, 10 );
			ImGui::SliderFloat( "Weight target", &scenario->swarm->SteeringWeight, 0, 10 );
			ImGui::SliderFloat( "Blindspotangle (degrees)", &scenario->swarm->BlindspotAngleDeg, 1, 180 );
			ImGui::SliderFloat( "Perception radius", &scenario->swarm->PerceptionRadius, 1, 100 );

			ImGui::NewLine();
			ImGui::Text( "Color controls" );

			ImVec4 expColor = ImGui::ColorConvertU32ToFloat4( scenario->expensiveColor );
			float *expColorPointer = (float *)&expColor;
			if ( ImGui::ColorEdit3( "Expensive color ", expColorPointer ) )
			{ scenario->expensiveColor = ImGui::ColorConvertFloat4ToU32( expColor ); }

			ImVec4 cheapColor = ImGui::ColorConvertU32ToFloat4( scenario->cheapColor );
			float *cheapColorPointer = (float *)&cheapColor;
			if ( ImGui::ColorEdit3( "Cheap color ", cheapColorPointer ) )
			{ scenario->cheapColor = ImGui::ColorConvertFloat4ToU32( cheapColor ); }

			ImVec4 accColor = ImGui::ColorConvertU32ToFloat4( scenario->boidAcceleration );
			float *accColorPointer = (float *)&accColor;
			if ( ImGui::ColorEdit3( "Acceleration color ", accColorPointer ) )
			{ scenario->boidAcceleration = ImGui::ColorConvertFloat4ToU32( accColor ); }

			ImVec4 velColor = ImGui::ColorConvertU32ToFloat4( scenario->boidVelocity );
			float *velColorPointer = (float *)&velColor;
			if ( ImGui::ColorEdit3( "Velocity color: ", velColorPointer ) )
			{ scenario->boidVelocity = ImGui::ColorConvertFloat4ToU32( velColor ); }

			ImGui::EndTabItem();
		}

		if ( ImGui::BeginTabItem( "Graphs" ) )
		{
			for ( int i = 0; i < graphs.size(); i++ )
			{
				char title[40];
				snprintf( title, sizeof( title ), "%s (scale = %.1fms)", graphs[i]->m_name, graphs[i]->m_scaleMax );
				if ( graphs[i]->m_showGraph && !graphs[i]->m_graphDestroyed )
				{
					ImGui::BeginGroup();
					ImGui::PlotHistogram( "", graphs[i]->m_graphData, graphs[i]->m_graphWidth, 0, title, graphs[i]->m_scaleMin, graphs[i]->m_scaleMax, ImVec2( 0, 80 ) );
					float min = minimum( graphs[i]->m_graphData, 100 );
					float max = maximum( graphs[i]->m_graphData, 100 );
					float avg = average( graphs[i]->m_graphData, 100 );
					float std = stdev( graphs[i]->m_graphData, 100 );

					char buffer[50];
					snprintf( buffer, sizeof( buffer ), "Statistics" );
					ImGui::Text( buffer );
					snprintf( buffer, sizeof( buffer ), "min: %7.4f", min );
					ImGui::Text( buffer );
					snprintf( buffer, sizeof( buffer ), "max: %7.4f", max );
					ImGui::Text( buffer );
					snprintf( buffer, sizeof( buffer ), "avg: %7.4f", avg );
					ImGui::Text( buffer );
					snprintf( buffer, sizeof( buffer ), "std: %7.4f", std );
					ImGui::Text( buffer );
					ImGui::EndGroup();
				}
			}

			ImGui::EndTabItem();
		}
		if ( ImGui::BeginTabItem( "Grid" ) )
		{
			ImGui::Text( "Todo" );

			ImGui::EndTabItem();
		}

		ImGui::EndTabBar();
	}
}

// -----------------------------------------------------------
// Initialize the application
// -----------------------------------------------------------
void Game::Init()
{
	union {
		float a[8];
		__m256 a4 = _mm256_set1_ps( 1.0f );
	};

	__m256 b = _mm256_set1_ps( 0.0f );
	__m256 d = _mm256_div_ps( a4, b );
	__m256 c = _mm256_blendv_ps( a4, b,
		_mm256_cmp_ps( d, _mm256_setzero_ps(), _CMP_EQ_OQ ) );

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

	graphTotal.Start();
	graphUpdate.Start();

	// if we're not paused
	if ( !paused )
	{
		graphUpdate.Start();
		scenario->Update( deltaTime * 0.01f );
		graphUpdate.StopAndStore();
	}

	// if we're paused, but we want to step
	if ( paused && step )
	{
		graphUpdate.Start();
		scenario->Update( deltaTime * 0.01f );
		step = false;
		graphUpdate.StopAndStore();
	}

	graphDraw.Start();
	scenario->Draw( screen );
	scenario->swarm->grid->DrawGrid( screen, 0xffffff );
	DrawGUI();
	graphDraw.StopAndStore();

	graphTotal.StopAndStore();
}