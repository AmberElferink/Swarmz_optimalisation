
#include "precomp.h"

// -------------------------------------------------------------------
// Non-static variables / functions
// -------------------------------------------------------------------

Graph::Graph( const char *name, int width, int color, float scaleMin, float scaleMax )
{
	// set all the variables.
	m_name = name;
	m_graphWidth = width;
	m_graphData = new float[width];
	m_color = color;
	m_current = 0;
	m_showGraph = true;
	m_scaleMin = scaleMin;
	m_scaleMax = scaleMax;

	// initialize the data array.
	memset( m_graphData, 0, width * sizeof( float ) );
	graphs.push_back( this );
}

Graph::~Graph()
{
	delete[] m_graphData;
	m_showGraph = false;
	m_graphDestroyed = true;
}

// starts the timer. Use 'pause' to stop the timer.
void Graph::Start()
{
	m_t.reset();
}

// pauses and stores the time between this pause and
// the last start.
void Graph::Stop()
{
	m_current += m_t.elapsed();
}

void Graph::StopAndStore()
{
	( *this ).Stop();
	( *this ).Store();
}

// Stores the timer value in the graph. Take note: always
// call 'pause' first.
void Graph::Store()
{
	// shift everything
	// memmove( m_graphData, m_graphData + 4, sizeof( float ) * ( (size_t)m_graphSize - 1 ) );
	for ( int j = 0; j < m_graphWidth - 1; j++ )
	{ m_graphData[j] = m_graphData[j + 1]; }

	// store the new value.
	m_graphData[m_graphWidth - 1] = m_current;
	m_current = 0;
}

