#pragma once

class Graph
{

  private:
	// the timer and the time already spend (due to pausing).
	timer m_t;

  public:
	const char *m_name;
	bool m_showGraph;			   // Whether or not the graph should be shown.
	float m_current;			   // The current (next) value of the graph. Use Store() to store.
	bool m_graphDestroyed = false; // Whether or not the destructor has been called.
	int m_graphWidth;
	float *m_graphData; // The (past) data of the graph.
	int m_color;		// The color of the graph.
	float m_scaleMin;   // visible scale in the Y direction
	float m_scaleMax;   // visible scale in the X direction

	//FLT_MAX is the default value in plots for imgui
	Graph( const char *name, int width, int color, float scaleMin = FLT_MAX, float scaleMax = FLT_MAX );
	~Graph();

	void Start();
	void Stop();
	void StopAndStore();
	void Store();
};
extern std::vector<Graph *> graphs;

