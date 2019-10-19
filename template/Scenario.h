#pragma once

class Scenario
{
  private:
	// what is privacy these days anyway?

  protected:
	// represents a few debugging colors
	Pixel boidPosition = 0xffffff;
	Pixel boidVelocity = 0xffff00;
	Pixel boidAcceleration = 0xff0000;

	// used for debugging purposes
	float draw_velocity_distance = 2.0f;
	float draw_acceleration_distance = 2.0f;

	// represents a basic camera
	float camera_x = 0.0f;
	float camera_y = 0.0f;
	float camera_speed = 5.0f;

	// represents the swarm
	Swarm *swarm;

  public:
	float camera_scale = 1.0f;
	float camera_scale_min = 0.1f;
	float camera_scale_max = 10.0f;
	float camera_scale_factor = 0.9f;

	Scenario()
	{
		swarm = new Swarm( &boids );
	}

	~Scenario()
	{
		targets.clear();
		boids.clear();
	}

	// represents the various targets of the boids.
	vector<Vec3> targets;

	// represents the boids.
	vector<Boid> boids;

	// prepares the scenario
	virtual void Init(int count) = 0;

	// makes the scenario interactive.
	void MouseDown( int key );
	void KeyDown( int key );

	// runs / draws the scenario. Is generally
	// the same, hence we implement it on this
	// level.
	void Update( float dt );
	void Draw( Surface *screen );
	void DrawVoxelDensity( Surface *screen );
};

class ScenarioRandom : public Scenario
{
  private:
	// represents the box in which boids can spawn.
	const int SPAWN_WIDTH = 800;
	const int SPAWN_HEIGHT = 500;
	const int SPAWN_ORIGIN_X = 400;
	const int SPAWN_ORIGIN_Y = 250;

  public:
	void Init( int count );
};

class ScenarioCircle : public Scenario
{
  private:
	// represents the box in which boids can spawn.
	const int SPAWN_RADIUS = 400;
	const int SPAWN_ORIGIN_X = 400;
	const int SPAWN_ORIGIN_Y = 250;

  public:
	void Init( int count );
};