#include "precomp.h"
#include "Bucket.h"

Bucket::Bucket( const int n )
{
	// the positions of the boids
	posX = new float[n];
	posY = new float[n];
	posZ = new float[n];

	// the velocities of the boids
	velX = new float[n];
	velY = new float[n];
	velZ = new float[n];

	// the indices of the boids
	indx = new int[n];

	maximum = n;
	count = 0;
}

Bucket::~Bucket()
{
	delete posX;
	delete posY;
	delete posZ;
	delete velX;
	delete velY;
	delete velZ;
	delete indx;
}

void Bucket::Clear()
{
	count = 0;
}

bool Bucket::CheckFull()
{
	return count >= maximum;
}