#pragma once

class Bucket
{
  public:

	// represents the information the bucket can store.
	float* posX;
	float* posY;
	float* posZ;
	float* velX;
	float* velY;
	float* velZ;
	int* indx;

	// represents the maximum number of
	// boids that this bucket can support.
	int maximum;
	
	// represents the number of boids
	// that are currently present.
	int count;

	// constructs a grid bucket with
	// the given size.
	Bucket( const int n );

	// hurrr, frees all the arrays.
	~Bucket();

	// clears the bucket. Take note,
	// does not actually erase the data.
	void Clear();

	// checks whether this bucket is full.
	bool CheckFull();

  private:
};
