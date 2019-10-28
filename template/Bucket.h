#pragma once

class Bucket
{
  public:

	// represents the information the bucket can store.
	__declspec( align( 64 ) ) float *posX;
	__declspec( align( 64 ) ) float *posY;
	__declspec( align( 64 ) ) float *posZ;
	__declspec( align( 64 ) ) float *velX;
	__declspec( align( 64 ) ) float *velY;
	__declspec( align( 64 ) ) float *velZ;
	__declspec( align( 64 ) ) int *indx;

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
