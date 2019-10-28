#pragma once

class Bucket
{
  public:

	// represents the information the bucket can store.
	__declspec( align( 64 ) ) union {
		float *posX;
		__m128 *posX4;
	};
	__declspec( align( 64 ) ) union {
		float *posY;
		__m128 *posY4;
	};
	__declspec( align( 64 ) ) union {
		float *posZ;
		__m128 *posZ4;
	};
	__declspec( align( 64 ) ) union {
		float *velX;
		__m128 *velX4;
	};
	__declspec( align( 64 ) ) union {
		float *velY;
		__m128 *velY4;
	};
	__declspec( align( 64 ) ) union {
		float *velZ;
		__m128 *velZ4;
	};
	__declspec( align( 64 ) ) union {
		int *indx;
		__m128 *indx4;
	};

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
