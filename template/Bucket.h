#pragma once

class Bucket
{
  public:

	// represents the information the bucket can store.
	__declspec( align( 64 ) ) union {
		float *posX;
		__m256 *posX4;
	};
	__declspec( align( 64 ) ) union {
		float *posY;
		__m256 *posY4;
	};
	__declspec( align( 64 ) ) union {
		float *posZ;
		__m256 *posZ4;
	};
	__declspec( align( 64 ) ) union {
		float *velX;
		__m256 *velX4;
	};
	__declspec( align( 64 ) ) union {
		float *velY;
		__m256 *velY4;
	};
	__declspec( align( 64 ) ) union {
		float *velZ;
		__m256 *velZ4;
	};
	__declspec( align( 64 ) ) union {
		int *indx;
		__m256i *indx4;
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
