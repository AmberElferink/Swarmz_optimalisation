#pragma once

class Bucket
{
  public:

	// represents the information the bucket can store.
	__declspec( align( 64 ) ) union {
		float posX[ELEMENTS_IN_BUCKET];
		__m256 posX4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};
	__declspec( align( 64 ) ) union {
		float posY[ELEMENTS_IN_BUCKET];
		__m256 posY4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};
	__declspec( align( 64 ) ) union {
		float posZ[ELEMENTS_IN_BUCKET];
		__m256 posZ4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};
	__declspec( align( 64 ) ) union {
		float velX[ELEMENTS_IN_BUCKET];
		__m256 velX4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};
	__declspec( align( 64 ) ) union {
		float velY[ELEMENTS_IN_BUCKET];
		__m256 velY4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};
	__declspec( align( 64 ) ) union {
		float velZ[ELEMENTS_IN_BUCKET];
		__m256 velZ4[ELEMENTS_IN_BUCKET / SIMDSIZE];
	};
	__declspec( align( 64 ) ) union {
		int indx[ELEMENTS_IN_BUCKET];
		__m256i indx4[ELEMENTS_IN_BUCKET / SIMDSIZE];
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
