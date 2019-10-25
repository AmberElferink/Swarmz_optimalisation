#pragma once

template <typename T>
class Bucket
{
  public:
	// represents the boids within this
	// bucket.
	T *ts;

	// represents the maximum number of
	// boids that this bucket can support.
	int maximum;

	// represents the number of boids
	// that are currently present.
	int count;

	// constructs a grid bucket with
	// the given size.
	Bucket( const int n );

	// clears the bucket. Take note,
	// does not actually erase the data.
	inline void Clear();

	// checks whether this bucket is full.
	inline bool CheckFull();

  private:
};
