#pragma once

class BucketPool
{

  public:
	// represents the bucket size and the
	// pre-allocated number of buckets.
	BucketPool( int n, int b );

	// Clears out all the buckets.
	~BucketPool();

	// Retreives the bucket from the bucketpool, so
	// that you can it.
	Bucket *GetBucket( int i );

	// reserves a bucket and returns the bucket,
	// if no bucket is available the function
	// will crash.
	int ReserveBucket();

	// reserves a bucket and returns the bucket,
	// if no bucket is available a new bucket
	// is allocated.
	int ReserveBucketSafe();

	// frees all buckets. Take note: does not
	// remove the data within the buckets.
	void Clear();

  private:
	// represents the number of entries
	// within a bucket.
	int entriesInBucket;

	// represents the number of buckets in this
	// factory
	int maximumBuckets;

	// represents the next free bucket.
	int nextBucket;

	// represents all of the buckets.
	vector<Bucket *> buckets;
};