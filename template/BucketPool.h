#pragma once

template <typename T>
class BucketPool
{

  public:
	// represents the bucket size and the
	// pre-allocated number of buckets.
	BucketPool( int b, int n );

	// Retreives the bucket from the bucketpool, so 
	// that you can it.
	Bucket<T> *GetBucket( int i );

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
	void FreeBuckets();

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
	vector<Bucket<T>*> buckets;
};

template <typename T>
inline BucketPool<T>::BucketPool( int b, int n )
{
	// initialize our state
	entriesInBucket = b;
	maximumBuckets = n;

	// construct the buckets
	buckets.reserve( n );
	for ( int j = 0; j < n; j++ )
		buckets.push_back( new Bucket<T>( b ) );
}

template <typename T>
inline Bucket<T> *BucketPool<T>::GetBucket( int i )
{
	// return the bucket
	// take note: the bucket may not exist.
	return buckets[i];
}

template <typename T>
inline int BucketPool<T>::ReserveBucket()
{
	if ( nextBucket >= maximumBuckets )
	{
		// the 'unsafe' version: crash. This is what happens on the GPU.
		throw new error( "Bucket factory is empty." );
	}

	// store the id, prepare for the next cycle
	int bucketid = nextBucket;
	nextBucket++;

	// return the id
	return bucketid;
}

template <typename T>
inline int BucketPool<T>::ReserveBucketSafe()
{
	if ( nextBucket >= maximumBuckets )
	{
		// the 'safe' version: add a new bucket.
		maximumBuckets++;
		buckets.push_back( new Bucket<T>() )
	}

	// store the id, prepare for the next cycle
	int bucketid = nextBucket;
	nextBucket++;

	// return the id
	return bucketid;
}

template <typename T>
inline void BucketPool<T>::FreeBuckets()
{
	// clear 'dem out, men.
	for ( Bucket<T> *bucket : buckets )
	{
		bucket->Clear();
	}
}
