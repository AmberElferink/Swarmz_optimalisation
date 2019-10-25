#pragma once

template <typename T>
class BucketPool
{

  public:
	// represents the bucket size and the
	// pre-allocated number of buckets.
	BucketPool( int b, int n );

	// 
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
	//vector<*Bucket<T>> buckets;
};

template <typename T>
inline BucketPool<T>::BucketPool( int b, int n )
{
	entriesInBucket = b;
	maximumBuckets = n;

	buckets.reserve( n );

	for ( int j = 0; j < n; j++ )
		buckets.push_back( new Bucket<T>( b ) );
}

template <typename T>
inline Bucket<T> *BucketPool<T>::GetBucket( int i )
{
	return NULL;
}

template <typename T>
inline int BucketPool<T>::ReserveBucket()
{
	if ( nextBucket >= maximumBuckets )
	{
		throw new error( "Bucket factory is empty." );
	}

	int bucketid = nextBucket;
	nextBucket++;
}

template <typename T>
inline int BucketPool<T>::ReserveBucketSafe()
{
	if ( nextBucket >= maximumBuckets )
	{
		maximumBuckets++;
		buckets.push_back( new Bucket<T>() )
	}
}

template <typename T>
inline void BucketPool<T>::FreeBuckets()
{
}
