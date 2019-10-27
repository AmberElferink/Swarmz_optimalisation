#include "precomp.h"

BucketPool::BucketPool( int n, int b )
{
	// initialize our state
	maximumBuckets = n;
	entriesInBucket = b;
	nextBucket = 0;

	// construct the buckets
	buckets.resize( n );
	for ( int j = 0; j < n; j++ )
		buckets[j] = ( new Bucket( b ) );
}

BucketPool::~BucketPool()
{
	for ( int j = 0; j < maximumBuckets; j++ )
		buckets[j]->~Bucket();

	buckets.clear();
}

Bucket *BucketPool::GetBucket( int i )
{
	// return the bucket
	// take note: the bucket may not exist.
	return buckets[i];
}

int BucketPool::ReserveBucket()
{
	if ( nextBucket >= maximumBuckets )
	{
		// the 'unsafe' version: crash. This is what happens on the GPU.
		throw "Bucket factory is empty.";
	}

	// store the id, prepare for the next cycle
	int bucketid = nextBucket;
	nextBucket++;

	// return the id
	return bucketid;
}

int BucketPool::ReserveBucketSafe()
{
	if ( nextBucket >= maximumBuckets )
	{
		// the 'safe' version: add a new bucket.
		maximumBuckets++;
		buckets.push_back( new Bucket( entriesInBucket ) );
	}

	// store the id, prepare for the next cycle
	int bucketid = nextBucket;
	nextBucket++;

	// return the id
	return bucketid;
}

void BucketPool::Clear()
{
	// clear 'dem out, men.
	for ( Bucket *bucket : buckets )
	{
		bucket->Clear();
	}

	printf( "Number of buckets in use: %i\r\n", nextBucket );
	nextBucket = 0;
}

