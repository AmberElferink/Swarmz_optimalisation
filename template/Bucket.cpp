#include "precomp.h"
#include "Bucket.h"

template <typename T>
inline Bucket<T>::Bucket( const int n )
{
	ts = new T[n];
	maximum = n;
	count = 0;
}

template <typename T>
inline void Bucket<T>::Clear()
{
	// todo: call destructors on T?
	count = 0;
}

template <typename T>
inline bool Bucket<T>::CheckFull()
{
	return count >= maximum;
}