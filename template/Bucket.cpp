#include "precomp.h"
#include "Bucket.h"

Bucket::Bucket( const int n )
{
	maximum = ELEMENTS_IN_BUCKET;
	count = 0;
}

Bucket::~Bucket()
{
	delete posX;
	delete posY;
	delete posZ;
	delete velX;
	delete velY;
	delete velZ;
	delete indx;
}

void Bucket::Clear()
{
	count = 0;
}

bool Bucket::CheckFull()
{
	return count >= maximum;
}