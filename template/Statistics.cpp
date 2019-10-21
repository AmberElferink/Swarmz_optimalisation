#pragma once

#include "precomp.h"

float average( float *i, int n )
{
	float factor = 1.0f / n;
	float average = 0.0f;
	for ( int j = 0; j < n; j++ )
	{
		average += factor * i[j];
	}

	return average;
}

float stdev( float *i, int n )
{
	float avg = average( i, n );

	float stdev = 0.0f;
	float factor = 1.0f / n;
	for ( int j = 0; j < n; j++ )
	{
		stdev += factor * ( i[j] - avg ) * ( i[j] - avg );
	}

	return sqrtf( stdev );
}

float minimum( float *i, int n )
{
	float min = i[1];
	for ( int j = 1; j < n; j++ )
	{
		if ( i[j] < min && i[j] > 0 )
		{
			min = i[j];
		}
	}

	return min;
}

float maximum( float *i, int n )
{
	float max = i[1];
	for ( int j = 1; j < n; j++ )
	{
		if ( i[j] > max )
		{
			max = i[j];
		}
	}

	return max;
}