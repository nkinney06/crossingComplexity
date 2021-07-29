#include <stdlib.h>
#include <string.h>
#include <stdio.h> 
#include "utilities.h" 
#include "engine.h"
#include "library.h"
#include "shapes.h"
#include "initialization.h"
#include "contacts.h"
#include "outputs.h"

/////////////////////////////////////////////////////////////////////////////////////////
// three ways to make a new shape for the library and function to retrieve a shape
/////////////////////////////////////////////////////////////////////////////////////////
void joinShapes( int a, int b, int c ) {
	int i;
	int x=0,y=0;
	fstPnt[numShapes+1] = fstPnt[numShapes];
	for ( i = 0 ; i < searchPath.L ; i++ ) {
		if ( searchPath.region[i] == shapes.region[a] )
			shapes.points[fstPnt[numShapes+1]++] = shapes.points[fstPnt[a]+x++];
		if ( searchPath.region[i] == shapes.region[b] )
			shapes.points[fstPnt[numShapes+1]++] = shapes.points[fstPnt[b]+y++];
		}
	shapes.region[numShapes] = c;
	numShapes++;
}

void getShape( int s ) {
	int i;
	for ( i = fstPnt[s] ; i < fstPnt[s+1] ; i++ )
		searchPath.points[i-fstPnt[s]] = shapes.points[i];
	searchPath.L = fstPnt[s+1] - fstPnt[s];
}

/////////////////////////////////////////////////////////////////////////////////////////
// delete stuff from the library
/////////////////////////////////////////////////////////////////////////////////////////
void deleteShape( int n ) {
	int i,j;
	j = fstPnt[n+1] - fstPnt[n];
	for ( i = fstPnt[n] ; i < fstPnt[numShapes]-j ; i++ )
		shapes.points[i] = shapes.points[i+j];
	numShapes--;
	for ( i = n ; i < numShapes ; i++ ) {
		shapes.region[i] = shapes.region[i+1];
		fstPnt[i+1] = fstPnt[i+2]-j;
	}
}

void deleteRegion( int region ) {
	int shape = firstShape( region );
	while ( shape >= 0 ) {
		deleteShape( shape );
		shape = firstShape( region );
	}
}

void deleteLibrary() {
	while ( numShapes > 0 )
		deleteShape( numShapes - 1 );
}

/////////////////////////////////////////////////////////////////////////////////////////
// funding stuff in the shapes library
/////////////////////////////////////////////////////////////////////////////////////////
int firstShape( int region ) {
	int i;
	for ( i = 0 ; i < numShapes ; i++ )
		if ( shapes.region[i] == region )
			return i;
	return -1;
}

int lastShape( int region ) {
	int i,j=-1;
	for ( i = 0 ; i < numShapes ; i++ )
		if ( shapes.region[i] == region )
			j = i;
	return j;
}

int countRegion( int region ) {
	int i,j=0;
	for ( i = 0 ; i < numShapes ; i++ )
		if ( shapes.region[i] == region )
			j++;
	return j;
}

int maxRegion() {
	int i,j=0;
	for ( i = 0 ; i < numShapes ; i++ )
		if ( shapes.region[i] > j )
			j = shapes.region[i];
	return j;
}

int nullRegion() {
	int i;
	for ( i = 1 ; i < maxRegion() ; i++ )
		if ( countRegion(i) == 0 )
			return i;
	return -1;
}

/////////////////////////////////////////////////////////////////////////////////////////
// reset the region numbering
/////////////////////////////////////////////////////////////////////////////////////////
void reSet( int regionOne, int regionTwo ) {
	int i,j,k;
	i = nullRegion();
	while ( i > 0 ) {
		j = maxRegion();
		for ( k = 0 ; k < numShapes ; k++ )
			if ( shapes.region[k] == j )
				shapes.region[k] = i;
		for ( k = 0 ; k < searchPath.L ; k++ )
			if ( searchPath.region[k] == j )
				searchPath.region[k] = i;
		i = nullRegion();
	}
}