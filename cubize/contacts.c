#include <stdio.h> 
#include <string.h>       
#include <stdlib.h>
#include "utilities.h"
#include "engine.h"

/////////////////////////////////////////////////////////////////////////////////////////
// functions related to various contact maps
/////////////////////////////////////////////////////////////////////////////////////////
void currentHiC() { 
	int i,j,k;     
	ply++;
	memset(searchPath.matrix, 0, sizeof searchPath.matrix);
	for ( i = 0; i < searchPath.L ; i++ ) {
		fM[ply+1] = fM[ply];
		pseudoLegalMoves(searchPath.points[i]);
		for ( j = fM[ply] ; j < fM[ply+1] ; j++ ) {
			k = arrayCon(searchPath.points,searchPath.L,mL[j]);
			if (k >= 0) {
				searchPath.matrix[i*searchPath.L+k] = 1;
				searchPath.matrix[k*searchPath.L+i] = 1;
			}
		}
	}
	ply--;
}

void getRegionPair ( int a, int b ) {
	int i,j;
	int bn = arrayCnt(searchPath.region,b,L[3]);
	memset(regionPath.matrix, 0, sizeof regionPath.matrix);
	int u = 0;
	for ( i = 0 ; i < L[3] ; i++ ) {
		if ( searchPath.region[i] == a ) {
			int v = 0;
			for ( j = 0 ; j < L[3] ; j++ ) {
				if ( searchPath.region[j] == b ) {
					regionPath.matrix[u*bn+v] = searchPath.matrix[i*L[3]+j];
					v++;
				}
			}
			u++;
		}
	}
}

void getRegion( int region ) {             
	int i,j,m,n;
	m = arrayFst(searchPath.region,region,L[3]);
	n = arrayLst(searchPath.region,region,L[3]);
	memset(regionPath.matrix, 0, sizeof regionPath.matrix);
	memset(regionPath.bounds, 0, sizeof regionPath.bounds);
	for ( i = m ; i <= n ; i++ )
		for ( j = m ; j <= n ; j++ )
			regionPath.matrix[ (i-m)*(n-m+1)+(j-m) ] = searchPath.matrix[ i*L[3]+j ];
	for ( i = m ; i <= n ; i++ ) {
		regionPath.bounds[ i - m ] = searchPath.bounds[i];
		regionPath.region[ i - m ] = region;
		}
	regionPath.L = n-m+1;
}

void bestRegionPair(int* row, int* col) {
	int i,j,k;
	regionPath.L = arrayMax(searchPath.region,0,searchPath.L);
	memset(regionPath.matrix, 0, sizeof regionPath.matrix);
	memset(regionPath.bounds, 0, sizeof regionPath.bounds);
	for ( i = 0 ; i < L[6] ; i++ )
		for ( j = 0; j < regionPath.L-1 ; j++ )
			for ( k = j + 1 ; k < regionPath.L ; k++ )
				if ( ( searchPath.region[i/L[3]] == j+1 ) & ( searchPath.region[i%L[3]] == k+1 ) )
					regionPath.matrix[regionPath.L*j+k]+=searchPath.matrix[i];
	k = arrayFst(regionPath.matrix,arrayMax(regionPath.matrix,0,regionPath.L*regionPath.L),regionPath.L*regionPath.L);
	*row = (k/regionPath.L)+1;
	*col = (k%regionPath.L)+1;
}

void currentBoundary() {
	int i;
	for ( i = 0; i < searchPath.L; i++ )
		searchPath.bounds[i] = boundType[searchPath.points[i]];
}