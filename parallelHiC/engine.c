#include <stdio.h> 
#include <string.h>       
#include <stdlib.h>
#include "utilities.h"
#include "engine.h"

//////////////////////////////////////////////////////////////////////////////////////
// coordinate functions
//////////////////////////////////////////////////////////////////////////////////////
int pwr( int x,int y ) { return y ? (y%2?x:1)*pwr(x*x,y>>1) : 1; }
int x  ( int m ) { return (m/L[0])%L[1]; }
int y  ( int m ) { return (m/L[1])%L[1]; }
int z  ( int m ) { return (m/L[2])%L[1]; }

/////////////////////////////////////////////////////////////////////////////////////////
// depth first search
/////////////////////////////////////////////////////////////////////////////////////////
void pseudoLegalMoves( int v ) {
	int i,j;
	for ( i = 0 ; i < 3 ; i++ ) {
		j = (v/L[i]) % L[1];
		if ( j == 0 )
			mL[fM[ply+1]++] = v+L[i];
		else if ( j == L[1]-1 )
			mL[fM[ply+1]++] = v-L[i];
		else {
			mL[fM[ply+1]++] = v+L[i];
			mL[fM[ply+1]++] = v-L[i];
			}
		}
	return;
}

void legalMoves( int v ) {
	int i,j,n;
	if ( v == L[3] ) {
		for ( i = 0 ; i < L[3] ; i++ )
			if ( ( searchPath.bounds[ply] == 0 ) | ( searchPath.bounds[ply] == boundType[i] ) )
				mL[fM[ply+1]++] = i;
		return;
		}
	for ( i = 0 ; i < 3 ; i++ ) {
		j = (v/L[i]) % L[1];
			if ( j == 0 ) {
				n = v+L[i];
				if ( ( searchPath.bounds[ply] == 0 ) | ( searchPath.bounds[ply] == boundType[n] ) )
					mL[fM[ply+1]++] = n;
			} else if ( j == L[1]-1 ) {
				n = v-L[i];
				if ( ( searchPath.bounds[ply] == 0 ) | ( searchPath.bounds[ply] == boundType[n] ) )
					mL[fM[ply+1]++] = n;
			} else {
				n = v+L[i];
				if ( ( searchPath.bounds[ply] == 0 ) | ( searchPath.bounds[ply] == boundType[n] ) )
					mL[fM[ply+1]++] = n;
				n = v-L[i];
				if ( ( searchPath.bounds[ply] == 0 ) | ( searchPath.bounds[ply] == boundType[n] ) )
					mL[fM[ply+1]++] = n;
			}
		}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////
// depth first search
/////////////////////////////////////////////////////////////////////////////////////////
void search ( int total ) {
	int i;
	fM[ply+1] = fM[ply];
	if (ply==searchPath.L) {
		(*savePtr)();
		return;
	}
	legalMoves(pos);
	for ( i = fM[ply] ; i < fM[ply+1] ; i++ ) {
		if ( arrayCon(searchPath.points,searchPath.L,mL[i]) >= 0 )
			continue;
		if ( cpn == total )
			return;
		searchPath.points[ply++] = mL[i];
		pos = mL[i];
		search(total);
		searchPath.points[--ply] = L[3];
		pos = L[3];
	}
}

void clearPly() {
	ply = 0;
	fM[0] = ply;
	fM[ply+1] = fM[ply];
}

/////////////////////////////////////////////////////////////////////////////////////////
// algorithms for shuffling an initial condition
/////////////////////////////////////////////////////////////////////////////////////////
void shuffle() { 
	int i,j,n = fM[ply+1] - fM[ply];
    if (n > 1) {
		for ( i = 0 ; i < n ; i++ ) {
          j = i + ( rand() / (RAND_MAX / (n - i) + 1) ) + fM[ply];
          int t = mL[j];
          mL[j] = mL[i+fM[ply]];
          mL[i+fM[ply]] = t;
        }
    }
}

void multiMansfield( int iterations ) {
	int i,j;
	int numPaths = arrayMax(searchPath.chains,0,searchPath.L) + 1;
	int *pathLength = (int*)malloc(sizeof(int)*numPaths);
	memset(pathLength, 0, sizeof(int) *numPaths);
	int **listedPath = (int **)malloc(sizeof(int*)*numPaths);
	for ( i = 0 ; i < numPaths ; i++ )
		listedPath[i]=(int*)malloc(sizeof(int)*searchPath.L);
	for ( i = 0 ; i < searchPath.L ; i++ )
		listedPath[searchPath.chains[i]][pathLength[searchPath.chains[i]]++] = searchPath.points[i];
	for ( i = 0 ; i < iterations ; i++ ) {
		clearPly();
		int pathGoto = -1;
		int pathFrom = rand() % numPaths;
		int branch = -1;
		if ( rand() % 2 )
			arrayRev(listedPath[pathFrom],0,pathLength[pathFrom]);
		legalMoves(listedPath[pathFrom][pathLength[pathFrom]-1]);
		shuffle();
		int choice = rand()%6;    
		if ( choice < fM[ply+1] ) {
			do { pathGoto++;
				 branch = arrayCon(listedPath[pathGoto],pathLength[pathGoto],mL[choice]);
			   } while ( branch == -1 );
			if ( pathFrom == pathGoto ) { // normal mansfield move
				 arrayRev(listedPath[pathFrom],0,pathLength[pathFrom]);
				 branch = arrayCon(listedPath[pathGoto],pathLength[pathGoto],mL[choice]);
				 arrayRev(listedPath[pathGoto],0,branch);
			} else {                      // multiple chain mansfield
				if ( rand() % 2 ) 
					arrayRev(listedPath[pathGoto],0,pathLength[pathGoto]);
				branch = arrayCon(listedPath[pathGoto],pathLength[pathGoto],mL[choice]);
				if ( ( branch > 0 ) & ( branch < pathLength[pathGoto]-1 ) ) {
					for ( j = branch ; j < pathLength[pathGoto] ; j++ )
						listedPath[pathFrom][pathLength[pathFrom]++] = listedPath[pathGoto][j];
					pathLength[pathGoto] = branch;
				}
			}
		}
	}
	int k = 0;
	for ( i = 0 ; i < numPaths; i++ )
		for ( j = 0 ; j < pathLength[i]; j++ ) {
			searchPath.points[k] = listedPath[i][j];
			searchPath.chains[k] = i;
			searchPath.bounds[k] = boundType[listedPath[i][j]];
			k++;
		}
	for( i = 0 ; i < numPaths; i++ )
		free(listedPath[i]);
	free(listedPath);
	free(pathLength);
}