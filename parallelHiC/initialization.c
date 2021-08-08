#include <stdio.h> 
#include <string.h>       
#include <stdlib.h>
#include "utilities.h"
#include "engine.h"
#include "initialization.h"
#include "contacts.h"
#include "library.h"
#include "shapes.h"

/////////////////////////////////////////////////////////////////////////////////////////
// save and unsave the searchPath
/////////////////////////////////////////////////////////////////////////////////////////
void saveToTarget() { 
	targetPath = searchPath; 
	cpn++; 
}
	
void unsaveTarget() { 
	searchPath = targetPath;
}
	
void saveToRegion() {
	regionPath = searchPath; 
}
	
void unsaveRegion() {
	searchPath = regionPath;
}

/////////////////////////////////////////////////////////////////////////////////////////
// sort chains in the search path by size
/////////////////////////////////////////////////////////////////////////////////////////
void sortChains() {
	int i,j;
	int numPaths = arrayMax(searchPath.chains,0,searchPath.L) + 1;
	int *pathLength = (int*)malloc(sizeof(int)*numPaths);
	int *sizeOrders = (int*)malloc(sizeof(int)*numPaths);
	memset(pathLength, 0, sizeof(int) *numPaths);
	memset(sizeOrders, 0, sizeof(int) *numPaths);
	int **listedPath = (int **)malloc(sizeof(int*)*numPaths);
	for ( i = 0 ; i < numPaths ; i++ )
		listedPath[i]=(int*)malloc(sizeof(int)*searchPath.L);
	for ( i = 0 ; i < searchPath.L ; i++ )
		listedPath[searchPath.chains[i]][pathLength[searchPath.chains[i]]++] = searchPath.points[i];
	for ( i = 0 ; i < numPaths ; i++ )
		sizeOrders[i] = i;
	
    // https://stackoverflow.com/questions/36714030/c-sort-float-array-while-keeping-track-of-indices
	int switched = 0;
	do {
		switched = 0;
		for(i = 1; i < numPaths; i++) {
			if( pathLength[sizeOrders[i-1]] < pathLength[sizeOrders[i]] ) {
				int temp = sizeOrders[i];
				sizeOrders[i] = sizeOrders[i-1];
				sizeOrders[i-1] = temp;
				switched = 1;
			}
		}
	} while(switched);
	
	int k = 0;
	for ( i = 0 ; i < numPaths; i++ )
		for ( j = 0 ; j < pathLength[sizeOrders[i]]; j++ ) {
			searchPath.points[k] = listedPath[sizeOrders[i]][j];
			searchPath.chains[k] = i;
			searchPath.bounds[k] = boundType[listedPath[sizeOrders[i]][j]];
			k++;
		}
	for( i = 0 ; i < numPaths; i++ )
		free(listedPath[i]);
	free(listedPath);
	free(pathLength);
	free(sizeOrders);
}

/////////////////////////////////////////////////////////////////////////////////////////
// divide the initial path and final path into seperate chains
/////////////////////////////////////////////////////////////////////////////////////////
void divide() {
	int i;
	for ( i = 0 ; i < searchPath.L ; i++ )
		searchPath.chains[i] = (i)/(((searchPath.L-1)+sumPaths)/sumPaths);
}

void redivide() {
	int i,j=0;
	searchPath.chains[0] = j;
	for ( i = 1 ; i < searchPath.L ; i++ ) {
		if ( neighbors[searchPath.points[i] * L[3] + searchPath.points[i-1]] == 0 )
			j++;
		searchPath.chains[i] = j;
	}
}
	
/////////////////////////////////////////////////////////////////////////////////////////
// initialization functions
/////////////////////////////////////////////////////////////////////////////////////////
void fillNeighbors() {
	int i,j;	
	memset(neighbors, 0, sizeof neighbors);
	for ( i = 0 ; i < L[3] ; i++ ) {
		clearPly();
		pseudoLegalMoves(i);
		for ( j = fM[ply] ; j < fM[ply+1] ; j++ )
			neighbors[ i * L[3] + mL[j] ] = 1;
	}
}

void fillLayerType() {
	int i,j,k;
	for ( i = 0 ; i < L[3] ; i++ ) {
		int layer = 1000;
		for ( j = 0 ; j <= 2 ; j++ ) {
			k = (i/L[j]) % L[1];
			if ( ( L[1] - 1 - k ) < layer )
				layer = ( L[1] - 1 - k );
			if ( k < layer )
				layer = k;
		}
		layerType[i] = layer;
	}
}

void fillBoundType() {
	int v,i,j,k;
	for ( v = 0 ; v < L[3] ; v++ ) {
		k = 0;
		for ( i = 0 ; i < 3 ; i++ ) {
			j = (v/L[i]) % L[1];
			if ( j == 0 )
				k += 1;
			else if ( j == L[1]-1 )
				k += 1;
			else
				k+= 2;
		}
		boundType[v] = k;
	}
}
	
void initialize() {
	int i;
	fillNeighbors();
	fillLayerType();
	fillBoundType();
	savePtr = &saveToTarget;
	cpn = 0;
	searchPath.L = L[3];
	pos = searchPath.L;
	clearPly();
	for ( i = 0 ; i < searchPath.L ; i++ )
		searchPath.points[i] = L[3];
	memset(searchPath.bounds, 0, sizeof searchPath.bounds);
}

int sizesMatched() {
	int i,sizesMatched = 0;
	for ( i = 0 ; i < sumPaths ; i++ )
		if ( mySizes[i] == arrayCnt(searchPath.chains, i, searchPath.L) )
			sizesMatched++;
	return sizesMatched;
}

void initConfig() {
	initialize();
	search(1);             
	unsaveTarget();
	divide();
	multiMansfield(30000);
	sortChains();
	int attempts = 0;
	if ( mySizes[0] > 0 ) { // user has specified sizes of each chain
		while ( sizesMatched() < sumPaths ) {
			multiMansfield(1);
			sortChains();
			attempts++;
			if ( attempts >= 50000 ) {
					initialize();
					search(1);             
					unsaveTarget();
					divide();
					multiMansfield(30000);
					sortChains();
					attempts = 0;
			}
		}
	}
	currentHiC();
	currentBoundary();
}

int scoreBoundary( int n ) {
	int i,score;
	score = n*L[3];
	for ( i = 0 ; i < L[3] ; i++ )
		score = score - abs( ( n-regionPath.matrix[i+n*L[3]] )-( n*(searchPath.bounds[i]/6) ) );
	return score/n;
}

void fitBoundary( FILE *out,int n,int max ) {
	int i=0;
	initialize();
	search(1);             
	unsaveTarget();
	divide();
	int bestScore = 0;
	printf("fitting boundary to score greater than %d\n",max);
	while ( bestScore < max ) {
		currentBoundary();
		if ( scoreBoundary( n ) > bestScore ) {
			bestScore = scoreBoundary( n );
			fprintf(out,"edge_score_%d\t%d\n",i++,bestScore);
			printf("edge_score_%d\t%d\n",i++,bestScore);
		}
		multiMansfield(1);
	}
}

void convolution( int convolve,char *matFile,char *bndFile, char *fldFile ) {
	int i,j;
	FILE *matr = fopen(matFile,"w");
	if (matr == NULL)
		exit(1);
	FILE *bndr = fopen(bndFile,"w");
	if (bndr == NULL)
		exit(1);
	FILE *fold = fopen(fldFile,"w");
	if (fold == NULL)
		exit(1);
	memset(regionPath.matrix, 0, sizeof regionPath.matrix);
	memset(regionPath.bounds, 0, sizeof regionPath.bounds);
	for ( i = 0 ; i < convolve ; i++ ) {
		initConfig();
		for ( j = 0 ; j < L[6] ; j++ )
			regionPath.matrix[j] = regionPath.matrix[j] + searchPath.matrix[j];
		for ( j = 0 ; j < L[3] ; j++ )
			regionPath.bounds[j] = regionPath.bounds[j] + ( searchPath.bounds[j] < 6 ? 1 : 0 ) ;
		for ( j = 0 ; j < L[3] ; j++ )
			fprintf(bndr,"%d%s",searchPath.bounds[j],( j%L[3] == L[3]-1 ) ? "\n" : " ");
	}
	for ( j = 0 ; j < L[6] ; j++ )
		fprintf(matr,"%d%s",regionPath.matrix[j],( j%L[3] == L[3]-1 ) ? "\n" : " ");
	for ( j = 0 ; j < L[3] ; j++ )
		fprintf(bndr,"%d%s",regionPath.bounds[j],( j%L[3] == L[3]-1 ) ? "\n" : " ");
	for ( i = 0 ; i < searchPath.L ; i++ )
		fprintf(fold,"%d %d %d %d\n",x(searchPath.points[i]),y(searchPath.points[i]),z(searchPath.points[i]),searchPath.chains[i]);
	fclose(bndr);
	fclose(matr);
	fclose(fold);
}

/////////////////////////////////////////////////////////////////////////////////////////
// boundary segmentation functions
/////////////////////////////////////////////////////////////////////////////////////////
int breakpoint( int v ) {
	int numPaths = arrayMax(searchPath.chains,0,searchPath.L) + 1;
	int i,j,sum;
	for ( i = 0 ; i < numPaths ; i++ ) {
		sum = 0;
		for ( j = 0 ; j <= i ; j++ )
			sum += mySizes[j];
		if ( v == sum )
			return 1;
	}
	return 0;
}

void sortSegments() {
	int i,j,k;
	int switched;
	do { // sort the segmented regions
		switched = 0;
		for ( i = 1 ; i <= arrayMax(searchPath.region,0,L[3])-1 ; i++ ) {
			for ( j = i+1 ; j <= arrayMax(searchPath.region,0,L[3]) ; j++ ) {
				if ( arrayFst(searchPath.region,i,L[3]) > arrayFst(searchPath.region,j,L[3]) ) {
					switched = 1;
					for ( k = 0 ; k < L[3] ; k++ )
						if ( searchPath.region[k] == i )
							searchPath.region[k] = -1;
					for ( k = 0 ; k < L[3] ; k++ )
						if ( searchPath.region[k] == j )
							searchPath.region[k] = i;
					for ( k = 0 ; k < L[3] ; k++ )
						if ( searchPath.region[k] == -1 )
							searchPath.region[k] = j;
				}
			}
		}
	} while ( switched == 1 );
}
		
void segmentBoundary( int maxSize ) {
	int i,j,k;
	int a,b,c;
	memset(searchPath.region, searchPath.bounds[0] == 6 ? -1 : 0, sizeof searchPath.region);
	for ( i = 0 ; i < L[3]-1 ; i++ ) {
		searchPath.region[i+1] = searchPath.region[i];
		if ( ( searchPath.bounds[i] == 6 ) & ( searchPath.bounds[i+1] == 6 ) )
			continue;
		else if ( ( searchPath.bounds[i] == 6 ) & ( searchPath.bounds[i+1] < 6 ) )	
			searchPath.region[i+1] = searchPath.region[i] + 1;
		else if ( ( searchPath.bounds[i+1] == 6 ) & ( searchPath.bounds[i] < 6 ) )	
			searchPath.region[i+1] = searchPath.region[i] + 1;
		else if ( searchPath.matrix[i*L[3]+(i+1)] == 0 )
			searchPath.region[i+1] = searchPath.region[i] + 2;
		else if ( breakpoint(i+1) == 1 )
			searchPath.region[i+1] = searchPath.region[i] + 2;
		else
			continue;
	}
	for ( i = 0 ; i < L[3] ; i++ )
		if ( searchPath.region[i]%2 == 0 )
			searchPath.region[i] = searchPath.region[i]/2 + 1;
		else
			searchPath.region[i] = 0;
	for ( i = 1 ; i < L[3] ; i++ ) {
		k = arrayCnt(searchPath.region,i,L[3]);
		if ( k > maxSize ) {
			a = arrayFst(searchPath.region,i,L[3]);
			b = arrayLst(searchPath.region,i,L[3]);
			c = arrayMax(searchPath.region,0,L[3]);
			int n_piece = (k + maxSize - 1)/maxSize;
			int newSize = (k + n_piece - 1)/n_piece;
			for ( j = a + newSize; j <= b ; j++ )
				searchPath.region[j] = c + 1;
		}
	}
	sortSegments();
}







