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
#include <omp.h>

/////////////////////////////////////////////////////////////////////////////////////////
// functions for evaluating shapes in the library
/////////////////////////////////////////////////////////////////////////////////////////
void computeScoreMatrixSimple(int numShapesA,    int numShapesB,
							  int *shapes_a,     int *fstPnt_a,
							  int *shapes_b,     int *fstPnt_b,
							  int *contacts,     int *scores,
							  int evalFun ) {
	int a,b,i,j;	
	for ( a = 0 ; a < numShapesA ; a++ ) {
		#pragma omp parallel for
		for ( b = 0 ; b < numShapesB ; b++ ) {
			int score_for_pairs=0;
			int pointOne,pointTwo;
			for ( i = fstPnt_a[a] ; i < fstPnt_a[a+1] ; i++ ) {
				if ( score_for_pairs >= 1000 )
					break;
				pointOne = shapes_a[i];
				for ( j = fstPnt_b[b] ; j < fstPnt_b[b+1] ; j++ ) {
					pointTwo = shapes_b[j];
					if ( pointTwo == pointOne ) {
						score_for_pairs = 1000;
						break;
					}
					if ( evalFun == 2 ) {
						if ( neighbors[ pointOne * L[3] + pointTwo ] == 1 )
							if ( contacts[(i-fstPnt_a[a])*fstPnt_b[1] + (j-fstPnt_b[b])] == 0 )
								score_for_pairs += 2;
					} else {
						if ( neighbors[ pointOne * L[3] + pointTwo ] != contacts[(i-fstPnt_a[a])*fstPnt_b[1] + (j-fstPnt_b[b])] )
							score_for_pairs += 2;
					}
				}
			}
			scores[a*numShapesB+b] += score_for_pairs;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// shape generation functions
/////////////////////////////////////////////////////////////////////////////////////////
void joinAllShapesSimple( int a, int b ) {                     
	int i,j;
	int aFrom = firstShape(a);
	int bFrom = firstShape(b);
	int aGoto = lastShape(a);
	int bGoto = lastShape(b);
	int aSize = arrayCnt(searchPath.region,a,L[3]);
	int bSize = arrayCnt(searchPath.region,b,L[3]);
	int aTotl = aGoto-aFrom+1;
	int bTotl = bGoto-bFrom+1;
	int *scores   = (int *)malloc(aTotl*bTotl*sizeof(int));
	int *shapes_a = (int *)malloc(aTotl*aSize*sizeof(int));
	int *shapes_b = (int *)malloc(bTotl*bSize*sizeof(int));
	int *fstPnt_a = (int *)malloc((aTotl + 1)*sizeof(int));
	int *fstPnt_b = (int *)malloc((bTotl + 1)*sizeof(int));
	int *contacts;
	
	memset(scores,0,aTotl*bTotl*sizeof(int));	
	memcpy(shapes_a,shapes.points+fstPnt[aFrom],(fstPnt[aGoto+1]-fstPnt[aFrom])*sizeof(int));
	memcpy(shapes_b,shapes.points+fstPnt[bFrom],(fstPnt[bGoto+1]-fstPnt[bFrom])*sizeof(int));
	for ( i = aFrom ; i <= aGoto+1 ; i++ )
		fstPnt_a[i-aFrom] = fstPnt[i]-fstPnt[aFrom];
	for ( i = bFrom ; i <= bGoto+1 ; i++ )
		fstPnt_b[i-bFrom] = fstPnt[i]-fstPnt[bFrom];

	contacts = (int *)malloc(aSize*bSize*sizeof(int));
	getRegionPair(a,b);  
	memcpy(contacts,regionPath.matrix,aSize*bSize*sizeof(int));
	computeScoreMatrixSimple(aTotl,bTotl,shapes_a,fstPnt_a,shapes_b,fstPnt_b,contacts,scores,myEval);
	free(contacts);
	
	int k = arrayMax(searchPath.region,0,searchPath.L) + 1;
	for ( i = aFrom; i <= aGoto ; i++ ) 
		for ( j = bFrom; j <= bGoto ; j++ )
			if ( scores[(i-aFrom)*(bTotl)+(j-bFrom)] <= maxMisMatch )
				joinShapes( i, j, k );
	for ( i = 0 ; i < searchPath.L ; i++ )
		if ( ( searchPath.region[i] == a ) | ( searchPath.region[i] == b ) )
			searchPath.region[i] = k;
	
	free(scores);
	free(shapes_a);
	free(shapes_b);
	free(fstPnt_a);
	free(fstPnt_b);
}

/////////////////////////////////////////////////////////////////////////////////////////
// functions to populate library
/////////////////////////////////////////////////////////////////////////////////////////
int evalSearchPath() {
	int i,j;
	int score_for_shape = 0;
	for ( i = 0 ; i < searchPath.L-1 ; i++ ) {
		for ( j = i+1 ; j < searchPath.L ; j++ ) {
			if ( searchPath.points[i] == searchPath.points[j] )
				return 1000;
			if ( myEval == 2 ) {
				if ( ( neighbors[ searchPath.points[i] * L[3] + searchPath.points[j] ] == 1 ) & ( regionPath.matrix[ i * searchPath.L + j ] == 0 ) )
					score_for_shape += 2;
			} else {
				if ( neighbors[ searchPath.points[i] * L[3] + searchPath.points[j] ] != regionPath.matrix[ i * searchPath.L + j ] )
					score_for_shape += 2;
			}
		}
	}
	return score_for_shape;
}

void saveShape() {
	int i,e;                                  
	e = evalSearchPath();                                      
	if ( e <= maxMisMatch ) {
		fstPnt[numShapes+1] = fstPnt[numShapes];
		for ( i = 0 ; i < searchPath.L ; i++ )
			shapes.points[fstPnt[numShapes+1]++] = searchPath.points[i];
		shapes.region[numShapes] = searchPath.region[0];
		numShapes++;
	}
	cpn++;
}

void populateLibrary() {         
	int i,j;
	for ( i = 1 ; i <= arrayMax(searchPath.region,0,L[3]) ; i++ ) {
		saveToTarget();                        
		clearPly();
		savePtr = &saveShape;
		getRegion( i );
		unsaveRegion();                          
		for ( j = 0 ; j < searchPath.L ; j++ )
			searchPath.points[j] = L[3];
		search(1000000);
		unsaveTarget();
	}
}

void rePopulateLibrary() { 
	int i,j,k;
	for ( i = 0 ; i < searchPath.L ; i++ ) {
		if ( searchPath.region[i] == 0 ) {
			j = maxRegion() + 1;
			searchPath.region[i] = j;
			for ( k = 0 ; k < L[3] ; k++ ) {
				if ( boundType[k] != 6 )
					continue;
				fstPnt[numShapes+1] = fstPnt[numShapes];
				shapes.points[fstPnt[numShapes+1]++] = k;
				shapes.region[numShapes] = j;
				numShapes++;
			}
		}
	}
}
