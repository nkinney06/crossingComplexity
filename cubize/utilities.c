#include <string.h>       
#include <stdlib.h>
#include <stdio.h> 
#include "utilities.h"

/////////////////////////////////////////////////////////////////////////////////////////
// list manipulation functions. these function all have the same signature
/////////////////////////////////////////////////////////////////////////////////////////
void arrayRev( int *arr, int s, int to ) { 
	int i,j,k = to-1;
	for ( i = s ; i < to/2 ; i++ ) {
		j = arr[i];
		arr[i] = arr[k];
		arr[k] = j;
		k--;
	}
	return;
}

int arrayCon( int *arr, int end, int v ) {
	int i;
	for ( i = 0 ; i < end ; i++ )
		if ( arr[i] == v )
			return i;
	return -1;
}
 
int arrayMax( int *arr, int n, int end ) { 
	int i;
	for ( i = 0 ; i < end ; i++ )   
		if ( arr[i] > n )
			n = arr[i];
	return n;
}

int arrayCnt( int *arr, int n, int end ) { 
	int i,j=0;
	for ( i = 0 ; i < end ; i++ )    
		if ( arr[i] == n )
			j++;
	return j;
}

int arrayFst( int *arr, int n, int end ) { 
	int i;
	for ( i = 0 ; i < end ; i++ )    
		if ( arr[i] == n )
			return i;
	return 0;
}

int arrayLst( int *arr, int n, int end ) { 
	int i,j=0;
	for ( i = 0 ; i < end ; i++ )    
		if ( arr[i] == n )
			j = i;
	return j;
}
