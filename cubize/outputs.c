#include <string.h>       
#include <stdlib.h>
#include <stdio.h> 
#include "utilities.h"
#include "engine.h"
#include "library.h"
#include "shapes.h"
#include "outputs.h"

int readHicMatrix( char *matFile ) {
	int lineNum = 0,fieldNum = 0;
	FILE * fp;
	char * lineArray;
	char * line = NULL;
	size_t size = 0;
	int field;
	ssize_t input;
	fp = fopen(matFile, "r");
		if (fp == NULL)
			exit(EXIT_FAILURE);
	while ((input = getline(&line, &size, fp)) != -1) {
		lineArray = strtok (line," ");
		while (lineArray != NULL) {
			field = (int) strtol(lineArray, (char **)NULL, 10);
			if ( field > 0 )
				searchPath.matrix[fieldNum++] = 1;
			else
				searchPath.matrix[fieldNum++] = 0;
			lineArray = strtok (NULL," ");
			}
			lineNum++;
		}
	fclose(fp);
	return lineNum;
}

int readBoundary( char *bndFile ) {
	int lineNum = 0,fieldNum = 0;
	FILE * fp;
	char * lineArray;
	char * line = NULL;
	size_t size = 0;
	ssize_t input;
	fp = fopen(bndFile, "r");
	if (fp == NULL)
		exit(EXIT_FAILURE);
	while ((input = getline(&line, &size, fp)) != -1) {
		lineArray = strtok (line," ");
		while (lineArray != NULL) {
			regionPath.matrix[fieldNum++] = (int) strtol(lineArray, (char **)NULL, 10);
			lineArray = strtok (NULL," ");
			}
			lineNum++;
		}		
	fclose(fp);
	return lineNum;
}

void summarizeLibrary() {
	int i,j,k,max=0;
	for ( i = 0 ; i < numShapes ; i++ )
		if ( shapes.region[i] > max )
			max = shapes.region[i];
	for ( i = 0 ; i <= max ; i++ ) {
		k = 0;
		for ( j = 0 ; j < numShapes ; j++ )
			if ( shapes.region[j] == i )
				k++;
		printf("shapes in region %d = %d\n",i,k);
	}
	printf("%d total matching shapes\n",numShapes);
	printf("\n");	
}

void printPath(FILE *out) {
	int i;
	for ( i = 0 ; i < searchPath.L ; i++ )
		fprintf(out,"%d %d %d %d\n",x(searchPath.points[i]),y(searchPath.points[i]),z(searchPath.points[i]),searchPath.chains[i]);
}

void printMatrix(FILE *out) {
	int i;
	for ( i = 0 ; i < pwr(searchPath.L,2) ; i++ ) {
		fprintf(out,"%s%s",searchPath.matrix[i] == 0 ? " " : "1"," ");
		if ( i%searchPath.L == searchPath.L-1 )
			fprintf(out," %d%d\n",searchPath.bounds[i/searchPath.L],searchPath.region[i/searchPath.L]);
	} fprintf(out,"\n");
	for ( i = 0 ; i < searchPath.L ; i++ )
		fprintf(out,"%d ",searchPath.bounds[i]);
	fprintf(out,"\n");
	for ( i = 0 ; i < searchPath.L ; i++ )
		fprintf(out,"%d ",searchPath.region[i]);
	fprintf(out,"\n");
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// split array functions
///////////////////////////////////////////////////////////////////////////////////////////////////////
int split(const char *str, char sep, split_fn fun, int *array) {
    unsigned int start = 0, stop, position = 0;
    for (stop = 0; str[stop]; stop++) {
        if (str[stop] == sep) {
            fun(str + start, stop - start, array, position++);
            start = stop + 1;
        }
    }
    fun(str + start, stop - start, array, position++);
	return position;
}

void saveopts(const char *str, size_t len, int *array, int position) {
	char string[8];
    sprintf(string,"%.*s", (int)len, str);
	array[position] = atoi(str);
}