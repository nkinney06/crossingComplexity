#include <stdio.h> 
#include <string.h>       
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <getopt.h>
#include "engine.h"
#include "initialization.h"
#include "outputs.h"
#include "contacts.h"
#include "library.h"
#include "shapes.h"
// gcc -Wall -Werror -O3 engine.c main.c utilities.c outputs.c initialization.c contacts.c shapes.c library.c -fopenmp
// scp -r -P 2222 nick@198.82.232.118:coolerinjupyter/data/*six ./
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820075_Cell19_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820068_Cell12_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820070_Cell14_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820059_Cell3_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820067_Cell11_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820057_Cell1_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820058_Cell2_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820066_Cell10_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820069_Cell13_model_six
// ./a.out -d -e2 -b -t96,81,39 -v -s6 -oGSM3820074_Cell18_model_six

int ftime_ok = 0;  /* does ftime return milliseconds? */
int get_ms() {
	struct timeb timebuffer;
	ftime(&timebuffer);
	if (timebuffer.millitm != 0)
		ftime_ok = 1;
	return (timebuffer.time * 1000) + timebuffer.millitm;
}

//////////////////////////////////////////////////////////////////////////////////////
// main programming section
//////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv){
	
	//////////////////////////////////////////////////////////
	// program parameters and input/output files
	////////////////////////////////////////////////////////// 
	pathSize         = 4;  // size of fundamental unit
	sumPaths         = 3;  // number of paths in the system
	segsize          = 8;  // boundary segmentation size
	verbose          = 0;  // flag for verbose
	myEval           = 0;  // evaluation function 1 or 2
	bndopt           = -1; // flag/score for boundary optimization
	convolve         = -1; // number of matricies to convolve
	maxMisMatch      = -1; // deconvolution mode
	numShapes        = 0;
	randSeed = time(NULL); // random seed
	memset(mySizes,-1,sizeof mySizes);
	char saveDir[127] = "./results";
	char matFile[127];
	char bndFile[127];
	char fldFile[127];
	char outMatr[127];
	char outFold[127];
	char outLogs[127];
	char outTarg[127];
	char outSols[127];
	
	//////////////////////////////////////////////////////////
	// https://azrael.digipen.edu/~mmead/www/Courses/CS180/getopt.html#OPTOPTARGS
	////////////////////////////////////////////////////////// 
	int c;
	while (1) 
	{
		int option_index = 0;
		static struct option long_options[] = 
		{
			{"pathSize",    required_argument, NULL,  's'},
			{"sizes",       required_argument, NULL,  't'},
			{"segsize",     required_argument, NULL,  'z'},
			{"numPaths",    required_argument, NULL,  'n'},
			{"eval",        required_argument, NULL,  'e'},
			{"seed",        required_argument, NULL,  'i'},
			{"saveDir",     required_argument, NULL,  'o'},
			{"help",        no_argument,       NULL,  'h'},
			{"verbose",     no_argument,       NULL,  'v'},
			{"convolve",    optional_argument, NULL,  'c'},
			{"deconvolve",  optional_argument, NULL,  'd'},
			{"bndopt",      optional_argument, NULL,  'b'},
			{NULL,          0,                 NULL,    0}
		};

		c = getopt_long(argc, argv, "-:s:t:e:z:n:i:o:hvc::d::b::", long_options, &option_index);
		if (c == -1)
		  break;

		switch (c) 
		{
		  case 0:
			break;
		  case 1:
			break;
		  case 's':
			pathSize = atoi(optarg);
			break;
		  case 't':
			split(optarg, ',', saveopts, mySizes);
			break;
		  case 'i':
			randSeed = atoi(optarg);
			break;
		  case 'o':
		    strcat(saveDir,"/");
		    strcat(saveDir,optarg);
			break;
		  case 'z':
			segsize = atoi(optarg);
			break;
		  case 'e':
			myEval = atoi(optarg);
			break;
		  case 'n':
			sumPaths = atoi(optarg);
			break;
		  case 'v':
			verbose = 1;
			break;
		  case 'c':
		    optarg ? convolve = atoi(optarg) : ( convolve = 1 );
            break;
		  case 'd':
		    optarg ? maxMisMatch = atoi(optarg) : ( maxMisMatch = 0 );
            break;
		  case 'b':
		    optarg ? bndopt = atoi(optarg) : ( bndopt = ( pwr(pathSize,3)/3 )*2 );
            break;			
		  case 'h':		  
			printf("\n\nSee the readme file for help\n\n");
			return 0;
			break;
		  case '?':
			printf("Unknown option %c\n", optopt);
			break;
		  case ':':
			printf("Missing option for %c\n", optopt);
			break;
		  default:
			printf("?? getopt returned character code 0%o ??\n", c);
		 }
	}

	//////////////////////////////////////////////////////////
	// setup the output files
	////////////////////////////////////////////////////////// 
	sprintf(matFile,"%s%s",saveDir,"/inputMatrix.txt");
	sprintf(bndFile,"%s%s",saveDir,"/inputBounds.txt");
	sprintf(fldFile,"%s%s",saveDir,"/inputFolded.txt");
	sprintf(outMatr,"%s%s",saveDir,"/outputMatrix.txt");
	sprintf(outFold,"%s%s",saveDir,"/outputFolded.txt");
	sprintf(outLogs,"%s%s",saveDir,"/outputLogged.txt");
	sprintf(outTarg,"%s%s",saveDir,"/outputTarget.txt");
	sprintf(outSols,"%s%s",saveDir,"/outputSolved.txt");
	FILE *sols = fopen(outSols, "w");
	FILE *fold = fopen(outFold, "w");
	FILE *matr = fopen(outMatr, "w");
	FILE *outf = fopen(outTarg, "w");
	FILE *logs = fopen(outLogs, "w");
	if (sols == NULL) {  exit(1);  }
	if (fold == NULL) {  exit(1);  }
	if (matr == NULL) {  exit(1);  }
	if (outf == NULL) {  exit(1);  }
	if (logs == NULL) {	 exit(1);  }
	
	//////////////////////////////////////////////////////////
	// initialization steps
	////////////////////////////////////////////////////////// 
	int i,j,k,starttime,stoptime;;
	srand(randSeed);
	printf("\nset random seed to %d\n",randSeed);
	for ( i = 0 ; i < 7 ; i++ )
		L[i] = pwr(pathSize,i);

	if ( convolve > 0 )    // convolve the user specified number of matricies
		convolution(convolve,matFile,bndFile,fldFile);
	else
		initialize();
	
	if ( maxMisMatch < 0 ) // if we are not deconvolving the matrix then exit
		return 0;
	
	int b = readBoundary(bndFile);	
	if ( bndopt >= 0 )     // read the boundary file and try to fit
		fitBoundary( logs,b-1,bndopt );
	else
		for ( i = 0 ; i < L[3] ; i++ )
			searchPath.bounds[i] = regionPath.matrix[i+0*L[3]];

	searchPath.L = readHicMatrix(matFile);	
	segmentBoundary(segsize);	
	printMatrix(outf);
		
	//////////////////////////////////////////////////////////
	// algorithm section 1
    //////////////////////////////////////////////////////////
	printf("\nstarting part one of the algorithm\n\n");	
	starttime = get_ms();	
	populateLibrary();
	stoptime = get_ms();
	fprintf(logs,"library_time_1\t%d\n", (stoptime - starttime)/1000);
	
	//////////////////////////////////////////////////////////
	// algorithm section 2
    //////////////////////////////////////////////////////////
	int totIterations = maxRegion() + (pathSize-2)*(pathSize-2)*(pathSize-2) - 1; 
	int numIterations = 0;
	starttime = get_ms();
	while ( maxRegion() > 1 ) {
		fprintf(logs,"library_size_%d\t%d\n",numIterations,numShapes);
		bestRegionPair(&j,&k);	
		printf("joining region %d and %d, %d/%d iterations\n",j,k,++numIterations,totIterations);
		if ( verbose )
			summarizeLibrary();		
		joinAllShapesSimple( j,k );
		deleteRegion( j );
		deleteRegion( k );
		reSet();
	}
	stoptime = get_ms();
	fprintf(logs,"algo_time_1\t%d\n", (stoptime - starttime)/1000);
	
	//////////////////////////////////////////////////////////
	// algorithm section 3
    //////////////////////////////////////////////////////////
	printf("\nstarting part two of the algorithm\n\n");
	starttime = get_ms();
	rePopulateLibrary();
	stoptime = get_ms();
	fprintf(logs,"library_time_2\t%d\n", (stoptime - starttime)/1000);
	
	//////////////////////////////////////////////////////////
	// algorithm section 4
    //////////////////////////////////////////////////////////
	starttime = get_ms();
	while ( maxRegion() > 1 ) {
		bestRegionPair(&j,&k);
		printf("joining region %d and %d, %d/%d iterations\n",j,k,++numIterations,totIterations);
		if ( verbose )
			summarizeLibrary();
		joinAllShapesSimple( j,k );
		deleteRegion( j );
		deleteRegion( k );
		reSet();
	}
	stoptime = get_ms();
	fprintf(logs,"algo_time_2\t%d\n", (stoptime - starttime)/1000);
	
	//////////////////////////////////////////////////////////
	// finish up the run
    ////////////////////////////////////////////////////////// 
	printf("found %d shapes matching the input matrix\n\n",numShapes);
	fprintf(logs,"number_paths\t%d\n",sumPaths);
	fprintf(logs,"square_size\t%d\n",L[1]);
	fprintf(logs,"max_mismatch\t%d\n",maxMisMatch);
	fprintf(logs,"seg_size\t%d\n",segsize);
	fprintf(logs,"random_seed\t%d\n",randSeed);
	fprintf(logs,"solutions\t%d\n",numShapes);
	
	getShape(0);           // final output
	redivide();
	currentHiC();
	printPath(fold);
	printMatrix(matr);

	if (numShapes <= 48) { // print the first 48 solutions
		for ( i = 0 ; i < numShapes ; i++ ) {
			getShape(i);
			redivide();
			for ( j = 0 ; j < searchPath.L ; j++ )
				fprintf(sols,"%d %d %d %d %d\n",x(searchPath.points[j]),y(searchPath.points[j]),z(searchPath.points[j]),searchPath.chains[j],i);
		}
	}

	fclose(outf);
	fclose(matr);
	fclose(fold);
	fclose(logs);
	fclose(sols);
	deleteLibrary();
    return 0;
}
