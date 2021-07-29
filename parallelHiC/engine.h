struct path {                  // results whill be stored here
	int points[4096];          // current working path
	int chains[4096];          // chain each point belongs to	
	int bounds[4096];          // boundary conditions
	int region[4096];          // region segmentation array
	int matrix[16777216];      // contact matrix
	int L;                     // length of the path
} searchPath, 
  targetPath,
  regionPath;

// table to quickly check if points are neighbors
int neighbors[16777216];
int layerType[4096];
int boundType[4096];
  
// powers of the fundamental length
int L[7];

// sizes specified by the user
int mySizes[24];
int myEval;

// number of paths specified by the user
int sumPaths;
int pathSize;
int segsize;
int verbose;
int evalFun;
int bndopt;
int convolve;
int deconvolve;
int randSeed;

// move lists
int mL[16384];    // used to store lists of moves used by the search function
int fM[16384];    // keeps track of where the moves lists start and stop

// search function counters
int cpn;
int ply;
int pos;

// defines where the search results are saved
void (*savePtr)();

// functions
int pwr( int x,int y );
int x( int m );
int y( int m );
int z( int m );
void pseudoLegalMoves( int v ); 
void legalMoves( int v );
void search ( int total );
void clearPly();
void multiMansfield( int iterations );