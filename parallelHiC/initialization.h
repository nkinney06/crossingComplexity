#define pass (void)0
void segmentBoundary( int maxSize );
void saveToTarget();
void unsaveTarget();
void saveToRegion();
void unsaveRegion();
void divide();
void redivide();
void initialize();
void initConfig();
void fitBoundary( FILE *out,int n,int max );
void convolution( int convolve,char *matFile,char *bndFile, char *fldFile );