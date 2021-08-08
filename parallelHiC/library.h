typedef struct {
    int points[3000000000];
	int region[300000000];
} shape;
shape shapes;
int numShapes;
int maxMisMatch;
int fstPnt[300000000];

// functions
void deleteShape( int n );
void joinShapes( int a, int b, int c );
void getShape( int s );
int firstShape( int region );
int lastShape( int region );
int countRegion( int region );
int maxRegion();
int nullRegion();
void reSet();
void deleteRegion( int region );
void deleteLibrary();
void getShape( int s );
