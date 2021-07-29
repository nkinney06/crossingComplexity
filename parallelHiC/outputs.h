int readHicMatrix( char *matFile );
int readBoundary( char *bndFile );
void summarizeLibrary();
void printPath(FILE *out);
void printShapeMatrix(int a);
void printMatrix(FILE *out);
typedef void(*split_fn)(const char *, size_t, int *array, int position);
int split(const char *str, char sep, split_fn fun, int *array);
void saveopts(const char *str, size_t len, int *array, int position);