#include <stdio.h> 
#include <string.h>       
#include <stdlib.h> 
#include <time.h>
#include <getopt.h>

//////////////////////////////////////////////////////////////////////////////////////
//
// README
//
// this is a minimalist implementation of the "multi-Mansfield" algorithm for generating
// two space filling curves on a simple cubic lattice described in the paper "Crossing 
// Complexity of Space-Filling Curves reveals Entanglement of S-Phase DNA". The code is 
// not optimized and THERE IS NOT DYNAMIC MEMORY ALLOCATION. Some parts of this code are
// rather difficult. Questions can be directed to nkinney06@gmail.com.
//
// I thought about doing a major revision including renaming variables and function but
// decided to keep the code identical to that used to generate data in our publication.
// updates will be made in the following git repo:
//
// comments in this version of the code are sparse. Again see the git repo for updates
//
// to compile this code:
// gcc -Wall -Werror -O3 S1code.c -o multiHamiltonian
//
// this code has only been tested in a linux environment (16.04)
// if issues arise contact nkinney06@gmail.com. 
//
// usage: 
// ./multiHamiltonian --msd <int> --numberPaths <int> --initialSize <int> --cycles <string>
//
// the options are as follows:
//
// "msd" - specifies the maximum mean squared deviation for the space filling curves
//         for equal size curves use --msd 0
//		
// "numberPaths" - specifies the number of curves to place on the same lattice
//                 CAUTION: the program as is does not use dynamic memory allocation
//                 to thing might break if you exceed ~10 paths on the same lattice.
//                 Again, this code is a minimalist starting point
//			
// "initialSize" - starting size of the cubic lattice; for example, 4 will create the
//                 initial configuration on a 4x4x4 lattice. CAUTION: there is no
//                 dynamic memory allocation. For very large lattices array sizes need
//                 to be increased.
//				
// "cycles" - this is a string represent the order of "nesting operations". Use intgers
//            2 or 3 separated by commas; i.e. "2,2" or "2,3,2". For example the string
//            "2,2" will nest 2x2x2 configurations on top of each point in the initial
//            configuration; then, will next 2x2x2 configurations again. the string "2x3x2"
//            will nest 2x2x2 configurations on the initial configuration, followed by 3x3x3
//            configurations, followed by 2x2x2 configurations.
//	   
//            CAUTION: nesting 3x3x3 configurations require up to 10 minutes of execution.
//            CAUTION: nesting causes the lattice size to increase VERY quickly. Nesting
//            depth above ~4 would require increase in array sizes.
//				
// the output is space delimited with the follow fields:
// <x> <y> <z> <n>
//
// Here x, y, and z are the integer coordinates. n is the curve number
//	            
//////////////////////////////////////////////////////////////////////////////////////
					
//////////////////////////////////////////////////////////////////////////////////////
// global variables and structures
//////////////////////////////////////////////////////////////////////////////////////

int cpn;               // current path number
int pos;               // current position
int ply;               // search depth used by the search function
int npp;               // number of precomputed paths

int moveLists[16384];  // used to store lists of moves used by the search function
int firstMove[16384];  // keeps track of where the moves lists start and stop

struct computed {
	int cwp[16384];    // current working path. 
	int l;
	int L;
	int n;
};

struct precomputed {
	int cwp[27];
	int cycl;
	int mark;
	int l;
	int L;
};

struct computed cp, tp;
struct computed np[12];
struct precomputed pp[4960752];

void (*savePtr)();
typedef int (*markPtr)( int i, int u, int v );
typedef void (*split_fn)(const char  *, size_t, int *array, int position);

//////////////////////////////////////////////////////////////////////////////////////
// coordinate functions
//////////////////////////////////////////////////////////////////////////////////////

int pwr( int x,int y ) { return y ? (y%2?x:1)*pwr(x*x,y>>1) : 1; } // to skip math.h
int x  ( int m,int d ) { return (m/pwr(d,0))%d; }
int y  ( int m,int d ) { return (m/pwr(d,1))%d; }
int z  ( int m,int d ) { return (m/pwr(d,2))%d; }

//////////////////////////////////////////////////////////////////////////////////////
// move generator functions
//////////////////////////////////////////////////////////////////////////////////////

void legalMoves ( int v, int d  ) {
	int i,j;
	if ( v == (d*d*d) ) {
		for ( i = 0 ; i < (d*d*d) ; i++ )
			moveLists[firstMove[ply+1]++] = i;
		return;
		}
	for ( i = 0 ; i < 3 ; i++ ) {
		j = (v/pwr(d,i))%d;
		if ( j == 0 )
			moveLists[firstMove[ply+1]++] = v+pwr(d,i);
		else if ( j == d-1 )
			moveLists[firstMove[ply+1]++] = v-pwr(d,i);
		else {
				moveLists[firstMove[ply+1]++] = v+pwr(d,i);
				moveLists[firstMove[ply+1]++] = v-pwr(d,i);
			}
		}
	return;
}

void shuffle() {
	int i,j,n = firstMove[ply+1] - firstMove[ply];
    if (n > 1) {
		for ( i = 0 ; i < n ; i++ ) {
          j = i + ( rand() / (RAND_MAX / (n - i) + 1) ) + firstMove[ply];
          int t = moveLists[j];
          moveLists[j] = moveLists[i+firstMove[ply]];
          moveLists[i+firstMove[ply]] = t;
        }
    }
}

int pathContains( int *path,int end,int v ) {
	int i;
	for ( i = 0 ; i < end ; i++ )
		if ( path[i] == v )
			return i;
	return -1;
}

void search ( int total ) {
	int i;
	firstMove[ply+1] = firstMove[ply];
	if (ply==cp.L) {
		(*savePtr)();
		return;
	}
	legalMoves(pos,cp.l);
	for ( i = firstMove[ply] ; i < firstMove[ply+1] ; i++ ) {
		if ( pathContains(cp.cwp,cp.L,moveLists[i]) >= 0 )
			continue;
		if ( cpn == total )
			return;
		cp.cwp[ply++] = moveLists[i];
		pos = moveLists[i];
		search(total);
		cp.cwp[--ply] = cp.L;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// save functions
/////////////////////////////////////////////////////////////////////////////////////////

void saveNP() {
	np[cpn] = cp;
	cpn++;
}

void savePP() {
	int i;
	for ( i = 0 ; i < cp.L ; i++ )
		pp[cpn].cwp[i] = cp.cwp[i];
	pp[cpn].cycl = 0;
	pp[cpn].mark = 0;
	pp[cpn].l = cp.l;
	pp[cpn].L = cp.L;
	legalMoves(cp.cwp[cp.L-1],cp.l);
	for ( i = firstMove[ply] ; i < firstMove[ply+1] ; i++ )
		if (moveLists[i] == cp.cwp[0])
			pp[cpn].cycl = 1;
	cpn++;
}

/////////////////////////////////////////////////////////////////////////////////////////
// input and output functions
/////////////////////////////////////////////////////////////////////////////////////////

typedef void(*split_fn)(const char *, size_t, int *array, int position);
 
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

void printUsage() {
	printf("usage: ./nestedPath \"3,2,3\" (for now only use comma separated 2s and 3s)\n");
	exit(1);
}

void checkopts(char *str) {
	int i, len = strlen(str);
	char warn = str[0];
	for (i = 1; i < len; i++) {
		if ( (str[i] == '2') | (str[i] == '3') | (str[i] == ',') ) {
			if ( str[i] == warn )
				printUsage();
			else
				warn = str[i];
		}
		else
			printUsage();
	}
}

void printPath() {
	int i;
	for ( i = 0 ; i < cp.L ; i++ )
		printf("%d %d %d %d\n",x(cp.cwp[i],cp.l),
							   y(cp.cwp[i],cp.l),
							   z(cp.cwp[i],cp.l),
							   cp.n);
}

/////////////////////////////////////////////////////////////////////////////////////////
// functions to fetch and move paths
/////////////////////////////////////////////////////////////////////////////////////////

void getNP( int n ) { cp = np[n]; }
void getTP()        { cp = tp;    }
void putNP( int n ) { np[n] = cp; }
void delTP()        { tp.L = 0;   }

void getPP( int n ) {
	int i;
	cp.l = pp[n].l;
	cp.L = pp[n].L;
	cp.n = n;
	for ( i = 0 ; i < pp[n].L ; i++ )
		cp.cwp[i] = pp[n].cwp[i];
}

/////////////////////////////////////////////////////////////////////////////////////////
// initialization functions
/////////////////////////////////////////////////////////////////////////////////////////

void clearMoves() {
	ply = 0;
	firstMove[0] = ply;
	firstMove[ply+1] = firstMove[ply];
}

void initialize( int d ) {
	int i;
	cp.l = d;
	cp.L = pwr(d,3);
	cp.n = 0;
	pos = cp.L;
	clearMoves();
	for ( i = 0 ; i < cp.L ; i++ )
		cp.cwp[i] = cp.L;
}

void initConfig( int d ) {
	cpn = 0;
	savePtr = &saveNP; 
	initialize(d);
	search(1);
	getNP(0);
}

void precomputePaths( int d ) {
	initialize(d);
	cpn = ( d == 2 ) ? 0 : 144;
	savePtr = &savePP;          
	search( d == 2 ? 144 : 4960752 );
	npp = d == 2 ? 144 : 4960752;
}

/////////////////////////////////////////////////////////////////////////////////////////
// counting and finding functions
/////////////////////////////////////////////////////////////////////////////////////////

int countMarked( int n ) {
	int i,j = 0;
	for ( i = 0; i < npp ; i++ )
		if ( pp[i].mark == n )
			j++;
	return j;
}

int chooseRandom( int n ) {
	int i,j=0,k;
	k = rand() % countMarked(n);
	for ( i = 0; i < npp ; i++ )
		if ( pp[i].mark == n )
			if ( j++ == k )
				return i;
	return 0;
}
   
int comparePathLength ( int n ) {
	int i,j,msd = 0;
	for ( i = 0 ; i < n ; i++ )
		for ( j = i+1 ; j < n ; j++ )
			msd += pwr(np[i].L - np[j].L,2);
	if ( n == 1 )
		return 0;
	else
		return msd / (((n*n)-n)/2);
}

int maxMark() {
	int i,max = 0;
	for ( i = 0 ; i < npp ; i++ )
		if ( pp[i].mark > max )
			max = pp[i].mark;
	return max;
}

/////////////////////////////////////////////////////////////////////////////////////////
// index manipulation functions
/////////////////////////////////////////////////////////////////////////////////////////

void setIndex( int n ) {
	int i;
	for ( i = 0 ; i < cp.L ; i++ )
		cp.cwp[i] = x(cp.cwp[i],cp.l)*pwr(n,0) + 
					y(cp.cwp[i],cp.l)*pwr(n,1) + 
					z(cp.cwp[i],cp.l)*pwr(n,2);
	cp.l = n;
}

void translate( int a, int b, int c ) {
	int i;
	for ( i = 0 ; i < cp.L ; i++ )
		cp.cwp[i] = (x(cp.cwp[i],cp.l)+a)*pwr(cp.l,0) +
			        (y(cp.cwp[i],cp.l)+b)*pwr(cp.l,1) +
				    (z(cp.cwp[i],cp.l)+c)*pwr(cp.l,2);
}

/////////////////////////////////////////////////////////////////////////////////////////
// list manipulation functions
/////////////////////////////////////////////////////////////////////////////////////////

void append() {
	int i;
	for ( i = 0 ; i < cp.L ; i++ )
		tp.cwp[tp.L++] = cp.cwp[i];
}

void reverse( int *path,int n ) { 
	int i,j,k = n-1;
	for ( i = 0 ; i < n/2 ; i++ ) {
		j = path[i];
		path[i] = path[k];
		path[k] = j;
		k--;
	}
	return;
}
 
void dividecwp( int pieces ) {
	int i;
	for ( i = 0 ; i < cp.L ; i++ ) {
		np[i/((cp.L+pieces)/pieces)].cwp[i%((cp.L+pieces)/pieces)] = cp.cwp[i];
		np[i/((cp.L+pieces)/pieces)].L = i%((cp.L+pieces)/pieces) + 1;
		np[i/((cp.L+pieces)/pieces)].l = cp.l;
		np[i/((cp.L+pieces)/pieces)].n = i/((cp.L+pieces)/pieces);
	}
	return;
}

void swapSubPaths( int pathFrom,int pathGoto,int branch) {
	int i;
	for ( i = branch ; i < np[pathGoto].L ; i++ )
		np[pathFrom].cwp[np[pathFrom].L++] = np[pathGoto].cwp[i];
	np[pathGoto].L = branch;
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////
// mark functions to help select from the precomputed path list
/////////////////////////////////////////////////////////////////////////////////////////

int markSize( int i, int u, int v ) {
	if ( ( pp[i].l == u ) && ( pp[i].L = v ) )
		return 1;
	return 0;
}

int markWall( int i, int u, int v ) {
	if ( (pp[i].cwp[pp[i].L-1]/pwr(pp[i].l,u)) % pp[i].l == v )
		return 2;
	return 0;
}

int markValu( int i, int u, int v ) {
	if ( pp[i].cwp[u] == v )
		return 4;
	return 0;
}

void mark( int u, int v, markPtr fun) {
	int i;
	for ( i = 0 ; i < npp ; i++ )
		pp[i].mark += fun(i,u,v);
}

void unmark( int n ) {
	int i;
	for ( i = 0 ; i < npp ; i++ )
		pp[i].mark &= n;
}

int direction( int i, int j, int d ) {
	return ( x(i,d) - x(j,d) ) ? 0 : ( y(i,d) - y(j,d) ) ? 1 : 2;
}

int nextValue( int i, int j, int n ) {
	return ( ( i - j ) > 0 ) ? n-1 : 0;
}

int mustBegin( int a, int b, int c ) {
	return c ? a - (cp.l-1)*pwr(cp.l,b) : a + (cp.l-1)*pwr(cp.l,b);	
}

/////////////////////////////////////////////////////////////////////////////////////////
// algorithms
/////////////////////////////////////////////////////////////////////////////////////////

void mansfield( int iterations ) {
	int i;
	for ( i = 0 ; i < iterations ; i++ ) {
		clearMoves();
		if ( rand() % 2 )
			reverse(cp.cwp,cp.L);       
		legalMoves(cp.cwp[0],cp.l);
		shuffle();
		int choice = rand()%6;    
		if ( choice < firstMove[ply+1] )
			reverse(cp.cwp,pathContains(cp.cwp,cp.L,moveLists[choice]));
	}
}

void multiMansfield( int iterations,int n ) {
	int i;
	for ( i = 0 ; i < iterations ; i++ ) {
		clearMoves();
		int pathGoto = -1;
		int pathFrom = rand() % n;
		int branch = -1;
		if ( rand() % 2 )
			reverse(np[pathFrom].cwp,np[pathFrom].L);
		legalMoves(np[pathFrom].cwp[np[pathFrom].L-1],np[pathFrom].l);
		shuffle();
		int choice = rand()%6;    
		if ( choice < firstMove[ply+1] ) {
			do { pathGoto++;
				 branch = pathContains(np[pathGoto].cwp,np[pathGoto].L,moveLists[choice]);
			   } while ( branch == -1 );
			if ( pathFrom == pathGoto ) {
				 reverse(np[pathFrom].cwp,np[pathFrom].L);
				 branch = pathContains(np[pathGoto].cwp,np[pathGoto].L,moveLists[choice]);
				 reverse(np[pathGoto].cwp,branch);
			} else {  
				 if ( rand() % 2 )
					 reverse(np[pathGoto].cwp,np[pathGoto].L);
				 branch = pathContains(np[pathGoto].cwp,np[pathGoto].L,moveLists[choice]);
				 if ( ( branch > 0 ) & ( branch < np[pathGoto].L-1 ) )
					swapSubPaths(pathFrom,pathGoto,branch);
			}
		}
	}
}

void nest(int n,int m) {
	int i;
	delTP();
	unmark(0);
	mark(n,pwr(n,3),&markSize);
	int a = -1, b = -1, c = -1;
	for ( i = 0 ; i < np[m].L ; i++ ) {
		int thisPoint = np[m].cwp[i];
		int nextPoint = np[m].cwp[i+1];
		unmark(1);
		if ( i < np[m].L-1 ) {
			b = direction(nextPoint,thisPoint,np[m].l);
			c = nextValue(nextPoint,thisPoint,n);
			mark(b,c,&markWall);
		}
		if ( a > -1 )
			mark(0,a,&markValu);
		getPP(chooseRandom(maxMark()));
		if ( b >= 0 )
			a = mustBegin(cp.cwp[cp.L-1],b,c);
		setIndex(np[m].l*n);
		translate(x(thisPoint,np[m].l)*n,
			      y(thisPoint,np[m].l)*n,
			      z(thisPoint,np[m].l)*n);
		append();
	}
	getTP();
	cp.l = np[m].l*n;
	cp.n = m;
	putNP(m);
}

//////////////////////////////////////////////////////////////////////////////////////
// main programming section
//////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{	

    //////////////////////////////////////////////////////////
    // variables
    ////////////////////////////////////////////////////////// 
	srand(time(NULL)); 
	int pathSize   = 4;                      
	int numPaths   = 2;                      
	int msd        = 0;
	int ls[10];
	int iterations = 0;
	int i,j;
	
    //////////////////////////////////////////////////////////
    // get the parameters
    //////////////////////////////////////////////////////////   
    int c;
    int option_index    = 0;
    struct option long_options[] = {
			   {"msd",        required_argument, 0, 'm' },
			   {"numberPaths",required_argument, 0, 'n' },
			   {"initialSize",required_argument, 0, 'i' },
			   {"cycles",     required_argument, 0, 's' },
			   {0,            0,                 0, 0   }
		   };

	while (1) {
		   c = getopt_long(argc, argv,"",long_options, &option_index);
		   if (c == -1) 
			   break;
		   switch (c) {
			    case 0:
					printf("option %s", long_options[option_index].name);
					if (optarg)
						printf(" with arg %s", optarg);
					printf("\n");
					break;
				case 'm':
					msd      = atoi(optarg);
					break;
				case 'n':
					numPaths = atoi(optarg);
					break;
				case 'i':
					pathSize = atoi(optarg);
					break;
				case 's':
					checkopts(optarg);
					iterations = split(optarg, ',', saveopts, ls);
					break;
		   		case '?':
					break;
				default:
					printf("?? getopt returned character code 0%o ??\n", c);
					break;
               }
           }

    //////////////////////////////////////////////////////////
    // precompute h-paths ( avoid 3x3x3 paths if possible )
    ////////////////////////////////////////////////////////// 
	precomputePaths(2);
	for ( i = 0 ; i < iterations ; i++ ) {
		if ( ls[i] == 3 ) {
			precomputePaths(3);
			break;
		}
	}

    //////////////////////////////////////////////////////////
    // construct the curves and print output
    ////////////////////////////////////////////////////////// 
	initConfig(pathSize);	
	getNP(0);
	dividecwp(numPaths);        
	multiMansfield(100000,numPaths);
	
	while( comparePathLength(numPaths) > msd )
		multiMansfield(1,numPaths);

	for ( i = 0; i < numPaths ; i++ )
		for ( j = 0; j < iterations ; j++ ) 
			nest(ls[j],i);

	for ( i = 0; i < numPaths ; i++ ) {
		getNP(i);
		printPath();
	}

	return 0;
}


