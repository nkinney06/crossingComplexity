#include <stdio.h> 
#include <string.h>       
#include <stdlib.h>
#include <time.h>
#include <sys/timeb.h>
#include <getopt.h>

//////////////////////////////////////////////////////////////////////////////////////
// gcc -Wall -Werror -O3 cubizeDraftThree.c -mcmodel=large -o cubize            
//////////////////////////////////////////////////////////////////////////////////////
					
//////////////////////////////////////////////////////////////////////////////////////
// global variables and structures
//////////////////////////////////////////////////////////////////////////////////////

struct computed {    // results whill be stored here
	int p[16384];    // current working path. 
	int L;           // length of the path
	int n;           // path number
};

// counters
int cpn,cpm,pos,ply;  // search depth used by the search function
int TPN,TOL;          // total pathNumber, tolerance
int L[7];             // powers of the fundamental length

int sL[1000001000];
int sR[100001000];
int sE[100001000];
int fS[100001000];
int mL[16384];        // used to store lists of moves used by the search function
int fM[16384];        // keeps track of where the moves lists start and stop
int bC[16384];        // b oundary conditions
int bS[16384];        // boundary segments
int lp[531441];       // for computing the hic matrix
int Mat[532899];      // current working path matrix
int HiC[532899];      // this is the target matrix that we want to match
int Reg[532899];      // space to hold the submatrix for reconstructions

struct computed cp, tp;     // single paths (current path and temp path)
struct computed np[1000];   // numbered paths

void (*savePtr)(); // pointer for how the search function saves paths
int (*evalPtr)();  // pointer for how the search function evals paths

//////////////////////////////////////////////////////////////////////////////////////
// timing functions
//////////////////////////////////////////////////////////////////////////////////////
int ftime_ok = 0;  /* does ftime return milliseconds? */
int get_ms() {
	struct timeb timebuffer;
	ftime(&timebuffer);
	if (timebuffer.millitm != 0)
		ftime_ok = 1;
	return (timebuffer.time * 1000) + timebuffer.millitm;
}

//////////////////////////////////////////////////////////////////////////////////////
// coordinate functions
//////////////////////////////////////////////////////////////////////////////////////
int pwr( int x,int y ) { return y ? (y%2?x:1)*pwr(x*x,y>>1) : 1; } // to skip math.h
int x  ( int m ) { return (m/L[0])%L[1]; }
int y  ( int m ) { return (m/L[1])%L[1]; }
int z  ( int m ) { return (m/L[2])%L[1]; }

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

int arraySum( int *arr, int n, int end ) {
	int i,j=0;
	for ( i = n ; i < end ; i++ )    
		j += arr[i];
	return j;
}

//////////////////////////////////////////////////////////////////////////////////////
// move generator functions and move shuffling
//////////////////////////////////////////////////////////////////////////////////////
void shuffle() {
	int i,j,n = fM[ply+1] - fM[ply];
    if (n > 1) {
		for ( i = 0 ; i < n ; i++ ) {
          j = i + ( rand() / (RAND_MAX / (n - i) + 1) ) + fM[ply];
          int t = mL[j];
          mL[j] = mL[i+fM[ply]];
          mL[i+fM[ply]] = t;
        }
    }
}

int boundaryType( int v ) {
	int i,j,k = 0;
	for ( i = 0 ; i < 3 ; i++ ) {
		j = (v/L[i]) % L[1];
		if ( j == 0 )
			k += 1;
		else if ( j == L[1]-1 )
			k += 1;
		else
			k+= 2;
	}
	return k;
}

void pseudoLegalMoves( int v ) { // for computing the hic Matrix
	int i,j;
	for ( i = 0 ; i < 3 ; i++ ) {
		j = (v/L[i]) % L[1];
		if ( j == 0 )
			mL[fM[ply+1]++] = v+L[i];
		else if ( j == L[1]-1 )
			mL[fM[ply+1]++] = v-L[i];
		else {
			mL[fM[ply+1]++] = v+L[i];
			mL[fM[ply+1]++] = v-L[i];
			}
		}
	return;
}

void legalMoves( int v, int s ) { // s is the start point
	int i,j,n;
	if ( v == L[3] ) {
		if ( s < 0 ) {
			for ( i = 0 ; i < L[3] ; i++ )
				mL[fM[ply+1]++] = i;
		} else
			mL[fM[ply+1]++] = s;
		return;
		}
	for ( i = 0 ; i < 3 ; i++ ) {
		j = (v/L[i]) % L[1];
		if ( j == 0 ) {
			n = v+L[i];
			if ( ( bC[ply] == 0 ) | ( bC[ply] == boundaryType(n) ) )
				mL[fM[ply+1]++] = n;
		}
		else if ( j == L[1]-1 ) {
			n = v-L[i];
			if ( ( bC[ply] == 0 ) | ( bC[ply] == boundaryType(n) ) )
				mL[fM[ply+1]++] = n;
		}
		else {
			n = v+L[i];
			if ( ( bC[ply] == 0 ) | ( bC[ply] == boundaryType(n) ) )
				mL[fM[ply+1]++] = n;
			n = v-L[i];
			if ( ( bC[ply] == 0 ) | ( bC[ply] == boundaryType(n) ) )
				mL[fM[ply+1]++] = n;
			}
		}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////
// depth first search
/////////////////////////////////////////////////////////////////////////////////////////
void search ( int total, int s ) {
	int i;
	fM[ply+1] = fM[ply];
	if (ply==cp.L) {
		(*savePtr)();
		return;
	}
	legalMoves(pos,s);
	for ( i = fM[ply] ; i < fM[ply+1] ; i++ ) {
		if ( arrayCon(cp.p,cp.L,mL[i]) >= 0 )
			continue;
		if ( cpn == total )
			return;
		cp.p[ply++] = mL[i];
		pos = mL[i];
		search(total,s);
		cp.p[--ply] = L[3];
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// output functions
/////////////////////////////////////////////////////////////////////////////////////////
void  printPath(FILE *out) {
	int i;
	for ( i = 0 ; i < cp.L ; i++ )
		fprintf(out,"%d %d %d %d\n",x(cp.p[i]),y(cp.p[i]),z(cp.p[i]),cp.n);
}

void printsCwpMatrix( FILE *out,int pathSize ) {  // prints the current working path matrix
	int j;
	for ( j = 0 ; j < pwr(pathSize,2) ; j++ ) {
		fprintf(out,"%s%s",Mat[j] == 0 ? " " : "1"," ");
		if ( j%pathSize == pathSize-1 )
			fprintf(out,"  %d\n",Mat[pwr(pathSize,2)+j/pathSize]);
	} fprintf(out,"\n");
	for ( j = pwr(pathSize,2) ; j < pwr(pathSize,2)+pathSize ; j++ )
		fprintf(out,"%d ",Mat[j]);
	fprintf(out,"\n");
}

void printscBoundary( FILE *out,int pathSize ) {
	int i;
	for ( i = 0 ; i < pathSize ; i++ )
		fprintf(out,"%d ",bC[i]);
	fprintf(out,"\n");
}

void printsbSegments() {
	int i;
	printf("\ncurrent boundary segmentation\n");
	for ( i = 0 ; i < L[3] ; i++ )
		printf("%d ",bS[i]);
	printf("\n");
}

void printRegion( int *pairs, int size ) {
	int i;
	printf("\nsummary of contacts between segmented regions\n");
	for ( i = 0 ; i < size*size ; i++ )
		printf("%d%s",pairs[i],i % size == ( size - 1 ) ? "\n" : " ");
	printf("\n");
}

/////////////////////////////////////////////////////////////////////////////////////////
// save functions, either save to numbered paths or subregions array
/////////////////////////////////////////////////////////////////////////////////////////
void saveNP() { 
	np[cpn++] = cp;
}

void saveSR( int cpeval ) {
	int i;
	fS[cpm+1] = fS[cpm];
	for ( i = 0 ; i < cp.L ; i++ )
		sL[fS[cpm+1]++] = cp.p[i];
	sE[cpm] = cpeval;
	cpm++;
}

/////////////////////////////////////////////////////////////////////////////////////////
// initialization functions
/////////////////////////////////////////////////////////////////////////////////////////
void clearPly() {
	ply = 0;
	fM[0] = ply;
	fM[ply+1] = fM[ply];
}

void initialize() {
	int i;
	cpn = 0;
	savePtr = &saveNP; 
	cp.L = L[3];
	cp.n = 0;
	pos = cp.L;
	clearPly();
	memset(bC, 0, sizeof bC);
	for ( i = 0 ; i < cp.L ; i++ )
		cp.p[i] = L[3];
}

void initConfig() {
	initialize();
	memset(bC, 0, sizeof bC);
	search(1,-1); // search for 1 path with any starting position
	cp = np[0];   // set current path as the first result
}

/////////////////////////////////////////////////////////////////////////////////////////
// list manipulation functions
/////////////////////////////////////////////////////////////////////////////////////////
void dividecwp() {
	int i;
	for ( i = 0 ; i < cp.L ; i++ ) {
		np[i/((cp.L+TPN)/TPN)].p[i%((cp.L+TPN)/TPN)] = cp.p[i];
		np[i/((cp.L+TPN)/TPN)].L = i%((cp.L+TPN)/TPN) + 1;
		np[i/((cp.L+TPN)/TPN)].n = i/((cp.L+TPN)/TPN);
	}
	return;
}

void swapSubPaths( int pathFrom,int pathGoto,int branch) {
	int i;
	for ( i = branch ; i < np[pathGoto].L ; i++ )
		np[pathFrom].p[np[pathFrom].L++] = np[pathGoto].p[i];
	np[pathGoto].L = branch;
	return;
}

//////////////////////////////////////////////////////////////////////////////////////
// HiC matrix functions
//////////////////////////////////////////////////////////////////////////////////////
void numPathHiC() { // computes HiC matrix for numbered paths
	int i,j,k;
	int end=0;
	ply++;
	memset(Mat, 0, sizeof Mat);
	for ( i = 0; i < TPN ; i++ )
		for ( j = 0; j < np[i].L ; j++ )
			lp[np[i].p[j]] = end++;
	for ( i = 0; i < TPN ; i++ ) {
		for ( j = 0; j < np[i].L ; j++ ) {
			fM[ply+1] = fM[ply];
			pseudoLegalMoves(np[i].p[j]);
			for ( k = fM[ply] ; k < fM[ply+1] ; k++ ) {
				Mat[lp[np[i].p[j]]*L[3]+lp[mL[k]]] = 1;
				Mat[lp[mL[k]]*L[3]+lp[np[i].p[j]]] = 1;
			}
		}
	}
	for ( i = 0; i < end; i++ ) {
		k = 0;                                    
		for ( j = 0; j < end; j++ )
			k += Mat[i*end+j];
		Mat[end*end+i] = k;
	}
	ply--;
}

void currentHiC() { // function that computes the hic matrix
	int j,k;        // and boundary for the cp (current path)
	int end=0;      // see line 602
	ply++;
	memset(Mat, 0, sizeof Mat);
	for ( j = 0; j < cp.L ; j++ )
		lp[cp.p[j]] = end++;
	for ( j = 0; j < cp.L ; j++ ) {
		fM[ply+1] = fM[ply];
		pseudoLegalMoves(cp.p[j]);
		for ( k = fM[ply] ; k < fM[ply+1] ; k++ ) {
			if (arrayCon(cp.p,cp.L,mL[k]) >= 0) {
				Mat[lp[cp.p[j]]*cp.L+lp[mL[k]]] = 1;
				Mat[lp[mL[k]]*cp.L+lp[cp.p[j]]] = 1;
			}
		}
	}
	for ( j = pwr(cp.L,2); j < pwr(cp.L,2)+cp.L; j++ )
		Mat[j] = boundaryType(cp.p[j-pwr(cp.L,2)]);
	ply--;
}

void saveSubrgMatrix( int a, int b ) { // get a region from HiC matrix
	int i,j;
	memset(Reg, 0, sizeof Reg);
	memset(bC, 0, sizeof bC);
	for ( i = a ; i <= b ; i++ )
		for ( j = a ; j <= b ; j++ )
			Reg[ (i-a)*(b-a+1)+(j-a) ] = HiC[ i*L[3]+j ];
	for ( i = a + L[6] ; i <= b + L[6] ; i++ ) 
		Reg[ pwr(( b - a + 1 ),2)+( i - a - L[6]) ] = HiC[i];
	for ( i = a + L[6] ; i <= b + L[6] ; i++ )
		bC[ i - a - L[6] ] = HiC[i];
}

//////////////////////////////////////////////////////////////////////////////////////
// functions for fetching and comparing matricies
//////////////////////////////////////////////////////////////////////////////////////
void saveNumPathHiC() {
	int i;
	for ( i = 0 ; i < L[6]+L[3] ; i++ )
		HiC[i] = Mat[i];
}

void getsnpHiC() {
	int i;
	for ( i = 0 ; i < L[6]+L[3] ; i++ )
		Mat[i] = HiC[i];
}

void getsSubrgMatrix( int size ) {
	int i;
	for ( i = 0 ; i < pwr(size,2)+size ; i++ )
		Mat[i] = Reg[i];
}

int eval() { // Mat is the current path (shape)
	int i,msd=0;
	for ( i = 0; i < pwr(cp.L,2) ; i++ )
		msd += (Reg[i]-Mat[i])*(Reg[i]-Mat[i]);
	return msd;
}

int comp() { // Mat is the current path (shape)
	int i,j,score=0;
	for ( i = 0; i < pwr(cp.L,2) ; i++ )
		if ( Mat[i] == 1 )
			if ( Reg[i] == 0 )
				score++;
	for ( i = 0; i < cp.L-1 ; i++ )
		for ( j = i+1; j < cp.L ; j++ )
			if ( cp.p[i] == cp.p[j] )
				return 1000;
	return score;
}

int norm() { // Mat is the current path (shape)
	int i,j,score=0;
	for ( i = 0; i < pwr(cp.L,2) ; i++ )
		if ( Mat[i] == 1 )
			if ( Reg[i] == 0 )
				score++;
	score = score * ( L[6] / pwr(cp.L,2));
	for ( i = 0; i < cp.L-1 ; i++ )
		for ( j = i+1; j < cp.L ; j++ )
			if ( cp.p[i] == cp.p[j] )
				return 1000;
	return score;
}

int nmsd() { // Mat is the current path (shape)
	int i,j,msd=0;
	for ( i = 0; i < cp.L-1 ; i++ )
		for ( j = i+1; j < cp.L ; j++ )
			if ( cp.p[i] == cp.p[j] )
				return 1000;
	for ( i = 0; i < pwr(cp.L,2) ; i++ )
		if ( Mat[i] != Reg[i] )
				msd++;
	return ( msd * 100 ) / pwr(cp.L,2);
}

int diag() { // Mat is the current path (shape)
	int diagMat[531441];
	int diagReg[531441];
	memset(diagMat, 0, sizeof diagMat);
	memset(diagReg, 0, sizeof diagReg);
	int i,j,m,n,score=0;
	for ( i = 0; i < cp.L; i++ )
		for ( j = 0; j < cp.L; j++ ) {
			diagMat[abs(i-j)] += Mat[i*cp.L+j];
			diagReg[abs(i-j)] += Reg[i*cp.L+j];
		}
	m = arraySum(diagMat,0,cp.L);
	n = arraySum(diagReg,0,cp.L);
	if ( m == 0 )
		m++;
	if ( n == 0 )
		n++;
	for ( i = 0; i < cp.L ; i++ )
		score += abs( ( (diagMat[i]*pwr(cp.L,2)) / m ) - ( (diagReg[i]*pwr(cp.L,2)) / n ) );
	for ( i = 0; i < cp.L-1 ; i++ )
		for ( j = i+1; j < cp.L ; j++ )
			if ( cp.p[i] == cp.p[j] )
				return 1000;
	return score;
}

int verify() {
	int i,msd=0;
	for ( i = 0; i < L[6] ; i++ )
		msd += (HiC[i]-Mat[i])*(HiC[i]-Mat[i]);
	return msd;
}

/////////////////////////////////////////////////////////////////////////////////////////
// algorithms for shuffling an initial condition
/////////////////////////////////////////////////////////////////////////////////////////
void multiMansfield( int iterations ) {
	int i;
	for ( i = 0 ; i < iterations ; i++ ) {
		clearPly();
		int pathGoto = -1;
		int pathFrom = rand() % TPN;
		int branch = -1;
		if ( rand() % 2 )
			arrayRev(np[pathFrom].p,0,np[pathFrom].L);
		legalMoves(np[pathFrom].p[np[pathFrom].L-1],-1);
		shuffle();
		int choice = rand()%6;    
		if ( choice < fM[ply+1] ) {
			do { pathGoto++;
				 branch = arrayCon(np[pathGoto].p,np[pathGoto].L,mL[choice]);
			   } while ( branch == -1 );
			if ( pathFrom == pathGoto ) {
				 arrayRev(np[pathFrom].p,0,np[pathFrom].L);
				 branch = arrayCon(np[pathGoto].p,np[pathGoto].L,mL[choice]);
				 arrayRev(np[pathGoto].p,0,branch);
			} else {  
				 if ( rand() % 2 ) 
					 arrayRev(np[pathGoto].p,0,np[pathGoto].L);
				 branch = arrayCon(np[pathGoto].p,np[pathGoto].L,mL[choice]);
				 if ( ( branch > 0 ) & ( branch < np[pathGoto].L-1 ) )
					swapSubPaths(pathFrom,pathGoto,branch);
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// gets random space filling curve
/////////////////////////////////////////////////////////////////////////////////////////
void randomHiC() {
	initConfig();	
	dividecwp();
	multiMansfield(30000);
	numPathHiC();
	saveNumPathHiC();
}

void startPosition() {
	initConfig();	
	dividecwp();
}

/////////////////////////////////////////////////////////////////////////////////////////
// boundary segmentation function i.e. voodoo
/////////////////////////////////////////////////////////////////////////////////////////
void segmentBoundary( int maxSize ) { // segments the boundary regions into pieces
	int i,j,k;	
	memset(bS, HiC[L[6]] == 6 ? -1 : 0, sizeof bS);
	for ( i = L[6] ; i < L[6]+L[3]-1 ; i++ ) {
		j = i-L[6];
		bS[j+1] = bS[j];
		clearPly(); // we need to know the boundary types of each neighbor
		pseudoLegalMoves( HiC[i] );
		for ( k = fM[ply] ; k < fM[ply+1] ; k++ )
			mL[k] = boundaryType(mL[k]);
		if ( ( HiC[i] == 6 ) & ( HiC[i+1] < 6 ) )	
			bS[j+1] = bS[j] + 1;
		else if ( ( HiC[i+1] == 6 ) & ( HiC[i] < 6 ) )	
			bS[j+1] = bS[j] + 1;
		else if ( HiC[j*L[3]+(j+1)] == 0 )
			bS[j+1] = bS[j] + 2;
		else if ( arrayCon(mL,fM[ply+1],boundaryType(HiC[i+1])) < 0 )
			bS[j+1] = bS[j] + 2;
		else
			continue;
	}
	for ( i = 0 ; i < L[3] ; i++ )
		if ( bS[i]%2 == 0 )
			bS[i] = bS[i]/2 + 1;
		else
			bS[i] = 0;
	j = arrayMax(bS,0,L[3]);
	for ( i = 1 ; i <= arrayMax(bS,0,L[3]) ; i++ )
		if ( arrayCnt(bS,i,L[3]) >= maxSize ) {
			for ( k = arrayFst(bS,i,L[3]) + arrayCnt(bS,i,L[3])/((arrayCnt(bS,i,L[3])/maxSize) + 1) ; k <= arrayLst(bS,i,L[3]) ; k++ )
				bS[k] = j + 1;
			j = arrayMax(bS,0,L[3]);
		}
}

int scoreBoundary( int n ) {
	int i,score;
	score = n*L[3];
	for ( i = 0 ; i < L[3] ; i++ )
		score = score - abs( ( n-Reg[i+n*L[3]] )-( n*(Mat[i+L[6]]/6) ) );
	return score/n;
}

/////////////////////////////////////////////////////////////////////////////////////////
// shape library management
/////////////////////////////////////////////////////////////////////////////////////////
void deleteShape( int n ) {
	int i,j;
	j = fS[n+1] - fS[n];
	for ( i = fS[n] ; i < fS[cpm]-j ; i++ )
		sL[i] = sL[i+j];
	cpm--;
	for ( i = n ; i < cpm ; i++ ) {
		sR[i] = sR[i+1];
		sE[i] = sE[i+1];
		fS[i+1] = fS[i+2]-j;
	}
}
		
int nextShape( int n ) {
	int i;
	for ( i = 0 ; i < cpm ; i++ )
		if ( sR[i] == n )
			return i;
	return -1;
}

int countRegions( int n ) {
	int i,total=0;
	for ( i = 0 ; i < cpm ; i++ )
		if ( sR[i] == n )
			total++;
	return total;
	}

void deleteRegion( int region ) {
	while ( countRegions(region) > 0 )
		deleteShape( nextShape(region) ); 
}
	
void resetsRegion( int i, int j, int k ) {
	int a;
	for ( a = 0 ; a < cpm ; a++ )
		if ( sR[a] == i )
			sR[a] = j;
	for ( a = 0 ; a < L[3] ; a++ )
		if ( bS[a] == k )
			bS[a] = j;
}

int worsteShape() {
	int i,j=0,k=-1;
	for ( i = 0 ; i < cpm ; i++ )
		if ( sE[i] > k ) {
			j = i;
			k = sE[i];
		}
	return j;
}

int bestShape() {
	int i,j=0,k=10000000;
	for ( i = 0 ; i < cpm ; i++ )
		if ( sE[i] < k ) {
			j = i;
			k = sE[i];
		}
	return j;
}

/////////////////////////////////////////////////////////////////////////////////////////
// subregion search functions
/////////////////////////////////////////////////////////////////////////////////////////
void evalSR() {
	int cpeval;
	currentHiC(); // stored in MAT line line 375
	cpeval = (*evalPtr)();
	if ( cpeval <= TOL )
		saveSR(cpeval);
	cpn++;
}

void pathsMatched( int a, int b, int c, int r, int s ) {
	int i;
	initialize();
	savePtr = &evalSR; 
	cp.L = b - a + 1;
	saveSubrgMatrix( a, b );
	search(1000000,c);
	for ( i = s ; i < cpm ; i++ )
		sR[i] = r;
}

void enumerateSR() {
	int i;
	printf("\nboundary shapes library (%d total)\n",cpm);
	for ( i = 1 ; i <= arrayMax(bS,0,L[3])+1 ; i++ )
		printf("\tnumber of paths in region %d is %d\n",i,countRegions(i));
	}

/////////////////////////////////////////////////////////////////////////////////////////
// functions for handleing pairs of shapes
/////////////////////////////////////////////////////////////////////////////////////////
void joinReg ( int a, int b ) {
	clearPly();
	mL[fM[ply+1]++] = a;
	if ( b > -1 )
		mL[fM[ply+1]++] = b;
	int i,j,k,sq=0,end=0;
	memset(Reg, 0, sizeof Reg);
	for ( i = fM[ply] ; i < fM[ply+1] ; i++ )
		sq += arrayCnt(bS,mL[i],L[3]); 
	for ( k = 0 ; k < L[6] ; k++ )
		for ( i = fM[ply] ; i < fM[ply+1] ; i++ )
			for ( j = fM[ply] ; j < fM[ply+1] ; j++ )
				if ( ( bS[k/L[3]] == mL[i] ) & ( bS[k%L[3]] == mL[j] ) )
					Reg[end++] = HiC[k];
	for ( k = L[6] ; k < L[6]+L[3] ; k++ )
		for ( i = fM[ply] ; i < fM[ply+1] ; i++ )
			if ( bS[k-L[6]] == mL[i] )
				Reg[end++] = HiC[k];
	getsSubrgMatrix(sq);
}

int bestRegionPair( int region,int verbose ) {
	int i,j,k;
	int* pairs;
	pairs = (int*)malloc((region * region) * sizeof(int));
	memset(pairs, 0, (region * region) * sizeof(int));
	for ( k = 0 ; k < L[6] ; k++ )
		for ( i = 0; i < region ; i++ )
			for ( j = i + 1 ; j < region ; j++ )
				if ( ( bS[k/L[3]] == i+1 ) & ( bS[k%L[3]] == j+1 ) )
					pairs[region*i+j]+=HiC[k];
	k = arrayFst(pairs,arrayMax(pairs,0,region*region),region*region);
	if ( verbose )
		printRegion( pairs,region );
	free(pairs);
	return k;
}

int mergeRegions( int a, int b ) {
	int k,cpeval,m=0,n=0;
	int region = arrayMax(bS,0,L[3]);
	cp.L = 0;
	for ( k = 0 ; k < L[3] ; k++ ) {
		if ( bS[k] == sR[a] )
			cp.p[cp.L++] = sL[fS[a] + n++];
		if ( bS[k] == sR[b] )
			cp.p[cp.L++] = sL[fS[b] + m++];
	}
	currentHiC();
	cpeval = (*evalPtr)();
	if ( cpeval <= TOL ) {
		saveSR(cpeval);
		sR[cpm-1] = region+1;
		cpn++;
		if ( cpm > 100000000 )
			deleteShape( worsteShape() );
		if( fS[cpm] > 1000000000 )
			deleteShape( worsteShape() );
	}
	return 0;
}
	
/////////////////////////////////////////////////////////////////////////////////////////
// input functions
/////////////////////////////////////////////////////////////////////////////////////////
int readHicMatrix( char *matFile ) {
	int lineNum = 0,fieldNum = 0;
	FILE * fp;
	char * lineArray;
	char * line = NULL;
	size_t size = 0;
	ssize_t input;
	fp = fopen(matFile, "r");
		if (fp == NULL)
			exit(EXIT_FAILURE);
	while ((input = getline(&line, &size, fp)) != -1) {
		lineArray = strtok (line," ");
		while (lineArray != NULL) {
			HiC[fieldNum++] = (int) strtol(lineArray, (char **)NULL, 10);
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
			Reg[fieldNum++] = (int) strtol(lineArray, (char **)NULL, 10);
			lineArray = strtok (NULL," ");
			}
			lineNum++;
		}		
	fclose(fp);
	return lineNum;
}

/////////////////////////////////////////////////////////////////////////////////////////
// split and merge paths
/////////////////////////////////////////////////////////////////////////////////////////
void getShape( int n ) {
	int i;
	for ( i = fS[n] ; i < fS[n+1] ; i++ )
		cp.p[i] = sL[i];
	cp.n = n;
}

void divideSR( int n ) {
	cpn = 0;
	int i,j = 0;
	for ( i = fS[n]+1 ; i < fS[n+1] ; i++ ) {
		clearPly();
		np[cpn].p[j++] = sL[i-1];
		pseudoLegalMoves(sL[i-1]);
		if ( arrayCon( mL+fM[ply],fM[ply+1]-fM[ply],sL[i] ) == -1 ) {
			np[cpn].L = j;
			np[cpn].n = cpn;
			j = 0;
			cpn++;
		}
	}
	np[cpn].p[j] = sL[fS[n+1]-1];
	np[cpn].L = ++j;
	np[cpn].n = cpn;
	cpn++;
}

int aggregate( int verbose ) {
	int a,b,c;
	int i,j,k;
	while (arrayMax(bS,0,L[3]) > 1) {
		c = arrayMax(bS,0,L[3]);
		k = bestRegionPair( c,verbose );
		i = ( k / c ) + 1;
		j = ( k % c ) + 1;
		if ( verbose )
			printf("max is pair %d,%d making submatrix...\n\n",i,j);
		joinReg( i,j ); // make the pair matrix
		cpn = 0;
		for ( a = 0 ; a < cpm ; a++ )
			if ( sR[a] == i )
				for ( b = 0 ; b < cpm ; b++ )
					if ( sR[b] == j )
						if ( mergeRegions( a, b ) )
							return 1;				
		if ( verbose ) {
			enumerateSR();
			printf("\nremoving pairs %d and %d\n",i,j);
			printsbSegments();
		}
		deleteRegion(i);
		deleteRegion(j);
		resetsRegion(arrayMax(bS,0,L[3])+1,i,j);
	}
	return 0;
}

void convolution( int convolve,char *matFile,char *bndFile  ) {
	int i,j;
	FILE *matr = fopen(matFile,"w");
	if (matr == NULL)
		exit(1);
	FILE *bndr = fopen(bndFile,"w");
	if (bndr == NULL)
		exit(1);
	memset(Reg, 0, sizeof Reg);
	for ( i = 0 ; i < convolve ; i++ ) {
		randomHiC();
		for ( j = 0 ; j < L[6] ; j++ )
			Reg[j] = Reg[j] + HiC[j];
		for ( j = L[6] ; j < L[6]+L[3] ; j++ )
			Reg[j] = Reg[j] + ( HiC[j] < 6 ? 1 : 0 ) ;
		for ( j = 0 ; j < L[3] ; j++ )
			fprintf(bndr,"%d%s",HiC[L[6]+j],( j%L[3] == L[3]-1 ) ? "\n" : " ");
	}
	for ( j = 0 ; j < L[6] ; j++ )
		fprintf(matr,"%d%s",Reg[j],( j%L[3] == L[3]-1 ) ? "\n" : " ");
	for ( j = 0 ; j < L[3] ; j++ )
		fprintf(bndr,"%d%s",Reg[L[6]+j],( j%L[3] == L[3]-1 ) ? "\n" : " ");
	fclose(bndr);
	fclose(matr);
}

void fitBoundary( FILE *out,int n,int max ) {
	int i=0;
	startPosition();
	int bestScore = 0;
	while ( bestScore < max ) {
		numPathHiC();
		if ( scoreBoundary( n ) > bestScore ) {
			bestScore = scoreBoundary( n );
			saveNumPathHiC();
			fprintf(out,"edge_score_%d\t%d\n",i++,bestScore);
		}
		multiMansfield(1);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
// main programming section
/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	//////////////////////////////////////////////////////////
	// variables with default values
	////////////////////////////////////////////////////////// 
	int c;               // variable for getopt.h
	int pathSize   = 4;  // total path length = 4*4*4 
	int numPaths   = 3;  // number of paths in the HiC matrix
	int verbose    = 0;  // flag for verbose
	int segsize    = 7;  // ends up being segsize + 1
	int bndopt     = -1; // flag/score for boundary optimization
	int convolve   = -1; // number of matricies to convolve
	int deconvolve = -1; // deconvolution mode
	evalPtr = &eval;
	
	// the file paths are a bit sloppy
	char matFile[96] = "";
	char bndFile[96] = "";
	char outMatr[96] = "";
	char outFold[96] = "";
	char outLogs[96] = "";
	char outTarg[96] = "";
	char saveDir[96] = "";
	
	// make sure the shape library is initialized to zeros
	memset(sL, 0, sizeof sL);
	memset(sE, 0, sizeof sE);
	memset(sR, 0, sizeof sR);
	memset(fS, 0, sizeof fS);
	
	//////////////////////////////////////////////////////////
	// https://azrael.digipen.edu/~mmead/www/Courses/CS180/getopt.html#OPTOPTARGS
	////////////////////////////////////////////////////////// 
	while (1) 
	{
		int option_index = 0;
		static struct option long_options[] = 
		{
			{"pathSize",    required_argument, NULL,  's'},
			{"segsize",     required_argument, NULL,  'z'},
			{"numPaths",    required_argument, NULL,  'n'},
			{"evaluation",  required_argument, NULL,  'e'},
			{"saveDir",     required_argument, NULL,  'o'},
			{"help",        no_argument,       NULL,  'h'},
			{"verbose",     no_argument,       NULL,  'v'},
			{"convolve",    optional_argument, NULL,  'c'},
			{"deconvolve",  optional_argument, NULL,  'd'},
			{"bndopt",      optional_argument, NULL,  'b'},
			{NULL,          0,                 NULL,    0}
		};

		c = getopt_long(argc, argv, "-:s:z:n:e:o:hvc::d::b::", long_options, &option_index);
		if (c == -1)
		  break;

		switch (c) 
		{
		  case 0:
			printf("long option %s", long_options[option_index].name);
			if (optarg)
			   printf(" with arg %s", optarg);
			printf("\n");
			break;

		  case 1:
			printf("regular argument '%s'\n", optarg);
			break;

		  case 's':
			printf("setting pathSize to '%s'\n", optarg);
			pathSize = atoi(optarg);
			break;
			
		  case 'z':
			printf("setting segmentation size to '%s'\n", optarg);
			segsize = atoi(optarg);
			break;
			
		  case 'n':
			printf("setting numPaths to '%s'\n", optarg);
			numPaths = atoi(optarg);
			break;
			
		  case 'o':
			printf("setting output directory to '%s'\n", optarg);
			strcpy(saveDir,"./data/");
			strcat(saveDir,optarg);
			strcpy(matFile,saveDir);
			strcpy(bndFile,saveDir);
			strcpy(outMatr,saveDir);
			strcpy(outFold,saveDir);
			strcpy(outLogs,saveDir);
			strcpy(outTarg,saveDir);
			strcat(matFile,"/inputMatrix.txt");
			strcat(bndFile,"/inputBounds.txt");
			strcat(outMatr,"/outputMatrix.txt");
			strcat(outFold,"/outputFolded.txt");
			strcat(outLogs,"/outputLogged.txt");
			strcat(outTarg,"/outputTarget.txt");
			break;
			
		  case 'e':
			if ( !strcmp(optarg,"eval") ) {
				printf("setting eval function to eval\n");
				evalPtr = &eval;
			} else if ( !strcmp(optarg,"comp") ) {
				printf("setting eval function to comp\n");
				evalPtr = &comp;
			} else if ( !strcmp(optarg,"diag") ) {
				printf("setting eval function to diag\n");
				evalPtr = &diag;
			} else if ( !strcmp(optarg,"norm") ) {
				printf("setting eval function to norm\n");
				evalPtr = &norm;
			}  else if ( !strcmp(optarg,"nmsd") ) {
				printf("setting eval function to nmsd\n");
				evalPtr = &nmsd;
			} else
				evalPtr = &eval;
			break;
						
		  case 'v':
			verbose = 1;
			printf("set verbose output\n");
			break;

		  case 'c':
		    optarg ? convolve = atoi(optarg) : ( convolve = 1 );
            printf("convolving %d random matricies\n",convolve);
            break;
			
		  case 'd':
		    optarg ? deconvolve = atoi(optarg) : ( deconvolve = 0 );
            printf("deconvolving with tolerance set to %d\n",deconvolve);
            break;

		  case 'b':
		    optarg ? bndopt = atoi(optarg) : ( bndopt = ( pwr(pathSize,3)/3 )*2 );
            printf("simulating boundary to max score %d\n",bndopt);
            break;			

		  case 'h':		  
			printf("\n\nCode for paper: My Paper Title Here\n\n");
		    printf(" --pathSize (s) size of the box\n");
			printf(" --segsize (z) size of boundary segmentation\n");
			printf(" --numPaths (n) num of paths withing the box\n");
			printf(" --savedir (o) where to save and read the input matrix\n");
			printf(" --evaluation (e) type of evaluation function\n");
			printf("   -e eval is quick and suitable for 1 matrix\n");
			printf("   -e comp is slow and used for deconvolution\n");
			printf(" --bndopt (b) turns on boundary simulation\n");
			printf(" --help (h) print help message\n");
	        printf(" --verbose (v) prints verbose output to stdout\n");
			printf(" --convolve (c) convolved random matricies\n");
			printf("   -c outputs 1 matrix and boundary\n");
			printf("   -c2 convolves 2 matricies and boundary\n");
			printf("   -cn convolves n matrix and boundary\n");
			printf(" --deconvolve (d) deconvolve random matricies\n");
			printf("   -d looks for an exact match withing the convolved matrix\n");
			printf("   -d2 allows up to to mismatches in the solution matrix\n");
			printf("   -dn allows up to n mismatches in the solution matrix\n\n\n");
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

	if ( saveDir[0] == '\0' ) {
		strcpy(matFile,"./data/results/inputMatrix.txt");
		strcpy(bndFile,"./data/results/inputBounds.txt");
		strcat(outMatr,"./data/results/outputMatrix.txt");
		strcat(outFold,"./data/results/outputFolded.txt");
		strcat(outLogs,"./data/results/outputLogged.txt");
		strcat(outTarg,"./data/results/outputTarget.txt");
	}
	
	//////////////////////////////////////////////////////////
    // update a few global variables
    ////////////////////////////////////////////////////////// 
	int i,j,k,a,b;
	int t = time(NULL);
	srand(t); 
	TOL = deconvolve;
	TPN = numPaths;
	for ( i = 0 ; i < 7 ; i++ )
		L[i] = pwr(pathSize,i);
	
	//////////////////////////////////////////////////////////
	// if we are generating a random matrix
    ////////////////////////////////////////////////////////// 
	if ( convolve > 0 )
		convolution(convolve,matFile,bndFile);
	if ( deconvolve < 0 )
		return 0;

	//////////////////////////////////////////////////////////
    // the rest of the code if for deconvolution, first open outfiles
    ////////////////////////////////////////////////////////// 
	FILE *fold = fopen(outFold, "w");
	if (fold == NULL)
		exit(1);
	FILE *matr = fopen(outMatr, "w");
	if (matr == NULL) 
		exit(1);
	FILE *outf = fopen(outTarg, "w");
	if (outf == NULL)
		exit(1);
	FILE *logs = fopen(outLogs, "w");
	if (outf == NULL)
		exit(1);
			
	//////////////////////////////////////////////////////////
    // read in the hic matrix we want to match or make one
    //////////////////////////////////////////////////////////
	b = readBoundary(bndFile);	
	if ( bndopt >= 0 )
		fitBoundary( logs,b-1,bndopt );
	else
		for ( i = 0 ; i < L[3] ; i++ )
			HiC[L[6]+i] = Reg[i+0*L[3]];
	readHicMatrix(matFile);
	getsnpHiC();
	printsCwpMatrix( outf,L[3] );
	segmentBoundary( segsize + 1 );
	if ( verbose ) {
		printsbSegments();
		printf("\nread HiC matrix from %s file\n",matFile);
		printf("read boundary from %s file\n",bndFile);
	}
	
	
	//////////////////////////////////////////////////////////
	// part one of the algorithm begins with library of boundary shapes
    ////////////////////////////////////////////////////////// 
	int starttime = get_ms();
	for ( i = 1 ; i <= arrayMax(bS,0,L[3]) ; i++ ) {
		a = arrayFst(bS,i,L[3]);
		b = arrayLst(bS,i,L[3]);
		for ( j = 0 ; j < L[3] ; j++ ) {
			if ( boundaryType(j) != HiC[L[6]+a] )
				continue;
			pathsMatched(a,b,j,i,cpm);
		}
	}
	if ( verbose )
		enumerateSR();

	//////////////////////////////////////////////////////////
	// part two of the algorithm find shapes that share the most contacts
    //////////////////////////////////////////////////////////
	if ( aggregate( verbose ) )
		return 0;
	if ( verbose )
		enumerateSR();

	//////////////////////////////////////////////////////////
	// part three of the algorithm, fill in missing verticies one by one
    //////////////////////////////////////////////////////////
	j = 2;
	for ( i = 0 ; i < L[3] ; i++ )
		if ( bS[i] == 0 )
			bS[i] = j++;
	cpn = 0;
	for ( k = 2 ; k < j ; k++ )
		for ( i = 0 ; i < L[3] ; i++ ) {
			cp.p[0] = i;
			cp.L = 1;
			saveSR(0);
			sR[cpm-1] = k;
			cpn++;
		}
		
	//////////////////////////////////////////////////////////
	// part four of the algorithm find the shapes that share the most contacts
    //////////////////////////////////////////////////////////	
	if ( aggregate( verbose ) )
		return 0;
	if ( verbose )
		enumerateSR( outf );
	
	//////////////////////////////////////////////////////////
	// cleanup and output the solution
    //////////////////////////////////////////////////////////
	fprintf(logs,"total_time\t%d\n", (get_ms() - starttime));
	i = bestShape();
	getShape(i);
	currentHiC();
	if ( verbose )
		printf("\nreconstructed %d boundary paths, eval is %d\n\n",cpm,verify());
	printsCwpMatrix( matr,L[3] );
	divideSR(i);
	for ( i = 0; i < cpn ; i++ ) {
		cp = np[i];
		printPath(fold);
	}
		
	//////////////////////////////////////////////////////////
	// finish up the logs
    //////////////////////////////////////////////////////////
	for ( i = 0; i <= TOL ; i++ ) {
		k = 0;
		for ( j = 0 ; j < cpm ; j++ )
			if ( sE[j] == i )
				k++;
		fprintf(logs,"num_eval_%d\t%d\n",i,k);
	}
	fprintf(logs,"random_seed\t%d\n",t);
	fprintf(logs,"tolerance\t%d\n",TOL);
	fprintf(logs,"number_paths\t%d\n",TPN);
	fprintf(logs,"square_size\t%d\n",L[1]);
	fprintf(logs,"seg_size\t%d\n",segsize);
	
	//////////////////////////////////////////////////////////
	// exit the program
    //////////////////////////////////////////////////////////
	fclose(outf );
	fclose(matr);
	fclose(fold);
	fclose(logs);
	return 0;
}


