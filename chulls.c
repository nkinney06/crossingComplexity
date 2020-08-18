#include <stdlib.h>   /* malloc */
#include <stdio.h>    /* printf */
#include <math.h>

/* Define Boolean type */
typedef	enum { FALSE, TRUE }	bool;

////////////////////////////////////////////////////////////////////////////////
/* Macros and definitions */
////////////////////////////////////////////////////////////////////////////////
#define X   0
#define Y   1
#define Z   2
#define ONHULL   	TRUE
#define REMOVED  	TRUE
#define VISIBLE  	TRUE
#define PROCESSED	TRUE
#define SAFE		1000000

#define SWAP(t,x,y)	{ t = x; x = y; y = t; }

#define NEW(p,type)	if ((p=(type *) malloc (sizeof(type))) == NULL) { \
				printf ("Out of Memory!\n"); \
				exit(0); \
			}
			
#define FREE(p)		if (p) { free ((char *) p); p = NULL; }
			
#define ADD( head, p )  if ( head )  { \
				p->next = head; \
				p->prev = head->prev; \
				head->prev = p; \
				p->prev->next = p; \
			} \
			else { \
				head = p; \
				head->next = head->prev = p; \
			}
			
#define DELETE( head, p ) if ( head )  { \
				if ( head == head->next ) \
					head = NULL;  \
				else if ( p == head ) \
					head = head->next; \
				p->next->prev = p->prev;  \
				p->prev->next = p->next;  \
				FREE( p ); \
			} 

////////////////////////////////////////////////////////////////////////////////
/* Define structures for vertices, edges and faces */
////////////////////////////////////////////////////////////////////////////////
typedef struct tVertexStructure tsVertex;
typedef tsVertex *tVertex;
typedef struct tEdgeStructure tsEdge;
typedef tsEdge *tEdge;
typedef struct tFaceStructure tsFace;
typedef tsFace *tFace;

struct vector {
    int x;     
	int y;
	int z;
	int invader;
	int inhull;
};

struct tVertexStructure {
   int      v[3];
   int	    vnum;
   tEdge    duplicate;	/* pointer to incident cone edge (or NULL) */
   bool     onhull;		/* T iff point on hull. */
   bool	    mark;		/* T iff point already processed. */
   tVertex  next, prev;
};

struct tEdgeStructure {
   tFace    adjface[2];
   tVertex  endpts[2];
   tFace    newface;    /* pointer to incident cone face. */
   bool     delete;		/* T iff edge should be delete. */
   tEdge    next, prev;
};

struct tFaceStructure {
   tEdge    edge[3];
   tVertex  vertex[3];
   bool	    visible;	/* T iff face visible from new point. */
   tFace    next, prev;
};

////////////////////////////////////////////////////////////////////////////////
/* Global variable definitions */
////////////////////////////////////////////////////////////////////////////////
tVertex vertices = NULL;
tEdge edges    	 = NULL;
tFace faces    	 = NULL;
bool debug = FALSE;
bool check = FALSE;

tVertex	MakeNullVertex( void )
{
   tVertex  v;
   NEW( v, tsVertex );
   v->duplicate = NULL;
   v->onhull = !ONHULL;
   v->mark = !PROCESSED;
   ADD( vertices, v );
   return v;
}

/*---------------------------------------------------------------------
Collinear checks to see if the three points given are collinear,
by checking to see if each element of the cross product is zero.
---------------------------------------------------------------------*/
bool	Collinear( tVertex a, tVertex b, tVertex c )
{
   return 
         ( c->v[Z] - a->v[Z] ) * ( b->v[Y] - a->v[Y] ) -
         ( b->v[Z] - a->v[Z] ) * ( c->v[Y] - a->v[Y] ) == 0
      && ( b->v[Z] - a->v[Z] ) * ( c->v[X] - a->v[X] ) -
         ( b->v[X] - a->v[X] ) * ( c->v[Z] - a->v[Z] ) == 0
      && ( b->v[X] - a->v[X] ) * ( c->v[Y] - a->v[Y] ) -
         ( b->v[Y] - a->v[Y] ) * ( c->v[X] - a->v[X] ) == 0  ;
}

/*---------------------------------------------------------------------
MakeNullEdge creates a new cell and initializes all pointers to NULL
and sets all flags to off.  It returns a pointer to the empty cell.
---------------------------------------------------------------------*/
tEdge 	MakeNullEdge( void )
{
   tEdge  e;

   NEW( e, tsEdge );
   e->adjface[0] = e->adjface[1] = e->newface = NULL;
   e->endpts[0] = e->endpts[1] = NULL;
   e->delete = !REMOVED;
   ADD( edges, e );
   return e;
}

/*--------------------------------------------------------------------
MakeNullFace creates a new face structure and initializes all of its
flags to NULL and sets all the flags to off.  It returns a pointer
to the empty cell.
---------------------------------------------------------------------*/
tFace 	MakeNullFace( void )
{
   tFace  f;
   int    i;

   NEW( f, tsFace);
   for ( i=0; i < 3; ++i ) {
      f->edge[i] = NULL;
      f->vertex[i] = NULL;
   }
   f->visible = !VISIBLE;
   ADD( faces, f );
   return f;
}

/*---------------------------------------------------------------------
MakeFace creates a new face structure from three vertices (in ccw
order).  It returns a pointer to the face.
---------------------------------------------------------------------*/
tFace   MakeFace( tVertex v0, tVertex v1, tVertex v2, tFace fold )
{
   tFace  f;
   tEdge  e0, e1, e2;

   /* Create edges of the initial triangle. */
   if( !fold ) {
     e0 = MakeNullEdge();
     e1 = MakeNullEdge();
     e2 = MakeNullEdge();
   }
   else { /* Copy from fold, in reverse order. */
     e0 = fold->edge[2];
     e1 = fold->edge[1];
     e2 = fold->edge[0];
   }
   e0->endpts[0] = v0;              e0->endpts[1] = v1;
   e1->endpts[0] = v1;              e1->endpts[1] = v2;
   e2->endpts[0] = v2;              e2->endpts[1] = v0;
	
   /* Create face for triangle. */
   f = MakeNullFace();
   f->edge[0]   = e0;  f->edge[1]   = e1; f->edge[2]   = e2;
   f->vertex[0] = v0;  f->vertex[1] = v1; f->vertex[2] = v2;
	
   /* Link edges to face. */
   e0->adjface[0] = e1->adjface[0] = e2->adjface[0] = f;
	
   return f;
}

/*---------------------------------------------------------------------
VolumeSign returns the sign of the volume of the tetrahedron determined by f
and p.  VolumeSign is +1 iff p is on the negative side of f,
where the positive side is determined by the rh-rule.  So the volume 
is positive if the ccw normal to f points outside the tetrahedron.
The final fewer-multiplications form is due to Bob Williamson.
---------------------------------------------------------------------*/
int  VolumeSign( tFace f, tVertex p )
{
   double  vol;
   int     voli;
   double  ax, ay, az, bx, by, bz, cx, cy, cz;

   ax = f->vertex[0]->v[X] - p->v[X];
   ay = f->vertex[0]->v[Y] - p->v[Y];
   az = f->vertex[0]->v[Z] - p->v[Z];
   bx = f->vertex[1]->v[X] - p->v[X];
   by = f->vertex[1]->v[Y] - p->v[Y];
   bz = f->vertex[1]->v[Z] - p->v[Z];
   cx = f->vertex[2]->v[X] - p->v[X];
   cy = f->vertex[2]->v[Y] - p->v[Y];
   cz = f->vertex[2]->v[Z] - p->v[Z];

   vol =   ax * (by*cz - bz*cy)
         + ay * (bz*cx - bx*cz)
         + az * (bx*cy - by*cx);

   /* The volume should be an integer. */
   if      ( vol >  0.5 )  return  1;
   else if ( vol < -0.5 )  return -1;
   else                    return  0;
}

/*---------------------------------------------------------------------
 DoubleTriangle builds the initial double triangle.  It first finds 3 
 noncollinear points and makes two faces out of them, in opposite order.
 It then finds a fourth point that is not coplanar with that face.  The  
 vertices are stored in the face structure in counterclockwise order so 
 that the volume between the face and the point is negative. Lastly, the
 3 newfaces to the fourth point are constructed and the data structures
 are cleaned up. 
---------------------------------------------------------------------*/
void    DoubleTriangle( void )
{
   tVertex  v0, v1, v2, v3, t;
   tFace    f0, f1 = NULL;
   tEdge    e0, e1, e2, s;
   int      vol;
	
   /* Find 3 noncollinear points. */
   v0 = vertices;
   while ( Collinear( v0, v0->next, v0->next->next ) )
      if ( ( v0 = v0->next ) == vertices )
         printf("DoubleTriangle:  All points are Collinear!\n"), exit(0);
   v1 = v0->next;
   v2 = v1->next;
	
   /* Mark the vertices as processed. */
   v0->mark = PROCESSED;
   v1->mark = PROCESSED;
   v2->mark = PROCESSED;
   
   /* Create the two "twin" faces. */
   f0 = MakeFace( v0, v1, v2, f1 );
   f1 = MakeFace( v2, v1, v0, f0 );
   
   /* Link adjacent face fields. */
   f0->edge[0]->adjface[1] = f1;
   f0->edge[1]->adjface[1] = f1;
   f0->edge[2]->adjface[1] = f1;
   f1->edge[0]->adjface[1] = f0;
   f1->edge[1]->adjface[1] = f0;
   f1->edge[2]->adjface[1] = f0;
   
   /* Find a fourth, noncoplanar point to form tetrahedron. */
   v3 = v2->next;
   vol = VolumeSign( f0, v3 );
   while ( !vol )   {
      if ( ( v3 = v3->next ) == v0 ) 
         printf("DoubleTriangle:  All points are coplanar!\n"), exit(0);
      vol = VolumeSign( f0, v3 );
   }
   
   vertices = v3;
}

/*---------------------------------------------------------------------
MakeCcw puts the vertices in the face structure in counterclock wise 
order.  We want to store the vertices in the same 
order as in the visible face.  The third vertex is always p.

Although no specific ordering of the edges of a face are used
by the code, the following condition is maintained for each face f:
one of the two endpoints of f->edge[i] matches f->vertex[i]. 
But note that this does not imply that f->edge[i] is between
f->vertex[i] and f->vertex[(i+1)%3].  (Thanks to Bob Williamson.)
---------------------------------------------------------------------*/
void	MakeCcw( tFace f, tEdge e, tVertex p )
{
   tFace  fv;   /* The visible face adjacent to e */
   int    i;    /* Index of e->endpoint[0] in fv. */
   tEdge  s;	/* Temporary, for swapping */
      
   if  ( e->adjface[0]->visible )      
        fv = e->adjface[0];
   else fv = e->adjface[1];
       
   /* Set vertex[0] & [1] of f to have the same orientation
      as do the corresponding vertices of fv. */ 
   for ( i=0; fv->vertex[i] != e->endpts[0]; ++i )
      ;
   /* Orient f the same as fv. */
   if ( fv->vertex[ (i+1) % 3 ] != e->endpts[1] ) {
      f->vertex[0] = e->endpts[1];  
      f->vertex[1] = e->endpts[0];    
   }
   else {                               
      f->vertex[0] = e->endpts[0];   
      f->vertex[1] = e->endpts[1];      
      SWAP( s, f->edge[1], f->edge[2] );
   }
   /* This swap is tricky. e is edge[0]. edge[1] is based on endpt[0],
      edge[2] on endpt[1].  So if e is oriented "forwards," we
      need to move edge[1] to follow [0], because it precedes. */
   
   f->vertex[2] = p;
}

/*---------------------------------------------------------------------
MakeConeFace makes a new face and two new edges between the 
edge and the point that are passed to it. It returns a pointer to
the new face.
---------------------------------------------------------------------*/
tFace	MakeConeFace( tEdge e, tVertex p )
{
   tEdge  new_edge[2];
   tFace  new_face;
   int 	  i, j;

   /* Make two new edges (if don't already exist). */
   for ( i=0; i < 2; ++i ) 
      /* If the edge exists, copy it into new_edge. */
      if ( !( new_edge[i] = e->endpts[i]->duplicate) ) {
	  /* Otherwise (duplicate is NULL), MakeNullEdge. */
	  new_edge[i] = MakeNullEdge();
	  new_edge[i]->endpts[0] = e->endpts[i];
	  new_edge[i]->endpts[1] = p;
	  e->endpts[i]->duplicate = new_edge[i];
      }

   /* Make the new face. */
   new_face = MakeNullFace();   
   new_face->edge[0] = e;
   new_face->edge[1] = new_edge[0];
   new_face->edge[2] = new_edge[1];
   MakeCcw( new_face, e, p ); 
        
   /* Set the adjacent face pointers. */
   for ( i=0; i < 2; ++i )
      for ( j=0; j < 2; ++j )  
	 /* Only one NULL link should be set to new_face. */
	 if ( !new_edge[i]->adjface[j] ) {
	    new_edge[i]->adjface[j] = new_face;
	    break;
	 }
        
   return new_face;
}

/*---------------------------------------------------------------------
AddOne is passed a vertex.  It first determines all faces visible from 
that point.  If none are visible then the point is marked as not 
onhull.  Next is a loop over edges.  If both faces adjacent to an edge
are visible, then the edge is marked for deletion.  If just one of the
adjacent faces is visible then a new face is constructed.
---------------------------------------------------------------------*/
bool 	AddOne( tVertex p )
{
   tFace  f; 
   tEdge  e, temp;
   int 	  vol;
   bool	  vis = FALSE;

   /* Mark faces visible from p. */
   f = faces;
   do {
      vol = VolumeSign( f, p );
      if ( vol < 0 ) {
	     f->visible = VISIBLE;  
	     vis = TRUE;                      
      }
      f = f->next;
   } while ( f != faces );

   /* If no faces are visible from p, then p is inside the hull. */
   if ( !vis ) {
      p->onhull = !ONHULL;  
      return FALSE; 
   }

   /* Mark edges in interior of visible region for deletion.
      Erect a newface based on each border edge. */
   e = edges;
   do {
      temp = e->next;
      if ( e->adjface[0]->visible && e->adjface[1]->visible )
		 e->delete = REMOVED;
      else if ( e->adjface[0]->visible || e->adjface[1]->visible ) 
		 e->newface = MakeConeFace( e, p );
      e = temp;
   } while ( e != edges );
   return TRUE;
}

/*---------------------------------------------------------------------
CleanEdges runs through the edge list and cleans up the structure.
If there is a newface then it will put that face in place of the 
visible face and NULL out newface. It also deletes so marked edges.
---------------------------------------------------------------------*/
void	CleanEdges( void )
{
   tEdge  e;	/* Primary index into edge list. */
   tEdge  t;	/* Temporary edge pointer. */
		
   /* Integrate the newface's into the data structure. */
   /* Check every edge. */
   e = edges;
   do {
      if ( e->newface ) { 
	 if ( e->adjface[0]->visible )
	    e->adjface[0] = e->newface; 
	 else	e->adjface[1] = e->newface;
	 e->newface = NULL;
      }
      e = e->next;
   } while ( e != edges );

   /* Delete any edges marked for deletion. */
   while ( edges && edges->delete ) { 
      e = edges;
      DELETE( edges, e );
   }
   e = edges->next;
   do {
      if ( e->delete ) {
	 t = e;
	 e = e->next;
	 DELETE( edges, t );
      }
      else e = e->next;
   } while ( e != edges );
}

/*---------------------------------------------------------------------
CleanFaces runs through the face list and deletes any face marked visible.
---------------------------------------------------------------------*/
void	CleanFaces( void )
{
   tFace  f;	/* Primary pointer into face list. */
   tFace  t;	/* Temporary pointer, for deleting. */
	

   while ( faces && faces->visible ) { 
      f = faces;
      DELETE( faces, f );
   }
   f = faces->next;
   do {
      if ( f->visible ) {
	 t = f;
	 f = f->next;
	 DELETE( faces, t );
      }
      else f = f->next;
   } while ( f != faces );
}

/*---------------------------------------------------------------------
CleanVertices runs through the vertex list and deletes the 
vertices that are marked as processed but are not incident to any 
undeleted edges. 
The pointer to vnext, pvnext, is used to alter vnext in
ConstructHull() if we are about to delete vnext.
---------------------------------------------------------------------*/
void	CleanVertices( tVertex *pvnext )
{
   tEdge    e;
   tVertex  v, t;
	
   /* Mark all vertices incident to some undeleted edge as on the hull. */
   e = edges;
   do {
      e->endpts[0]->onhull = e->endpts[1]->onhull = ONHULL;
      e = e->next;
   } while (e != edges);
	
   /* Delete all vertices that have been processed but
      are not on the hull. */
   while ( vertices && vertices->mark && !vertices->onhull ) { 
      /* If about to delete vnext, advance it first. */
      v = vertices;
      if( v == *pvnext )
         *pvnext = v->next;
      DELETE( vertices, v );
   }
   v = vertices->next;
   do {
      if ( v->mark && !v->onhull ) {    
	 t = v; 
	 v = v->next;
	 if( t == *pvnext )
         *pvnext = t->next;
	 DELETE( vertices, t );
      }
      else v = v->next;
   } while ( v != vertices );
	
   /* Reset flags. */
   v = vertices;
   do {
      v->duplicate = NULL; 
      v->onhull = !ONHULL; 
      v = v->next;
   } while ( v != vertices );
}

/*---------------------------------------------------------------------
CleanUp goes through each data structure list and clears all
flags and NULLs out some pointers.  The order of processing
(edges, faces, vertices) is important.
---------------------------------------------------------------------*/
void	CleanUp( tVertex *pvnext )
{
   CleanEdges();
   CleanFaces();
   CleanVertices( pvnext );
}

/*---------------------------------------------------------------------
ConstructHull adds the vertices to the hull one at a time.  The hull
vertices are those in the list marked as onhull.
---------------------------------------------------------------------*/
void	ConstructHull( void )
{
   tVertex  v, vnext;
   int 	    vol;
   bool	    changed;	/* T if addition changes hull; not used. */

   v = vertices;
   do {
      vnext = v->next;
      if ( !v->mark ) {
         v->mark = PROCESSED;
	     changed = AddOne( v );
	     CleanUp( &vnext ); /* Pass down vnext in case it gets deleted. */
      }
      v = vnext;
   } while ( v != vertices );
}

void DeleteAll( void ) 
{
   tVertex  v, vnext;
   tEdge    e, enext;
   tFace    f, fnext;

   v = vertices->next;
   do {                    
	 vnext = v->next;
     DELETE( vertices, v );
	 v = vnext;
   } while ( v != vertices );
   FREE(vertices);
   
   e = edges->next;
   do {                    
	 enext = e->next;
     DELETE( edges, e );
	 e = enext;
   } while ( e != edges );
   FREE(edges);
   
   f = faces->next;
   do {                    
	 fnext = f->next;
     DELETE( faces, f );
	 f = fnext;
   } while ( f != faces );
   FREE(faces);
   
   tVertex vertices  = NULL;
   tEdge edges    	 = NULL;
   tFace faces    	 = NULL;
}

void chromosome(struct vector *p, int length)
{
   tVertex  v; 
   int i;
   for ( i = 0 ; i < length ; i++ ) {
	 	 v = MakeNullVertex();
		 v->v[X] = p[i].x;
		 v->v[Y] = p[i].y;
		 v->v[Z] = p[i].z;
		 v->mark = 0;
		 v->vnum = i;
   }
}

int HullPoints(struct vector *p, int length)
{
   tVertex  v;
   v = vertices;
   struct vector q;
   do {
      q.x = v->v[X];
	  q.y = v->v[Y];
	  q.z = v->v[Z];
	  p[length++] = q;
      v = v->next;
   } while ( v != vertices );
   return length;
}

int InHull(struct vector p)
{
   int onHull = 1;
   tVertex  v;
   v = vertices;
   do {
	  if ( (p.x == v->v[X]) & (p.y == v->v[Y]) & (p.z == v->v[Z]) )
		  onHull = 0;
      v = v->next;
   } while ( v != vertices );
   return onHull;
}

int readFile ( struct vector *p, int length, char* argv ) 
{
	FILE *fp;
	fp = fopen(argv, "r");
	struct vector q;
	int line = 1;
    while ( line > 0 ) {
        line = fscanf(fp, "%d %d %d", &q.x, &q.y, &q.z);
		q.inhull = 0;
		q.invader = 0;
		if ( line > 0 )
			p[length++] = q;
		}
	fclose(fp);
	return length;
}

// gcc -O3 terrIndex.c -o calcTerritories -lm
int main(int argc, char** argv)
{
	int i,j;
	int p_length = 0, q_length = 0;
	struct vector p[1000], q[1000];
	int hp_length = 0, hq_length = 0;
	struct vector hp[1000], hq[1000];

	if (argc < 3) {
        printf("Usage: ./calcTerritories <pointList> <pointList>\n");
        exit(0);
    }
	
	p_length = readFile(p,p_length,argv[1]); // read the first list of points 
	chromosome(p,p_length);                  // load the first chromosome
	DoubleTriangle();
	ConstructHull();
	hp_length = HullPoints(hp,hp_length);    // save the convex hull for curve p
	DeleteAll();
	
	// loop over points in p and check if they are on hull surface
    for ( i = 0 ; i < p_length ; i++ )
		for ( j = 0 ; j < hp_length ; j++ )
			if ( (p[i].x == hp[j].x) & (p[i].y == hp[j].y) & (p[i].z == hp[j].z) )
				p[i].inhull = 1;
	
	q_length = readFile(q,q_length,argv[2]); // read the secon list of points
	chromosome(q,q_length);                  // load the first chromosome
	DoubleTriangle();
	ConstructHull();
	hq_length = HullPoints(hq,hq_length);    // save the convex hull for curve q
	DeleteAll();
	
	// loop over points in p and check if they are on hull surface
    for ( i = 0 ; i < q_length ; i++ )
		for ( j = 0 ; j < hq_length ; j++ )
			if ( (q[i].x == hq[j].x) & (q[i].y == hq[j].y) & (q[i].z == hq[j].z) )
				q[i].inhull = 1;
		
	// loop over points in p and check if they are in convex hull of q
    for ( i = 0 ; i < p_length ; i++ ) {
		hq[hq_length].x = p[i].x;
		hq[hq_length].y = p[i].y;
		hq[hq_length].z = p[i].z;
		chromosome(hq,(hq_length+1));
		DoubleTriangle();
		ConstructHull();		
		p[i].invader = InHull(hq[hq_length]);
	    DeleteAll();
	}

	// loop over points in q and check if they are in convex hull of p
    for ( i = 0 ; i < q_length ; i++ ) {
		hp[hp_length].x = q[i].x;
		hp[hp_length].y = q[i].y;
		hp[hp_length].z = q[i].z;
		chromosome(hp,(hp_length+1));
		DoubleTriangle();
		ConstructHull();
		q[i].invader = InHull(hp[hp_length]);
		DeleteAll();
	}
				
	// printf("x y z n othType selfType\n");
				
	// print results
	for ( i = 0 ; i < p_length ; i++ )
		printf("%d %d %d 0 %s %s\n",p[i].x, p[i].y, p[i].z, p[i].invader ? "1" : "0", p[i].inhull ? "0" : "1");
		
	for ( i = 0 ; i < q_length ; i++ )
		printf("%d %d %d 1 %s %s\n",q[i].x, q[i].y, q[i].z, q[i].invader ? "1" : "0", q[i].inhull ? "0" : "1");

	return 0;
}

















