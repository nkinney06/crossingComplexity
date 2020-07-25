#include <stdlib.h>   /* malloc */
#include <stdio.h>    /* printf */
#include <math.h>
#include <omp.h>

//////////////////////////////////////////////////////////////////////////
// 
// README
//
// this is a minimalist implementation of the "crossing complexity" for two
// curves described in the paper "Crossing Complexity of Space-Filling
// Curves reveals Entanglement of S-Phase DNA". The code is not optimized 
// and there is no dynamic memory allocation. It should be relatively easy
// to study and adapt. Questions can be directed to nkinney06@gmail.com
// 
// to compile this code:
// gcc -Wall -Werror -O3 -fopenmp S2code.c -o calcEntanglement -lm
//
// this code has only been tested in a linux environment (16.04)
// if issues arise contact nkinney06@gmail.com. 
//
// Usage: ./calcEntanglement <verticies> <chr 1> <chr 2>
//
// here <verticies> <chr 1> and <chr 2> are text files. formatting for each
// file is as follows:
//
// 0.051 -0.131 0.992
// -0.218 0.107 0.973
// 0.299 0.087 0.954
// -0.190 -0.314 0.935
// -0.080 0.406 0.911
// 0.367 -0.269 0.892
// -0.490 -0.053 0.873
// 0.3471 0.395 0.858
//
// i.e. each file is a list of x, y, and z coordinates (floats) separated
// by spaces WITHOUT a header line. The order does not matter for the test
// directions. The order DOES matter for <chr 1> and <chr 2>; for these we
// assume that the point on line 1 is bonded to the point on line 2; line 2
// is bonded to line 3; etc.
//
// The program translates curve 1 along each test direction while enumerating
// the crossings with respect to curve 2. Note on lines 325 and 326 that the
// test directions are scaled to a magnitude of 42; this is arbitraty; as long
// as the two curves are fully separated you get the same answer. You just
// need to make sure the magnitude of direction vectors exceed the dimensions
// of your two curves.
//
// Also note that the program uses no dynamic memory allocation. If you are
// dealing with extra large structure youll need to increate the array sizes.
//
// the output is tab delimited with the follow fields:
// thread     just a sanity check for parallelization that I never removed
// phi        polar coordinate of the direction vector
// psi        polar coordinate of the direction vector
// x          x coordinate of the direction vector
// y          y coordinate of the direction vector
// z          z coordinate of the direction vector
// cx         "cubized" x coordinate for the direction vector
// cy         "cubized" y coordinate for the direction vector
// cz         "cubized" z coordinate for the direction vector
// crossings  number of times curve 1 crossed curve 2 in this direction
//
//////////////////////////////////////////////////////////////////////////

struct vector {
    float x;                    
	float y;
	float z;
} buffer[1000000];

float det(struct vector a, struct vector b, struct vector c) {
	float d = 0.0;
	d = d + (a.x)*(b.y)*(c.z)+(a.y)*(b.z)*(c.x);
	d = d + (a.z)*(b.x)*(c.y)-(a.z)*(b.y)*(c.x);
	d = d - (a.y)*(b.x)*(c.z)-(a.x)*(b.z)*(c.y);
	return d;
}

float detX(struct vector a, struct vector b, struct vector c) {
	float d = 0.0;
	d = d + (1.)*(b.y)*(c.z)+(a.y)*(b.z)*(1.);
	d = d + (a.z)*(1.)*(c.y)-(a.z)*(b.y)*(1.);
	d = d - (a.y)*(1.)*(c.z)-(1.)*(b.z)*(c.y);
	return d;
}

float detY(struct vector a, struct vector b, struct vector c) {
	float d = 0.0;
	d = d + (a.x)*(1.)*(c.z)+(1.)*(b.z)*(c.x);
	d = d + (a.z)*(b.x)*(1.)-(a.z)*(1.)*(c.x);
	d = d - (1.)*(b.x)*(c.z)-(a.x)*(b.z)*(1.);
	return d;
}

float detZ(struct vector a, struct vector b, struct vector c) {
	float d = 0.0;
	d = d + (a.x)*(b.y)*(1.)+(a.y)*(1.)*(c.x);
	d = d + (1.)*(b.x)*(c.y)-(1.)*(b.y)*(c.x);
	d = d - (a.y)*(b.x)*(1.)-(a.x)*(1.)*(c.y);
	return d;
}

int readFile ( struct vector *p, char* argv ) {
	FILE *fp;
	fp = fopen(argv, "r");
	struct vector q;
	int line = 1;
	int lineNumber = 0;
    	while ( line > 0 ) {
        line = fscanf(fp, "%f %f %f", &q.x, &q.y, &q.z);
		if ( line > 0 )
			buffer[lineNumber++] = q;
		}
	for ( line = 0 ; line < lineNumber ; line++ )
		p[line] = buffer[line];
	fclose(fp);
	return lineNumber;
}

struct vector add_vectors(struct vector u, struct vector v) {
	struct vector sum;
	sum.x = u.x + v.x;
	sum.y = u.y + v.y;
	sum.z = u.z + v.z;
	return sum;
}

struct vector multiply_scalar(struct vector u, int n) {
	struct vector product;
	product.x = u.x * n;
	product.y = u.y * n;
	product.z = u.z * n;
	return product;
}

struct vector subtract(struct vector a, struct vector b) {
	struct vector difference;
	difference.x = a.x - b.x;
	difference.y = a.y - b.y;
	difference.z = a.z - b.z;
	return difference;
}

struct vector normalize(struct vector s) {
	struct vector norm;
	float size = sqrt( s.x*s.x + s.y*s.y + s.z*s.z );
	norm.x = s.x / size;
	norm.y = s.y / size;
	norm.z = s.z / size;
	return norm;
}

float dot_product(struct vector a, struct vector b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}


//////////////////////////////////////////////////////////////////////////
// I did not write this function and take no credit.
// see the answers to the following stackoverflow question:
// https://stackoverflow.com/questions/2656899/mapping-a-sphere-to-a-cube
//////////////////////////////////////////////////////////////////////////
struct vector cubizePoint2(struct vector sphere)
{
	struct vector position;
	position = sphere;
	
    float x,y,z;
    x = position.x;
    y = position.y;
    z = position.z;

    float fx, fy, fz;
    fx = fabsf(x);
    fy = fabsf(y);
    fz = fabsf(z);

    const double inverseSqrt2 = 0.70710676908493042;

    if (fy >= fx && fy >= fz) {
        float a2 = x * x * 2.0;
        float b2 = z * z * 2.0;
        float inner = -a2 + b2 -3;
        float innersqrt = -sqrtf((inner * inner) - 12.0 * a2);

        if(x == 0.0 || x == -0.0) { 
            position.x = 0.0; 
        }
        else {
            position.x = sqrtf(innersqrt + a2 - b2 + 3.0) * inverseSqrt2;
        }

        if(z == 0.0 || z == -0.0) {
            position.z = 0.0;
        }
        else {
            position.z = sqrtf(innersqrt - a2 + b2 + 3.0) * inverseSqrt2;
        }

        if(position.x > 1.0) position.x = 1.0;
        if(position.z > 1.0) position.z = 1.0;

        if(x < 0) position.x = -position.x;
        if(z < 0) position.z = -position.z;

        if (y > 0) {
            // top face
            position.y = 1.0;
        }
        else {
            // bottom face
            position.y = -1.0;
        }
    }
    else if (fx >= fy && fx >= fz) {
        float a2 = y * y * 2.0;
        float b2 = z * z * 2.0;
        float inner = -a2 + b2 -3;
        float innersqrt = -sqrtf((inner * inner) - 12.0 * a2);

        if(y == 0.0 || y == -0.0) { 
            position.y = 0.0; 
        }
        else {
            position.y = sqrtf(innersqrt + a2 - b2 + 3.0) * inverseSqrt2;
        }

        if(z == 0.0 || z == -0.0) {
            position.z = 0.0;
        }
        else {
            position.z = sqrtf(innersqrt - a2 + b2 + 3.0) * inverseSqrt2;
        }

        if(position.y > 1.0) position.y = 1.0;
        if(position.z > 1.0) position.z = 1.0;

        if(y < 0) position.y = -position.y;
        if(z < 0) position.z = -position.z;

        if (x > 0) {
            // right face
            position.x = 1.0;
        }
        else {
            // left face
            position.x = -1.0;
        }
    }
    else {
        float a2 = x * x * 2.0;
        float b2 = y * y * 2.0;
        float inner = -a2 + b2 -3;
        float innersqrt = -sqrtf((inner * inner) - 12.0 * a2);

        if(x == 0.0 || x == -0.0) { 
            position.x = 0.0; 
        }
        else {
            position.x = sqrtf(innersqrt + a2 - b2 + 3.0) * inverseSqrt2;
        }

        if(y == 0.0 || y == -0.0) {
            position.y = 0.0;
        }
        else {
            position.y = sqrtf(innersqrt - a2 + b2 + 3.0) * inverseSqrt2;
        }

        if(position.x > 1.0) position.x = 1.0;
        if(position.y > 1.0) position.y = 1.0;

        if(x < 0) position.x = -position.x;
        if(y < 0) position.y = -position.y;

        if (z > 0) {
            // front face
            position.z = 1.0;
        }
        else {
            // back face
            position.z = -1.0;
        }
    }
	return position;
}

int main(int argc, char** argv)
{
    if (argc != 4) {
        printf("Usage: ./calcEntanglement <verticies> <chr 1> <chr 2>\n");
        exit(0);
    }

    int a,b,c,i,thread_id;
	struct vector directions[100000];
	struct vector dirOnCube[100000];
	struct vector curveOne[100000];
	struct vector curveTwo[100000];
	
	a = readFile(directions,argv[1]);
	b = readFile(curveOne,argv[2]);
	c = readFile(curveTwo,argv[3]);
	for ( i = 0 ; i < a ; i++ )
		dirOnCube[i] = cubizePoint2(directions[i]);
	
	printf("thread\tphi\tpsi\tx\ty\tz\tcx\tcy\tcz\tcrossings\n");
	#pragma omp parallel private(thread_id)
	{
		
	#pragma omp for
	for ( i = 0 ; i < a ; i++ ) { // loop over the direction vectors
	
		int j, crossings = 0;
		
		for ( j = 0 ; j < b-1 ; j++ ) { // loop over segments in curve 1
		
			struct vector m,n,p,q;
			float W,X,Y,Z;
			int k;
			
			m = curveOne[j];
			n = curveOne[j+1];
			p = add_vectors(curveOne[j],multiply_scalar(directions[i],42));
			q = add_vectors(curveOne[j+1],multiply_scalar(directions[i],42));

			W = -1*det(m,n,p); // a few different types of determinants are used...
			X = detX(m,n,p);   
			Y = detY(m,n,p);
			Z = detZ(m,n,p);

			for ( k = 0 ; k < c-1 ; k++ ) { // loop over segments in curve 2
			
				struct vector u,v;
				float numerator,denominator,divide;
				
				u = curveTwo[k];
				v = curveTwo[k+1];
				numerator = X*u.x+Y*u.y+Z*u.z+W;
				denominator = X*(u.x-v.x)+Y*(u.y-v.y)+Z*(u.z-v.z);
				
				if (denominator!=0.0) {
					struct vector s,e,f,g,h;
					divide = numerator/denominator;
					s.x = u.x+divide*(v.x-u.x); // the intesection of the line and plane is store
					s.y = u.y+divide*(v.y-u.y); // in the point s. we find s using the parameteric
					s.z = u.z+divide*(v.z-u.z); // equation of the line
					
					e = normalize(subtract(m,s)); // the next 9 lines find if the intersection point s
					f = normalize(subtract(n,s)); // is in the boundary of a,b,c,and d
					g = normalize(subtract(p,s));
					h = normalize(subtract(q,s));

					float sum_angles = acos(dot_product(e,g)) + acos(dot_product(g,h)) + acos(dot_product(f,h)) + acos(dot_product(e,f));
					if ((sum_angles > 2*3.1415) && (divide>=0 && divide<=1))
						crossings++; // tally the number of crossings for this direction
				}
			}
		}
		thread_id = omp_get_thread_num();
		if ( directions[i].x == 0 )
			printf("thread->%d\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%d\n",thread_id,1.5707963,acos(directions[i].z),directions[i].x,directions[i].y,directions[i].z,dirOnCube[i].x,dirOnCube[i].y,dirOnCube[i].z,crossings);
		else
			printf("thread->%d\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t%d\n",thread_id,atan2(directions[i].y,directions[i].x),acos(directions[i].z),directions[i].x,directions[i].y,directions[i].z,dirOnCube[i].x,dirOnCube[i].y,dirOnCube[i].z,crossings);
		}
		
	}
    return 0;
}
