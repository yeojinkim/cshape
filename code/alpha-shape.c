#include "alpha-shape.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
//headers from qhull
#include "qhull.h"
#include "poly.h"
#include "qset.h"
#include "geom.h"


/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

/*-------------------------------------------------------------------*/
void AlphaShape( unsigned int alpha )
{
	tVertex ptr_v;
	tVertex *all_v = NULL;
	int vsize = 0;
	int id = 0;

	static char* options = (char*)"delaunay C-4 QJ Pp";
	coordT *pt = NULL;
	int curlong, totlong;
	tTetra tetra;
	int vid = 0;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	tFace faceTetra = NULL;
	tVertex vertexTetra = NULL;
	int signVolTetra = -1;
	coordT* center = NULL;
	double a=.0, b=.0, c=.0, d=.0;
	double radius=0.0;

	NEW(faceTetra, tsFace);

	//
	// create points in 4D (x,y,z,x^2+y^2+z^2)
	//

	//Count the number of points
	ptr_v = vertices;
	do {
		vsize++;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//Allocate memory
	pt = (coordT*)calloc(vsize * 4, sizeof(coordT));
	all_v = (tVertex*)calloc(vsize, sizeof(tVertex));
	assert(pt && all_v);

	//Copy points and compute 4th element
	ptr_v = vertices;
	do {
		pt[id++] = ptr_v->v[0];
		pt[id++] = ptr_v->v[1];
		pt[id++] = ptr_v->v[2];
		pt[id++] = ptr_v->v[0] * ptr_v->v[0] + ptr_v->v[1] * ptr_v->v[1] + ptr_v->v[2] * ptr_v->v[2];
		all_v[ptr_v->vnum] = ptr_v;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	//
	// compute convex hull in 4D by calling qhull
	//

	qh_init_A(stdin, stdout, stderr, 0, NULL);
	qh_initflags(options);
	qh_init_B(pt, vsize, 4, false);
	qh_qhull();
	qh_check_output();

	//loop through all faces
	FORALLfacets
	{
		tetra = MakeNullTetra();

		//get vertices of facet
		//loop through each vertex
		vid = 0;
		FOREACHvertex_(facet->vertices)
		{
			//get the id of the vertex (4 points)
			tetra->vertex[vid++] = all_v[qh_pointid(vertex->point)];
		}
		
		//Compute Radius
		a = facet->center[0] - tetra->vertex[0]->v[0];
		b = facet->center[1] - tetra->vertex[0]->v[1];
		c = facet->center[2] - tetra->vertex[0]->v[2];
		d = facet->center[3] - tetra->vertex[0]->v[3];
		radius = sqrt(a * a + b * b + c * c + d * d);
		
		if (radius > alpha)
		{
			DELETE(tetras, tetra);
			continue;
		}

		//Compute the sign volume of tetrahedron
		faceTetra->vertex[0] = tetra->vertex[0];
		faceTetra->vertex[1] = tetra->vertex[1];
		faceTetra->vertex[2] = tetra->vertex[2];
		vertexTetra = tetra->vertex[3];
		signVolTetra = VolumeSign(faceTetra, vertexTetra);

		//if the normal vector of tetrahedron points downward(=lower convex hull) and the volume is not zero, 
		//generate faces. (I didn't care about duplications of faces)
		if (facet->normal[3] < 0 && signVolTetra != 0)
		{
			MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[2], NULL);
			MakeFace(tetra->vertex[1], tetra->vertex[2], tetra->vertex[3], NULL);
			MakeFace(tetra->vertex[2], tetra->vertex[3], tetra->vertex[0], NULL);
			MakeFace(tetra->vertex[3], tetra->vertex[0], tetra->vertex[1], NULL);
		}
	}

	free(pt);
	free(all_v);
	free(faceTetra);
	pt = NULL;
	all_v = NULL;
	faceTetra = NULL;
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
}

