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

	static char* options = (char*)"delaunay QJ Pp";
	coordT *pt = NULL;
	int curlong, totlong;
	tTetra tetra;
	int vid = 0;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	double radius=0.0;
	double* center = NULL;
	double volume = 0.0;

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
		
		//if the normal vector of tetrahedron points downward(=lower convex hull)
		if (facet->normal[3] < 0)
		{
			//Compute the volume of tetrahedron
			volume = qh_facetarea(facet);

			//if the volume is not zero, 
			//generate faces. (I didn't care about duplications of faces)
			if (volume > 0.5)
			{
				//Compute Radius
				center = qh_facetcenter(facet->vertices);
				radius = qh_pointdist(center, tetra->vertex[0]->v, 3);

				//If radius is bigger than alpha, remove tetra. 
				if (SQR(radius) > alpha)
				{
					DELETE(tetras, tetra);
					continue;
				}

				MakeFace(tetra->vertex[0], tetra->vertex[1], tetra->vertex[2], NULL);
				MakeFace(tetra->vertex[1], tetra->vertex[2], tetra->vertex[3], NULL);
				MakeFace(tetra->vertex[2], tetra->vertex[3], tetra->vertex[0], NULL);
				MakeFace(tetra->vertex[3], tetra->vertex[0], tetra->vertex[1], NULL);
			}
		}
	}

	free(pt);
	free(all_v);
	pt = NULL;
	all_v = NULL;
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
}

