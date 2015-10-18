#include "crust.h"
#include "delauany.h"
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
tVertex MakeCenterVertex(double* center)
{
	tVertex centerVertex;
	//copy the center of delaunay triangle
	NEW(centerVertex, tsVertex);
	centerVertex->v[0] = center[0];
	centerVertex->v[1] = center[1];
	centerVertex->v[2] = center[2];
	centerVertex->next = centerVertex->prev = NULL;
	centerVertex->ispole = false;
	centerVertex->vvlist = NULL;

	return centerVertex;
}

tList MakeVoronoiVertex(tVertex vertex)
{
	tList VornoiVertexList;
	
	NEW(VornoiVertexList, tsList);
	VornoiVertexList->p = (void*)vertex;
	VornoiVertexList->next = VornoiVertexList->prev = NULL;

	return VornoiVertexList;
}

bool HasVornoiVertex(tList voronoiVertices, double* vertex)
{
	tList ptr_vorvs = voronoiVertices;
	tVertex vorv;

	if (voronoiVertices == NULL) return false;

	do{
		vorv = (tVertex)ptr_vorvs->p;
		if (SQR(vorv->v[0] - vertex[0]) + SQR(vorv->v[1] - vertex[1]) + SQR(vorv->v[2] - vertex[2])<0.5)
			return true;
		ptr_vorvs = ptr_vorvs->next;
	} while (ptr_vorvs != voronoiVertices);
	return false;
}

bool HasVertex(tVertex vertex)
{
	tVertex ptr_v = vertices;

	if (vertices == NULL) return false;
	do{
		if (SQR(ptr_v->v[0] - vertex->v[0]) + SQR(ptr_v->v[1] - vertex->v[1]) + SQR(ptr_v->v[2] - vertex->v[2])<0.5)
			return true;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);
	return false;
}

bool isBoundedVornoicell(tList voronoiVertices)
{
	int num = 0;
	tList ptr_vorvs = voronoiVertices;

	do{
		num++;
		ptr_vorvs = ptr_vorvs->next;
	} while (ptr_vorvs != voronoiVertices);

	if (num < 4) return false;
	return true;
}
void FindPoleAntipole(int vsize)
{
	tVertex ptr_v;
	tList ptr_vorv = NULL;
	tVertex vorv = NULL;
	double vorDistance = 0.0;
	double maxDistance = 0.0;
	double minDistance = 1e+16;
	tVertex pole = NULL;
	tVertex antipole = NULL;
	tVertex normal = NULL;

	NEW(normal, tsVertex);

	int i = 0;

	//compute a pole and an antipole
	ptr_v = vertices;
	do {

		ptr_vorv = ptr_v->vvlist;
		//if there is no voronoi verticies, it means this point is added pole/antipole.
		if (ptr_vorv == NULL) break;

		//count the number of vornoi vertices. if voronoi verticies < 4, it is unbounded cell
		if (isBoundedVornoicell(ptr_vorv))
		{
			//if voronoi cell is bounded, find a pole and a normal
			do{
				//Compute the distance between a voronoi cite and a voronoi vertex
				vorv = (tVertex)(ptr_vorv->p);
				vorDistance = qh_pointdist(vorv->v, ptr_v->v, 3);

				//If current voronoi vertex is farther than other candidates, store it as current and best candidates for pole
				if (vorDistance > maxDistance)
				{
					pole = vorv;
					maxDistance = vorDistance;
				}

				ptr_vorv = ptr_vorv->next;

			} while (ptr_vorv != ptr_v->vvlist);
			if (pole != NULL)
			{
				pole->ispole = true;
				//compute normal  vector sp (s:voronoi cite, p:pole)
				normal->v[0] = pole->v[0] - ptr_v->v[0];
				normal->v[1] = pole->v[1] - ptr_v->v[1];
				normal->v[2] = pole->v[2] - ptr_v->v[2];
			}
		}
		else
		{
			//if voronoi cell is unbounded, find a normal 
			

			//----- adjacent triangles + if there is any variable for bounded/unbounded?
		}
		
		//find a antipole
		ptr_vorv = ptr_v->vvlist;
		do{
			vorv = (tVertex)(ptr_vorv->p);
			
			if (!(vorv->ispole))
			{
				//inner product between 
				vorDistance = (vorv->v[0] - ptr_v->v[0]) * (pole->v[0] - ptr_v->v[0]) + (vorv->v[1] - ptr_v->v[1]) * (pole->v[1] - ptr_v->v[1]) + (vorv->v[2] - ptr_v->v[2]) * (pole->v[2] - ptr_v->v[2]);
				if (vorDistance < minDistance)
				{
					antipole = vorv;
					minDistance = vorDistance;
				}
			}
			ptr_vorv = ptr_vorv->next;

		} while (ptr_vorv != ptr_v->vvlist);
		
		if (antipole != NULL) antipole->ispole = true;
		
		//add a pole and an antipole to point set
		if (pole != NULL && !HasVertex(pole))
		{
			pole->vnum = vsize + i; i++;
			ADD(vertices, pole);
		}
		if (antipole != NULL && !HasVertex(antipole))
		{
			antipole->vnum = vsize + i; i++;
			ADD(vertices, antipole);
			
		}
	
		pole = NULL;
		antipole = NULL;
		maxDistance = 0.0;
		minDistance = 1e+16;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);
	free(normal);
}

void	Crust( void )
{
	tVertex ptr_v;
	tVertex *all_v = NULL;
	int vsize = 0;
	int id = 0;

	//static char* options = (char*)"delaunay QJ Pp";
	static char* options = (char*)"voronoi v n FA QJ Pp";
	coordT *pt = NULL;
	int curlong, totlong;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	double volume = 0.0;

	tVertex center = NULL;
	tList voronoiVertex = NULL;

	tTetra testt = NULL;
	tEdge teste = NULL;
	tFace testf = NULL;
	
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

	//Copy points 
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
	// compute convex hull in 3D by calling qhull
	//

	qh_init_A(stdin, stdout, stderr, 0, NULL);
	qh_initflags(options);
	qh_init_B(pt, vsize, 4, false);
	qh_qhull();
	qh_check_output();

	qh_setvoronoi_all();

	//loop through all faces, compute vornoi vertices.
	FORALLfacets
	{
		//get the vertex of voronoi diagram
		//loop through each vertex

		//Compute the volume of tetrahedron
		volume = qh_facetarea(facet);
		
		FOREACHvertex_(facet->vertices)
		{	
		
			//if this facet is lower hull and volume is bigger than 0, this tetrahderon consists of delaunay triangulation 
			if (!(facet->upperdelaunay) && (volume > 0.5))
			{
				//for degenerate case (circumesphere has 4 points), create the voronoi vertex once. 
				if (HasVornoiVertex(all_v[qh_pointid(vertex->point)]->vvlist, facet->center)) 
					continue;
				//copy the center of delaunay triangle		
				center = MakeCenterVertex(facet->center);
				
				//create the vornoi vertex using the center of delanauy triangle
				voronoiVertex = MakeVoronoiVertex(center);

				//Add voronoi vertex to vornoi cite
				ADD(all_v[qh_pointid(vertex->point)]->vvlist, voronoiVertex);
			}		
		}
	}
	
	//compute a pole and an antipole among vornoi verticies
	FindPoleAntipole(vsize);
	
	//compute delaunay triangulation of new points set
	//report faces of tetrahedron which ahve end points from only original point set
	Delaunay();
	
	free(pt);
	free(all_v);
	all_v = NULL;
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
	
}

