#include "crust.h"
#include "delauany.h"
#include "shape.h"

/*-------------------------------------------------------------------*/
//headers from qhull
#include "qhull.h"
#include "poly.h"
#include "qset.h"
#include "geom.h"
#include "io.h"

/*-------------------------------------------------------------------*/
//defined in shape.h/.c
extern tVertex vertices;
extern tEdge edges;
extern tFace faces;
extern tTetra tetras;

/*-------------------------------------------------------------------*/

#define POLE_MAX_THRESHOLD 1000

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

tList MakeBoundedInformation(vertexT* neigborCites)
{
	tList VornoiVertexList;

	NEW(VornoiVertexList, tsList);
	VornoiVertexList->p = (void*)neigborCites;
	VornoiVertexList->next = VornoiVertexList->prev = NULL;

	return VornoiVertexList;
}

bool HasVornoiVertex(tList voronoiVertices, double* vertex)
{
	tList ptr_vorvs = voronoiVertices;
	tVertex vorv;

	if (voronoiVertices == NULL) return false;
	if (ptr_vorvs->next == ptr_vorvs) return false;
	ptr_vorvs = ptr_vorvs->next;	//skip bounded/unbounded information
	do{
		vorv = (tVertex)ptr_vorvs->p;
		if (SQR(vorv->v[0] - vertex[0]) + SQR(vorv->v[1] - vertex[1]) + SQR(vorv->v[2] - vertex[2]) < 0.5)
			return true;
		ptr_vorvs = ptr_vorvs->next;
	} while (ptr_vorvs != voronoiVertices);
	return false;
}

bool HasBoundedInformation(tList voronoiVertices)
{
	tList ptr_vorvs = voronoiVertices;
	if (ptr_vorvs == NULL) return false;
	return true;
}

bool HasVertex(tVertex vertex)
{
	tVertex ptr_v = vertices;

	if (vertices == NULL) return false;
	do{
		if (SQR(ptr_v->v[0] - vertex->v[0]) + SQR(ptr_v->v[1] - vertex->v[1]) + SQR(ptr_v->v[2] - vertex->v[2]) < 0.5)
			return true;
		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);
	return false;
}

void FindPoleAntipole(int vsize)
{
	tVertex cite;
	tList voronoiVertices = NULL;
	tVertex vorv = NULL;
	double vorDistance = 0.0;
	double maxDistance = 0.5;
	double minProduct = 1e+16;
	double minDistance = 0.5;
	tVertex pole = NULL;
	tVertex antipole = NULL;
	tVertex normal = NULL;
	vertexT* vertex = NULL;
	double faceNormal[3] = { 0.0, 0.0, 0.0 };
	double dist = 0.0;
	facetT *neighbor, **neighborp;
	setT* neigborcites = NULL;
	vertexT* ncite;
	int id = 0;
	int i = 0;

	NEW(normal, tsVertex);

	//loop through the valid cites, compute a pole and an antipole
	cite = vertices;
	do {

		voronoiVertices = cite->vvlist;

		//if there is no voronoi verticies, skip the cites. 
		if (voronoiVertices == NULL)
		{
			cite = cite->next;
			continue;
		}

		//check if this cell is bounded/unbounded
		//if vertex==NULL this cell is bounded
		//otherwise, this cell is unbounded
		vertex = (vertexT*)(voronoiVertices->p);
		voronoiVertices = voronoiVertices->next;

		//if voronoi cell is unbounded, find an average normal of adjacent triangles 
		if (vertex != NULL)
		{
			normal->v[0] = 0.0; normal->v[1] = 0.0; normal->v[2] = 0.0;
			FOREACHneighbor_(vertex)
			{
				neigborcites = neighbor->vertices;
				for (i = 0; i < neigborcites->maxsize; i++)
				{
					ncite = (vertexT*)neigborcites->e[i].p;
					//Test if ncite is unbounded or bounded. if it is bounded, skip the neighbor cite
					if (ncite->seen2 == false) continue;


					//compute the outer normal of adjacent triangle with cite
					if (neighbor->toporient){
						normal->v[0] -= neighbor->normal[0];
						normal->v[1] -= neighbor->normal[1];
						normal->v[2] -= neighbor->normal[2];
					}
					else{
						normal->v[0] += neighbor->normal[0];
						normal->v[1] += neighbor->normal[1];
						normal->v[2] += neighbor->normal[2];
					}
				}
			}

			dist = sqrt(normal->v[0] * normal->v[0] + normal->v[1] * normal->v[1] + normal->v[2] * normal->v[2]);

			normal->v[0] /= dist;
			normal->v[1] /= dist;
			normal->v[2] /= dist;
		}
		//if voronoi cell is bounded, find a pole and a normal
		else
		{
			do{
				//compute the distance between a voronoi cite and a voronoi vertex
				vorv = (tVertex)(voronoiVertices->p);
				vorDistance = qh_pointdist(vorv->v, cite->v, 3);

				//If current voronoi vertex is farther than other candidates, store it as current and best candidates for pole
				if (vorDistance > maxDistance)
				{
					pole = vorv;
					maxDistance = vorDistance;

				}

				voronoiVertices = voronoiVertices->next;

			} while (voronoiVertices != cite->vvlist);

			//Test if pole is too far away from the model  
			if (maxDistance > POLE_MAX_THRESHOLD) pole = NULL;
			else{
				//compute normal vector cp (c:voronoi cite, p:pole)
				pole->ispole = true;
				normal->v[0] = pole->v[0] - cite->v[0];
				normal->v[1] = pole->v[1] - cite->v[1];
				normal->v[2] = pole->v[2] - cite->v[2];
				dist = sqrt(normal->v[0] * normal->v[0] + normal->v[1] * normal->v[1] + normal->v[2] * normal->v[2]);
				normal->v[0] /= dist;
				normal->v[1] /= dist;
				normal->v[2] /= dist;
			}
		}

		//find a antipole
		voronoiVertices = cite->vvlist;
		voronoiVertices = voronoiVertices->next;		//skip bounded/unbounded information
		do{

			vorv = (tVertex)(voronoiVertices->p);

			//if a pole exists and this voronoi vertex is a pole, skip the loop
			if (pole != NULL && SQR(vorv->v[0] - pole->v[0]) + SQR(vorv->v[1] - pole->v[1]) + SQR(vorv->v[2] - pole->v[2]) < 0.5)
			{
				voronoiVertices = voronoiVertices->next;
				continue;
			}

			//inner product normal with vector cv (c:voronoi cite, v:candidate antipole)
			vorDistance = (vorv->v[0] - cite->v[0]) * normal->v[0] + (vorv->v[1] - cite->v[1]) * normal->v[1] + (vorv->v[2] - cite->v[2]) * normal->v[2];

			if (vorDistance < minProduct)
			{
				antipole = vorv;
				minProduct = vorDistance;
				minDistance = qh_pointdist(vorv->v, cite->v, 3);
			}

			voronoiVertices = voronoiVertices->next;

		} while (voronoiVertices != cite->vvlist);

		if (antipole != NULL) {
			if (minDistance > POLE_MAX_THRESHOLD) antipole = NULL;
			else  antipole->ispole = true;
		}

		//add a pole and an antipole to point set
		if (pole != NULL && !HasVertex(pole))
		{
			pole->vnum = vsize + id; id++;
			ADD(vertices, pole);
		}
		if (antipole != NULL && !HasVertex(antipole))
		{
			antipole->vnum = vsize + id; id++;
			ADD(vertices, antipole);

		}

		pole = NULL;
		antipole = NULL;
		maxDistance = 0.5;
		minProduct = 1e+16;
		minDistance = 0.5;
		normal->v[0] = 0.0;
		normal->v[1] = 0.0;
		normal->v[2] = 0.0;
		cite = cite->next;
	} while (cite != vertices);

	free(normal);
}

void	Crust(void)
{
	tVertex ptr_v;
	tVertex *all_v = NULL;
	int vsize = 0;
	int id = 0;

	static char* options = (char*)"voronoi v Qbb Qz QJ Pp";
	coordT *pt = NULL;
	int curlong, totlong;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	facetT *neighbor, **neighborp;

	double volume = 0.0;
	tVertex center = NULL;
	tList voronoiVertex = NULL;
	tList boundedInformation = NULL;
	tVertex voronoiCite = NULL;
	vertexT* neighborCites = NULL;

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
		//Compute the volume of tetrahedron
		volume = qh_facetarea(facet);

		FOREACHvertex_(facet->vertices)
		{
			//if this facet is lower hull and volume is bigger than 0, this tetrahderon consists of delaunay triangulation 
			if (!(facet->upperdelaunay) && (volume > 0.5))
			{
				//obtain current voronoi cite
				voronoiCite = all_v[qh_pointid(vertex->point)];

				if (!HasBoundedInformation(voronoiCite->vvlist))
				{
					//check if this cell is bounded or unbounded and add unbounded/bounded information to voronoi cell
					FOREACHneighbor_(vertex)
					{
						if (neighbor->upperdelaunay)	//this is unbounded cell
						{
							//mark this is unbounded cell
							vertex->seen2 = true;

							//keep neighbor cites
							neighborCites = vertex;
						}
						else neighborCites = NULL;
					}

					boundedInformation = MakeBoundedInformation(neighborCites);
					ADD(voronoiCite->vvlist, boundedInformation);

				}
				//for degenerate case (circumesphere has 4 points), create the voronoi vertex once. 
				if (HasVornoiVertex(voronoiCite->vvlist, facet->center))
					continue;

				//copy the center of delaunay triangle and add voronoi vertex to vornoi cite		
				center = MakeCenterVertex(facet->center);
				voronoiVertex = MakeVoronoiVertex(center);
				ADD(voronoiCite->vvlist, voronoiVertex);

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
	pt = NULL;
	all_v = NULL;
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);

}

