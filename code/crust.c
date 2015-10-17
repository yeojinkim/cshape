#include "crust.h"
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

void	Crust( void )
{
	tVertex ptr_v;
	tVertex *all_v = NULL;
	int vsize = 0;
	int id = 0;

	//static char* options = (char*)"delaunay QJ Pp";
	static char* options = (char*)"voronoi v FA QJ Pp";
	coordT *pt = NULL;
	int curlong, totlong;
	facetT *facet = NULL;
	vertexT *vertex = NULL;
	vertexT **vertexp = NULL;
	double volume = 0.0;
	tVertex center = NULL;
	tList voronoiVertex = NULL;
	tList ptr_vorv = NULL;
	tVertex currVertex = NULL;
	double testRadius = 0.0;
/*	
	
	double maxDistance = -1.0;
	double vorDistance = 0.0;
	tList ptr_vorv = NULL;
	tVertex poleVertex = NULL;
	tVertex currVertex = NULL;*/
	
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

	//loop through all faces
	FORALLfacets
	{
		//get the vertex of voronoi diagram
		//loop through each vertex

		//Compute the volume of tetrahedron
		volume = qh_facetarea(facet);
		printf("volume : %lf\n", volume);

		FOREACHvertex_(facet->vertices)
		{	
		
			//If this facet is lower hull and volume is bigger than 0, this tetrahderon consists of delaunay triangulation 
			if (!(facet->upperdelaunay))// && (volume > 0.5))
			{
				printf("point : %lf %lf %lf\n", vertex->point[0], vertex->point[1], vertex->point[2]);
				
				//copy the center of delaunay triangle		
				center = MakeCenterVertex(facet->center);
				printf("add center : %lf %lf %lf\n", center->v[0], center->v[1], center->v[2]);
//				testRadius = qh_pointdist(center->v, vertex->point, 3);
//				printf("center Radius : %lf \n", testRadius);

				//printf("Vertex: %lf %lf %lf\n", vertex->point[0], vertex->point[1], vertex->point[2]);
				//printf("fac center : %lf %lf %lf\n", facet->center[0], facet->center[1], facet->center[2]);
				//printf("add center : %lf %lf %lf\n", center->v[0], center->v[1], center->v[2]);
				//printf("test center: %lf %lf %lf\n", testCenter[0], testCenter[1], testCenter[2]);
				
				//create the vornoi vertex using the center of delanauy triangle
				voronoiVertex = MakeVoronoiVertex(center);

				//Add voronoi vertex to vornoi cite
				ADD(all_v[qh_pointid(vertex->point)]->vvlist, voronoiVertex);
			}	
			
		}
		printf("\n");
	}
		/*
	ptr_v = vertices;
	do {
		printf("vertex : %lf %lf %lf \n", ptr_v->v[0], ptr_v->v[1], ptr_v->v[2]);
		ptr_vorv = ptr_v->vvlist;
		do{
			currVertex = (tVertex)(ptr_vorv->p);
			printf("vv : %lf %lf %lf \n", currVertex->v[0], currVertex->v[1], currVertex->v[2]);
			ptr_vorv = ptr_vorv->next;
		} while (ptr_vorv != ptr_v->vvlist);

		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);
	*/
	/*
	//Compute a pole and an antipole
	ptr_v = vertices;
	do {		
		//Find a pole
		ptr_vorv = ptr_v->vvlist;
		do{
			//Compute the distance between a voronoi cite and a voronoi vertex
			currVertex = (tVertex) (ptr_vorv->p);
			vorDistance = qh_pointdist(currVertex->v, ptr_v->v, 3);
			
			//If current voronoi vertex is farther than other candidates, store it as current and best candidates for pole
			if (vorDistance > maxDistance)
			{
				//Erase the last pole vertex
				if(poleVertex!=NULL) poleVertex->ispole = false;

				//Save the new pole vertex
				poleVertex = currVertex;
				poleVertex->ispole = true;
				maxDistance = vorDistance;
			}

			ptr_vorv = ptr_vorv->next;
		} while (ptr_vorv != ptr_v->vvlist);
		
		//Find an antipole 
		ptr_vorv = ptr_v->vvlist;
		do{
			ptr_vorv = ptr_vorv->next;
		} while (ptr_vorv != ptr_v->vvlist);

		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);

	ptr_v = vertices;
	do {
		printf("vertex : %lf %lf %lf\n", ptr_v->v[0], ptr_v->v[1], ptr_v->v[2]);
		//Find a pole
		ptr_vorv = ptr_v->vvlist;
		do{
			currVertex = (tVertex)ptr_vorv->p;
			if (currVertex->ispole)
				printf("have a pole:%lf %lf %lf\n\n",currVertex->v[0], currVertex->v[1], currVertex->v[2]);
			ptr_vorv = ptr_vorv->next;
		} while (ptr_vorv != ptr_v->vvlist);

		ptr_v = ptr_v->next;
	} while (ptr_v != vertices);
	*/
	free(pt);
	free(all_v);
	all_v = NULL;
	qh_freeqhull(!qh_ALL);
	qh_memfreeshort(&curlong, &totlong);
}

//Find a pole

//Find an antipole


//if the normal vector of tetrahedron points downward(=lower convex hull) and the volume is not zero, 
//generate faces. (I didn't care about duplications of faces)


