/*
===========================================================================

Doom 3 GPL Source Code
Copyright (C) 1999-2011 id Software LLC, a ZeniMax Media company. 

This file is part of the Doom 3 GPL Source Code (?Doom 3 Source Code?).  

Doom 3 Source Code is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Doom 3 Source Code is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Doom 3 Source Code.  If not, see <http://www.gnu.org/licenses/>.

In addition, the Doom 3 Source Code is also subject to certain additional terms. You should have received a copy of these additional terms immediately following the terms and conditions of the GNU General Public License which accompanied the Doom 3 Source Code.  If not, please request a copy in writing from id Software at the address below.

If you have questions concerning this license or the applicable additional terms, you may contact in writing id Software LLC, c/o ZeniMax Media Inc., Suite 120, Rockville, Maryland 20850 USA.

===========================================================================
*/

#ifndef __CM_LOCAL_H__
#define __CM_LOCAL_H__

#include "../idlib/q_shared.h"
#include "../idlib/surfaceflags.h"
#include "../idlib/idlib_public.h"
#include "../idlib/q_containers.h"
#include "../idlib/mapfile.h"
#include "../idlib/bsp/bsp.h"
#include "cm_public.h"

#ifndef ON_EPSILON
#define ON_EPSILON 0.1f
#endif

// NOTE: needs to be same as in bg_public.h
#define	MAX_SUBMODELS		1024	// max bsp models, q3map2 limits to 1024 via MAX_MAP_MODELS
#define MAX_MATERIALS       1024

typedef struct cm_material_s {
    char    name[MAX_QPATH];
    int     surfaceFlags;
    int     contentFlags;
} cm_material_t;

/*
===============================================================================

	Trace model vs. polygonal model collision detection.

===============================================================================
*/

#define MIN_NODE_SIZE						64.0f
#define MAX_NODE_POLYGONS					128
#define CM_MAX_POLYGON_EDGES				64
#define CIRCLE_APPROXIMATION_LENGTH			64.0f

//#define	MAX_SUBMODELS						2048
#define	TRACE_MODEL_HANDLE					MAX_SUBMODELS

#define VERTEX_HASH_BOXSIZE					(1<<6)	// must be power of 2
#define VERTEX_HASH_SIZE					(VERTEX_HASH_BOXSIZE*VERTEX_HASH_BOXSIZE)
#define EDGE_HASH_SIZE						(1<<14)

#define NODE_BLOCK_SIZE_SMALL				8
#define NODE_BLOCK_SIZE_LARGE				256
#define REFERENCE_BLOCK_SIZE_SMALL			8
#define REFERENCE_BLOCK_SIZE_LARGE			256

#define MAX_WINDING_LIST					128		// quite a few are generated at times
#define INTEGRAL_EPSILON					0.01f
#define VERTEX_EPSILON						0.1f
#define CHOP_EPSILON						0.1f


typedef struct cm_windingList_s {
	int					numWindings;			// number of windings
	fixedWinding_t		w[MAX_WINDING_LIST];	// windings
	vec3_t				normal;					// normal for all windings
	vec3_t			    bounds[2];					// bounds of all windings in list
	vec3_t				origin;					// origin for radius
	float				radius;					// radius relative to origin for all windings
	int					contents;				// winding surface contents
	int					primitiveNum;			// number of primitive the windings came from
} cm_windingList_t;

/*
===============================================================================

Collision model

===============================================================================
*/

typedef struct cm_vertex_s {
	vec3_t					p;					// vertex point
	int						checkcount;			// for multi-check avoidance
	unsigned long			side;				// each bit tells at which side this vertex passes one of the trace model edges
	unsigned long			sideSet;			// each bit tells if sidedness for the trace model edge has been calculated yet
} cm_vertex_t;

typedef struct cm_edge_s {
	int						checkcount;			// for multi-check avoidance
	unsigned short			internal;			// a trace model can never collide with internal edges
	unsigned short			numUsers;			// number of polygons using this edge
	unsigned long			side;				// each bit tells at which side of this edge one of the trace model vertices passes
	unsigned long			sideSet;			// each bit tells if sidedness for the trace model vertex has been calculated yet
	int						vertexNum[2];		// start and end point of edge
	vec3_t					normal;				// edge normal
} cm_edge_t;

typedef struct cm_polygonBlock_s {
	int						bytesRemaining;
	byte *					next;
} cm_polygonBlock_t;

typedef struct cm_polygon_s {
	vec3_t				    bounds[2];				// polygon bounds
	int						checkcount;			// for multi-check avoidance
	int						contents;			// contents behind polygon
	qhandle_t		 		material;			// material
	plane_t 				plane;				// polygon plane
	int						numEdges;			// number of edges
	int						edges[1];			// variable sized, indexes into cm_edge_t list
} cm_polygon_t;

typedef struct cm_polygonRef_s {
	cm_polygon_t *			p;					// pointer to polygon
	struct cm_polygonRef_s *next;				// next polygon in chain
} cm_polygonRef_t;

typedef struct cm_polygonRefBlock_s {
	cm_polygonRef_t *		nextRef;			// next polygon reference in block
	struct cm_polygonRefBlock_s *next;			// next block with polygon references
} cm_polygonRefBlock_t;

typedef struct cm_brushBlock_s {
	int						bytesRemaining;
	byte *					next;
} cm_brushBlock_t;

typedef struct cm_brush_s {
	int						checkcount;			// for multi-check avoidance
	vec3_t				    bounds[2];				// brush bounds
	int						contents;			// contents of brush
	qhandle_t		 		material;			// material
	int						primitiveNum;		// number of brush primitive
	int						numPlanes;			// number of bounding planes
	plane_t	    			planes[1];			// variable sized
} cm_brush_t;

typedef struct cm_brushRef_s {
	cm_brush_t *			b;					// pointer to brush
	struct cm_brushRef_s *	next;				// next brush in chain
} cm_brushRef_t;

typedef struct cm_brushRefBlock_s {
	cm_brushRef_t *			nextRef;			// next brush reference in block
	struct cm_brushRefBlock_s *next;			// next block with brush references
} cm_brushRefBlock_t;

typedef struct cm_node_s {
	int						planeType;			// node axial plane type
	float					planeDist;			// node plane distance
	cm_polygonRef_t *		polygons;			// polygons in node
	cm_brushRef_t *			brushes;			// brushes in node
	struct cm_node_s *		parent;				// parent of this node
	struct cm_node_s *		children[2];		// node children
} cm_node_t;

typedef struct cm_nodeBlock_s {
	cm_node_t *				nextNode;			// next node in block
	struct cm_nodeBlock_s *next;				// next block with nodes
} cm_nodeBlock_t;

typedef struct cm_model_s {
	char					name[MAX_QPATH];	// model name
	vec3_t				    bounds[2];				// model bounds
	int						contents;			// all contents of the model ored together
	qboolean				isConvex;			// set if model is convex
	// model geometry
	int						maxVertices;		// size of vertex array
	int						numVertices;		// number of vertices
	cm_vertex_t *			vertices;			// array with all vertices used by the model
	int						maxEdges;			// size of edge array
	int						numEdges;			// number of edges
	cm_edge_t *				edges;				// array with all edges used by the model
	cm_node_t *				node;				// first node of spatial subdivision
	// blocks with allocated memory
	cm_nodeBlock_t *		nodeBlocks;			// list with blocks of nodes
	cm_polygonRefBlock_t *	polygonRefBlocks;	// list with blocks of polygon references
	cm_brushRefBlock_t *	brushRefBlocks;		// list with blocks of brush references
	cm_polygonBlock_t *		polygonBlock;		// memory block with all polygons
	cm_brushBlock_t *		brushBlock;			// memory block with all brushes
	// statistics
	int						numPolygons;
	int						polygonMemory;
	int						numBrushes;
	int						brushMemory;
	int						numNodes;
	int						numBrushRefs;
	int						numPolygonRefs;
	int						numInternalEdges;
	int						numSharpEdges;
	int						numRemovedPolys;
	int						numMergedPolys;
	int						usedMemory;
} cm_model_t;

/*
===============================================================================

Data used during collision detection calculations

===============================================================================
*/

typedef struct cm_trmVertex_s {
	int used;										// true if this vertex is used for collision detection
	vec3_t p;										// vertex position
	vec3_t endp;									// end point of vertex after movement
	int polygonSide;								// side of polygon this vertex is on (rotational collision)
	vec6_t pl;				    					// pluecker coordinate for vertex movement
	vec3_t rotationOrigin;							// rotation origin for this vertex
	vec3_t rotationBounds[2];						// rotation bounds for this vertex
} cm_trmVertex_t;

typedef struct cm_trmEdge_s {
	int used;										// true when vertex is used for collision detection
	vec3_t start;									// start of edge
	vec3_t end;										// end of edge
	int vertexNum[2];								// indexes into cm_traceWork_t->vertices
	vec6_t pl;  									// pluecker coordinate for edge
	vec3_t cross;									// (z,-y,x) of cross product between edge dir and movement dir
	vec3_t rotationBounds[2];						// rotation bounds for this edge
	vec6_t plzaxis;								// pluecker coordinate for rotation about the z-axis
	unsigned short bitNum;							// vertex bit number
} cm_trmEdge_t;

typedef struct cm_trmPolygon_s {
	int used;
	plane_t plane;									// polygon plane
	int numEdges;									// number of edges
	int edges[MAX_TRACEMODEL_POLYEDGES];			// index into cm_traceWork_t->edges
	vec3_t rotationBounds[2];						// rotation bounds for this polygon
} cm_trmPolygon_t;

typedef struct cm_traceWork_s {
	int numVerts;
	cm_trmVertex_t vertices[MAX_TRACEMODEL_VERTS];	// trm vertices
	int numEdges;
	cm_trmEdge_t edges[MAX_TRACEMODEL_EDGES+1];		// trm edges
	int numPolys;
	cm_trmPolygon_t polys[MAX_TRACEMODEL_POLYS];	// trm polygons
	cm_model_t *model;								// model colliding with
	vec3_t start;									// start of trace
	vec3_t end;										// end of trace
	vec3_t dir;										// trace direction
	vec3_t bounds[2];								// bounds of full trace
	vec3_t size[2];									// bounds of transformed trm relative to start
	vec3_t extents; 								// largest of abs(size[0]) and abs(size[1]) for BSP trace
	int contents;									// ignore polygons that do not have any of these contents flags
	cm_trace_t trace;									// collision detection result

	qboolean rotation;									// true if calculating rotational collision
	qboolean pointTrace;								// true if only tracing a point
	qboolean positionTest;								// true if not tracing but doing a position test
	qboolean isConvex;									// true if the trace model is convex
	qboolean axisIntersectsTrm;							// true if the rotation axis intersects the trace model
	qboolean getContacts;								// true if retrieving contacts
	qboolean quickExit;									// set to quickly stop the collision detection calculations

	vec3_t origin;									// origin of rotation in model space
	vec3_t axis;									// rotation axis in model space
	vec3_t matrix[3];								// rotates axis of rotation to the z-axis
	float angle;									// angle for rotational collision
	float maxTan;									// max tangent of half the positive angle used instead of fraction
	float radius;									// rotation radius of trm start
	rotation_t modelVertexRotation;					// inverse rotation for model vertices

	contactInfo_t *contacts;						// array with contacts
	int maxContacts;								// max size of contact array
	int numContacts;								// number of contacts found

	plane_t heartPlane1;							// polygons should be near anough the trace heart planes
	float maxDistFromHeartPlane1;
	plane_t heartPlane2;
	float maxDistFromHeartPlane2;
	vec6_t polygonEdgePlueckerCache[CM_MAX_POLYGON_EDGES];
	vec6_t polygonVertexPlueckerCache[CM_MAX_POLYGON_EDGES];
	vec3_t polygonRotationOriginCache[CM_MAX_POLYGON_EDGES];
} cm_traceWork_t;

/*
===============================================================================

Collision Map

===============================================================================
*/

typedef struct cm_procNode_s {
	plane_t plane;
	int children[2];				// negative numbers are (-1 - areaNumber), 0 = solid
} cm_procNode_t;

// cm_translate.c
int				CM_TranslateEdgeThroughEdge( vec3_t cross, vec6_t l1, vec6_t l2, float *fraction );
void			CM_TranslateTrmEdgeThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmEdge_t *trmEdge );
void			CM_TranslateTrmVertexThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmVertex_t *v, int bitNum );
void			CM_TranslatePointThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmVertex_t *v );
void			CM_TranslateVertexThroughTrmPolygon( cm_traceWork_t *tw, cm_trmPolygon_t *trmpoly, cm_polygon_t *poly, cm_vertex_t *v, vec3_t endp, vec6_t pl );
qboolean		CM_TranslateTrmThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *p );
void			CM_SetupTranslationHeartPlanes( cm_traceWork_t *tw );
void			CM_SetupTrm( cm_traceWork_t *tw, const traceModel_t *trm );

// cm_rotate.c
int CM_CollisionBetweenEdgeBounds( cm_traceWork_t *tw, const vec3_t va, const vec3_t vb,
                                    const vec3_t vc, const vec3_t vd, float tanHalfAngle,
                                    vec3_t collisionPoint, vec3_t collisionNormal );
int				CM_RotateEdgeThroughEdge( cm_traceWork_t *tw, const vec6_t pl1,
                                    const vec3_t vc, const vec3_t vd,
                                    const float minTan, float *tanHalfAngle );
int				CM_EdgeFurthestFromEdge( cm_traceWork_t *tw, const vec6_t pl1,
											const vec3_t vc, const vec3_t vd,
											float *tanHalfAngle, float *dir );
void			CM_RotateTrmEdgeThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmEdge_t *trmEdge );
int				CM_RotatePointThroughPlane( const cm_traceWork_t *tw, const vec3_t point, const plane_t plane,
											const float angle, const float minTan, float *tanHalfAngle );
int				CM_PointFurthestFromPlane( const cm_traceWork_t *tw, const vec3_t point, const plane_t plane,
											const float angle, float *tanHalfAngle, float *dir );
int				CM_RotatePointThroughEpsilonPlane( const cm_traceWork_t *tw, const vec3_t point, const vec3_t endPoint,
											const plane_t plane, const float angle, const vec3_t origin,
											float *tanHalfAngle, vec3_t collisionPoint, vec3_t endDir );
void			CM_RotateTrmVertexThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmVertex_t *v, int vertexNum);
void			CM_RotateVertexThroughTrmPolygon( cm_traceWork_t *tw, cm_trmPolygon_t *trmpoly, cm_polygon_t *poly,
											cm_vertex_t *v, vec3_t rotationOrigin );
qboolean		CM_RotateTrmThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *p );
void			CM_BoundsForRotation( const vec3_t origin, const vec3_t axis, const vec3_t start, const vec3_t end, vec3_t bounds[2] );
void			CM_Rotation180( cm_trace_t *results, const vec3_t rorg, const vec3_t axis,
									const float startAngle, const float endAngle, const vec3_t start,
									const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
									cmHandle_t model, const vec3_t origin, const vec3_t modelAxis[3] );

// CollisionMap_contents.cpp
qboolean		CM_TestTrmVertsInBrush( cm_traceWork_t *tw, cm_brush_t *b );
qboolean		CM_TestTrmInPolygon( cm_traceWork_t *tw, cm_polygon_t *p );
cm_node_t *		CM_PointNode( const vec3_t p, cm_model_t *model );
int				CM_PointContents( const vec3_t p, cmHandle_t model );
int				CM_TransformedPointContents( const vec3_t p, cmHandle_t model, const vec3_t origin, const vec3_t modelAxis[3] );
int				CM_ContentsTrm( cm_trace_t *results, const vec3_t start,
                                const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
                                cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] );

// CollisionMap_trace.cpp
void			CM_TraceTrmThroughNode( cm_traceWork_t *tw, cm_node_t *node );
void			CM_TraceThroughAxialBSPTree_r( cm_traceWork_t *tw, cm_node_t *node, float p1f, float p2f, vec3_t p1, vec3_t p2);
void			CM_TraceThroughModel( cm_traceWork_t *tw );
void			CM_RecurseProcBSP_r( cm_trace_t *results, int parentNodeNum, int nodeNum, float p1f, float p2f, const vec3_t p1, const vec3_t p2 );

// CollisionMap_load.cpp
void			CM_Clear( void );
void			CM_FreeTrmModelStructure( void );
                // model deallocation
void			CM_RemovePolygonReferences_r( cm_node_t *node, cm_polygon_t *p );
void			CM_RemoveBrushReferences_r( cm_node_t *node, cm_brush_t *b );
void			CM_FreeNode( cm_node_t *node );
void			CM_FreePolygonReference( cm_polygonRef_t *pref );
void			CM_FreeBrushReference( cm_brushRef_t *bref );
void			CM_FreePolygon( cm_model_t *model, cm_polygon_t *poly );
void			CM_FreeBrush( cm_model_t *model, cm_brush_t *brush );
void			CM_FreeTree_r( cm_model_t *model, cm_node_t *headNode, cm_node_t *node );
void            CM_FreeMaterial( cm_material_t *material );
void			CM_FreeModel( cm_model_t *model );
                // merging polygons
void			CM_ReplacePolygons( cm_model_t *model, cm_node_t *node, cm_polygon_t *p1, cm_polygon_t *p2, cm_polygon_t *newp );
cm_polygon_t *	CM_TryMergePolygons( cm_model_t *model, cm_polygon_t *p1, cm_polygon_t *p2 );
qboolean		CM_MergePolygonWithTreePolygons( cm_model_t *model, cm_node_t *node, cm_polygon_t *polygon );
void			CM_MergeTreePolygons( cm_model_t *model, cm_node_t *node );
                // finding internal edges
qboolean		CM_PointInsidePolygon( cm_model_t *model, cm_polygon_t *p, vec3_t v );
void			CM_FindInternalEdgesOnPolygon( cm_model_t *model, cm_polygon_t *p1, cm_polygon_t *p2 );
void			CM_FindInternalPolygonEdges( cm_model_t *model, cm_node_t *node, cm_polygon_t *polygon );
void			CM_FindInternalEdges( cm_model_t *model, cm_node_t *node );
void			CM_FindContainedEdges( cm_model_t *model, cm_polygon_t *p );
                // loading of proc BSP tree
void			CM_ParseProcNodes( source_t *src );
void			CM_LoadProcBSP( const char *name );
                // removal of contained polygons
int				R_ChoppedAwayByProcBSP( int nodeNum, fixedWinding_t *w, const vec3_t normal, const vec3_t origin, const float radius );
int				ChoppedAwayByProcBSP( const fixedWinding_t *w, const plane_t plane, int contents );
void			ChopWindingListWithBrush( cm_windingList_t *list, cm_brush_t *b );
void			R_ChopWindingListWithTreeBrushes( cm_windingList_t *list, cm_node_t *node );
fixedWinding_t *WindingOutsideBrushes( fixedWinding_t *w, const plane_t plane, int contents, int patch, cm_node_t *headNode );
                // creation of axial BSP tree
cm_material_t*  CM_AllocMaterial( void );
cm_model_t *	CM_AllocModel( void );
cm_node_t *		AllocNode( cm_model_t *model, int blockSize );
cm_polygonRef_t*AllocPolygonReference( cm_model_t *model, int blockSize );
cm_brushRef_t *	AllocBrushReference( cm_model_t *model, int blockSize );
cm_polygon_t *	AllocPolygon( cm_model_t *model, int numEdges );
cm_brush_t *	AllocBrush( cm_model_t *model, int numPlanes );
void			AddPolygonToNode( cm_model_t *model, cm_node_t *node, cm_polygon_t *p );
void			AddBrushToNode( cm_model_t *model, cm_node_t *node, cm_brush_t *b );
void			SetupTrmModelStructure( void );
void			R_FilterPolygonIntoTree( cm_model_t *model, cm_node_t *node, cm_polygonRef_t *pref, cm_polygon_t *p );
void			R_FilterBrushIntoTree( cm_model_t *model, cm_node_t *node, cm_brushRef_t *pref, cm_brush_t *b );
cm_node_t *		R_CreateAxialBSPTree( cm_model_t *model, cm_node_t *node, const vec3_t bounds[2] );
cm_node_t *		CreateAxialBSPTree( cm_model_t *model, cm_node_t *node );
                // creation of raw polygons
void			SetupHash(void);
void			ShutdownHash(void);
void			ClearHash( vec3_t bounds[2] );
int				GetVertex( cm_model_t *model, const vec3_t v, int *vertexNum );
int				GetEdge( cm_model_t *model, const vec3_t v1, const vec3_t v2, int *edgeNum, int v1num );
void			CreatePolygon( cm_model_t *model, fixedWinding_t *w, const plane_t plane, const qhandle_t material, int primitiveNum );
void			PolygonFromWinding( cm_model_t *model, fixedWinding_t *w, const plane_t plane, const qhandle_t material, int primitiveNum );
void			CM_CalculateEdgeNormals( cm_model_t *model, cm_node_t *node );
//void			CreatePatchPolygons( cm_model_t *model, idSurface_Patch &mesh, const idMaterial *material, int primitiveNum );
//void			ConvertPatch( cm_model_t *model, const idMapPatch *patch, int primitiveNum );
void			ConvertBrushSides( cm_model_t *model, const mapPrimitive_t *mapBrush, int primitiveNum );
void			ConvertBrush( cm_model_t *model, const mapPrimitive_t *mapBrush, int primitiveNum );
void			PrintModelInfo( const cm_model_t *model );
void			AccumulateModelInfo( cm_model_t *model );
void			RemapEdges( cm_node_t *node, int *edgeRemap );
void			OptimizeArrays( cm_model_t *model );
void			FinishModel( cm_model_t *model );
void			BuildModels( const mapFile_t *mapFile );
cmHandle_t		CM_FindModel( const char *name );
cm_model_t *	CM_CollisionModelForMapEntity( const mapEntity_t *mapEnt );	// brush/patch model from .map
cm_model_t *	CM_LoadRenderModel( const char *fileName );					// ASE/LWO models
qboolean		CM_TrmFromModel_r( traceModel_t *trm, cm_node_t *node );
qboolean		CM_TrmFromModel2( const cm_model_t *model, traceModel_t *trm );

// cm_bsp.c
extern cvar_t *cm_noAreas;
extern cvar_t *cm_noCurves;
extern cvar_t *cm_playerCurveClip;

void		CM_DrawDebugSurface( void (*drawPoly)(int color, int numPoints, float *points) );

int			CM_NumClusters (void);
int			CM_NumInlineModels( void );
char		*CM_EntityString (void);
qboolean	CM_GetEntityToken( int *parseOffset, char *token, int size );

byte		*CM_ClusterPVS (int cluster);

int			CM_PointLeafnum( const vec3_t p );

// only returns non-solid leafs
// overflow if return listsize and if *lastLeaf != list[listsize-1]
int			CM_BoxLeafnums( const vec3_t mins, const vec3_t maxs, int *list,
		 					int listsize, int *lastLeaf );

int			CM_LeafCluster (int leafnum);
int			CM_LeafArea (int leafnum);

void		CM_AdjustAreaPortalState( int area1, int area2, qboolean open );
qboolean	CM_AreasConnected( int area1, int area2 );

int			CM_WriteAreaBits( byte *buffer, int area );

qhandle_t   CM_RegisterMaterial( const char *name );
void		CM_GetMaterialName( qhandle_t hShader, char *buffer, int bufferSize );
int			CM_GetMaterialContentFlags( qhandle_t material );
int         CM_GetMaterialSurfaceFlags( qhandle_t material );
qboolean	CM_MaterialNeedsBackSide( qhandle_t material );

// cm_files.c
                // writing
void			CM_WriteNodes( fileHandle_t fp, cm_node_t *node );
int				CM_CountPolygonMemory( cm_node_t *node );
void			CM_WritePolygons( fileHandle_t fp, cm_node_t *node );
int				CM_CountBrushMemory( cm_node_t *node );
void			CM_WriteBrushes( fileHandle_t fp, cm_node_t *node );
void			CM_WriteCollisionModel( fileHandle_t fp, cm_model_t *model );
void			CM_WriteCollisionModelsToFile( const char *filename, int firstModel, int lastModel, unsigned int mapFileCRC );
                // loading
cm_node_t *		CM_ParseNodes( source_t *src, cm_model_t *model, cm_node_t *parent );
void			CM_ParseVertices( source_t *src, cm_model_t *model );
void			CM_ParseEdges( source_t *src, cm_model_t *model );
void			CM_ParsePolygons( source_t *src, cm_model_t *model );
void			CM_ParseBrushes( source_t *src, cm_model_t *model );
qboolean		CM_ParseCollisionModel( source_t *src );
qboolean		CM_LoadCollisionModelFile( const char *name, unsigned int mapFileCRC );

// cm_debug.c

extern cvar_t	*cm_drawMask;
extern cvar_t	*cm_drawColor;
extern cvar_t	*cm_drawFilled;
extern cvar_t	*cm_drawInternal;
extern cvar_t	*cm_drawNormals;
extern cvar_t	*cm_backFaceCull;
extern cvar_t	*cm_debugCollision;

int				CM_ContentsFromString( const char *string );
const char *	CM_StringFromContents( const int contents );
void			CM_DrawEdge( cm_model_t *model, int edgeNum, const vec3_t origin, const vec3_t axis[3] );
void			CM_DrawPolygon( cm_model_t *model, cm_polygon_t *p, const vec3_t origin, const vec3_t axis[3],
                            const vec3_t viewOrigin );
void			CM_DrawNodePolygons( cm_model_t *model, cm_node_t *node, const vec3_t origin, const vec3_t axis[3],
                            const vec3_t viewOrigin, const float radius );

typedef struct {
	cplane_t	*plane;
	int			children[2];		// negative numbers are leafs
} cNode_t;

typedef struct {
	int			cluster;
	int			area;

	int			firstLeafBrush;
	int			numLeafBrushes;

	int			firstLeafSurface;
	int			numLeafSurfaces;
} cLeaf_t;

typedef struct {
	int			floodnum;
	int			floodvalid;
} cmArea_t;

typedef struct leafList_s {
	int		count;
	int		maxcount;
	qboolean	overflowed;
	int		*list;
	vec3_t	bounds[2];
	int		lastLeaf;		// for overflows where each leaf can't be stored individually
	void	(*storeLeafs)( struct leafList_s *ll, int nodenum );
} leafList_t;

// collision map data
typedef struct cmLocal_s {
	char			mapName[MAX_QPATH];
	int			    mapFileTime;
	int				loaded;
					// for multi-check avoidance
	int				checkCount;
                    // materials
    int             maxMaterials;
    int             numMaterials;
    cm_material_t **materials;
    hashIndex_t     materialsHash;
					// models
	int				maxModels;
	int				numModels;
	cm_model_t **	models;
					// polygons and brush for trm model
	cm_polygonRef_t*trmPolygons[MAX_TRACEMODEL_POLYS];
	cm_brushRef_t *	trmBrushes[1];
	qhandle_t		trmMaterial;
					// for data pruning
	int				numProcNodes;
	cm_procNode_t *	procNodes;
					// for retrieving contact points
	qboolean		getContacts;
	contactInfo_t *	contacts;
	int				maxContacts;
	int				numContacts;

	//
	// cm1 compatibility (PVS & linking/unlinking entities)
	//
	int			numPlanes;
	cplane_t	*planes;

	int			numNodes;
	cNode_t		*nodes;

	int			numLeafs;
	cLeaf_t		*leafs;

	int			numClusters;
	int			clusterBytes;
	byte		*visibility;
	qboolean	vised;			// if false, visibility is just a single cluster of ffs

	int			numEntityChars;
	char		*entityString;

	int			numAreas;
	cmArea_t	*areas;
	int			*areaPortals;	// [ numAreas*numAreas ] reference counts

	int			floodvalid;
} cmLocal_t;

extern cmLocal_t    cmLocal;
extern const bspFile_t	*cm_bsp;
extern cmimport_t	cmi;

void CM_Init( void );
void CM_Shutdown( void );

// Loads collision models from a map file.
void CM_LoadMap2( const mapFile_t *mapFile, qboolean clear );
void CM_LoadMap( const mapFile_t *mapFile );
void CM_LoadBSP( const char *fileName, int *checksum );
// Frees all the collision models.
void CM_FreeMap( void );

// Gets the clip handle for a model.
cmHandle_t CM_LoadModel( const char *modelName, const qboolean precache );
// Sets up a trace model for collision with other trace models.
cmHandle_t CM_SetupTrmModel( const traceModel_t *trm, qhandle_t material );
// Creates a trace model from a collision model, returns true if succesfull.
qboolean CM_TrmFromModel( const char *modelName, traceModel_t *trm );

// Gets the name of a model.
const char *CM_GetModelName( cmHandle_t model );
// Gets the bounds of a model.
qboolean CM_GetModelBounds( cmHandle_t model, vec3_t bounds[2] );
// Gets all contents flags of brushes and polygons of a model ored together.
qboolean CM_GetModelContents( cmHandle_t model, int *contents );
// Gets a vertex of a model.
qboolean CM_GetModelVertex( cmHandle_t model, int vertexNum, vec3_t vertex );
// Gets an edge of a model.
qboolean CM_GetModelEdge( cmHandle_t model, int edgeNum, vec3_t start, vec3_t end );
// Gets a polygon of a model.
qboolean CM_GetModelPolygon( cmHandle_t model, int polygonNum, fixedWinding_t *winding );

// Translates a trace model and reports the first collision if any.
void CM_Translation( cm_trace_t *results, const vec3_t start, const vec3_t end,
								const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
								cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] );
// Rotates a trace model and reports the first collision if any.
void CM_Rotation( cm_trace_t *results, const vec3_t start, const rotation_t *rotation,
								const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
								cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] );
// Returns the contents touched by the trace model or 0 if the trace model is in free space.
int CM_Contents( const vec3_t start,
                            const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
                            cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] );
// Stores all contact points of the trace model with the model, returns the number of contacts.
int CM_Contacts( contactInfo_t *contacts, const int maxContacts, const vec3_t start, const vec6_t dir, const float depth,
                            const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
                            cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] );

// Tests collision detection.
void CM_DebugOutput( const vec3_t origin );
// Draws a model.
void CM_DrawModel( cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3],
												const vec3_t viewOrigin, const float radius );
// Prints model information, use -1 handle for accumulated model info.
void CM_ModelInfo( cmHandle_t model );
// Lists all loaded models.
void CM_ListModels( void );
// Writes a collision model file for the given map entity.
qboolean CM_WriteCollisionModelForMapEntity( const mapEntity_t *mapEnt, const char *filename, const qboolean testTraceModel );

#endif /* !__CM_LOCAL_H__ */
