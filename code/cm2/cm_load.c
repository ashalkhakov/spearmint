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

/*
===============================================================================

	Trace model vs. polygonal model collision detection.

	It is more important to minimize the number of collision polygons
	than it is to minimize the number of edges used for collision
	detection (total edges - internal edges).

	Stitching the world tends to minimize the number of edges used
	for collision detection (more internal edges). However stitching
	also results in more collision polygons which usually makes a
	stitched world slower.

	In an average map over 30% of all edges is internal.

===============================================================================
*/

#include "cm_local.h"
#include "../idlib/l_script.h"
#include "../idlib/l_precomp.h"

#define PROC_FILE_EXT				"proc"
#define	PROC_FILE_ID				"mapProcFile003"



cm_windingList_t *				cm_windingList;
cm_windingList_t *				cm_outList;
cm_windingList_t *				cm_tmpList;

hashIndex_t * 					cm_vertexHash;
hashIndex_t * 					cm_edgeHash;

vec3_t						    cm_modelBounds[2];
int								cm_vertexShift;


/*
===============================================================================

Proc BSP tree for data pruning

===============================================================================
*/

/*
================
HashIndexGenerateKeyForPair
================
*/
static ID_INLINE int HashIndexGenerateKeyForPair( hashIndex_t *h, const int n1, const int n2 ) {
	return ( ( n1 + n2 ) & h->hashMask );
}

/*
================
CM_ParseProcNodes
================
*/
void CM_ParseProcNodes( source_t *src ) {
	int i;

	PC_ExpectTokenString( src, "{" );

	cmLocal.numProcNodes = PC_ParseInt( src );
	if ( cmLocal.numProcNodes < 0 ) {
		PC_SourceError( src, "CM_ParseProcNodes: bad numProcNodes" );
	}
	cmLocal.procNodes = (cm_procNode_t *)ii.GetMemory( cmLocal.numProcNodes * sizeof( cm_procNode_t ) );
    memset( cmLocal.procNodes, 0, cmLocal.numProcNodes * sizeof( cm_procNode_t ) );

	for ( i = 0; i < cmLocal.numProcNodes; i++ ) {
		cm_procNode_t *node;

		node = &cmLocal.procNodes[i];

		PC_Parse1DMatrix( src, 4, node->plane );
		node->children[0] = PC_ParseInt( src );
		node->children[1] = PC_ParseInt( src );
	}

	PC_ExpectTokenString( src, "}" );
}

/*
================
CM_LoadProcBSP

  FIXME: if the nodes would be at the start of the .proc file it would speed things up considerably
================
*/
void CM_LoadProcBSP( const char *name ) {
	char filename[MAX_QPATH];
	token_t token;
	source_t *src;

	// load it
    Q_strncpyz( filename, name, sizeof( filename ) );
    COM_SetExtension( filename, sizeof( filename ), "." PROC_FILE_EXT );
	src = PC_LoadSource( filename, NULL, NULL );
    if ( !src ) {
		ii.Com_Printf( "CM_LoadProcBSP: couldn't load %s", filename );
		return;
	}
	//src = new idLexer( filename, LEXFL_NOSTRINGCONCAT | LEXFL_NODOLLARPRECOMPILE );
	SetScriptFlags( src->scriptstack, SCFL_NOSTRINGWHITESPACES );

	if ( !PC_ReadToken( src, &token ) || Q_stricmp( token.string, PROC_FILE_ID ) ) {
		ii.Com_Printf( "CM_LoadProcBSP: bad id '%s' instead of '%s'", token.string, PROC_FILE_ID );
		PC_FreeSource( src );
		return;
	}

	// parse the file
	while ( 1 ) {
		if ( !PC_ReadToken( src, &token ) ) {
			break;
		}

		if ( !strcmp( token.string, "model" ) ) {
			PC_SkipBracedSection( src, qtrue );
			continue;
		}

		if ( !strcmp( token.string, "shadowModel" ) ) {
			PC_SkipBracedSection( src, qtrue );
			continue;
		}

		if ( !strcmp( token.string, "interAreaPortals" ) ) {
			PC_SkipBracedSection( src, qtrue );
			continue;
		}

		if ( !strcmp( token.string, "nodes" ) ) {
			CM_ParseProcNodes( src );
			break;
		}

		ii.Com_Error( ERR_DROP, "CM_LoadProcBSP: bad token \"%s\"", token.string );
	}

	PC_FreeSource( src );
}

/*
===============================================================================

Free map

===============================================================================
*/

/*
================
CM_Clear
================
*/
void CM_Clear( void ) {
	cmLocal.mapName[0] = 0;
	cmLocal.mapFileTime = 0;
	cmLocal.loaded = 0;
	cmLocal.checkCount = 0;
    cmLocal.maxMaterials = 0;
    cmLocal.numMaterials = 0;
    cmLocal.materials = NULL;
    HashIndexFree( &cmLocal.materialsHash );
	cmLocal.maxModels = 0;
	cmLocal.numModels = 0;
	cmLocal.models = NULL;
	memset( cmLocal.trmPolygons, 0, sizeof( cmLocal.trmPolygons ) );
	cmLocal.trmBrushes[0] = NULL;
	cmLocal.trmMaterial = 0;
	cmLocal.numProcNodes = 0;
	cmLocal.procNodes = NULL;
	cmLocal.getContacts = qfalse;
	cmLocal.contacts = NULL;
	cmLocal.maxContacts = 0;
	cmLocal.numContacts = 0;

	BSP_Free( cm_bsp );
	cm_bsp = NULL;

	cmLocal.numPlanes = 0;
	cmLocal.planes = NULL;
	cmLocal.numNodes = 0;
	cmLocal.nodes = NULL;
	cmLocal.numLeafs = 0;
	cmLocal.leafs = NULL;
	cmLocal.numClusters = 0;
	cmLocal.clusterBytes = 0;
	cmLocal.visibility = NULL;
	cmLocal.vised = qfalse;
	cmLocal.numEntityChars = 0;
	cmLocal.entityString = NULL;
	cmLocal.numAreas = 0;
	cmLocal.areas = NULL;
	cmLocal.areaPortals = NULL;
	cmLocal.floodvalid = 0;
}

/*
================
CM_RemovePolygonReferences_r
================
*/
void CM_RemovePolygonReferences_r( cm_node_t *node, cm_polygon_t *p ) {
	cm_polygonRef_t *pref;

	while( node ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			if ( pref->p == p ) {
				pref->p = NULL;
				// cannot return here because we can have links down the tree due to polygon merging
				//return;
			}
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		if ( p->bounds[0][node->planeType] > node->planeDist ) {
			node = node->children[0];
		}
		else if ( p->bounds[1][node->planeType] < node->planeDist ) {
			node = node->children[1];
		}
		else {
			CM_RemovePolygonReferences_r( node->children[1], p );
			node = node->children[0];
		}
	}
}

/*
================
CM_RemoveBrushReferences_r
================
*/
void CM_RemoveBrushReferences_r( cm_node_t *node, cm_brush_t *b ) {
	cm_brushRef_t *bref;

	while( node ) {
		for ( bref = node->brushes; bref; bref = bref->next ) {
			if ( bref->b == b ) {
				bref->b = NULL;
				return;
			}
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		if ( b->bounds[0][node->planeType] > node->planeDist ) {
			node = node->children[0];
		}
		else if ( b->bounds[1][node->planeType] < node->planeDist ) {
			node = node->children[1];
		}
		else {
			CM_RemoveBrushReferences_r( node->children[1], b );
			node = node->children[0];
		}
	}
}

/*
================
CM_FreeNode
================
*/
void CM_FreeNode( cm_node_t *node ) {
	// don't free the node here
	// the nodes are allocated in blocks which are freed when the model is freed
}

/*
================
CM_FreePolygonReference
================
*/
void CM_FreePolygonReference( cm_polygonRef_t *pref ) {
	// don't free the polygon reference here
	// the polygon references are allocated in blocks which are freed when the model is freed
}

/*
================
CM_FreeBrushReference
================
*/
void CM_FreeBrushReference( cm_brushRef_t *bref ) {
	// don't free the brush reference here
	// the brush references are allocated in blocks which are freed when the model is freed
}

/*
================
CM_FreePolygon
================
*/
void CM_FreePolygon( cm_model_t *model, cm_polygon_t *poly ) {
	model->numPolygons--;
	model->polygonMemory -= sizeof( cm_polygon_t ) + ( poly->numEdges - 1 ) * sizeof( poly->edges[0] );
	if ( model->polygonBlock == NULL ) {
		ii.FreeMemory( poly );
	}
}

/*
================
CM_FreeBrush
================
*/
void CM_FreeBrush( cm_model_t *model, cm_brush_t *brush ) {
	model->numBrushes--;
	model->brushMemory -= sizeof( cm_brush_t ) + ( brush->numPlanes - 1 ) * sizeof( brush->planes[0] );
	if ( model->brushBlock == NULL ) {
		ii.FreeMemory( brush );
	}
}

/*
================
CM_FreeTree_r
================
*/
void CM_FreeTree_r( cm_model_t *model, cm_node_t *headNode, cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;
	cm_brushRef_t *bref;
	cm_brush_t *b;

	// free all polygons at this node
	for ( pref = node->polygons; pref; pref = node->polygons ) {
		p = pref->p;
		if ( p ) {
			// remove all other references to this polygon
			CM_RemovePolygonReferences_r( headNode, p );
			CM_FreePolygon( model, p );
		}
		node->polygons = pref->next;
		CM_FreePolygonReference( pref );
	}
	// free all brushes at this node
	for ( bref = node->brushes; bref; bref = node->brushes ) {
		b = bref->b;
		if ( b ) {
			// remove all other references to this brush
			CM_RemoveBrushReferences_r( headNode, b );
			CM_FreeBrush( model, b );
		}
		node->brushes = bref->next;
		CM_FreeBrushReference( bref );
	}
	// recurse down the tree
	if ( node->planeType != -1 ) {
		CM_FreeTree_r( model, headNode, node->children[0] );
		node->children[0] = NULL;
		CM_FreeTree_r( model, headNode, node->children[1] );
		node->children[1] = NULL;
	}
	CM_FreeNode( node );
}

/*
================
CM_FreeMaterial
================
*/
void CM_FreeMaterial( cm_material_t *material ) {
    // FIXME: remove from the hash?
	// free the material
	ii.FreeMemory( material );
}

/*
================
CM_FreeModel
================
*/
void CM_FreeModel( cm_model_t *model ) {
	cm_polygonRefBlock_t *polygonRefBlock, *nextPolygonRefBlock;
	cm_brushRefBlock_t *brushRefBlock, *nextBrushRefBlock;
	cm_nodeBlock_t *nodeBlock, *nextNodeBlock;

	// free the tree structure
	if ( model->node ) {
		CM_FreeTree_r( model, model->node, model->node );
	}
	// free blocks with polygon references
	for ( polygonRefBlock = model->polygonRefBlocks; polygonRefBlock; polygonRefBlock = nextPolygonRefBlock ) {
		nextPolygonRefBlock = polygonRefBlock->next;
		ii.FreeMemory( polygonRefBlock );
	}
	// free blocks with brush references
	for ( brushRefBlock = model->brushRefBlocks; brushRefBlock; brushRefBlock = nextBrushRefBlock ) {
		nextBrushRefBlock = brushRefBlock->next;
		ii.FreeMemory( brushRefBlock );
	}
	// free blocks with nodes
	for ( nodeBlock = model->nodeBlocks; nodeBlock; nodeBlock = nextNodeBlock ) {
		nextNodeBlock = nodeBlock->next;
		ii.FreeMemory( nodeBlock );
	}
	// free block allocated polygons
    if ( model->polygonBlock ) {
	    ii.FreeMemory( model->polygonBlock );
    }
	// free block allocated brushes
    if ( model->brushBlock ) {
	    ii.FreeMemory( model->brushBlock );
    }
	// free edges
	ii.FreeMemory( model->edges );
	// free vertices
	ii.FreeMemory( model->vertices );
	// free the model
	ii.FreeMemory( model );
}

/*
================
CM_FreeMap
================
*/
void CM_FreeMap( void ) {
	int i;

	if ( !cmLocal.loaded ) {
		CM_Clear();
		return;
	}

	for ( i = 0; i < cmLocal.maxModels; i++ ) {
		if ( !cmLocal.models[i] ) {
			continue;
		}
		CM_FreeModel( cmLocal.models[i] );
	}

	CM_FreeTrmModelStructure();

	ii.FreeMemory( cmLocal.models );

    for ( i = 0; i < cmLocal.maxMaterials; i++ ) {
        if ( !cmLocal.materials[i] ) {
            continue;
        }
        CM_FreeMaterial( cmLocal.materials[i] );
    }

    ii.FreeMemory( cmLocal.materials );

	CM_Clear();

	ShutdownHash();
}

/*
================
FreeTrmModelStructure
================
*/
void CM_FreeTrmModelStructure( void ) {
	int i;

	assert( cmLocal.models );
	if ( !cmLocal.models[MAX_SUBMODELS] ) {
		return;
	}

	for ( i = 0; i < MAX_TRACEMODEL_POLYS; i++ ) {
		CM_FreePolygon( cmLocal.models[MAX_SUBMODELS], cmLocal.trmPolygons[i]->p );
	}
	CM_FreeBrush( cmLocal.models[MAX_SUBMODELS], cmLocal.trmBrushes[0]->b );

	cmLocal.models[MAX_SUBMODELS]->node->polygons = NULL;
	cmLocal.models[MAX_SUBMODELS]->node->brushes = NULL;
	CM_FreeModel( cmLocal.models[MAX_SUBMODELS] );
}


/*
===============================================================================

Edge normals

===============================================================================
*/

/*
================
CM_CalculateEdgeNormals
================
*/
#define SHARP_EDGE_DOT	-0.7f

void CM_CalculateEdgeNormals( cm_model_t *model, cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;
	cm_edge_t *edge;
	float dot, s;
	int i, edgeNum;
	vec3_t dir, tmp1, tmp2, tmp3;

	while( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			p = pref->p;
			// if we checked this polygon already
			if ( p->checkcount == cmLocal.checkCount ) {
				continue;
			}
			p->checkcount = cmLocal.checkCount;

			for ( i = 0; i < p->numEdges; i++ ) {
				edgeNum = p->edges[i];
				edge = model->edges + abs( edgeNum );
				if ( edge->normal[0] == 0.0f && edge->normal[1] == 0.0f && edge->normal[2] == 0.0f ) {
					// if the edge is only used by this polygon
					if ( edge->numUsers == 1 ) {
                        VectorSubtract( model->vertices[ edge->vertexNum[edgeNum < 0]].p, model->vertices[ edge->vertexNum[edgeNum > 0]].p, dir );
                        CrossProduct( p->plane, dir, edge->normal );
						VectorNormalize( edge->normal );
					} else {
						// the edge is used by more than one polygon
                        VectorCopy( p->plane, edge->normal );
					}
				} else {
					dot = DotProduct( edge->normal, p->plane );
					// if the two planes make a very sharp edge
					if ( dot < SHARP_EDGE_DOT ) {
						// max length normal pointing outside both polygons
                        VectorSubtract( model->vertices[ edge->vertexNum[edgeNum > 0]].p, model->vertices[ edge->vertexNum[edgeNum < 0]].p, dir );
                        CrossProduct( edge->normal, dir, tmp1 );
                        VectorNegate( dir, tmp2 );
                        CrossProduct( p->plane, tmp2, tmp3 );
						VectorAdd( tmp1, tmp3, edge->normal );

                        VectorScale( edge->normal, ( 0.5f / ( 0.5f + 0.5f * SHARP_EDGE_DOT ) ) / VectorLength( edge->normal ), edge->normal );

						model->numSharpEdges++;
					} else {
						s = 0.5f / ( 0.5f + 0.5f * dot );
                        VectorAdd( edge->normal, p->plane, tmp1 );
                        VectorScale( tmp1, s, edge->normal );
					}
				}
			}
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		CM_CalculateEdgeNormals( model, node->children[1] );
		node = node->children[0];
	}
}

/*
===============================================================================

Trace model to general collision model

===============================================================================
*/

/*
================
CM_AllocMaterial
================
*/
cm_material_t *CM_AllocMaterial( void ) {
    cm_material_t *material;

    material = ii.GetMemory( sizeof( *material ) );
    material->contentFlags = 0;
    material->surfaceFlags = 0;
    material->name[0] = 0;

    return material;
}

/*
================
CM_AllocModel
================
*/
cm_model_t *CM_AllocModel( void ) {
	cm_model_t *model;

	model = ii.GetMemory( sizeof( *model ) );
	model->contents = 0;
	model->isConvex = qfalse;
	model->maxVertices = 0;
	model->numVertices = 0;
	model->vertices = NULL;
	model->maxEdges = 0;
	model->numEdges = 0;
	model->edges= NULL;
	model->node = NULL;
	model->nodeBlocks = NULL;
	model->polygonRefBlocks = NULL;
	model->brushRefBlocks = NULL;
	model->polygonBlock = NULL;
	model->brushBlock = NULL;
	model->numPolygons = model->polygonMemory =
	model->numBrushes = model->brushMemory =
	model->numNodes = model->numBrushRefs =
	model->numPolygonRefs = model->numInternalEdges =
	model->numSharpEdges = model->numRemovedPolys =
	model->numMergedPolys = model->usedMemory = 0;

	return model;
}

/*
================
AllocNode
================
*/
cm_node_t *AllocNode( cm_model_t *model, int blockSize ) {
	int i;
	cm_node_t *node;
	cm_nodeBlock_t *nodeBlock;

	if ( !model->nodeBlocks || !model->nodeBlocks->nextNode ) {
		nodeBlock = (cm_nodeBlock_t *) ii.GetMemory( sizeof( cm_nodeBlock_t ) + blockSize * sizeof(cm_node_t) );
        memset( nodeBlock, 0, sizeof( cm_nodeBlock_t ) + blockSize * sizeof(cm_node_t) );
		nodeBlock->nextNode = (cm_node_t *) ( ( (byte *) nodeBlock ) + sizeof( cm_nodeBlock_t ) );
		nodeBlock->next = model->nodeBlocks;
		model->nodeBlocks = nodeBlock;
		node = nodeBlock->nextNode;
		for ( i = 0; i < blockSize - 1; i++ ) {
			node->parent = node + 1;
			node = node->parent;
		}
		node->parent = NULL;
	}

	node = model->nodeBlocks->nextNode;
	model->nodeBlocks->nextNode = node->parent;
	node->parent = NULL;

	return node;
}

/*
================
AllocPolygonReference
================
*/
cm_polygonRef_t *AllocPolygonReference( cm_model_t *model, int blockSize ) {
	int i;
	cm_polygonRef_t *pref;
	cm_polygonRefBlock_t *prefBlock;

	if ( !model->polygonRefBlocks || !model->polygonRefBlocks->nextRef ) {
		prefBlock = (cm_polygonRefBlock_t *) ii.GetMemory( sizeof( cm_polygonRefBlock_t ) + blockSize * sizeof(cm_polygonRef_t) );
        memset( prefBlock, 0, sizeof( cm_polygonRefBlock_t ) + blockSize * sizeof(cm_polygonRef_t) );
		prefBlock->nextRef = (cm_polygonRef_t *) ( ( (byte *) prefBlock ) + sizeof( cm_polygonRefBlock_t ) );
		prefBlock->next = model->polygonRefBlocks;
		model->polygonRefBlocks = prefBlock;
		pref = prefBlock->nextRef;
		for ( i = 0; i < blockSize - 1; i++ ) {
			pref->next = pref + 1;
			pref = pref->next;
		}
		pref->next = NULL;
	}

	pref = model->polygonRefBlocks->nextRef;
	model->polygonRefBlocks->nextRef = pref->next;

	return pref;
}

/*
================
AllocBrushReference
================
*/
cm_brushRef_t *AllocBrushReference( cm_model_t *model, int blockSize ) {
	int i;
	cm_brushRef_t *bref;
	cm_brushRefBlock_t *brefBlock;

	if ( !model->brushRefBlocks || !model->brushRefBlocks->nextRef ) {
		brefBlock = (cm_brushRefBlock_t *) ii.GetMemory( sizeof(cm_brushRefBlock_t) + blockSize * sizeof(cm_brushRef_t) );
        memset( brefBlock, 0, sizeof(cm_brushRefBlock_t) + blockSize * sizeof(cm_brushRef_t) );
		brefBlock->nextRef = (cm_brushRef_t *) ( ( (byte *) brefBlock ) + sizeof(cm_brushRefBlock_t) );
		brefBlock->next = model->brushRefBlocks;
		model->brushRefBlocks = brefBlock;
		bref = brefBlock->nextRef;
		for ( i = 0; i < blockSize - 1; i++ ) {
			bref->next = bref + 1;
			bref = bref->next;
		}
		bref->next = NULL;
	}

	bref = model->brushRefBlocks->nextRef;
	model->brushRefBlocks->nextRef = bref->next;

	return bref;
}

/*
================
AllocPolygon
================
*/
cm_polygon_t *AllocPolygon( cm_model_t *model, int numEdges ) {
	cm_polygon_t *poly;
	int size;

	size = sizeof( cm_polygon_t ) + ( numEdges - 1 ) * sizeof( poly->edges[0] );
	model->numPolygons++;
	model->polygonMemory += size;
	if ( model->polygonBlock && model->polygonBlock->bytesRemaining >= size ) {
		poly = (cm_polygon_t *) model->polygonBlock->next;
		model->polygonBlock->next += size;
		model->polygonBlock->bytesRemaining -= size;
	} else {
		poly = (cm_polygon_t *) ii.GetMemory( size );
	}
	return poly;
}

/*
================
AllocBrush
================
*/
cm_brush_t *AllocBrush( cm_model_t *model, int numPlanes ) {
	cm_brush_t *brush;
	int size;

	size = sizeof( cm_brush_t ) + ( numPlanes - 1 ) * sizeof( brush->planes[0] );
	model->numBrushes++;
	model->brushMemory += size;
	if ( model->brushBlock && model->brushBlock->bytesRemaining >= size ) {
		brush = (cm_brush_t *) model->brushBlock->next;
		model->brushBlock->next += size;
		model->brushBlock->bytesRemaining -= size;
	} else {
		brush = (cm_brush_t *) ii.GetMemory( size );
	}
	return brush;
}

/*
================
AddPolygonToNode
================
*/
void AddPolygonToNode( cm_model_t *model, cm_node_t *node, cm_polygon_t *p ) {
	cm_polygonRef_t *pref;

	pref = AllocPolygonReference( model, model->numPolygonRefs < REFERENCE_BLOCK_SIZE_SMALL ? REFERENCE_BLOCK_SIZE_SMALL : REFERENCE_BLOCK_SIZE_LARGE );
	pref->p = p;
	pref->next = node->polygons;
	node->polygons = pref;
	model->numPolygonRefs++;
}

/*
================
AddBrushToNode
================
*/
void AddBrushToNode( cm_model_t *model, cm_node_t *node, cm_brush_t *b ) {
	cm_brushRef_t *bref;

	bref = AllocBrushReference( model, model->numBrushRefs < REFERENCE_BLOCK_SIZE_SMALL ? REFERENCE_BLOCK_SIZE_SMALL : REFERENCE_BLOCK_SIZE_LARGE );
	bref->b = b;
	bref->next = node->brushes;
	node->brushes = bref;
	model->numBrushRefs++;
}

/*
================
SetupTrmModelStructure
================
*/
void SetupTrmModelStructure( void ) {
	int i;
	cm_node_t *node;
	cm_model_t *model;

	// setup model
	model = CM_AllocModel();

	assert( cmLocal.models );
	cmLocal.models[MAX_SUBMODELS] = model;
	// create node to hold the collision data
	node = (cm_node_t *) AllocNode( model, 1 );
	node->planeType = -1;
	model->node = node;
	// allocate vertex and edge arrays
	model->numVertices = 0;
	model->maxVertices = MAX_TRACEMODEL_VERTS;
	model->vertices = (cm_vertex_t *) ii.GetMemory( model->maxVertices * sizeof(cm_vertex_t) );
    memset( model->vertices, 0, model->maxVertices * sizeof(cm_vertex_t) );
	model->numEdges = 0;
	model->maxEdges = MAX_TRACEMODEL_EDGES+1;
	model->edges = (cm_edge_t *) ii.GetMemory( model->maxEdges * sizeof(cm_edge_t) );
    memset(  model->edges, 0, model->maxEdges * sizeof(cm_edge_t) );
	// create a material for the trace model polygons
	cmLocal.trmMaterial = CM_RegisterMaterial( "_tracemodel" ); // TODO: makeDefault = false?
	if ( !cmLocal.trmMaterial ) {
		ii.Com_Error( ERR_DROP, "_tracemodel material not found" );
	}

	// allocate polygons
	for ( i = 0; i < MAX_TRACEMODEL_POLYS; i++ ) {
		cmLocal.trmPolygons[i] = AllocPolygonReference( model, MAX_TRACEMODEL_POLYS );
		cmLocal.trmPolygons[i]->p = AllocPolygon( model, MAX_TRACEMODEL_POLYEDGES );
        ClearBounds( cmLocal.trmPolygons[i]->p->bounds[0], cmLocal.trmPolygons[i]->p->bounds[1] );
        PlaneZero( cmLocal.trmPolygons[i]->p->plane );
		cmLocal.trmPolygons[i]->p->checkcount = 0;
		cmLocal.trmPolygons[i]->p->contents = -1;		// all contents
		cmLocal.trmPolygons[i]->p->material = cmLocal.trmMaterial;
		cmLocal.trmPolygons[i]->p->numEdges = 0;
	}
	// allocate brush for position test
	cmLocal.trmBrushes[0] = AllocBrushReference( model, 1 );
	cmLocal.trmBrushes[0]->b = AllocBrush( model, MAX_TRACEMODEL_POLYS );
	cmLocal.trmBrushes[0]->b->primitiveNum = 0;
	ClearBounds( cmLocal.trmBrushes[0]->b->bounds[0], cmLocal.trmBrushes[0]->b->bounds[1] );
	cmLocal.trmBrushes[0]->b->checkcount = 0;
	cmLocal.trmBrushes[0]->b->contents = -1;		// all contents
	cmLocal.trmBrushes[0]->b->numPlanes = 0;
}

/*
================
CM_SetupTrmModel

Trace models (item boxes, etc) are converted to collision models on the fly, using the last model slot
as a reusable temporary buffer
================
*/
cmHandle_t CM_SetupTrmModel( const traceModel_t *trm, qhandle_t material ) {
	int i, j;
	cm_vertex_t *vertex;
	cm_edge_t *edge;
	cm_polygon_t *poly;
	cm_model_t *model;
	const traceModelVert_t *trmVert;
	const traceModelEdge_t *trmEdge;
	const traceModelPoly_t *trmPoly;

	assert( cmLocal.models );

	if ( material == 0 ) {
		material = cmLocal.trmMaterial;
	}

	model = cmLocal.models[MAX_SUBMODELS];
	model->node->brushes = NULL;
	model->node->polygons = NULL;
	// if not a valid trace model
	if ( trm->type == TRM_INVALID || !trm->numPolys ) {
		return TRACE_MODEL_HANDLE;
	}
	// vertices
	model->numVertices = trm->numVerts;
	vertex = model->vertices;
	trmVert = trm->verts;
	for ( i = 0; i < trm->numVerts; i++, vertex++, trmVert++ ) {
        VectorCopy( *trmVert, vertex->p );
		vertex->sideSet = 0;
	}
	// edges
	model->numEdges = trm->numEdges;
	edge = model->edges + 1;
	trmEdge = trm->edges + 1;
	for ( i = 0; i < trm->numEdges; i++, edge++, trmEdge++ ) {
		edge->vertexNum[0] = trmEdge->v[0];
		edge->vertexNum[1] = trmEdge->v[1];
        VectorCopy( trmEdge->normal, edge->normal );
		edge->internal = qfalse;
		edge->sideSet = 0;
	}
	// polygons
	model->numPolygons = trm->numPolys;
	trmPoly = trm->polys;
	for ( i = 0; i < trm->numPolys; i++, trmPoly++ ) {
		poly = cmLocal.trmPolygons[i]->p;
		poly->numEdges = trmPoly->numEdges;
		for ( j = 0; j < trmPoly->numEdges; j++ ) {
			poly->edges[j] = trmPoly->edges[j];
		}
		VectorCopy( trmPoly->normal, poly->plane );
        PlaneSetDist( poly->plane, trmPoly->dist );
        VectorCopy( trmPoly->bounds[0], poly->bounds[0] );
        VectorCopy( trmPoly->bounds[1], poly->bounds[1] );
		poly->material = material;
		// link polygon at node
		cmLocal.trmPolygons[i]->next = model->node->polygons;
		model->node->polygons = cmLocal.trmPolygons[i];
	}
	// if the trace model is convex
	if ( trm->isConvex ) {
		// setup brush for position test
		cmLocal.trmBrushes[0]->b->numPlanes = trm->numPolys;
		for ( i = 0; i < trm->numPolys; i++ ) {
            PlaneCopy( cmLocal.trmPolygons[i]->p->plane, cmLocal.trmBrushes[0]->b->planes[i] );
		}
        VectorCopy( trm->bounds[0], cmLocal.trmBrushes[0]->b->bounds[0] );
        VectorCopy( trm->bounds[1], cmLocal.trmBrushes[0]->b->bounds[1] );
		// link brush at node
		cmLocal.trmBrushes[0]->next = model->node->brushes;
		model->node->brushes = cmLocal.trmBrushes[0];
	}
	// model bounds
    VectorCopy( trm->bounds[0], model->bounds[0] );
    VectorCopy( trm->bounds[1], model->bounds[1] );
	// convex
	model->isConvex = trm->isConvex;

	return TRACE_MODEL_HANDLE;
}

/*
===============================================================================

Optimisation, removal of polygons contained within brushes or solid

===============================================================================
*/

/*
============
R_ChoppedAwayByProcBSP
============
*/
int R_ChoppedAwayByProcBSP( int nodeNum, fixedWinding_t *w, const vec3_t normal, const vec3_t origin, const float radius ) {
	int res;
	fixedWinding_t back;
	cm_procNode_t *node;
	float dist;

	do {
		node = cmLocal.procNodes + nodeNum;
		dist = PlaneDistance( node->plane, origin );
		if ( dist > radius ) {
			res = SIDE_FRONT;
		}
		else if ( dist < -radius ) {
			res = SIDE_BACK;
		}
		else {
			res = SplitFixedWinding( w, &back, node->plane, CHOP_EPSILON );
		}
		if ( res == SIDE_FRONT ) {
			nodeNum = node->children[0];
		}
		else if ( res == SIDE_BACK ) {
			nodeNum = node->children[1];
		}
		else if ( res == SIDE_ON ) {
			// continue with the side the winding faces
			if ( DotProduct( node->plane, normal ) > 0.0f ) {
				nodeNum = node->children[0];
			}
			else {
				nodeNum = node->children[1];
			}
		}
		else {
			// if either node is not solid
			if ( node->children[0] < 0 || node->children[1] < 0 ) {
				return qfalse;
			}
			// only recurse if the node is not solid
			if ( node->children[1] > 0 ) {
				if ( !R_ChoppedAwayByProcBSP( node->children[1], &back, normal, origin, radius ) ) {
					return qfalse;
				}
			}
			nodeNum = node->children[0];
		}
	} while ( nodeNum > 0 );
	if ( nodeNum < 0 ) {
		return qfalse;
	}
	return qtrue;
}

/*
============
ChoppedAwayByProcBSP
============
*/
int ChoppedAwayByProcBSP( const fixedWinding_t *w, const plane_t plane, int contents ) {
	fixedWinding_t neww;
	vec3_t bounds[2];
	float radius;
	vec3_t origin, tmp;

	// if the .proc file has no BSP tree
	if ( cmLocal.procNodes == NULL ) {
		return qfalse;
	}
	// don't chop if the polygon is not solid
	if ( !(contents & CONTENTS_SOLID) ) {
		return qfalse;
	}
	// make a local copy of the winding
    CopyFixedWinding( w, &neww );
    FixedWindingBounds( &neww, bounds );
    VectorSubtract( bounds[1], bounds[0], tmp );
    VectorScale( tmp, 0.5f, origin );
	radius = VectorLength( origin ) + CHOP_EPSILON;
	VectorAdd( origin, bounds[0], origin );
	//
	return R_ChoppedAwayByProcBSP( 0, &neww, plane, origin, radius );
}

/*
=============
ChopWindingWithBrush

  returns the least number of winding fragments outside the brush
=============
*/
void ChopWindingListWithBrush( cm_windingList_t *list, cm_brush_t *b ) {
	int i, k, res, startPlane, planeNum, bestNumWindings;
	fixedWinding_t back, front;
	plane_t plane;
	qboolean chopped;
	int sidedness[MAX_POINTS_ON_WINDING];
	float dist;

	if ( b->numPlanes > MAX_POINTS_ON_WINDING ) {
		return;
	}

	// get sidedness for the list of windings
	for ( i = 0; i < b->numPlanes; i++ ) {
        PlaneNegate( b->planes[i], plane );

		dist = PlaneDistance( plane, list->origin );
		if ( dist > list->radius ) {
			sidedness[i] = SIDE_FRONT;
		}
		else if ( dist < -list->radius ) {
			sidedness[i] = SIDE_BACK;
		}
		else {
			sidedness[i] = BoxOnPlaneSideSlow( list->bounds[0], list->bounds[1], plane, ON_EPSILON );
			if ( sidedness[i] == PLANESIDE_FRONT ) {
				sidedness[i] = SIDE_FRONT;
			}
			else if ( sidedness[i] == PLANESIDE_BACK ) {
				sidedness[i] = SIDE_BACK;
			}
			else {
				sidedness[i] = SIDE_CROSS;
			}
		}
	}

	cm_outList->numWindings = 0;
	for ( k = 0; k < list->numWindings; k++ ) {
		//
		startPlane = 0;
		bestNumWindings = 1 + b->numPlanes;
		chopped = qfalse;
		do {
			front = list->w[k];
			cm_tmpList->numWindings = 0;
			for ( planeNum = startPlane, i = 0; i < b->numPlanes; i++, planeNum++ ) {

				if ( planeNum >= b->numPlanes ) {
					planeNum = 0;
				}

				res = sidedness[planeNum];

				if ( res == SIDE_CROSS ) {
                    PlaneNegate( b->planes[planeNum], plane );
					res = SplitFixedWinding( &front, &back, plane, CHOP_EPSILON );
				}

				// NOTE:	disabling this can create gaps at places where Z-fighting occurs
				//			Z-fighting should not occur but what if there is a decal brush side
				//			with exactly the same size as another brush side ?
				// only leave windings on a brush if the winding plane and brush side plane face the same direction
				if ( res == SIDE_ON && list->primitiveNum >= 0 && DotProduct( list->normal, b->planes[planeNum] ) > 0 ) {
					// return because all windings in the list will be on this brush side plane
					return;
				}

				if ( res == SIDE_BACK ) {
					if ( cm_outList->numWindings >= MAX_WINDING_LIST ) {
						ii.Com_DPrintf( "ChopWindingWithBrush: primitive %d more than %d windings", list->primitiveNum, MAX_WINDING_LIST );
						return;
					}
					// winding and brush didn't intersect, store the original winding
					cm_outList->w[cm_outList->numWindings] = list->w[k];
					cm_outList->numWindings++;
					chopped = qfalse;
					break;
				}

				if ( res == SIDE_CROSS ) {
					if ( cm_tmpList->numWindings >= MAX_WINDING_LIST ) {
						ii.Com_Printf( "ChopWindingWithBrush: primitive %d more than %d windings", list->primitiveNum, MAX_WINDING_LIST );
						return;
					}
					// store the front winding in the temporary list
					cm_tmpList->w[cm_tmpList->numWindings] = back;
					cm_tmpList->numWindings++;
					chopped = qtrue;
				}

				// if already found a start plane which generates less fragments
				if ( cm_tmpList->numWindings >= bestNumWindings ) {
					break;
				}
			}

			// find the best start plane to get the least number of fragments outside the brush
			if ( cm_tmpList->numWindings < bestNumWindings ) {
				bestNumWindings = cm_tmpList->numWindings;
				// store windings from temporary list in the out list
				for ( i = 0; i < cm_tmpList->numWindings; i++ ) {
					if ( cm_outList->numWindings + i >= MAX_WINDING_LIST ) {
						ii.Com_DPrintf( "ChopWindingWithBrush: primitive %d more than %d windings", list->primitiveNum, MAX_WINDING_LIST );
						return;
					}
					cm_outList->w[cm_outList->numWindings+i] = cm_tmpList->w[i];
				}
				// if only one winding left then we can't do any better
				if ( bestNumWindings == 1 ) {
					break;
				}
			}

			// try the next start plane
			startPlane++;

		} while ( chopped && startPlane < b->numPlanes );
		//
		if ( chopped ) {
			cm_outList->numWindings += bestNumWindings;
		}
	}
	for ( k = 0; k < cm_outList->numWindings; k++ ) {
		list->w[k] = cm_outList->w[k];
	}
	list->numWindings = cm_outList->numWindings;
}

/*
============
R_ChopWindingListWithTreeBrushes
============
*/
void R_ChopWindingListWithTreeBrushes( cm_windingList_t *list, cm_node_t *node ) {
	int i;
	cm_brushRef_t *bref;
	cm_brush_t *b;

	while( 1 ) {
		for ( bref = node->brushes; bref; bref = bref->next ) {
			b = bref->b;
			// if we checked this brush already
			if ( b->checkcount == cmLocal.checkCount ) {
				continue;
			}
			b->checkcount = cmLocal.checkCount;
			// if the windings in the list originate from this brush
			if ( b->primitiveNum == list->primitiveNum ) {
				continue;
			}
			// if brush has a different contents
			if ( b->contents != list->contents ) {
				continue;
			}
			// brush bounds and winding list bounds should overlap
			for ( i = 0; i < 3; i++ ) {
				if ( list->bounds[0][i] > b->bounds[1][i] ) {
					break;
				}
				if ( list->bounds[1][i] < b->bounds[0][i] ) {
					break;
				}
			}
			if ( i < 3 ) {
				continue;
			}
			// chop windings in the list with brush
			ChopWindingListWithBrush( list, b );
			// if all windings are chopped away we're done
			if ( !list->numWindings ) {
				return;
			}
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		if ( list->bounds[0][node->planeType] > node->planeDist ) {
			node = node->children[0];
		}
		else if ( list->bounds[1][node->planeType] < node->planeDist ) {
			node = node->children[1];
		}
		else {
			R_ChopWindingListWithTreeBrushes( list, node->children[1] );
			if ( !list->numWindings ) {
				return;
			}
			node = node->children[0];
		}
	}
}

/*
============
WindingOutsideBrushes

  Returns one winding which is not fully contained in brushes.
  We always favor less polygons over a stitched world.
  If the winding is partly contained and the contained pieces can be chopped off
  without creating multiple winding fragments then the chopped winding is returned.
============
*/
fixedWinding_t *WindingOutsideBrushes( fixedWinding_t *w, const plane_t plane, int contents, int primitiveNum, cm_node_t *headNode ) {
	int i, windingLeft;
    vec3_t tmp;

	ClearBounds( cm_windingList->bounds[0], cm_windingList->bounds[1] );
	for ( i = 0; i < w->numPoints; i++ ) {
		AddPointToBounds( w->points[i], cm_windingList->bounds[0], cm_windingList->bounds[1] );
	}

    VectorSubtract( cm_windingList->bounds[1], cm_windingList->bounds[0], tmp );
    VectorScale( tmp, 0.5f, cm_windingList->origin );
	cm_windingList->radius = VectorLength( cm_windingList->origin ) + CHOP_EPSILON;
    VectorAdd( cm_windingList->origin, cm_windingList->bounds[0], cm_windingList->origin );
    for ( i = 0; i < 3; i++ ) {
	    cm_windingList->bounds[0][i] -= CHOP_EPSILON;
	    cm_windingList->bounds[1][i] += CHOP_EPSILON;
    }

    CopyFixedWinding( w, &cm_windingList->w[0] );
	cm_windingList->numWindings = 1;
    VectorCopy( plane, cm_windingList->normal );
	cm_windingList->contents = contents;
	cm_windingList->primitiveNum = primitiveNum;
	//
	cmLocal.checkCount++;
	R_ChopWindingListWithTreeBrushes( cm_windingList, headNode );
	//
	if ( !cm_windingList->numWindings ) {
		return NULL;
	}
	if ( cm_windingList->numWindings == 1 ) {
		return &cm_windingList->w[0];
	}
	// if not the world model
	if ( cmLocal.numModels != 0 ) {
		return w;
	}
	// check if winding fragments would be chopped away by the proc BSP tree
	windingLeft = -1;
	for ( i = 0; i < cm_windingList->numWindings; i++ ) {
		if ( !ChoppedAwayByProcBSP( &cm_windingList->w[i], plane, contents ) ) {
			if ( windingLeft >= 0 ) {
				return w;
			}
			windingLeft = i;
		}
	}
	if ( windingLeft >= 0 ) {
		return &cm_windingList->w[windingLeft];
	}
	return NULL;
}

/*
===============================================================================

Merging polygons

===============================================================================
*/

/*
=============
CM_ReplacePolygons

  does not allow for a node to have multiple references to the same polygon
=============
*/
void CM_ReplacePolygons( cm_model_t *model, cm_node_t *node, cm_polygon_t *p1, cm_polygon_t *p2, cm_polygon_t *newp ) {
	cm_polygonRef_t *pref, *lastpref, *nextpref;
	cm_polygon_t *p;
	qboolean linked;

	while( 1 ) {
		linked = qfalse;
		lastpref = NULL;
		for ( pref = node->polygons; pref; pref = nextpref ) {
			nextpref = pref->next;
			//
			p = pref->p;
			// if this polygon reference should change
			if ( p == p1 || p == p2 ) {
				// if the new polygon is already linked at this node
				if ( linked ) {
					if ( lastpref ) {
						lastpref->next = nextpref;
					}
					else {
						node->polygons = nextpref;
					}
					CM_FreePolygonReference( pref );
					model->numPolygonRefs--;
				}
				else {
					pref->p = newp;
					linked = qtrue;
					lastpref = pref;
				}
			}
			else {
				lastpref = pref;
			}
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		if ( p1->bounds[0][node->planeType] > node->planeDist && p2->bounds[0][node->planeType] > node->planeDist ) {
			node = node->children[0];
		}
		else if ( p1->bounds[1][node->planeType] < node->planeDist && p2->bounds[1][node->planeType] < node->planeDist ) {
			node = node->children[1];
		}
		else {
			CM_ReplacePolygons( model, node->children[1], p1, p2, newp );
			node = node->children[0];
		}
	}
}

/*
=============
CM_TryMergePolygons
=============
*/
#define	CONTINUOUS_EPSILON	0.005f
#define NORMAL_EPSILON		0.01f

cm_polygon_t *CM_TryMergePolygons( cm_model_t *model, cm_polygon_t *p1, cm_polygon_t *p2 ) {
	int i, j, nexti, prevj;
	int p1BeforeShare, p1AfterShare, p2BeforeShare, p2AfterShare;
	int newEdges[CM_MAX_POLYGON_EDGES], newNumEdges;
	int edgeNum, edgeNum1, edgeNum2, newEdgeNum1, newEdgeNum2;
	cm_edge_t *edge;
	cm_polygon_t *newp;
	vec3_t delta, normal;
	float dot;
	qboolean keep1, keep2;

	if ( p1->material != p2->material ) {
		return NULL;
	}
	if ( fabs( PlaneGetDist( p1->plane ) - PlaneGetDist( p2->plane ) ) > NORMAL_EPSILON ) {
		return NULL;
	}
	for ( i = 0; i < 3; i++ ) {
		if ( fabs( p1->plane[i] - p2->plane[i] ) > NORMAL_EPSILON ) {
			return NULL;
		}
		if ( p1->bounds[0][i] > p2->bounds[1][i] ) {
			return NULL;
		}
		if ( p1->bounds[1][i] < p2->bounds[0][i] ) {
			return NULL;
		}
	}
	// this allows for merging polygons with multiple shared edges
	// polygons with multiple shared edges probably never occur tho ;)
	p1BeforeShare = p1AfterShare = p2BeforeShare = p2AfterShare = -1;
	for ( i = 0; i < p1->numEdges; i++ ) {
		nexti = (i+1)%p1->numEdges;
		for ( j = 0; j < p2->numEdges; j++ ) {
			prevj = (j+p2->numEdges-1)%p2->numEdges;
			//
			if ( abs(p1->edges[i]) != abs(p2->edges[j]) ) {
				// if the next edge of p1 and the previous edge of p2 are the same
				if ( abs(p1->edges[nexti]) == abs(p2->edges[prevj]) ) {
					// if both polygons don't use the edge in the same direction
					if ( p1->edges[nexti] != p2->edges[prevj] ) {
						p1BeforeShare = i;
						p2AfterShare = j;
					}
					break;
				}
			}
			// if both polygons don't use the edge in the same direction
			else if ( p1->edges[i] != p2->edges[j] ) {
				// if the next edge of p1 and the previous edge of p2 are not the same
				if ( abs(p1->edges[nexti]) != abs(p2->edges[prevj]) ) {
					p1AfterShare = nexti;
					p2BeforeShare = prevj;
					break;
				}
			}
		}
	}
	if ( p1BeforeShare < 0 || p1AfterShare < 0 || p2BeforeShare < 0 || p2AfterShare < 0 ) {
		return NULL;
	}

	// check if the new polygon would still be convex
	edgeNum = p1->edges[p1BeforeShare];
	edge = model->edges + abs(edgeNum);
    VectorSubtract(
            model->vertices[edge->vertexNum[INTSIGNBITNOTSET(edgeNum)]].p,
            model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p,
            delta );
    CrossProduct( p1->plane, delta, normal );
	VectorNormalize( normal );

	edgeNum = p2->edges[p2AfterShare];
	edge = model->edges + abs(edgeNum);
    VectorSubtract(
            model->vertices[edge->vertexNum[INTSIGNBITNOTSET(edgeNum)]].p,
            model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p,
            delta );

	dot = DotProduct( delta, normal );
	if (dot < -CONTINUOUS_EPSILON)
		return NULL;			// not a convex polygon
	keep1 = (qboolean)(dot > CONTINUOUS_EPSILON);

	edgeNum = p2->edges[p2BeforeShare];
	edge = model->edges + abs(edgeNum);
    VectorSubtract(
            model->vertices[edge->vertexNum[INTSIGNBITNOTSET(edgeNum)]].p,
            model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p,
            delta );
    CrossProduct( p1->plane, delta, normal );
	VectorNormalize( normal );

	edgeNum = p1->edges[p1AfterShare];
	edge = model->edges + abs(edgeNum);
    VectorSubtract(
            model->vertices[edge->vertexNum[INTSIGNBITNOTSET(edgeNum)]].p,
            model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p,
            delta );

	dot = DotProduct( delta, normal );
	if (dot < -CONTINUOUS_EPSILON)
		return NULL;			// not a convex polygon
	keep2 = (qboolean)(dot > CONTINUOUS_EPSILON);

	newEdgeNum1 = newEdgeNum2 = 0;
	// get new edges if we need to replace colinear ones
	if ( !keep1 ) {
		edgeNum1 = p1->edges[p1BeforeShare];
		edgeNum2 = p2->edges[p2AfterShare];
		GetEdge( model, model->vertices[model->edges[abs(edgeNum1)].vertexNum[INTSIGNBITSET(edgeNum1)]].p,
					model->vertices[model->edges[abs(edgeNum2)].vertexNum[INTSIGNBITNOTSET(edgeNum2)]].p,
						&newEdgeNum1, -1 );
		if ( newEdgeNum1 == 0 ) {
			keep1 = qtrue;
		}
	}
	if ( !keep2 ) {
		edgeNum1 = p2->edges[p2BeforeShare];
		edgeNum2 = p1->edges[p1AfterShare];
		GetEdge( model, model->vertices[model->edges[abs(edgeNum1)].vertexNum[INTSIGNBITSET(edgeNum1)]].p,
					model->vertices[model->edges[abs(edgeNum2)].vertexNum[INTSIGNBITNOTSET(edgeNum2)]].p,
						&newEdgeNum2, -1 );
		if ( newEdgeNum2 == 0 ) {
			keep2 = qtrue;
		}
	}
	// set edges for new polygon
	newNumEdges = 0;
	if ( !keep2 ) {
		newEdges[newNumEdges++] = newEdgeNum2;
	}
	if ( p1AfterShare < p1BeforeShare ) {
		for ( i = p1AfterShare + (!keep2); i <= p1BeforeShare - (!keep1); i++ ) {
			newEdges[newNumEdges++] = p1->edges[i];
		}
	}
	else {
		for ( i = p1AfterShare + (!keep2); i < p1->numEdges; i++ ) {
			newEdges[newNumEdges++] = p1->edges[i];
		}
		for ( i = 0; i <= p1BeforeShare - (!keep1); i++ ) {
			newEdges[newNumEdges++] = p1->edges[i];
		}
	}
	if ( !keep1 ) {
		newEdges[newNumEdges++] = newEdgeNum1;
	}
	if ( p2AfterShare < p2BeforeShare ) {
		for ( i = p2AfterShare + (!keep1); i <= p2BeforeShare - (!keep2); i++ ) {
			newEdges[newNumEdges++] = p2->edges[i];
		}
	}
	else {
		for ( i = p2AfterShare + (!keep1); i < p2->numEdges; i++ ) {
			newEdges[newNumEdges++] = p2->edges[i];
		}
		for ( i = 0; i <= p2BeforeShare - (!keep2); i++ ) {
			newEdges[newNumEdges++] = p2->edges[i];
		}
	}

	newp = AllocPolygon( model, newNumEdges );
	memcpy( newp, p1, sizeof(cm_polygon_t) );
	memcpy( newp->edges, newEdges, newNumEdges * sizeof(int) );
	newp->numEdges = newNumEdges;
	newp->checkcount = 0;
	// increase usage count for the edges of this polygon
	for ( i = 0; i < newp->numEdges; i++ ) {
		if ( !keep1 && newp->edges[i] == newEdgeNum1 ) {
			continue;
		}
		if ( !keep2 && newp->edges[i] == newEdgeNum2 ) {
			continue;
		}
		model->edges[abs(newp->edges[i])].numUsers++;
	}
	// create new bounds from the merged polygons
    ClearBounds( newp->bounds[0], newp->bounds[1] );
    AddBoundsToBounds( p1->bounds, newp->bounds );
    AddBoundsToBounds( p2->bounds, newp->bounds );

	return newp;
}

/*
=============
CM_MergePolygonWithTreePolygons
=============
*/
qboolean CM_MergePolygonWithTreePolygons( cm_model_t *model, cm_node_t *node, cm_polygon_t *polygon ) {
	int i;
	cm_polygonRef_t *pref;
	cm_polygon_t *p, *newp;

	while( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			p = pref->p;
			//
			if ( p == polygon ) {
				continue;
			}
			//
			newp = CM_TryMergePolygons( model, polygon, p );
			// if polygons were merged
			if ( newp ) {
				model->numMergedPolys++;
				// replace links to the merged polygons with links to the new polygon
				CM_ReplacePolygons( model, model->node, polygon, p, newp );
				// decrease usage count for edges of both merged polygons
				for ( i = 0; i < polygon->numEdges; i++ ) {
					model->edges[abs(polygon->edges[i])].numUsers--;
				}
				for ( i = 0; i < p->numEdges; i++ ) {
					model->edges[abs(p->edges[i])].numUsers--;
				}
				// free merged polygons
				CM_FreePolygon( model, polygon );
				CM_FreePolygon( model, p );

				return qtrue;
			}
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		if ( polygon->bounds[0][node->planeType] > node->planeDist ) {
			node = node->children[0];
		}
		else if ( polygon->bounds[1][node->planeType] < node->planeDist ) {
			node = node->children[1];
		}
		else {
			if ( CM_MergePolygonWithTreePolygons( model, node->children[1], polygon ) ) {
				return qtrue;
			}
			node = node->children[0];
		}
	}
	return qfalse;
}

/*
=============
CM_MergeTreePolygons

  try to merge any two polygons with the same surface flags and the same contents
=============
*/
void CM_MergeTreePolygons( cm_model_t *model, cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;
	qboolean merge;

	while( 1 ) {
		do {
			merge = qfalse;
			for ( pref = node->polygons; pref; pref = pref->next ) {
				p = pref->p;
				// if we checked this polygon already
				if ( p->checkcount == cmLocal.checkCount ) {
					continue;
				}
				p->checkcount = cmLocal.checkCount;
				// try to merge this polygon with other polygons in the tree
				if ( CM_MergePolygonWithTreePolygons( model, model->node, p ) ) {
					merge = qtrue;
					break;
				}
			}
		} while (merge);
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		CM_MergeTreePolygons( model, node->children[1] );
		node = node->children[0];
	}
}

/*
===============================================================================

Find internal edges

===============================================================================
*/

/*

	if (two polygons have the same contents)
		if (the normals of the two polygon planes face towards each other)
			if (an edge is shared between the polygons)
				if (the edge is not shared in the same direction)
					then this is an internal edge
			else
				if (this edge is on the plane of the other polygon)
					if (this edge if fully inside the winding of the other polygon)
						then this edge is an internal edge

*/

/*
=============
CM_PointInsidePolygon
=============
*/
qboolean CM_PointInsidePolygon( cm_model_t *model, cm_polygon_t *p, vec3_t v ) {
	int i, edgeNum;
	vec3_t *v1, *v2, dir1, dir2, vec;
	cm_edge_t *edge;

	for ( i = 0; i < p->numEdges; i++ ) {
		edgeNum = p->edges[i];
		edge = model->edges + abs(edgeNum);
		//
		v1 = &model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p;
		v2 = &model->vertices[edge->vertexNum[INTSIGNBITNOTSET(edgeNum)]].p;
        VectorSubtract( *v2, *v1, dir1 );
        VectorSubtract( v, *v1, vec );
        CrossProduct( dir1, p->plane, dir2 );
		if ( DotProduct( vec, dir2 ) > VERTEX_EPSILON ) {
			return qfalse;
		}
	}
	return qtrue;
}

/*
=============
CM_FindInternalEdgesOnPolygon
=============
*/
void CM_FindInternalEdgesOnPolygon( cm_model_t *model, cm_polygon_t *p1, cm_polygon_t *p2 ) {
	int i, j, k, edgeNum;
	cm_edge_t *edge;
	vec3_t *v1, *v2, dir1, dir2;
	float d;

	// bounds of polygons should overlap or touch
	for ( i = 0; i < 3; i++ ) {
		if ( p1->bounds[0][i] > p2->bounds[1][i] ) {
			return;
		}
		if ( p1->bounds[1][i] < p2->bounds[0][i] ) {
			return;
		}
	}
	//
	// FIXME: doubled geometry causes problems
	//
	for ( i = 0; i < p1->numEdges; i++ ) {
		edgeNum = p1->edges[i];
		edge = model->edges + abs(edgeNum);
		// if already an internal edge
		if ( edge->internal ) {
			continue;
		}
		//
		v1 = &model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p;
		v2 = &model->vertices[edge->vertexNum[INTSIGNBITNOTSET(edgeNum)]].p;
		// if either of the two vertices is outside the bounds of the other polygon
		for ( k = 0; k < 3; k++ ) {
			d = p2->bounds[1][k] + VERTEX_EPSILON;
			if ( (*v1)[k] > d || (*v2)[k] > d ) {
				break;
			}
			d = p2->bounds[0][k] - VERTEX_EPSILON;
			if ( (*v1)[k] < d || (*v2)[k] < d ) {
				break;
			}
		}
		if ( k < 3 ) {
			continue;
		}
		//
		k = abs(edgeNum);
		for ( j = 0; j < p2->numEdges; j++ ) {
			if ( k == abs(p2->edges[j]) ) {
				break;
			}
		}
		// if the edge is shared between the two polygons
		if ( j < p2->numEdges ) {
			// if the edge is used by more than 2 polygons
			if ( edge->numUsers > 2 ) {
				// could still be internal but we'd have to test all polygons using the edge
				continue;
			}
			// if the edge goes in the same direction for both polygons
			if ( edgeNum == p2->edges[j] ) {
				// the polygons can lay ontop of each other or one can obscure the other
				continue;
			}
		}
		// the edge was not shared
		else {
			// both vertices should be on the plane of the other polygon
			d = PlaneDistance( p2->plane, *v1 );
			if ( fabs(d) > VERTEX_EPSILON ) {
				continue;
			}
			d = PlaneDistance( p2->plane, *v2 );
			if ( fabs(d) > VERTEX_EPSILON ) {
				continue;
			}
		}
		// the two polygon plane normals should face towards each other
        VectorSubtract( *v2, *v1, dir1 );
        VectorSubtract( p1->plane, dir1, dir2 );
		if ( DotProduct( p2->plane, dir2 ) < 0 ) {
			//continue;
			break;
		}
		// if the edge was not shared
		if ( j >= p2->numEdges ) {
			// both vertices of the edge should be inside the winding of the other polygon
			if ( !CM_PointInsidePolygon( model, p2, *v1 ) ) {
				continue;
			}
			if ( !CM_PointInsidePolygon( model, p2, *v2 ) ) {
				continue;
			}
		}
		// we got another internal edge
		edge->internal = qtrue;
		model->numInternalEdges++;
	}
}

/*
=============
CM_FindInternalPolygonEdges
=============
*/
void CM_FindInternalPolygonEdges( cm_model_t *model, cm_node_t *node, cm_polygon_t *polygon ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;

	if ( CM_MaterialNeedsBackSide( polygon->material ) ) {
		return;
	}

	while( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			p = pref->p;
			//
			// FIXME: use some sort of additional checkcount because currently
			//			polygons can be checked multiple times
			//
			// if the polygons don't have the same contents
			if ( p->contents != polygon->contents ) {
				continue;
			}
			if ( p == polygon ) {
				continue;
			}
			CM_FindInternalEdgesOnPolygon( model, polygon, p );
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		if ( polygon->bounds[0][node->planeType] > node->planeDist ) {
			node = node->children[0];
		}
		else if ( polygon->bounds[1][node->planeType] < node->planeDist ) {
			node = node->children[1];
		}
		else {
			CM_FindInternalPolygonEdges( model, node->children[1], polygon );
			node = node->children[0];
		}
	}
}

/*
=============
CM_FindContainedEdges
=============
*/
void CM_FindContainedEdges( cm_model_t *model, cm_polygon_t *p ) {
	int i, edgeNum;
	cm_edge_t *edge;
	fixedWinding_t w;
	vec5_t tmp5;

	tmp5[3] = tmp5[4] = 0.0f;
	for ( i = 0; i < p->numEdges; i++ ) {
		edgeNum = p->edges[i];
		edge = model->edges + abs(edgeNum);
		if ( edge->internal ) {
			continue;
		}
        ClearFixedWinding( &w );
		VectorCopy( model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p, tmp5 );
        AddPointToFixedWinding( &w, tmp5 );
		VectorCopy( model->vertices[edge->vertexNum[INTSIGNBITNOTSET(edgeNum)]].p, tmp5 );
        AddPointToFixedWinding( &w, tmp5 );
		if ( ChoppedAwayByProcBSP( &w, p->plane, p->contents ) ) {
			edge->internal = qtrue;
		}
	}
}

/*
=============
CM_FindInternalEdges
=============
*/
void CM_FindInternalEdges( cm_model_t *model, cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;

	while( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			p = pref->p;
			// if we checked this polygon already
			if ( p->checkcount == cmLocal.checkCount ) {
				continue;
			}
			p->checkcount = cmLocal.checkCount;

			CM_FindInternalPolygonEdges( model, model->node, p );

			//CM_FindContainedEdges( model, p );
		}
		// if leaf node
		if ( node->planeType == -1 ) {
			break;
		}
		CM_FindInternalEdges( model, node->children[1] );
		node = node->children[0];
	}
}

/*
===============================================================================

Spatial subdivision

===============================================================================
*/

/*
================
CM_FindSplitter
================
*/
static int CM_FindSplitter( const cm_node_t *node, const vec3_t bounds[2], int *planeType, float *planeDist ) {
	int i, j, type, axis[3], polyCount;
	float dist, t, bestt, size[3];
	cm_brushRef_t *bref;
	cm_polygonRef_t *pref;
	const cm_node_t *n;
	qboolean forceSplit = qfalse;

	for ( i = 0; i < 3; i++ ) {
		size[i] = bounds[1][i] - bounds[0][i];
		axis[i] = i;
	}
	// sort on largest axis
	for ( i = 0; i < 2; i++ ) {
		if ( size[i] < size[i+1] ) {
			t = size[i];
			size[i] = size[i+1];
			size[i+1] = t;
			j = axis[i];
			axis[i] = axis[i+1];
			axis[i+1] = j;
			i = -1;
		}
	}
	// if the node is too small for further splits
	if ( size[0] < MIN_NODE_SIZE ) {
		polyCount = 0;
		for ( pref = node->polygons; pref; pref = pref->next) {
			polyCount++;
		}
		if ( polyCount > MAX_NODE_POLYGONS ) {
			forceSplit = qtrue;
		}
	}
	// find an axial aligned splitter
	for ( i = 0; i < 3; i++ ) {
		// start with the largest axis first
		type = axis[i];
		bestt = size[i];
		// if the node is small anough in this axis direction
		if ( !forceSplit && bestt < MIN_NODE_SIZE ) {
			break;
		}
		// find an axial splitter from the brush bounding boxes
		// also try brushes from parent nodes
		for ( n = node; n; n = n->parent ) {
			for ( bref = n->brushes; bref; bref = bref->next) {
				for ( j = 0; j < 2; j++ ) {
					dist = bref->b->bounds[j][type];
					// if the splitter is already used or outside node bounds
					if ( dist >= bounds[1][type] || dist <= bounds[0][type] ) {
						continue;
					}
					// find the most centered splitter
					t = abs((bounds[1][type] - dist) - (dist - bounds[0][type]));
					if ( t < bestt ) {
						bestt = t;
						*planeType = type;
						*planeDist = dist;
					}
				}
			}
		}
		// find an axial splitter from the polygon bounding boxes
		// also try brushes from parent nodes
		for ( n = node; n; n = n->parent ) {
			for ( pref = n->polygons; pref; pref = pref->next) {
				for ( j = 0; j < 2; j++ ) {
					dist = pref->p->bounds[j][type];
					// if the splitter is already used or outside node bounds
					if ( dist >= bounds[1][type] || dist <= bounds[0][type] ) {
						continue;
					}
					// find the most centered splitter
					t = abs((bounds[1][type] - dist) - (dist - bounds[0][type]));
					if ( t < bestt ) {
						bestt = t;
						*planeType = type;
						*planeDist = dist;
					}
				}
			}
		}
		// if we found a splitter on the largest axis
		if ( bestt < size[i] ) {
			// if forced split due to lots of polygons
			if ( forceSplit ) {
				return qtrue;
			}
			// don't create splitters real close to the bounds
			if ( bounds[1][type] - *planeDist > (MIN_NODE_SIZE*0.5f) &&
				*planeDist - bounds[0][type] > (MIN_NODE_SIZE*0.5f) ) {
				return qtrue;
			}
		}
	}
	return qfalse;
}

/*
================
CM_R_InsideAllChildren
================
*/
static int CM_R_InsideAllChildren( cm_node_t *node, const vec3_t bounds[2] ) {
	assert(node != NULL);
	if ( node->planeType != -1 ) {
		if ( bounds[0][node->planeType] >= node->planeDist ) {
			return qfalse;
		}
		if ( bounds[1][node->planeType] <= node->planeDist ) {
			return qfalse;
		}
		if ( !CM_R_InsideAllChildren( node->children[0], bounds ) ) {
			return qfalse;
		}
		if ( !CM_R_InsideAllChildren( node->children[1], bounds ) ) {
			return qfalse;
		}
	}
	return qtrue;
}

/*
================
R_FilterPolygonIntoTree
================
*/
void R_FilterPolygonIntoTree( cm_model_t *model, cm_node_t *node, cm_polygonRef_t *pref, cm_polygon_t *p ) {
	assert(node != NULL);
	while ( node->planeType != -1 ) {
		if ( CM_R_InsideAllChildren( node, p->bounds ) ) {
			break;
		}
		if ( p->bounds[0][node->planeType] >= node->planeDist ) {
			node = node->children[0];
		}
		else if ( p->bounds[1][node->planeType] <= node->planeDist ) {
			node = node->children[1];
		}
		else {
			R_FilterPolygonIntoTree( model, node->children[1], NULL, p );
			node = node->children[0];
		}
	}
	if ( pref ) {
		pref->next = node->polygons;
		node->polygons = pref;
	}
	else {
		AddPolygonToNode( model, node, p );
	}
}

/*
================
R_FilterBrushIntoTree
================
*/
void R_FilterBrushIntoTree( cm_model_t *model, cm_node_t *node, cm_brushRef_t *pref, cm_brush_t *b ) {
	assert(node != NULL);
	while ( node->planeType != -1 ) {
		if ( CM_R_InsideAllChildren( node, b->bounds ) ) {
			break;
		}
		if ( b->bounds[0][node->planeType] >= node->planeDist ) {
			node = node->children[0];
		}
		else if ( b->bounds[1][node->planeType] <= node->planeDist ) {
			node = node->children[1];
		}
		else {
			R_FilterBrushIntoTree( model, node->children[1], NULL, b );
			node = node->children[0];
		}
	}
	if ( pref ) {
		pref->next = node->brushes;
		node->brushes = pref;
	}
	else {
		AddBrushToNode( model, node, b );
	}
}

/*
================
R_CreateAxialBSPTree

  a brush or polygon is linked in the node closest to the root where
  the brush or polygon is inside all children
================
*/
cm_node_t *R_CreateAxialBSPTree( cm_model_t *model, cm_node_t *node, const vec3_t bounds[2] ) {
	int planeType;
	float planeDist;
	cm_polygonRef_t *pref, *nextpref, *prevpref;
	cm_brushRef_t *bref, *nextbref, *prevbref;
	cm_node_t *frontNode, *backNode, *n;
	vec3_t frontBounds[2], backBounds[2];

	if ( !CM_FindSplitter( node, bounds, &planeType, &planeDist ) ) {
		node->planeType = -1;
		return node;
	}
	// create two child nodes
	frontNode = AllocNode( model, NODE_BLOCK_SIZE_LARGE );
	memset( frontNode, 0, sizeof(cm_node_t) );
	frontNode->parent = node;
	frontNode->planeType = -1;
	//
	backNode = AllocNode( model, NODE_BLOCK_SIZE_LARGE );
	memset( backNode, 0, sizeof(cm_node_t) );
	backNode->parent = node;
	backNode->planeType = -1;
	//
	model->numNodes += 2;
	// set front node bounds
    VectorCopy( bounds[0], frontBounds[0] );
    VectorCopy( bounds[1], frontBounds[1] );
	frontBounds[0][planeType] = planeDist;
	// set back node bounds
    VectorCopy( bounds[0], backBounds[0] );
    VectorCopy( bounds[1], backBounds[1] );
	backBounds[1][planeType] = planeDist;
	//
	node->planeType = planeType;
	node->planeDist = planeDist;
	node->children[0] = frontNode;
	node->children[1] = backNode;
	// filter polygons and brushes down the tree if necesary
	for ( n = node; n; n = n->parent ) {
		prevpref = NULL;
		for ( pref = n->polygons; pref; pref = nextpref) {
			nextpref = pref->next;
			// if polygon is not inside all children
			if ( !CM_R_InsideAllChildren( n, pref->p->bounds ) ) {
				// filter polygon down the tree
				R_FilterPolygonIntoTree( model, n, pref, pref->p );
				if ( prevpref ) {
					prevpref->next = nextpref;
				}
				else {
					n->polygons = nextpref;
				}
			}
			else {
				prevpref = pref;
			}
		}
		prevbref = NULL;
		for ( bref = n->brushes; bref; bref = nextbref) {
			nextbref = bref->next;
			// if brush is not inside all children
			if ( !CM_R_InsideAllChildren( n, bref->b->bounds ) ) {
				// filter brush down the tree
				R_FilterBrushIntoTree( model, n, bref, bref->b );
				if ( prevbref ) {
					prevbref->next = nextbref;
				}
				else {
					n->brushes = nextbref;
				}
			}
			else {
				prevbref = bref;
			}
		}
	}
	R_CreateAxialBSPTree( model, frontNode, frontBounds );
	R_CreateAxialBSPTree( model, backNode, backBounds );
	return node;
}

/*
int cm_numSavedPolygonLinks;
int cm_numSavedBrushLinks;

int CM_R_CountChildren( cm_node_t *node ) {
	if ( node->planeType == -1 ) {
		return 0;
	}
	return 2 + CM_R_CountChildren(node->children[0]) + CM_R_CountChildren(node->children[1]);
}

void CM_R_TestOptimisation( cm_node_t *node ) {
	int polyCount, brushCount, numChildren;
	cm_polygonRef_t *pref;
	cm_brushRef_t *bref;

	if ( node->planeType == -1 ) {
		return;
	}
	polyCount = 0;
	for ( pref = node->polygons; pref; pref = pref->next) {
		polyCount++;
	}
	brushCount = 0;
	for ( bref = node->brushes; bref; bref = bref->next) {
		brushCount++;
	}
	if ( polyCount || brushCount ) {
		numChildren = CM_R_CountChildren( node );
		cm_numSavedPolygonLinks += (numChildren - 1) * polyCount;
		cm_numSavedBrushLinks += (numChildren - 1) * brushCount;
	}
	CM_R_TestOptimisation( node->children[0] );
	CM_R_TestOptimisation( node->children[1] );
}
*/

/*
================
CreateAxialBSPTree
================
*/
cm_node_t *CreateAxialBSPTree( cm_model_t *model, cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_brushRef_t *bref;
	vec3_t bounds[2];

	// get head node bounds
    ClearBounds( bounds[0], bounds[1] );
	for ( pref = node->polygons; pref; pref = pref->next) {
        AddBoundsToBounds( pref->p->bounds, bounds );
	}
	for ( bref = node->brushes; bref; bref = bref->next) {
        AddBoundsToBounds( bref->b->bounds, bounds );
	}

	// create axial BSP tree from head node
	node = R_CreateAxialBSPTree( model, node, bounds );

	return node;
}

/*
===============================================================================

Raw polygon and brush data

===============================================================================
*/

/*
================
SetupHash
================
*/
void SetupHash( void ) {
	if ( !cm_vertexHash ) {
        cm_vertexHash = ii.GetMemory( sizeof( *cm_vertexHash ) );
        HashIndexInit( cm_vertexHash, VERTEX_HASH_SIZE, 1024 );
	}
    if ( !cm_edgeHash ) {
        cm_edgeHash = ii.GetMemory( sizeof( *cm_edgeHash ) );
		HashIndexInit( cm_edgeHash, EDGE_HASH_SIZE, 1024 );
    }
	// init variables used during loading and optimization
	if ( !cm_windingList ) {
		cm_windingList = ii.GetMemory( sizeof( *cm_windingList ) );
        memset( cm_windingList, 0, sizeof( *cm_windingList ) );
	}
	if ( !cm_outList ) {
        cm_outList = ii.GetMemory( sizeof( *cm_outList ) );
        memset( cm_outList, 0, sizeof( *cm_outList ) );
	}
	if ( !cm_tmpList ) {
        cm_tmpList = ii.GetMemory( sizeof( *cm_tmpList ) );
        memset( cm_tmpList, 0, sizeof( *cm_tmpList ) );
 	}
}

/*
================
ShutdownHash
================
*/
void ShutdownHash( void ) {
    if ( cm_vertexHash ) {
        HashIndexFree( cm_vertexHash );
        ii.FreeMemory( cm_vertexHash );
        cm_vertexHash = NULL;
    }
    if ( cm_edgeHash ) {
        HashIndexFree( cm_edgeHash );
        ii.FreeMemory( cm_edgeHash );
        cm_edgeHash = NULL;
    }
    if ( cm_tmpList ) {
        ii.FreeMemory( cm_tmpList );
	    cm_tmpList = NULL;
    }
    if ( cm_outList ) {
        ii.FreeMemory( cm_outList );
	    cm_outList = NULL;
    }
    if ( cm_windingList ) {
        ii.FreeMemory( cm_windingList );
	    cm_windingList = NULL;
    }
}

/*
================
ClearHash
================
*/
void ClearHash( vec3_t bounds[2] ) {
	int i;
	float f, max;

	HashIndexClear( cm_vertexHash );
	HashIndexClear( cm_edgeHash );

    VectorCopy( bounds[0], cm_modelBounds[0] );
    VectorCopy( bounds[1], cm_modelBounds[1] );
	max = bounds[1][0] - bounds[0][0];
	f = bounds[1][1] - bounds[0][1];
	if ( f > max ) {
		max = f;
	}
	cm_vertexShift = (float) max / VERTEX_HASH_BOXSIZE;
	for ( i = 0; (1<<i) < cm_vertexShift; i++ ) {
	}
	if ( i == 0 ) {
		cm_vertexShift = 1;
	}
	else {
		cm_vertexShift = i;
	}
}

/*
================
:HashVec
================
*/
static ID_INLINE int HashVec(const vec3_t vec) {
	/*
	int x, y;

	x = (((int)(vec[0] - cm_modelBounds[0].x + 0.5 )) >> cm_vertexShift) & (VERTEX_HASH_BOXSIZE-1);
	y = (((int)(vec[1] - cm_modelBounds[0].y + 0.5 )) >> cm_vertexShift) & (VERTEX_HASH_BOXSIZE-1);

	assert (x >= 0 && x < VERTEX_HASH_BOXSIZE && y >= 0 && y < VERTEX_HASH_BOXSIZE);

	return y * VERTEX_HASH_BOXSIZE + x;
	*/
	int x, y, z;

	x = (((int) (vec[0] - cm_modelBounds[0][0] + 0.5)) + 2) >> 2;
	y = (((int) (vec[1] - cm_modelBounds[0][1] + 0.5)) + 2) >> 2;
	z = (((int) (vec[2] - cm_modelBounds[0][2] + 0.5)) + 2) >> 2;
	return (x + y * VERTEX_HASH_BOXSIZE + z) & (VERTEX_HASH_SIZE-1);
}

#define rint(f) (floorf( (f) + 0.5f ))

/*
================
GetVertex
================
*/
int GetVertex( cm_model_t *model, const vec3_t v, int *vertexNum ) {
	int i, hashKey, vn;
	vec3_t vert, *p;
	
	for (i = 0; i < 3; i++) {
		if ( fabs(v[i] - rint(v[i])) < INTEGRAL_EPSILON )
			vert[i] = rint(v[i]);
		else
			vert[i] = v[i];
	}

	hashKey = HashVec( vert );

	for (vn = HashIndexFirst( cm_vertexHash, hashKey ); vn >= 0; vn = HashIndexNext( cm_vertexHash, vn ) ) {
		p = &model->vertices[vn].p;
		// first compare z-axis because hash is based on x-y plane
		if (fabs(vert[2] - (*p)[2]) < VERTEX_EPSILON &&
			fabs(vert[0] - (*p)[0]) < VERTEX_EPSILON &&
			fabs(vert[1] - (*p)[1]) < VERTEX_EPSILON )
		{
			*vertexNum = vn;
			return qtrue;
		}
	}

	if ( model->numVertices >= model->maxVertices ) {
		cm_vertex_t *oldVertices;

		// resize vertex array
		model->maxVertices = (float) model->maxVertices * 1.5f + 1;
		oldVertices = model->vertices;
		model->vertices = (cm_vertex_t *) ii.GetMemory( model->maxVertices * sizeof(cm_vertex_t) );
        memset( model->vertices, 0, model->maxVertices * sizeof(cm_vertex_t) );
		memcpy( model->vertices, oldVertices, model->numVertices * sizeof(cm_vertex_t) );
		ii.FreeMemory( oldVertices );

		HashIndexResizeIndex( cm_vertexHash, model->maxVertices );
	}
    VectorCopy( vert, model->vertices[model->numVertices].p );
	model->vertices[model->numVertices].checkcount = 0;
	*vertexNum = model->numVertices;
	// add vertice to hash
	HashIndexAdd( cm_vertexHash, hashKey, model->numVertices );
	//
	model->numVertices++;
	return qfalse;
}

/*
================
GetEdge
================
*/
int GetEdge( cm_model_t *model, const vec3_t v1, const vec3_t v2, int *edgeNum, int v1num ) {
	int v2num, hashKey, e;
	int found, *vertexNum;

	// the first edge is a dummy
	if ( model->numEdges == 0 ) {
		model->numEdges = 1;
	}

	if ( v1num != -1 ) {
		found = 1;
	}
	else {
		found = GetVertex( model, v1, &v1num );
	}
	found &= GetVertex( model, v2, &v2num );
	// if both vertices are the same or snapped onto each other
	if ( v1num == v2num ) {
		*edgeNum = 0;
		return qtrue;
	}
	hashKey = HashIndexGenerateKeyForPair( cm_edgeHash, v1num, v2num );
	// if both vertices where already stored
	if (found) {
		for (e = HashIndexFirst( cm_edgeHash, hashKey ); e >= 0; e = HashIndexNext( cm_edgeHash, e ) )
		{
			// NOTE: only allow at most two users that use the edge in opposite direction
			if ( model->edges[e].numUsers != 1 ) {
				continue;
			}

			vertexNum = model->edges[e].vertexNum;
			if ( vertexNum[0] == v2num ) {
				if ( vertexNum[1] == v1num ) {
					// negative for a reversed edge
					*edgeNum = -e;
					break;
				}
			}
			/*
			else if ( vertexNum[0] == v1num ) {
				if ( vertexNum[1] == v2num ) {
					*edgeNum = e;
					break;
				}
			}
			*/
		}
		// if edge found in hash
		if ( e >= 0 ) {
			model->edges[e].numUsers++;
			return qtrue;
		}
	}
	if ( model->numEdges >= model->maxEdges ) {
		cm_edge_t *oldEdges;

		// resize edge array
		model->maxEdges = (float) model->maxEdges * 1.5f + 1;
		oldEdges = model->edges;
		model->edges = (cm_edge_t *) ii.GetMemory( model->maxEdges * sizeof(cm_edge_t) );
        memset( model->edges, 0, model->maxEdges * sizeof(cm_edge_t) );
		memcpy( model->edges, oldEdges, model->numEdges * sizeof(cm_edge_t) );
		ii.FreeMemory( oldEdges );

		HashIndexResizeIndex( cm_edgeHash, model->maxEdges );
	}
	// setup edge
	model->edges[model->numEdges].vertexNum[0] = v1num;
	model->edges[model->numEdges].vertexNum[1] = v2num;
	model->edges[model->numEdges].internal = qfalse;
	model->edges[model->numEdges].checkcount = 0;
	model->edges[model->numEdges].numUsers = 1; // used by one polygon atm
	VectorClear( model->edges[model->numEdges].normal );
	//
	*edgeNum = model->numEdges;
	// add edge to hash
	HashIndexAdd( cm_edgeHash, hashKey, model->numEdges );

	model->numEdges++;

	return qfalse;
}

/*
================
CreatePolygon
================
*/
void CreatePolygon( cm_model_t *model, fixedWinding_t *w, const plane_t plane, const qhandle_t material, int primitiveNum ) {
	int i, j, edgeNum, v1num;
	int numPolyEdges, polyEdges[MAX_POINTS_ON_WINDING];
	vec3_t bounds[2];
	cm_polygon_t *p;

	// turn the winding into a sequence of edges
	numPolyEdges = 0;
	v1num = -1;		// first vertex unknown
	for ( i = 0, j = 1; i < w->numPoints; i++, j++ ) {
		if ( j >= w->numPoints ) {
			j = 0;
		}
		GetEdge( model, w->points[i], w->points[j], &polyEdges[numPolyEdges], v1num );
		if ( polyEdges[numPolyEdges] ) {
			// last vertex of this edge is the first vertex of the next edge
			v1num = model->edges[ abs(polyEdges[numPolyEdges]) ].vertexNum[ INTSIGNBITNOTSET(polyEdges[numPolyEdges]) ];
			// this edge is valid so keep it
			numPolyEdges++;
		}
	}
	// should have at least 3 edges
	if ( numPolyEdges < 3 ) {
		return;
	}
	// the polygon is invalid if some edge is found twice
	for ( i = 0; i < numPolyEdges; i++ ) {
		for ( j = i+1; j < numPolyEdges; j++ ) {
			if ( abs(polyEdges[i]) == abs(polyEdges[j]) ) {
				return;
			}
		}
	}
	// don't overflow max edges
	if ( numPolyEdges > CM_MAX_POLYGON_EDGES ) {
		ii.Com_DPrintf( "CreatePolygon: polygon has more than %d edges", numPolyEdges );
		numPolyEdges = CM_MAX_POLYGON_EDGES;
	}

	FixedWindingBounds( w, bounds );

	p = AllocPolygon( model, numPolyEdges );
	p->numEdges = numPolyEdges;
	p->contents = CM_GetMaterialContentFlags( material );
	p->material = material;
	p->checkcount = 0;
    PlaneCopy( plane, p->plane );
	VectorCopy( bounds[0], p->bounds[0] );
	VectorCopy( bounds[1], p->bounds[1] );
	for ( i = 0; i < numPolyEdges; i++ ) {
		edgeNum = polyEdges[i];
		p->edges[i] = edgeNum;
	}
	R_FilterPolygonIntoTree( model, model->node, NULL, p );
}

/*
================
PolygonFromWinding

  NOTE: for patches primitiveNum < 0 and abs(primitiveNum) is the real number
================
*/
void PolygonFromWinding( cm_model_t *model, fixedWinding_t *w, const plane_t plane, const qhandle_t material, int primitiveNum ) {
	int contents;
	plane_t backplane;

	contents = CM_GetMaterialContentFlags( material );

	// if this polygon is part of the world model
	if ( cmLocal.numModels == 0 ) {
		// if the polygon is fully chopped away by the proc bsp tree
		if ( ChoppedAwayByProcBSP( w, plane, contents ) ) {
			model->numRemovedPolys++;
			return;
		}
	}

	// get one winding that is not or only partly contained in brushes
	w = WindingOutsideBrushes( w, plane, contents, primitiveNum, model->node );

	// if the polygon is fully contained within a brush
	if ( !w ) {
		model->numRemovedPolys++;
		return;
	}

	if ( FixedWindingIsHuge( w ) ) {
		ii.Com_DPrintf( "PolygonFromWinding: model %s primitive %d is degenerate", model->name, abs(primitiveNum) );
		return;
	}

	CreatePolygon( model, w, plane, material, primitiveNum );

	if ( CM_MaterialNeedsBackSide( material ) ) {
		ReverseFixedWindingSelf( w );
        PlaneNegate( plane, backplane );
		CreatePolygon( model, w, backplane, material, primitiveNum );
	}
}

/*
=================
CreatePatchPolygons
=================
*/
void CreatePatchPolygons( cm_model_t *model, surfacePatch_t *mesh, const qhandle_t material, int primitiveNum ) {
	int i, j;
	float dot;
	int v1, v2, v3, v4;
    surfVert_t *sv1, *sv2, *sv3, *sv4;
	fixedWinding_t w;
	plane_t plane;
	vec3_t d1, d2;
    vec5_t tmp;

    tmp[3] = tmp[4] = 0.0f;

	for ( i = 0; i < mesh->width - 1; i++ ) {
		for ( j = 0; j < mesh->height - 1; j++ ) {

			v1 = j * mesh->width + i;
            sv1 = SurfaceGetVertex( &mesh->surf, v1 );
			v2 = v1 + 1;
            sv2 = SurfaceGetVertex( &mesh->surf, v2 );
			v3 = v1 + mesh->width + 1;
            sv3 = SurfaceGetVertex( &mesh->surf, v3 );
			v4 = v1 + mesh->width;
            sv4 = SurfaceGetVertex( &mesh->surf, v4 );

			VectorSubtract( sv2->xyz, sv1->xyz, d1 );
			VectorSubtract( sv3->xyz, sv1->xyz, d2 );
            CrossProduct( d1, d2, plane );
			if ( VectorNormalize( plane ) != 0.0f ) {
				PlaneFitThroughPoint( plane, sv1->xyz );
				dot = PlaneDistance( plane, sv4->xyz );
				// if we can turn it into a quad
				if ( fabs(dot) < 0.1f ) {
					ClearFixedWinding( &w );
                    VectorCopy( sv1->xyz, tmp );
                    AddPointToFixedWinding( &w, tmp );
                    VectorCopy( sv2->xyz, tmp );
                    AddPointToFixedWinding( &w, tmp );
                    VectorCopy( sv3->xyz, tmp );
                    AddPointToFixedWinding( &w, tmp );
                    VectorCopy( sv4->xyz, tmp );
                    AddPointToFixedWinding( &w, tmp );

					PolygonFromWinding( model, &w, plane, material, -primitiveNum );
					continue;
				}
				else {
					// create one of the triangles
					ClearFixedWinding( &w );
                    VectorCopy( sv1->xyz, tmp );
                    AddPointToFixedWinding( &w, tmp );
                    VectorCopy( sv2->xyz, tmp );
                    AddPointToFixedWinding( &w, tmp );
                    VectorCopy( sv3->xyz, tmp );
                    AddPointToFixedWinding( &w, tmp );

					PolygonFromWinding( model, &w, plane, material, -primitiveNum );
				}
			}
			// create the other triangle
            VectorSubtract( sv3->xyz, sv1->xyz, d1 );
            VectorSubtract( sv4->xyz, sv1->xyz, d2 );
            CrossProduct( d1, d2, plane );
			if ( VectorNormalize( plane ) != 0.0f ) {
				PlaneFitThroughPoint( plane, sv1->xyz );

				ClearFixedWinding( &w );
                VectorCopy( sv1->xyz, tmp );
                AddPointToFixedWinding( &w, tmp );
                VectorCopy( sv3->xyz, tmp );
                AddPointToFixedWinding( &w, tmp );
                VectorCopy( sv4->xyz, tmp );
                AddPointToFixedWinding( &w, tmp );

				PolygonFromWinding( model, &w, plane, material, -primitiveNum );
			}
		}
	}
}

/*
=================
CM_EstimateVertsAndEdges
=================
*/
static void CM_EstimateVertsAndEdges( const mapEntity_t *mapEnt, int *numVerts, int *numEdges ) {
	int j, width, height;
	const mapPrimitive_t *mapPrim;

	*numVerts = *numEdges = 0;
	mapPrim = mapEnt->primitives;
	j = 0;
	while ( mapPrim ) {
		if ( mapPrim->type == PRIMTYPE_PATCH ) {
			// assume maximum tesselation without adding verts
			width = mapPrim->patch.width;
			height = mapPrim->patch.height;
			*numVerts += width * height;
			*numEdges += (width-1) * height + width * (height-1) + (width-1) * (height-1);
		}
		if ( mapPrim->type == PRIMTYPE_BRUSH ) {
			// assume cylinder with a polygon with (numSides - 2) edges ontop and on the bottom
			*numVerts += (GetNumSides( mapPrim ) - 2) * 2;
			*numEdges += (GetNumSides( mapPrim ) - 2) * 3;
		}
		j++;
		mapPrim = mapPrim->next;
	}
}

/*
=================
ConvertPatch
=================
*/
void ConvertPatch( cm_model_t *model, const mapPrimitive_t *patch, int primitiveNum ) {
	qhandle_t material;
	surfacePatch_t *cp;

	material = CM_FindMaterial( patch->material );
    // TODO: implement this?
	/*if ( !( material->GetContentFlags() & CONTENTS_REMOVE_UTIL ) ) {
		return;
	}*/

	// copy the patch
    cp = ( surfacePatch_t * )ii.GetMemory( sizeof( *cp ) );
    SurfacePatchInitFromPatch( cp, &patch->patch );

	// if the patch has an explicit number of subdivisions use it to avoid cracks
	if ( patch->explicitSubdivisions ) {
        SurfacePatchSubdivideExplicit( cp, patch->horzSubdivisions, patch->vertSubdivisions, qfalse, qtrue );
	} else {
        SurfacePatchSubdivide( cp, DEFAULT_CURVE_MAX_ERROR_CD, DEFAULT_CURVE_MAX_ERROR_CD, DEFAULT_CURVE_MAX_LENGTH_CD, qfalse );
	}

	// create collision polygons for the patch
	CreatePatchPolygons( model, cp, material, primitiveNum );

	SurfacePatchFree( cp );
    ii.FreeMemory( cp );
}

/*
================
ConvertBrushSides
================
*/
#define DEGENERATE_DIST_EPSILON		1e-4f

void ConvertBrushSides( cm_model_t *model, const mapPrimitive_t *mapBrush, int primitiveNum ) {
	int i, j;
	mapBrushSide_t *mapSide1, *mapSide;
	fixedWinding_t w;
	plane_t *planes;
	qhandle_t material;
	vec3_t tmp;
	plane_t pl;

	assert( mapBrush->type == PRIMTYPE_BRUSH );

	// fix degenerate planes
	planes = (plane_t *) ii.GetMemory( GetNumSides( mapBrush ) * sizeof( *planes ) );
	for ( i = 0, mapSide = mapBrush->sides; mapSide; i++, mapSide = mapSide->next ) {
        PlaneCopy( mapSide->plane, planes[i] );
		PlaneFixDegeneracies( planes[i], DEGENERATE_DIST_EPSILON );
	}

	// create a collision polygon for each brush side
	for ( i = 0, mapSide1 = mapBrush->sides; mapSide1; i++, mapSide1 = mapSide1->next ) {
		material = CM_RegisterMaterial( mapSide1->material );
		// TODO: what's this?
		//if ( !( CMI_GetMaterialContentFlags( material ) & CONTENTS_REMOVE_UTIL ) ) {
		//	continue;
		//}
		PlaneNegate( planes[i], pl );
		FixedWindingBaseForPlane( &w, pl );
		for ( j = 0, mapSide = mapBrush->sides; mapSide && w.numPoints; j++, mapSide = mapSide->next ) {
			if ( i == j ) {
				continue;
			}

			PlaneNegate( planes[j], pl );
			ClipFixedWindingInPlace( &w, pl, 0, qfalse );
		}

		if ( w.numPoints ) {
			PolygonFromWinding( model, &w, planes[i], material, primitiveNum );
		}
	}

	ii.FreeMemory( planes );
}

/*
================
ConvertBrush
================
*/
void ConvertBrush( cm_model_t *model, const mapPrimitive_t *mapBrush, int primitiveNum ) {
	int i, j, contents;
	vec3_t bounds[2];
	mapBrushSide_t *mapSide, *mapSide1;
	cm_brush_t *brush;
	plane_t *planes;
	fixedWinding_t w;
	qhandle_t material = 0;
	plane_t pl;

	contents = 0;
	ClearBounds( bounds[0], bounds[1] );

	// fix degenerate planes
	planes = (plane_t *) ii.GetMemory( GetNumSides( mapBrush ) * sizeof( *planes ) );
	for ( i = 0, mapSide = mapBrush->sides; mapSide; i++, mapSide = mapSide->next ) {
        PlaneCopy( mapSide->plane, planes[i] );
		PlaneFixDegeneracies( planes[i], DEGENERATE_DIST_EPSILON );
	}

	// we are only getting the bounds for the brush so there's no need
	// to create a winding for the last brush side
	for ( i = 0, mapSide = mapBrush->sides; mapSide && mapSide->next; i++, mapSide = mapSide->next ) {
		material = CM_RegisterMaterial( mapSide->material );
		// TODO: what's this?
		contents |= ( CM_GetMaterialContentFlags( material ) /*& CONTENTS_REMOVE_UTIL*/ );
        PlaneNegate( planes[i], pl );
		FixedWindingBaseForPlane( &w, pl );
		for ( j = 0, mapSide1 = mapBrush->sides; mapSide1 && w.numPoints; j++, mapSide1 = mapSide1->next ) {
			if ( i == j ) {
				continue;
			}
            PlaneNegate( planes[j], pl );
			ClipFixedWindingInPlace( &w, pl, 0, qfalse );
		}

		for ( j = 0; j < w.numPoints; j++ ) {
			AddPointToBounds( w.points[j], bounds[0], bounds[1] );
		}
	}
	if ( !contents ) {
		ii.FreeMemory( planes );
		return;
	}
	// create brush for position test
	brush = AllocBrush( model, GetNumSides( mapBrush ) );
	brush->checkcount = 0;
	brush->contents = contents;
	brush->material = material;
	brush->primitiveNum = primitiveNum;
	VectorCopy( bounds[0], brush->bounds[0] );
	VectorCopy( bounds[1], brush->bounds[1] );
	brush->numPlanes = GetNumSides( mapBrush );
	for (i = 0; i < GetNumSides( mapBrush ); i++) {
		PlaneCopy( planes[i], brush->planes[i] );
	}
	AddBrushToNode( model, model->node, brush );
	ii.FreeMemory( planes );
}

/*
================
CM_CountNodeBrushes
================
*/
static int CM_CountNodeBrushes( const cm_node_t *node ) {
	int count;
	cm_brushRef_t *bref;

	count = 0;
	for ( bref = node->brushes; bref; bref = bref->next ) {
		count++;
	}
	return count;
}

/*
================
CM_R_GetModelBounds
================
*/
static void CM_R_GetNodeBounds( vec3_t bounds[2], cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_brushRef_t *bref;

	while ( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
            AddPointToBounds( pref->p->bounds[0], bounds[0], bounds[1] );
            AddPointToBounds( pref->p->bounds[1], bounds[0], bounds[1] );
		}
		for ( bref = node->brushes; bref; bref = bref->next ) {
            AddPointToBounds( bref->b->bounds[0], bounds[0], bounds[1] );
            AddPointToBounds( bref->b->bounds[1], bounds[0], bounds[1] );
		}
		if ( node->planeType == -1 ) {
			break;
		}
		CM_R_GetNodeBounds( bounds, node->children[1] );
		node = node->children[0];
	}
}

/*
================
CM_GetNodeBounds
================
*/
void CM_GetNodeBounds( vec3_t bounds[2], cm_node_t *node ) {
	ClearBounds( bounds[0], bounds[1] );
	CM_R_GetNodeBounds( bounds, node );
	if ( bounds[0][0] > bounds[1][0] ) {
		VectorClear( bounds[0] );
		VectorClear( bounds[1] );
	}
}

/*
================
CM_GetNodeContents
================
*/
int CM_GetNodeContents( cm_node_t *node ) {
	int contents;
	cm_polygonRef_t *pref;
	cm_brushRef_t *bref;

	contents = 0;
	while ( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			contents |= pref->p->contents;
		}
		for ( bref = node->brushes; bref; bref = bref->next ) {
			contents |= bref->b->contents;
		}
		if ( node->planeType == -1 ) {
			break;
		}
		contents |= CM_GetNodeContents( node->children[1] );
		node = node->children[0];
	}
	return contents;
}

/*
==================
RemapEdges
==================
*/
void RemapEdges( cm_node_t *node, int *edgeRemap ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;
	int i;

	while ( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			p = pref->p;
			// if we checked this polygon already
			if ( p->checkcount == cmLocal.checkCount ) {
				continue;
			}
			p->checkcount = cmLocal.checkCount;
			for ( i = 0; i < p->numEdges; i++ ) {
				if ( p->edges[i] < 0 ) {
					p->edges[i] = -edgeRemap[ abs(p->edges[i]) ];
				}
				else {
					p->edges[i] = edgeRemap[ p->edges[i] ];
				}
			}
		}
		if ( node->planeType == -1 ) {
			break;
		}

		RemapEdges( node->children[1], edgeRemap );
		node = node->children[0];
	}
}

/*
==================
OptimizeArrays

  due to polygon merging and polygon removal the vertex and edge array
  can have a lot of unused entries.
==================
*/
void OptimizeArrays( cm_model_t *model ) {
	int i, newNumVertices, newNumEdges, *v;
	int *remap;
	cm_edge_t *oldEdges;
	cm_vertex_t *oldVertices;
	int remapSize;

	remapSize = ( model->numVertices > model->numEdges ? model->numVertices : model->numEdges ) * sizeof( int );
	remap = (int *) ii.GetMemory( remapSize );
	memset( remap, 0, remapSize );
	// get all used vertices
	for ( i = 0; i < model->numEdges; i++ ) {
		remap[ model->edges[i].vertexNum[0] ] = qtrue;
		remap[ model->edges[i].vertexNum[1] ] = qtrue;
	}
	// create remap index and move vertices
	newNumVertices = 0;
	for ( i = 0; i < model->numVertices; i++ ) {
		if ( remap[ i ] ) {
			remap[ i ] = newNumVertices;
			model->vertices[ newNumVertices ] = model->vertices[ i ];
			newNumVertices++;
		}
	}
	model->numVertices = newNumVertices;
	// change edge vertex indexes
	for ( i = 1; i < model->numEdges; i++ ) {
		v = model->edges[i].vertexNum;
		v[0] = remap[ v[0] ];
		v[1] = remap[ v[1] ];
	}

	// create remap index and move edges
	newNumEdges = 1;
	for ( i = 1; i < model->numEdges; i++ ) {
		// if the edge is used
		if ( model->edges[ i ].numUsers ) {
			remap[ i ] = newNumEdges;
			model->edges[ newNumEdges ] = model->edges[ i ];
			newNumEdges++;
		}
	}
	// change polygon edge indexes
	cmLocal.checkCount++;
	RemapEdges( model->node, remap );
	model->numEdges = newNumEdges;

	ii.FreeMemory( remap );

	// realloc vertices
	oldVertices = model->vertices;
	if ( oldVertices ) {
		model->vertices = (cm_vertex_t *) ii.GetMemory( model->numVertices * sizeof(cm_vertex_t) );
		memset( model->vertices, 0, model->numVertices * sizeof(cm_vertex_t) );
		memcpy( model->vertices, oldVertices, model->numVertices * sizeof(cm_vertex_t) );
		ii.FreeMemory( oldVertices );
	}

	// realloc edges
	oldEdges = model->edges;
	if ( oldEdges ) {
		model->edges = (cm_edge_t *) ii.GetMemory( model->numEdges * sizeof(cm_edge_t) );
		memset( model->edges, 0, model->numEdges * sizeof(cm_edge_t) );
		memcpy( model->edges, oldEdges, model->numEdges * sizeof(cm_edge_t) );
		ii.FreeMemory( oldEdges );
	}
}

/*
================
FinishModel
================
*/
void FinishModel( cm_model_t *model ) {
	// try to merge polygons
	cmLocal.checkCount++;
	CM_MergeTreePolygons( model, model->node );
	// find internal edges (no mesh can ever collide with internal edges)
	cmLocal.checkCount++;
	CM_FindInternalEdges( model, model->node );
	// calculate edge normals
	cmLocal.checkCount++;
	CM_CalculateEdgeNormals( model, model->node );

	ii.Com_Printf( "%s vertex hash spread is %d\n", model->name, HashIndexGetSpread( cm_vertexHash ) );
	ii.Com_Printf( "%s edge hash spread is %d\n", model->name, HashIndexGetSpread( cm_edgeHash ) );

	// remove all unused vertices and edges
	OptimizeArrays( model );
	// get model bounds from brush and polygon bounds
	CM_GetNodeBounds( model->bounds, model->node );
	// get model contents
	model->contents = CM_GetNodeContents( model->node );
	// total memory used by this model
	model->usedMemory = model->numVertices * sizeof(cm_vertex_t) +
						model->numEdges * sizeof(cm_edge_t) +
						model->polygonMemory +
						model->brushMemory +
						model->numNodes * sizeof(cm_node_t) +
						model->numPolygonRefs * sizeof(cm_polygonRef_t) +
						model->numBrushRefs * sizeof(cm_brushRef_t);
}

/*
================
CM_LoadRenderModel
================
*/
cm_model_t *CM_LoadRenderModel( const char *fileName ) {
#if 1
	return NULL;
#else
	int i, j;
	idRenderModel *renderModel;
	const modelSurface_t *surf;
	idFixedWinding w;
	cm_node_t *node;
	cm_model_t *model;
	idPlane plane;
	idBounds bounds;
	bool collisionSurface;
	idStr extension;

	// only load ASE and LWO models
	idStr( fileName ).ExtractFileExtension( extension );
	if ( ( extension.Icmp( "ase" ) != 0 ) && ( extension.Icmp( "lwo" ) != 0 ) && ( extension.Icmp( "ma" ) != 0 ) ) {
		return NULL;
	}

	if ( !renderModelManager->CheckModel( fileName ) ) {
		return NULL;
	}

	renderModel = renderModelManager->FindModel( fileName );

	model = CM_AllocModel();
	model->name = fileName;
	node = AllocNode( model, NODE_BLOCK_SIZE_SMALL );
	node->planeType = -1;
	model->node = node;

	model->maxVertices = 0;
	model->numVertices = 0;
	model->maxEdges = 0;
	model->numEdges = 0;

	bounds = renderModel->Bounds( NULL );

	collisionSurface = false;
	for ( i = 0; i < renderModel->NumSurfaces(); i++ ) {
		surf = renderModel->Surface( i );
		if ( surf->shader->GetSurfaceFlags() & SURF_COLLISION ) {
			collisionSurface = true;
		}
	}

	for ( i = 0; i < renderModel->NumSurfaces(); i++ ) {
		surf = renderModel->Surface( i );
		// if this surface has no contents
		if ( ! ( surf->shader->GetContentFlags() & CONTENTS_REMOVE_UTIL ) ) {
			continue;
		}
		// if the model has a collision surface and this surface is not a collision surface
		if ( collisionSurface && !( surf->shader->GetSurfaceFlags() & SURF_COLLISION ) ) {
			continue;
		}
		// get max verts and edges
		model->maxVertices += surf->geometry->numVerts;
		model->maxEdges += surf->geometry->numIndexes;
	}

	model->vertices = (cm_vertex_t *) Mem_ClearedAlloc( model->maxVertices * sizeof(cm_vertex_t) );
	model->edges = (cm_edge_t *) Mem_ClearedAlloc( model->maxEdges * sizeof(cm_edge_t) );

	// setup hash to speed up finding shared vertices and edges
	SetupHash();

	cm_vertexHash->ResizeIndex( model->maxVertices );
	cm_edgeHash->ResizeIndex( model->maxEdges );

	ClearHash( bounds );

	for ( i = 0; i < renderModel->NumSurfaces(); i++ ) {
		surf = renderModel->Surface( i );
		// if this surface has no contents
		if ( ! ( surf->shader->GetContentFlags() & CONTENTS_REMOVE_UTIL ) ) {
			continue;
		}
		// if the model has a collision surface and this surface is not a collision surface
		if ( collisionSurface && !( surf->shader->GetSurfaceFlags() & SURF_COLLISION ) ) {
			continue;
		}

		for ( j = 0; j < surf->geometry->numIndexes; j += 3 ) {
			w.Clear();
			w += surf->geometry->verts[ surf->geometry->indexes[ j + 2 ] ].xyz;
			w += surf->geometry->verts[ surf->geometry->indexes[ j + 1 ] ].xyz;
			w += surf->geometry->verts[ surf->geometry->indexes[ j + 0 ] ].xyz;
			w.GetPlane( plane );
			plane = -plane;
			PolygonFromWinding( model, &w, plane, surf->shader, 1 );
		}
	}

	// create a BSP tree for the model
	model->node = CreateAxialBSPTree( model, model->node );

	model->isConvex = false;

	FinishModel( model );

	// shutdown the hash
	ShutdownHash();

	common->Printf( "loaded collision model %s\n", model->name.c_str() );

	return model;
#endif
}

/*
================
CM_CollisionModelForMapEntity
================
*/
cm_model_t *CM_CollisionModelForMapEntity( const mapEntity_t *mapEnt ) {
	cm_model_t *model;
	vec3_t bounds[2];
	const char *name;
	int i, brushCount;
	mapPrimitive_t	*mapPrim;

	// if the entity has no primitives
	if ( mapEnt->numPrimitives < 1 ) {
		return NULL;
	}

	// get a name for the collision model
	DictGetString( &mapEnt->epairs, "model", "", &name );
	if ( !name[0] ) {
		DictGetString( &mapEnt->epairs, "name", "", &name );
		if ( !name[0] ) {
			if ( !cmLocal.numModels ) {
				// first model is always the world
				name = "worldMap";
			}
			else {
				name = "unnamed inline model";
			}
		}
	}

	model = CM_AllocModel();
	model->node = AllocNode( model, NODE_BLOCK_SIZE_SMALL );

	CM_EstimateVertsAndEdges( mapEnt, &model->maxVertices, &model->maxEdges );
	model->numVertices = 0;
	model->numEdges = 0;
	model->vertices = (cm_vertex_t *) ii.GetMemory( model->maxVertices * sizeof(cm_vertex_t) );
	Com_Memset( model->vertices, 0, model->maxVertices * sizeof(cm_vertex_t) );
	model->edges = (cm_edge_t *) ii.GetMemory( model->maxEdges * sizeof(cm_edge_t) );
	Com_Memset( model->edges, 0, model->maxEdges * sizeof(cm_edge_t) );

	HashIndexResizeIndex( cm_vertexHash, model->maxVertices );
	HashIndexResizeIndex( cm_edgeHash, model->maxEdges );

	Q_strncpyz( model->name, name, sizeof( model->name ) );
	model->isConvex = qfalse;

	// convert brushes
	for ( i = 0, mapPrim = mapEnt->primitives; mapPrim; i++, mapPrim = mapPrim->next ) {
		if ( mapPrim->type == PRIMTYPE_BRUSH ) {
			ConvertBrush( model, mapPrim, i );
		}
	}

	// create an axial bsp tree for the model if it has more than just a bunch brushes
	brushCount = CM_CountNodeBrushes( model->node );
	if ( brushCount > 4 ) {
		model->node = CreateAxialBSPTree( model, model->node );
	} else {
		model->node->planeType = -1;
	}

	// get bounds for hash
	if ( brushCount ) {
		CM_GetNodeBounds( bounds, model->node );
	} else {
		VectorSet( bounds[0], -256, -256, -256 );
		VectorSet( bounds[1], 256, 256, 256 );
	}

	// different models do not share edges and vertices with each other, so clear the hash
	ClearHash( bounds );

	// create polygons from patches and brushes
	for ( i = 0, mapPrim = mapEnt->primitives; mapPrim; i++, mapPrim = mapPrim->next ) {
		if ( mapPrim->type == PRIMTYPE_PATCH ) {
			ConvertPatch( model, mapPrim, i );
		}
		if ( mapPrim->type == PRIMTYPE_BRUSH ) {
			ConvertBrushSides( model, mapPrim, i );
		}
	}

	FinishModel( model );

	return model;
}

/*
================
CM_FindMaterial
================
*/
qhandle_t CM_FindMaterial( const char *name ) {
    int i;

    // check if this material is already loaded
    for ( i = 0; i < cmLocal.numMaterials; i++ ) {
        if ( !Q_stricmp( cmLocal.materials[i]->name, name ) ) {
            break;
        }
    }
    // if the material is already loaded
    if ( i < cmLocal.numMaterials ) {
        return i;
    }

    return 0;
}

/*
================
CM_FindModel
================
*/
cmHandle_t CM_FindModel( const char *name ) {
	int i;

	// check if this model is already loaded
	for ( i = 0; i < cmLocal.numModels; i++ ) {
		if ( !Q_stricmp( cmLocal.models[i]->name, name ) ) {
			break;
		}
	}
	// if the model is already loaded
	if ( i < cmLocal.numModels ) {
		return i;
	}
	return -1;
}

/*
==================
PrintModelInfo
==================
*/
void PrintModelInfo( const cm_model_t *model ) {
	ii.Com_Printf( "%6i vertices (%i KB)\n", model->numVertices, (model->numVertices * sizeof(cm_vertex_t))>>10 );
	ii.Com_Printf( "%6i edges (%i KB)\n", model->numEdges, (model->numEdges * sizeof(cm_edge_t))>>10 );
	ii.Com_Printf( "%6i polygons (%i KB)\n", model->numPolygons, model->polygonMemory>>10 );
	ii.Com_Printf( "%6i brushes (%i KB)\n", model->numBrushes, model->brushMemory>>10 );
	ii.Com_Printf( "%6i nodes (%i KB)\n", model->numNodes, (model->numNodes * sizeof(cm_node_t))>>10 );
	ii.Com_Printf( "%6i polygon refs (%i KB)\n", model->numPolygonRefs, (model->numPolygonRefs * sizeof(cm_polygonRef_t))>>10 );
	ii.Com_Printf( "%6i brush refs (%i KB)\n", model->numBrushRefs, (model->numBrushRefs * sizeof(cm_brushRef_t))>>10 );
	ii.Com_Printf( "%6i internal edges\n", model->numInternalEdges );
	ii.Com_Printf( "%6i sharp edges\n", model->numSharpEdges );
	ii.Com_Printf( "%6i contained polygons removed\n", model->numRemovedPolys );
	ii.Com_Printf( "%6i polygons merged\n", model->numMergedPolys );
	ii.Com_Printf( "%6i KB total memory used\n", model->usedMemory>>10 );
}

/*
================
AccumulateModelInfo
================
*/
void AccumulateModelInfo( cm_model_t *model ) {
	int i;

	memset( model, 0, sizeof( *model ) );
	// accumulate statistics of all loaded models
	for ( i = 0; i < cmLocal.numModels; i++ ) {
		model->numVertices += cmLocal.models[i]->numVertices;
		model->numEdges += cmLocal.models[i]->numEdges;
		model->numPolygons += cmLocal.models[i]->numPolygons;
		model->polygonMemory += cmLocal.models[i]->polygonMemory;
		model->numBrushes += cmLocal.models[i]->numBrushes;
		model->brushMemory += cmLocal.models[i]->brushMemory;
		model->numNodes += cmLocal.models[i]->numNodes;
		model->numBrushRefs += cmLocal.models[i]->numBrushRefs;
		model->numPolygonRefs += cmLocal.models[i]->numPolygonRefs;
		model->numInternalEdges += cmLocal.models[i]->numInternalEdges;
		model->numSharpEdges += cmLocal.models[i]->numSharpEdges;
		model->numRemovedPolys += cmLocal.models[i]->numRemovedPolys;
		model->numMergedPolys += cmLocal.models[i]->numMergedPolys;
		model->usedMemory += cmLocal.models[i]->usedMemory;
	}
}

/*
================
CM_ModelInfo
================
*/
void CM_ModelInfo( cmHandle_t model ) {
	cm_model_t modelInfo;

	if ( model == -1 ) {
		AccumulateModelInfo( &modelInfo );
		PrintModelInfo( &modelInfo );
		return;
	}
	if ( model < 0 || model > MAX_SUBMODELS || model > cmLocal.maxModels ) {
		ii.Com_Printf( "ModelInfo: invalid model handle\n" );
		return;
	}
	if ( !cmLocal.models[model] ) {
		ii.Com_Printf( "ModelInfo: invalid model\n" );
		return;
	}

	PrintModelInfo( cmLocal.models[model] );
}

/*
================
CM_ListModels
================
*/
void CM_ListModels( void ) {
	int i, totalMemory;

	totalMemory = 0;
	for ( i = 0; i < cmLocal.numModels; i++ ) {
		ii.Com_Printf( "%4d: %5d KB   %s\n", i, (cmLocal.models[i]->usedMemory>>10), cmLocal.models[i]->name );
		totalMemory += cmLocal.models[i]->usedMemory;
	}
	ii.Com_Printf( "%4d KB in %d models\n", (totalMemory>>10), cmLocal.numModels );
}

/*
================
BuildModels
================
*/
void BuildModels( const mapFile_t *mapFile ) {
	int i;
	const mapEntity_t *mapEnt;

	//idTimer timer;
	//timer.Start();

	if ( !CM_LoadCollisionModelFile( mapFile->name, mapFile->geometryCRC ) ) {

		if ( !mapFile->numEntities ) {
			return;
		}

		// load the .proc file bsp for data optimisation
		CM_LoadProcBSP( mapFile->name );

		// convert brushes and patches to collision data
		mapEnt = mapFile->entities;
		i = 0;
		while ( mapEnt ) {
			if ( cmLocal.numModels >= MAX_SUBMODELS ) {
				ii.Com_Error( ERR_DROP, "BuildModels: more than %d collision models", MAX_SUBMODELS );
				break;
			}
			cmLocal.models[cmLocal.numModels] = CM_CollisionModelForMapEntity( mapEnt );
			if ( cmLocal.models[ cmLocal.numModels] ) {
				cmLocal.numModels++;
			}
			mapEnt = mapEnt->next;
			i++;
		}

		// free the proc bsp which is only used for data optimization
        if ( cmLocal.procNodes ) {
            ii.FreeMemory( cmLocal.procNodes );
            cmLocal.procNodes = NULL;
        }

		// write the collision models to a file
		CM_WriteCollisionModelsToFile( mapFile->name, 0, cmLocal.numModels, mapFile->geometryCRC );
	}

	//timer.Stop();

	// print statistics on collision data
	cm_model_t model;
	AccumulateModelInfo( &model );
	ii.Com_Printf( "collision data:\n" );
	ii.Com_Printf( "%6i models\n", cmLocal.numModels );
	PrintModelInfo( &model );
	//ii.Com_Printf( "%.0f msec to load collision data.\n", timer.Milliseconds() );
}


/*
================
CM_LoadMap2
================
*/
void CM_LoadMap2( const mapFile_t *mapFile, qboolean clear ) {

	if ( mapFile == NULL ) {
		ii.Com_Error( ERR_DROP, "CM_LoadMap: NULL mapFile" );
	}

	// check whether we can keep the current collision map based on the mapName and mapFileTime
	if ( cmLocal.loaded ) {
		if ( Q_stricmp( cmLocal.mapName, mapFile->name ) == 0 ) {
			if ( mapFile->fileTime == cmLocal.mapFileTime ) {
				ii.Com_DPrintf( "Using loaded version\n" );
				return;
			}
			ii.Com_DPrintf( "Reloading modified map\n" );
		}
		CM_FreeMap();
	}

	// clear the collision map
    if ( clear ) {
	    CM_Clear();
    }

	// models
	cmLocal.maxModels = MAX_SUBMODELS;
	cmLocal.numModels = 0;
	cmLocal.models = (cm_model_t **) ii.GetMemory( (cmLocal.maxModels+1) * sizeof(cm_model_t *) );
	memset( cmLocal.models, 0, (cmLocal.maxModels+1) * sizeof(cm_model_t *) );

	// setup hash to speed up finding shared vertices and edges
	SetupHash();

	// setup trace model structure
	SetupTrmModelStructure();

	// build collision models
	BuildModels( mapFile );

	// save name and time stamp
	Q_strncpyz( cmLocal.mapName, mapFile->name, sizeof( cmLocal.mapName ) );
	cmLocal.mapFileTime = mapFile->fileTime;
	cmLocal.loaded = qtrue;

	// shutdown the hash
	ShutdownHash();
}

/*
================
CM_LoadMap
================
*/
void CM_LoadMap( const mapFile_t *mapFile ) {
    CM_LoadMap2( mapFile, qtrue );
}

/*
===================
CM_GetModelName
===================
*/
const char *CM_GetModelName( cmHandle_t model ) {
	if ( model < 0 || model > MAX_SUBMODELS || model >= cmLocal.numModels || !cmLocal.models[model] ) {
		ii.Com_Printf( "CM_GetModelName: invalid model handle\n" );
		return "";
	}
	return cmLocal.models[model]->name;
}

/*
===================
CM_GetModelBounds
===================
*/
qboolean CM_GetModelBounds( cmHandle_t model, vec3_t bounds[2] ) {

	if ( model < 0 || model > MAX_SUBMODELS || model >= cmLocal.numModels || !cmLocal.models[model] ) {
		ii.Com_Printf( "CM_GetModelBounds: invalid model handle\n" );
		return qfalse;
	}

	VectorCopy( cmLocal.models[model]->bounds[0], bounds[0] );
	VectorCopy( cmLocal.models[model]->bounds[1], bounds[1] );
	return qtrue;
}

/*
===================
CM_GetModelContents
===================
*/
qboolean CM_GetModelContents( cmHandle_t model, int *contents ) {
	if ( model < 0 || model > MAX_SUBMODELS || model >= cmLocal.numModels || !cmLocal.models[model] ) {
		ii.Com_Printf( PRINT_ALL, "GetModelContents: invalid model handle\n" );
		return qfalse;
	}

	*contents = cmLocal.models[model]->contents;

	return qtrue;
}

/*
===================
CM_GetModelVertex
===================
*/
qboolean CM_GetModelVertex( cmHandle_t model, int vertexNum, vec3_t vertex ) {
	if ( model < 0 || model > MAX_SUBMODELS || model >= cmLocal.numModels || !cmLocal.models[model] ) {
		ii.Com_Printf( "CM_GetModelVertex: invalid model handle\n" );
		return qfalse;
	}

	if ( vertexNum < 0 || vertexNum >= cmLocal.models[model]->numVertices ) {
		ii.Com_Printf( "CM_GetModelVertex: invalid vertex number\n" );
		return qfalse;
	}

	VectorCopy( cmLocal.models[model]->vertices[vertexNum].p, vertex );

	return qtrue;
}

/*
===================
CM_GetModelEdge
===================
*/
qboolean CM_GetModelEdge( cmHandle_t model, int edgeNum, vec3_t start, vec3_t end ) {
	if ( model < 0 || model > MAX_SUBMODELS || model >= cmLocal.numModels || !cmLocal.models[model] ) {
		ii.Com_Printf( "GetModelEdge: invalid model handle\n" );
		return qfalse;
	}

	edgeNum = abs( edgeNum );
	if ( edgeNum >= cmLocal.models[model]->numEdges ) {
		ii.Com_Printf( "GetModelEdge: invalid edge number\n" );
		return qfalse;
	}

	start = cmLocal.models[model]->vertices[cmLocal.models[model]->edges[edgeNum].vertexNum[0]].p;
	end = cmLocal.models[model]->vertices[cmLocal.models[model]->edges[edgeNum].vertexNum[1]].p;

	return qtrue;
}

/*
===================
CM_GetModelPolygon
===================
*/
qboolean CM_GetModelPolygon( cmHandle_t model, int polygonNum, fixedWinding_t *winding ) {
	int i, edgeNum;
	cm_polygon_t *poly;
	vec5_t tmp;

	if ( model < 0 || model > MAX_SUBMODELS || model >= cmLocal.numModels || !cmLocal.models[model] ) {
		ii.Com_Printf( "CM_GetModelPolygon: invalid model handle\n" );
		return qfalse;
	}

	poly = *(cm_polygon_t **)(&polygonNum);
	ClearFixedWinding( winding );
	for ( i = 0; i < poly->numEdges; i++ ) {
		edgeNum = poly->edges[i];
		VectorCopy( cmLocal.models[model]->vertices[ cmLocal.models[model]->edges[abs(edgeNum)].vertexNum[INTSIGNBITSET(edgeNum)] ].p, tmp );
		tmp[3] = 0.0f;
		tmp[4] = 0.0f;
		AddPointToFixedWinding( winding, tmp );
	}

	return qtrue;
}

/*
==================
CM_LoadModel
==================
*/
cmHandle_t CM_LoadModel( const char *modelName, const qboolean precache ) {
	int handle;

	handle = CM_FindModel( modelName );
	if ( handle >= 0 ) {
		return handle;
	}

	if ( cmLocal.numModels >= MAX_SUBMODELS ) {
		ii.Com_Error( ERR_DROP, "LoadModel: no free slots\n" );
		return 0;
	}

	// try to load a .cm file
	if ( CM_LoadCollisionModelFile( modelName, 0 ) ) {
		handle = CM_FindModel( modelName );
		if ( handle >= 0 ) {
			return handle;
		} else {
			ii.Com_DPrintf( "LoadModel: collision file for '%s' contains different model", modelName );
		}
	}

	// if only precaching .cm files do not waste memory converting render models
	if ( precache ) {
		return 0;
	}

	// try to load a .ASE or .LWO model and convert it to a collision model
	cmLocal.models[cmLocal.numModels] = CM_LoadRenderModel( modelName );
	if ( cmLocal.models[cmLocal.numModels] != NULL ) {
		cmLocal.numModels++;
		return ( cmLocal.numModels - 1 );
	}

	return 0;
}

/*
==================
TrmFromModel_r
==================
*/
qboolean TrmFromModel_r( traceModel_t *trm, cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;
	int i;

	while ( 1 ) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			p = pref->p;

			if ( p->checkcount == cmLocal.checkCount ) {
				continue;
			}

			p->checkcount = cmLocal.checkCount;

			if ( trm->numPolys >= MAX_TRACEMODEL_POLYS ) {
				return qfalse;
			}
			// copy polygon properties
			VectorCopy( p->bounds[0], trm->polys[ trm->numPolys ].bounds[0] );
			VectorCopy( p->bounds[1], trm->polys[ trm->numPolys ].bounds[1] );
			VectorCopy( p->plane, trm->polys[ trm->numPolys ].normal );
			trm->polys[ trm->numPolys ].dist = PlaneGetDist( p->plane );
			trm->polys[ trm->numPolys ].numEdges = p->numEdges;
			// copy edge index
			for ( i = 0; i < p->numEdges; i++ ) {
				trm->polys[ trm->numPolys ].edges[ i ] = p->edges[ i ];
			}
			trm->numPolys++;
		}
		if ( node->planeType == -1 ) {
			break;
		}
		if ( !TrmFromModel_r( trm, node->children[1] ) ) {
			return qfalse;
		}
		node = node->children[0];
	}
	return qtrue;
}

/*
==================
TrmFromModel2

  NOTE: polygon merging can merge colinear edges and as such might cause dangling edges.
==================
*/
qboolean TrmFromModel2( const cm_model_t *model, traceModel_t *trm ) {
	int i, j, numEdgeUsers[MAX_TRACEMODEL_EDGES+1];

	// if the model has too many vertices to fit in a trace model
	if ( model->numVertices > MAX_TRACEMODEL_VERTS ) {
		ii.Com_Printf( "TrmFromModel2: model %s has too many vertices.\n", model->name );
		PrintModelInfo( model );
		return qfalse;
	}

	// plus one because the collision model accounts for the first unused edge
	if ( model->numEdges > MAX_TRACEMODEL_EDGES+1 ) {
		ii.Com_Printf( "TrmFromModel2: model %s has too many edges.\n", model->name );
		PrintModelInfo( model );
		return qfalse;
	}

	trm->type = TRM_CUSTOM;
	trm->numVerts = 0;
	trm->numEdges = 1;
	trm->numPolys = 0;
    ClearBounds( trm->bounds[0], trm->bounds[1] );

	// copy polygons
	cmLocal.checkCount++;
	if ( !TrmFromModel_r( trm, model->node ) ) {
		ii.Com_Printf( "TrmFromModel2: model %s has too many polygons.\n", model->name );
		PrintModelInfo( model );
		return qfalse;
	}

	// copy vertices
	for ( i = 0; i < model->numVertices; i++ ) {
		VectorCopy( model->vertices[ i ].p, trm->verts[ i ] );
		AddPointToBounds( trm->verts[ i ], trm->bounds[0], trm->bounds[1] );
	}
	trm->numVerts = model->numVertices;

	// copy edges
	for ( i = 0; i < model->numEdges; i++ ) {
		trm->edges[ i ].v[0] = model->edges[ i ].vertexNum[0];
		trm->edges[ i ].v[1] = model->edges[ i ].vertexNum[1];
	}
	// minus one because the collision model accounts for the first unused edge
	trm->numEdges = model->numEdges - 1;

	// each edge should be used exactly twice
	memset( numEdgeUsers, 0, sizeof(numEdgeUsers) );
	for ( i = 0; i < trm->numPolys; i++ ) {
		for ( j = 0; j < trm->polys[i].numEdges; j++ ) {
			numEdgeUsers[ abs( trm->polys[i].edges[j] ) ]++;
		}
	}
	for ( i = 1; i <= trm->numEdges; i++ ) {
		if ( numEdgeUsers[i] != 2 ) {
			ii.Com_Printf( "TrmFromModel2: model %s has dangling edges, the model has to be an enclosed hull.\n", model->name );
			PrintModelInfo( model );
			return qfalse;
		}
	}

	// assume convex
	trm->isConvex = qtrue;
	// check if really convex
	for ( i = 0; i < trm->numPolys; i++ ) {
		// to be convex no vertices should be in front of any polygon plane
		for ( j = 0; j < trm->numVerts; j++ ) {
			if ( DotProduct( trm->polys[ i ].normal, trm->verts[ j ] ) - trm->polys[ i ].dist > 0.01f ) {
				trm->isConvex = qfalse;
				break;
			}
		}
		if ( j < trm->numVerts ) {
			break;
		}
	}

	// offset to center of model
	VectorAdd( trm->bounds[0], trm->bounds[1], trm->offset );
	VectorScale( trm->offset, 0.5f, trm->offset );

	TraceModelGenerateEdgeNormals( trm );

	return qtrue;
}

/*
==================
CM_TrmFromModel
==================
*/
qboolean CM_TrmFromModel( const char *modelName, traceModel_t *trm ) {
	cmHandle_t handle;

	handle = CM_LoadModel( modelName, qfalse );
	if ( !handle ) {
		ii.Com_Printf( "CM_TrmFromModel: model %s not found.\n", modelName );
		return qfalse;
	}

	return TrmFromModel2( cmLocal.models[ handle ], trm );
}
