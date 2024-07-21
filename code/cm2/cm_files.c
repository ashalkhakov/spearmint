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

===============================================================================
*/

#include "cm_local.h"
#include "../idlib/l_script.h"
#include "../idlib/l_precomp.h"

#define CM_FILE_EXT			"cm"
#define CM_FILEID			"CM"
#define CM_FILEVERSION		"1.00"


/*
===============================================================================

Writing of collision model file

===============================================================================
*/

void CM_GetNodeBounds( vec3_t bounds[2], cm_node_t *node );
int CM_GetNodeContents( cm_node_t *node );


/*
================
CM_WriteNodes
================
*/
void CM_WriteNodes( fileHandle_t fp, cm_node_t *node ) {
	FS_WriteFloatString( fp, "\t( %d %f )\n", node->planeType, node->planeDist );
	if ( node->planeType != -1 ) {
		CM_WriteNodes( fp, node->children[0] );
		CM_WriteNodes( fp, node->children[1] );
	}
}

/*
================
CM_CountPolygonMemory
================
*/
int CM_CountPolygonMemory( cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;
	int memory;

	memory = 0;
	for ( pref = node->polygons; pref; pref = pref->next ) {
		p = pref->p;
		if ( p->checkcount == cmLocal.checkCount ) {
			continue;
		}
		p->checkcount = cmLocal.checkCount;

		memory += sizeof( cm_polygon_t ) + ( p->numEdges - 1 ) * sizeof( p->edges[0] );
	}
	if ( node->planeType != -1 ) {
		memory += CM_CountPolygonMemory( node->children[0] );
		memory += CM_CountPolygonMemory( node->children[1] );
	}
	return memory;
}

/*
================
CM_WritePolygons
================
*/
void CM_WritePolygons( fileHandle_t fp, cm_node_t *node ) {
	cm_polygonRef_t *pref;
	cm_polygon_t *p;
	int i;
    char        materialName[MAX_QPATH];

	for ( pref = node->polygons; pref; pref = pref->next ) {
		p = pref->p;
		if ( p->checkcount == cmLocal.checkCount ) {
			continue;
		}
		p->checkcount = cmLocal.checkCount;
		FS_WriteFloatString( fp, "\t%d (", p->numEdges );
		for ( i = 0; i < p->numEdges; i++ ) {
			FS_WriteFloatString( fp, " %d", p->edges[i] );
		}
		FS_WriteFloatString( fp, " ) ( %f %f %f ) %f", p->plane[0], p->plane[1], p->plane[2], p->plane[3] );
		FS_WriteFloatString( fp, " ( %f %f %f )", p->bounds[0][0], p->bounds[0][1], p->bounds[0][2] );
		FS_WriteFloatString( fp, " ( %f %f %f )", p->bounds[1][0], p->bounds[1][1], p->bounds[1][2] );
        CM_GetMaterialName( p->material, materialName, sizeof( materialName ) );
		FS_WriteFloatString( fp, " \"%s\"\n", materialName );
	}
	if ( node->planeType != -1 ) {
		CM_WritePolygons( fp, node->children[0] );
		CM_WritePolygons( fp, node->children[1] );
	}
}

/*
================
CM_CountBrushMemory
================
*/
int CM_CountBrushMemory( cm_node_t *node ) {
	cm_brushRef_t *bref;
	cm_brush_t *b;
	int memory;

	memory = 0;
	for ( bref = node->brushes; bref; bref = bref->next ) {
		b = bref->b;
		if ( b->checkcount == cmLocal.checkCount ) {
			continue;
		}
		b->checkcount = cmLocal.checkCount;

		memory += sizeof( cm_brush_t ) + ( b->numPlanes - 1 ) * sizeof( b->planes[0] );
	}
	if ( node->planeType != -1 ) {
		memory += CM_CountBrushMemory( node->children[0] );
		memory += CM_CountBrushMemory( node->children[1] );
	}
	return memory;
}

/*
================
CM_WriteBrushes
================
*/
void CM_WriteBrushes( fileHandle_t fp, cm_node_t *node ) {
	cm_brushRef_t *bref;
	cm_brush_t *b;
	int i;

	for ( bref = node->brushes; bref; bref = bref->next ) {
		b = bref->b;
		if ( b->checkcount == cmLocal.checkCount ) {
			continue;
		}
		b->checkcount = cmLocal.checkCount;
		FS_WriteFloatString( fp, "\t%d {\n", b->numPlanes );
		for ( i = 0; i < b->numPlanes; i++ ) {
			FS_WriteFloatString( fp, "\t\t( %f %f %f ) %f\n", b->planes[i][0], b->planes[i][1], b->planes[i][2], b->planes[i][3] );
		}
		FS_WriteFloatString( fp, "\t} ( %f %f %f )", b->bounds[0][0], b->bounds[0][1], b->bounds[0][2] );
		FS_WriteFloatString( fp, " ( %f %f %f ) \"%s\"\n", b->bounds[1][0], b->bounds[1][1], b->bounds[1][2], CM_StringFromContents( b->contents ) );
	}
	if ( node->planeType != -1 ) {
		CM_WriteBrushes( fp, node->children[0] );
		CM_WriteBrushes( fp, node->children[1] );
	}
}

/*
================
CM_WriteCollisionModel
================
*/
void CM_WriteCollisionModel( fileHandle_t fp, cm_model_t *model ) {
	int i, polygonMemory, brushMemory;

	FS_WriteFloatString( fp, "collisionModel \"%s\" {\n", model->name );
	// vertices
	FS_WriteFloatString( fp, "\tvertices { /* numVertices = */ %d\n", model->numVertices );
	for ( i = 0; i < model->numVertices; i++ ) {
		FS_WriteFloatString( fp, "\t/* %d */ ( %f %f %f )\n", i, model->vertices[i].p[0], model->vertices[i].p[1], model->vertices[i].p[2] );
	}
	FS_WriteFloatString( fp, "\t}\n" );
	// edges
	FS_WriteFloatString( fp, "\tedges { /* numEdges = */ %d\n", model->numEdges );
	for ( i = 0; i < model->numEdges; i++ ) {
		FS_WriteFloatString( fp, "\t/* %d */ ( %d %d ) %d %d\n", i, model->edges[i].vertexNum[0], model->edges[i].vertexNum[1], model->edges[i].internal, model->edges[i].numUsers );
	}
	FS_WriteFloatString( fp, "\t}\n" );
	// nodes
	FS_WriteFloatString( fp, "\tnodes {\n" );
	CM_WriteNodes( fp, model->node );
	FS_WriteFloatString( fp, "\t}\n" );
	// polygons
	cmLocal.checkCount++;
	polygonMemory = CM_CountPolygonMemory( model->node );
	FS_WriteFloatString( fp, "\tpolygons /* polygonMemory = */ %d {\n", polygonMemory );
	cmLocal.checkCount++;
	CM_WritePolygons( fp, model->node );
	FS_WriteFloatString( fp, "\t}\n" );
	// brushes
	cmLocal.checkCount++;
	brushMemory = CM_CountBrushMemory( model->node );
	FS_WriteFloatString( fp, "\tbrushes /* brushMemory = */ %d {\n", brushMemory );
	cmLocal.checkCount++;
	CM_WriteBrushes( fp, model->node );
	FS_WriteFloatString( fp, "\t}\n" );
	// closing brace
	FS_WriteFloatString( fp, "}\n" );
}

/*
================
CM_WriteCollisionModelsToFile
================
*/
void CM_WriteCollisionModelsToFile( const char *filename, int firstModel, int lastModel, unsigned int mapFileCRC ) {
	int i;
	fileHandle_t fp;
	char name[MAX_QPATH];

    Q_strncpyz( name, filename, sizeof( name ) );
    COM_SetExtension( name, sizeof( name ), "." CM_FILE_EXT );

    ii.Com_Printf( PRINT_ALL, "writing %s\n", name );
	// _D3XP was saving to fs_cdpath
	if ( ii.FS_FOpenFile( name, &fp, FS_WRITE ) ) {
		ii.Com_DPrintf( "CM_WriteCollisionModelsToFile: Error opening file %s\n", name );
		return;
	}

	// write file id and version
	FS_WriteFloatString( fp, "%s \"%s\"\n\n", CM_FILEID, CM_FILEVERSION );
	// write the map file crc
	FS_WriteFloatString( fp, "%u\n\n", mapFileCRC );

	// write the collision models
	for ( i = firstModel; i < lastModel; i++ ) {
		CM_WriteCollisionModel( fp, cmLocal.models[ i ] );
	}

	ii.FS_FCloseFile( fp );
}

/*
================
CM_WriteCollisionModelForMapEntity
================
*/
qboolean CM_WriteCollisionModelForMapEntity( const mapEntity_t *mapEnt, const char *filename, const qboolean testTraceModel ) {
	fileHandle_t fp;
	char name[MAX_QPATH];
	cm_model_t *model;

	SetupHash();
	model = CM_CollisionModelForMapEntity( mapEnt );
    Q_strncpyz( model->name, filename, sizeof( model->name ) );

    Q_strncpyz( name, filename, sizeof( name ) );
    COM_SetExtension( name, sizeof( name ), CM_FILE_EXT );


	ii.Com_Printf( "writing %s\n", name);
	// TODO: write to fs_devpath
	if ( ii.FS_FOpenFile( name, &fp, FS_WRITE ) ) {
		ii.Com_DPrintf( "CM_WriteCollisionModelForMapEntity: Error opening file %s\n", name );
		CM_FreeModel( model );
		return qfalse;
	}

	// write file id and version
	FS_WriteFloatString( fp, "%s \"%s\"\n\n", CM_FILEID, CM_FILEVERSION );
	// write the map file crc
	FS_WriteFloatString( fp, "%u\n\n", 0 );

	// write the collision model
	CM_WriteCollisionModel( fp, model );

	ii.FS_FCloseFile( fp );

	if ( testTraceModel ) {
		traceModel_t trm;
		CM_TrmFromModel( model->name, &trm );
	}

	CM_FreeModel( model );

	return qtrue;
}


/*
===============================================================================

Loading of collision model file

===============================================================================
*/

/*
================
CM_ParseVertices
================
*/
void CM_ParseVertices( source_t *src, cm_model_t *model ) {
	int i;

    PC_ExpectTokenString( src, "{" );
	model->numVertices = PC_ParseInt( src );
	model->maxVertices = model->numVertices;
	model->vertices = (cm_vertex_t *) ii.GetMemory( model->maxVertices * sizeof( cm_vertex_t ) );
	for ( i = 0; i < model->numVertices; i++ ) {
		PC_Parse1DMatrix( src, 3, model->vertices[i].p );
		model->vertices[i].side = 0;
		model->vertices[i].sideSet = 0;
		model->vertices[i].checkcount = 0;
	}
	PC_ExpectTokenString( src, "}" );
}

/*
================
CM_ParseEdges
================
*/
void CM_ParseEdges( source_t *src, cm_model_t *model ) {
	int i;

	PC_ExpectTokenString( src, "{" );
	model->numEdges = PC_ParseInt( src );
	model->maxEdges = model->numEdges;
	model->edges = (cm_edge_t *) ii.GetMemory( model->maxEdges * sizeof( cm_edge_t ) );
	for ( i = 0; i < model->numEdges; i++ ) {
		PC_ExpectTokenString( src, "(" );
		model->edges[i].vertexNum[0] = PC_ParseInt( src );
		model->edges[i].vertexNum[1] = PC_ParseInt( src );
		PC_ExpectTokenString( src, ")" );
		model->edges[i].side = 0;
		model->edges[i].sideSet = 0;
		model->edges[i].internal = PC_ParseInt( src );
		model->edges[i].numUsers = PC_ParseInt( src );
		VectorCopy( vec3_origin, model->edges[i].normal );
		model->edges[i].checkcount = 0;
		model->numInternalEdges += model->edges[i].internal;
	}
	PC_ExpectTokenString( src, "}" );
}

/*
================
CM_ParseNodes
================
*/
cm_node_t *CM_ParseNodes( source_t *src, cm_model_t *model, cm_node_t *parent ) {
	cm_node_t *node;

	model->numNodes++;
	node = AllocNode( model, model->numNodes < NODE_BLOCK_SIZE_SMALL ? NODE_BLOCK_SIZE_SMALL : NODE_BLOCK_SIZE_LARGE );
	node->brushes = NULL;
	node->polygons = NULL;
	node->parent = parent;
	PC_ExpectTokenString( src, "(" );
	node->planeType = PC_ParseInt( src );
	node->planeDist = PC_ParseFloat( src, NULL );
	PC_ExpectTokenString( src, ")" );
	if ( node->planeType != -1 ) {
		node->children[0] = CM_ParseNodes( src, model, node );
		node->children[1] = CM_ParseNodes( src, model, node );
	}
	return node;
}

/*
================
CM_ParsePolygons
================
*/
void CM_ParsePolygons( source_t *src, cm_model_t *model ) {
	cm_polygon_t *p;
	int i, numEdges;
	vec3_t normal;
	token_t token;

	if ( PC_ExpectTokenType( src, TT_NUMBER, 0, &token ) ) {
		model->polygonBlock = (cm_polygonBlock_t *) ii.GetMemory( sizeof( cm_polygonBlock_t ) + token.intvalue );
		model->polygonBlock->bytesRemaining = token.intvalue;
		model->polygonBlock->next = ( (byte *) model->polygonBlock ) + sizeof( cm_polygonBlock_t );
	}

	PC_ExpectTokenString( src, "{" );
	while ( !PC_CheckTokenString( src, "}" ) ) {
		// parse polygon
		numEdges = PC_ParseInt( src );
		p = AllocPolygon( model, numEdges );
		p->numEdges = numEdges;
		PC_ExpectTokenString( src, "(" );
		for ( i = 0; i < p->numEdges; i++ ) {
			p->edges[i] = PC_ParseInt( src );
		}
		PC_ExpectTokenString( src, ")" );
		PC_Parse1DMatrix( src, 3, normal );
        VectorCopy( normal, p->plane );
		p->plane[3] = PC_ParseFloat( src, NULL );
		PC_Parse1DMatrix( src, 3, p->bounds[0] );
		PC_Parse1DMatrix( src, 3, p->bounds[1] );
		PC_ExpectTokenType( src, TT_STRING, 0, &token );
		// get material
		p->material = CM_RegisterMaterial( token.string );
		p->contents = CM_GetMaterialContentFlags( p->material );
		p->checkcount = 0;
		// filter polygon into tree
		R_FilterPolygonIntoTree( model, model->node, NULL, p );
	}
}

/*
================
CM_ParseBrushes
================
*/
void CM_ParseBrushes( source_t *src, cm_model_t *model ) {
	cm_brush_t *b;
	int i, numPlanes;
	vec3_t normal;
	token_t token;

	if ( PC_CheckTokenType( src, TT_NUMBER, 0, &token ) ) {
		model->brushBlock = (cm_brushBlock_t *) ii.GetMemory( sizeof( cm_brushBlock_t ) + token.intvalue );
		model->brushBlock->bytesRemaining = token.intvalue;
		model->brushBlock->next = ( (byte *) model->brushBlock ) + sizeof( cm_brushBlock_t );
	}

	PC_ExpectTokenString( src, "{" );
	while ( !PC_CheckTokenString( src, "}" ) ) {
		// parse brush
		numPlanes = PC_ParseInt( src );
		b = AllocBrush( model, numPlanes );
		b->numPlanes = numPlanes;
		PC_ExpectTokenString( src, "{" );
		for ( i = 0; i < b->numPlanes; i++ ) {
			PC_Parse1DMatrix( src, 3, normal );
			VectorCopy( normal, b->planes[i] );
			b->planes[i][3] = PC_ParseFloat( src, NULL );
		}
		PC_ExpectTokenString( src, "}" );
		PC_Parse1DMatrix( src, 3, b->bounds[0] );
		PC_Parse1DMatrix( src, 3, b->bounds[1] );
		PC_ReadToken( src, &token );
		if ( token.type == TT_NUMBER ) {
			b->contents = token.intvalue;		// old .cm files use a single integer
		} else {
			b->contents = CM_ContentsFromString( token.string );
		}
		b->checkcount = 0;
		b->primitiveNum = 0;
		// filter brush into tree
		R_FilterBrushIntoTree( model, model->node, NULL, b );
	}
}

/*
================
CM_ParseCollisionModel
================
*/
qboolean CM_ParseCollisionModel( source_t *src ) {
	cm_model_t *model;
	token_t token;

	if ( cmLocal.numModels >= MAX_SUBMODELS ) {
        ii.Com_Error( ERR_DROP, "LoadModel: no free slots" );
		return qfalse;
	}
	model = CM_AllocModel();
	cmLocal.models[cmLocal.numModels] = model;
	cmLocal.numModels++;
	// parse the file
	PC_ExpectTokenType( src, TT_STRING, 0, &token );
	Q_strncpyz( model->name, token.string, sizeof( model->name ) );
	PC_ExpectTokenString( src, "{" );
	while ( !PC_CheckTokenString( src, "}" ) ) {

		PC_ReadToken( src, &token );

		if ( !strcmp( token.string, "vertices" ) ) {
			CM_ParseVertices( src, model );
			continue;
		}

		if ( !strcmp( token.string, "edges" ) ) {
			CM_ParseEdges( src, model );
			continue;
		}

		if ( !strcmp( token.string, "nodes" ) ) {
			PC_ExpectTokenString( src, "{" );
			model->node = CM_ParseNodes( src, model, NULL );
			PC_ExpectTokenString( src, "}" );
			continue;
		}

		if ( !strcmp( token.string, "polygons" ) ) {
			CM_ParsePolygons( src, model );
			continue;
		}

		if ( !strcmp( token.string, "brushes" ) ) {
			CM_ParseBrushes( src, model );
			continue;
		}

		ii.Com_Error( ERR_DROP, "ParseCollisionModel: bad token \"%s\"", token.string );
	}
	// calculate edge normals
	cmLocal.checkCount++;
	CM_CalculateEdgeNormals( model, model->node );
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

	return qtrue;
}

/*
================
TokenGetUnsignedLongValue
================
*/
unsigned long TokenGetUnsignedLongValue( token_t *token ) {
    unsigned long intvalue;
	const char *p;

	assert( token->type == TT_NUMBER );
	p = token->string;
	intvalue = 0;
	if ( token->subtype & TT_DECIMAL ) {
		while( *p ) {
			intvalue = intvalue * 10 + (*p - '0');
			p++;
		}
	}
	else if ( token->subtype & TT_OCTAL ) {
		// step over the first zero
		p += 1;
		while( *p ) {
			intvalue = (intvalue << 3) + (*p - '0');
			p++;
		}
	}
	else if ( token->subtype & TT_HEX ) {
		// step over the leading 0x or 0X
		p += 2;
		while( *p ) {
			intvalue <<= 4;
			if (*p >= 'a' && *p <= 'f')
				intvalue += *p - 'a' + 10;
			else if (*p >= 'A' && *p <= 'F')
				intvalue += *p - 'A' + 10;
			else
				intvalue += *p - '0';
			p++;
		}
	}
	else if ( token->subtype & TT_BINARY ) {
		// step over the leading 0b or 0B
		p += 2;
		while( *p ) {
			intvalue = (intvalue << 1) + (*p - '0');
			p++;
		}
	}
    return intvalue;
}

/*
================
CM_LoadCollisionModelFile
================
*/
qboolean CM_LoadCollisionModelFile( const char *name, unsigned int mapFileCRC ) {
	char fileName[MAX_QPATH];
	token_t token;
	source_t *src;
	unsigned int crc;

	// load it
    Q_strncpyz( fileName, name, sizeof( fileName ) );
    COM_SetExtension( fileName, sizeof( fileName ), "." CM_FILE_EXT );
	src = PC_LoadSource( fileName, NULL, NULL );
	if ( !src ) {
		return qfalse;
	}

	//src->SetFlags( LEXFL_NOSTRINGCONCAT | LEXFL_NODOLLARPRECOMPILE );
	SetScriptFlags( src->scriptstack, SCFL_NOSTRINGWHITESPACES );

	if ( !PC_ExpectTokenString( src, CM_FILEID ) ) {
		ii.Com_DPrintf( "%s is not an CM file.", fileName );
        PC_FreeSource( src );
		return qfalse;
	}

	if ( !PC_ReadToken( src, &token ) || strcmp( token.string, CM_FILEVERSION ) ) {
		ii.Com_DPrintf( "%s has version %s instead of %s", fileName, token.string, CM_FILEVERSION );
		PC_FreeSource( src );
		return qfalse;
	}

	if ( !PC_ExpectTokenType( src, TT_NUMBER, TT_INTEGER, &token ) ) {
		ii.Com_Printf( "%s has no map file CRC", fileName );
		PC_FreeSource( src );
		return qfalse;
	}

	crc = TokenGetUnsignedLongValue( &token );
	if ( mapFileCRC && crc != mapFileCRC ) {
		ii.Com_Printf( "%s is out of date\n", fileName );
		PC_FreeSource( src );
		return qfalse;
	}

	// parse the file
	while ( 1 ) {
		if ( !PC_ReadToken( src, &token ) ) {
			break;
		}

		if ( !strcmp( token.string, "collisionModel" ) ) {
			if ( !CM_ParseCollisionModel( src ) ) {
				PC_FreeSource( src );
				return qfalse;
			}
			continue;
		}

		ii.Com_Error( ERR_DROP, "CM_LoadCollisionModelFile: bad token \"%s\"", token.string );
	}

	PC_FreeSource( src );

	return qtrue;
}
