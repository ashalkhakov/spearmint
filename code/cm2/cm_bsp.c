/*
===========================================================================
Copyright (C) 1999-2010 id Software LLC, a ZeniMax Media company.

This file is part of Spearmint Source Code.

Spearmint Source Code is free software; you can redistribute it
and/or modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 3 of the License,
or (at your option) any later version.

Spearmint Source Code is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Spearmint Source Code.  If not, see <http://www.gnu.org/licenses/>.

In addition, Spearmint Source Code is also subject to certain additional terms.
You should have received a copy of these additional terms immediately following
the terms and conditions of the GNU General Public License.  If not, please
request a copy in writing from id Software at the address below.

If you have questions concerning this license or the applicable additional
terms, you may contact in writing id Software LLC, c/o ZeniMax Media Inc.,
Suite 120, Rockville, Maryland 20850 USA.
===========================================================================
*/

#include "cm_local.h"
#include "../idlib/l_script.h"
#include "../idlib/l_precomp.h"
#include "../idlib/bsp/bsp.h"

// to allow boxes to be treated as brush models, we allocate
// some extra indexes along with those needed by the map
#define BOX_LEAF_BRUSHES	1
#define	BOX_BRUSHES		1
#define	BOX_SIDES		6
#define	BOX_LEAFS		2
#define	BOX_PLANES		12

#define MAX_PATCH_VERTS 1024

#if 1 // ZTM: FIXME: BSP is already swapped by BSP_Load, but removing these probably makes merging ioq3 changes harder...
#undef LittleShort
#define LittleShort
#undef LittleLong
#define LittleLong
#undef LittleFloat
#define LittleFloat
#define LL(x) /* nothing */
#else
#define	LL(x) x=LittleLong(x)
#endif


int			c_pointcontents;
int			c_traces, c_brush_traces, c_patch_traces;

int			capsule_contents;

void CM_FloodAreaConnections( void );

/*
===============================================================================

					MAP LOADING

===============================================================================
*/

/*
=================
CM_RegisterMaterial
=================
*/
qhandle_t CM_RegisterMaterial( const char *name ) {
    int i, hash;

    if ( !name || !name[0] ) {
        ii.Com_DPrintf( "CM_RegisterMaterial: NULL name" );
        return 0;
    }

    hash = HashIndexGenerateKeyForString( &cmLocal.materialsHash, name, qfalse );
    for ( i = HashIndexFirst( &cmLocal.materialsHash, hash ); i != -1; i = HashIndexNext( &cmLocal.materialsHash, i ) ) {
        if ( Q_stricmp( cmLocal.materials[i]->name, name ) == 0 ) {
            return i + 1;
        }
    }

    return 0;
}

/*
=================
CM_GetMaterialName
=================
*/
void CM_GetMaterialName( qhandle_t hShader, char *buffer, int bufferSize ) {
    if ( hShader < 1 || hShader >= cmLocal.numMaterials + 1 || !cmLocal.materials[hShader - 1] ) {
        buffer[0] = 0;
        return;
    }

    Q_strncpyz( buffer, cmLocal.materials[hShader - 1]->name, bufferSize );
}

/*
=================
CM_GetMaterialContentFlags
=================
*/
int CM_GetMaterialContentFlags( qhandle_t material ) {
    if ( material < 1 || material >= cmLocal.numMaterials + 1 || !cmLocal.materials[material - 1] ) {
        return 0;
    }

    return cmLocal.materials[material - 1]->contentFlags;
}

/*
=================
CM_GetMaterialSurfaceFlags
=================
*/
int CM_GetMaterialSurfaceFlags( qhandle_t material ) {
    if ( material < 1 || material >= cmLocal.numMaterials + 1 || !cmLocal.materials[material - 1] ) {
        return 0;
    }

    return cmLocal.materials[material - 1]->surfaceFlags;
}

/*
=================
CM_MaterialNeedsBackSide
=================
*/
qboolean CM_MaterialNeedsBackSide( qhandle_t material ) {
    //if ( material->GetCullType() == CT_TWO_SIDED || material->ShouldCreateBackSides() ) {
    return qfalse;
}

/*
=================
CMod_LoadShaders
=================
*/
void CMod_LoadShaders( void ) {
    int             i, key;
    cm_material_t   *material;
    dshader_t       *in;

	if (cm_bsp->numShaders < 1) {
		ii.Com_Error (ERR_DROP, "Map with no shaders");
	}

	cmLocal.maxMaterials = MAX_MATERIALS;
	cmLocal.numMaterials = 0;
	cmLocal.materials = (cm_material_t **) ii.GetMemory( (cmLocal.maxMaterials+1) * sizeof(cm_material_t *) );
	memset( cmLocal.materials, 0, (cmLocal.maxMaterials+1) * sizeof(cm_material_t *) );

    material = CM_AllocMaterial();

    Q_strncpyz( material->name, "_tracemodel", sizeof( material->name ) );
    material->contentFlags = 0;
    material->surfaceFlags = SURF_COLLISION;
    key = HashIndexGenerateKeyForString( &cmLocal.materialsHash, material->name, qfalse );
    HashIndexAdd( &cmLocal.materialsHash, key, cmLocal.numMaterials );

    cmLocal.materials[cmLocal.numMaterials] = material;
    cmLocal.numMaterials++;

    for ( i = 0; i < cm_bsp->numShaders; i++ ) {
        if ( cmLocal.numMaterials >= MAX_MATERIALS ) {
            ii.Com_Error( ERR_DROP, "CMod_LoadShaders: more than %d collision materials", MAX_MATERIALS );
            break;
        }

        in = &cm_bsp->shaders[i];
        material = CM_AllocMaterial();

        Q_strncpyz( material->name, in->shader, sizeof( material->name ) );
        material->contentFlags = in->contentFlags;
        material->surfaceFlags = in->surfaceFlags;
        key = HashIndexGenerateKeyForString( &cmLocal.materialsHash, material->name, qfalse );
        HashIndexAdd( &cmLocal.materialsHash, key, cmLocal.numMaterials );

        cmLocal.materials[cmLocal.numMaterials] = material;
        cmLocal.numMaterials++;
    }
}

// TODO: Remove?
#if 0
/*
=================
CMod_LoadSubmodels
=================
*/
void CMod_LoadSubmodels( void ) {
	dmodel_t	*in;
	cmodel_t	*out;
	int			i, j, count;
	int			*indexes;

	in = cm_bsp->submodels;
	count = cm_bsp->numSubmodels;

	if (count < 1)
		Com_Error (ERR_DROP, "Map with no models");
	cm.cmodels = Hunk_Alloc( count * sizeof( *cm.cmodels ), h_high );
	cm.numSubModels = count;

	for ( i=0 ; i<count ; i++, in++)
	{
		out = &cm.cmodels[i];

		for (j=0 ; j<3 ; j++)
		{	// spread the mins / maxs by a pixel
			out->mins[j] = LittleFloat (in->mins[j]) - 1;
			out->maxs[j] = LittleFloat (in->maxs[j]) + 1;
		}

		if ( i == 0 ) {
			continue;	// world model doesn't need other info
		}

		// make a "leaf" just to hold the model's brushes and surfaces
		out->leaf.numLeafBrushes = LittleLong( in->numBrushes );
		indexes = Hunk_Alloc( out->leaf.numLeafBrushes * 4, h_high );
		out->leaf.firstLeafBrush = indexes - cm.leafbrushes;
		for ( j = 0 ; j < out->leaf.numLeafBrushes ; j++ ) {
			indexes[j] = LittleLong( in->firstBrush ) + j;
		}

		out->leaf.numLeafSurfaces = LittleLong( in->numSurfaces );
		indexes = Hunk_Alloc( out->leaf.numLeafSurfaces * 4, h_high );
		out->leaf.firstLeafSurface = indexes - cm.leafsurfaces;
		for ( j = 0 ; j < out->leaf.numLeafSurfaces ; j++ ) {
			indexes[j] = LittleLong( in->firstSurface ) + j;
		}
	}
}
#endif


/*
=================
CMod_LoadNodes

=================
*/
void CMod_LoadNodes( void ) {
	dnode_t		*in;
	int			child;
	cNode_t		*out;
	int			i, j, count;
	
	in = cm_bsp->nodes;
	count = cm_bsp->numNodes;

	if (count < 1)
		ii.Com_Error (ERR_DROP, "Map has no nodes");
	cmLocal.nodes = ii.HunkAlloc( count * sizeof( *cmLocal.nodes ) );
	cmLocal.numNodes = count;

	out = cmLocal.nodes;

	for (i=0 ; i<count ; i++, out++, in++)
	{
		out->plane = cmLocal.planes + LittleLong( in->planeNum );
		for (j=0 ; j<2 ; j++)
		{
			child = LittleLong (in->children[j]);
			out->children[j] = child;
		}
	}

}

#if 0 // TODO: Remove?
/*
=================
CM_BoundBrush

=================
*/
void CM_BoundBrush( cbrush_t *b ) {
	b->bounds[0][0] = -b->sides[0].plane->dist;
	b->bounds[1][0] = b->sides[1].plane->dist;

	b->bounds[0][1] = -b->sides[2].plane->dist;
	b->bounds[1][1] = b->sides[3].plane->dist;

	b->bounds[0][2] = -b->sides[4].plane->dist;
	b->bounds[1][2] = b->sides[5].plane->dist;
}

/*
=================
CMod_LoadBrushes

=================
*/
void CMod_LoadBrushes( void ) {
	dbrush_t	*in;
	cbrush_t	*out;
	int			i, count;

	in = cm_bsp->brushes;
	count = cm_bsp->numBrushes;

	cm.brushes = Hunk_Alloc( ( BOX_BRUSHES + count ) * sizeof( *cm.brushes ), h_high );
	cm.numBrushes = count;

	out = cm.brushes;

	for ( i=0 ; i<count ; i++, out++, in++ ) {
		out->sides = cm.brushsides + LittleLong(in->firstSide);
		out->numsides = LittleLong(in->numSides);

		out->shaderNum = LittleLong( in->shaderNum );
		if ( out->shaderNum < 0 || out->shaderNum >= cm.numShaders ) {
			Com_Error( ERR_DROP, "CMod_LoadBrushes: bad shaderNum: %i", out->shaderNum );
		}
		out->contents = cm.shaders[out->shaderNum].contentFlags;

		CM_BoundBrush( out );
	}

}
#endif

/*
=================
CMod_LoadLeafs
=================
*/
void CMod_LoadLeafs( void )
{
	int			i;
	cLeaf_t		*out;
	dleaf_t 	*in;
	int			count;
	
	in = cm_bsp->leafs;
	count = cm_bsp->numLeafs;

	if (count < 1)
		ii.Com_Error (ERR_DROP, "Map with no leafs");

	cmLocal.leafs = ii.HunkAlloc( ( BOX_LEAFS + count ) * sizeof( *cmLocal.leafs ) );
	cmLocal.numLeafs = count;

	out = cmLocal.leafs;	
	for ( i=0 ; i<count ; i++, in++, out++)
	{
		out->cluster = LittleLong (in->cluster);
		out->area = LittleLong (in->area);
		out->firstLeafBrush = LittleLong (in->firstLeafBrush);
		out->numLeafBrushes = LittleLong (in->numLeafBrushes);
		out->firstLeafSurface = LittleLong (in->firstLeafSurface);
		out->numLeafSurfaces = LittleLong (in->numLeafSurfaces);

		if (out->cluster >= cmLocal.numClusters)
			cmLocal.numClusters = out->cluster + 1;
		if (out->area >= cmLocal.numAreas)
			cmLocal.numAreas = out->area + 1;
	}

	cmLocal.areas = ii.HunkAlloc( cmLocal.numAreas * sizeof( *cmLocal.areas ) );
	Com_Memset( cmLocal.areas, 0, cmLocal.numAreas * sizeof( *cmLocal.areas ) );
	cmLocal.areaPortals = ii.HunkAlloc( cmLocal.numAreas * cmLocal.numAreas * sizeof( *cmLocal.areaPortals ) );
	Com_Memset( cmLocal.areaPortals, 0, cmLocal.numAreas * cmLocal.numAreas * sizeof( *cmLocal.areaPortals ) );
}

/*
=================
CMod_LoadPlanes
=================
*/
void CMod_LoadPlanes( void )
{
	int			i, j;
	cplane_t	*out;
	dplane_t 	*in;
	int			count;
	int			bits;
	
	in = cm_bsp->planes;
	count = cm_bsp->numPlanes;

	if (count < 1)
		Com_Error (ERR_DROP, "Map with no planes");
	cmLocal.planes = ii.HunkAlloc( ( BOX_PLANES + count ) * sizeof( *cmLocal.planes ) );
	cmLocal.numPlanes = count;

	out = cmLocal.planes;	

	for ( i=0 ; i<count ; i++, in++, out++)
	{
		bits = 0;
		for (j=0 ; j<3 ; j++)
		{
			out->normal[j] = LittleFloat (in->normal[j]);
			if (out->normal[j] < 0)
				bits |= 1<<j;
		}

		out->dist = LittleFloat (in->dist);
		out->type = PlaneTypeForNormal( out->normal );
		out->signbits = bits;
	}
}

#if 0 // TODO: REMOVE?
/*
=================
CMod_LoadLeafBrushes
=================
*/
void CMod_LoadLeafBrushes( void )
{
	int			i;
	int			*out;
	int		 	*in;
	int			count;
	
	in = cm_bsp->leafBrushes;
	count = cm_bsp->numLeafBrushes;

	cm.leafbrushes = Hunk_Alloc( (BOX_LEAF_BRUSHES + count) * sizeof( *cm.leafbrushes ), h_high );
	cm.numLeafBrushes = count;

	out = cm.leafbrushes;

	for ( i=0 ; i<count ; i++, in++, out++) {
		*out = LittleLong (*in);
	}
}

/*
=================
CMod_LoadLeafSurfaces
=================
*/
void CMod_LoadLeafSurfaces( void )
{
	int			i;
	int			*out;
	int		 	*in;
	int			count;
	
	in = cm_bsp->leafSurfaces;
	count = cm_bsp->numLeafSurfaces;

	cm.leafsurfaces = Hunk_Alloc( count * sizeof( *cm.leafsurfaces ), h_high );
	cm.numLeafSurfaces = count;

	out = cm.leafsurfaces;

	for ( i=0 ; i<count ; i++, in++, out++) {
		*out = LittleLong (*in);
	}
}

/*
=================
CMod_LoadBrushSides
=================
*/
void CMod_LoadBrushSides ( void )
{
	int				i;
	cbrushside_t	*out;
	dbrushside_t 	*in;
	int				count;
	int				num;

	in = cm_bsp->brushSides;
	count = cm_bsp->numBrushSides;

	cm.brushsides = Hunk_Alloc( ( BOX_SIDES + count ) * sizeof( *cm.brushsides ), h_high );
	cm.numBrushSides = count;

	out = cm.brushsides;	

	for ( i=0 ; i<count ; i++, in++, out++) {
		num = LittleLong( in->planeNum );
		out->planeNum = num;
		out->plane = &cm.planes[num];
		out->shaderNum = LittleLong( in->shaderNum );
		if ( out->shaderNum < 0 || out->shaderNum >= cm.numShaders ) {
			Com_Error( ERR_DROP, "CMod_LoadBrushSides: bad shaderNum: %i", out->shaderNum );
		}
		out->surfaceFlags = cm.shaders[out->shaderNum].surfaceFlags;
		out->surfaceNum = LittleLong( in->surfaceNum );
	}
}
#endif

/*
=================
CMod_LoadEntityString
=================
*/
void CMod_LoadEntityString( void ) {
	cmLocal.entityString = cm_bsp->entityString;
	cmLocal.numEntityChars = cm_bsp->entityStringLength;
}

/*
=================
CMod_LoadVisibility
=================
*/
void CMod_LoadVisibility( void ) {
	if ( !cm_bsp->visibilityLength ) {
		cmLocal.clusterBytes = ( cmLocal.numClusters + 31 ) & ~31;
		cmLocal.visibility = ii.HunkAlloc( cmLocal.clusterBytes );
		Com_Memset( cmLocal.visibility, 255, cmLocal.clusterBytes );
		return;
	}

	cmLocal.vised = qtrue;
	cmLocal.visibility = cm_bsp->visibility;
	cmLocal.numClusters = cm_bsp->numClusters;
	cmLocal.clusterBytes = cm_bsp->clusterBytes;
}

//==================================================================

int		CM_NumClusters( void ) {
	return cmLocal.numClusters;
}

int		CM_NumInlineModels( void ) {
	return cmLocal.numModels; // TODO: isn't this slightly wrong?
}

char	*CM_EntityString( void ) {
	return cmLocal.entityString;
}

qboolean CM_GetEntityToken( int *parseOffset, char *token, int size ) {
	const char	*s;
	char	*parsePoint = cmLocal.entityString;

	if ( !cmLocal.entityString || *parseOffset < 0 || *parseOffset >= cmLocal.numEntityChars ) {
		return qfalse;
	}

	parsePoint = cmLocal.entityString + *parseOffset;

	s = COM_Parse( &parsePoint );
	Q_strncpyz( token, s, size );

	if ( !parsePoint && !s[0] ) {
		*parseOffset = 0;
		return qfalse;
	} else {
		*parseOffset = parsePoint - cmLocal.entityString;
		return qtrue;
	}
}

int		CM_LeafCluster( int leafnum ) {
	if (leafnum < 0 || leafnum >= cmLocal.numLeafs) {
		ii.Com_Error (ERR_DROP, "CM_LeafCluster: bad number");
	}
	return cmLocal.leafs[leafnum].cluster;
}

int		CM_LeafArea( int leafnum ) {
	if ( leafnum < 0 || leafnum >= cmLocal.numLeafs ) {
		Com_Error (ERR_DROP, "CM_LeafArea: bad number");
	}
	return cmLocal.leafs[leafnum].area;
}

static mapBrushSide_t *MapBrushSideInitFromBspBrushSide( const bspFile_t *bsp, dbrushside_t *in ) {
    mapBrushSide_t *out;
    dplane_t *plane;
    dshader_t *shader;

    out = ii.GetMemory( sizeof( *out ) );

    MapBrushSideClear( out );

    plane = &bsp->planes[in->planeNum];

    VectorCopy( plane->normal, out->plane );
    PlaneSetDist( out->plane, plane->dist );

    VectorScale( plane->normal, plane->dist, out->origin );

    shader = &bsp->shaders[in->shaderNum];
    Q_strncpyz( out->material, shader->shader, sizeof( out->material ) );

    return out;
}

static mapPrimitive_t *MapBrushInitFromBspBrush( const bspFile_t *bsp, dbrush_t *in ) {
    mapPrimitive_t *out;
    mapBrushSide_t *outBrushSide;
    dbrushside_t *brushSide;
    dshader_t *shader;
    int j;

    out = ii.GetMemory( sizeof( *out ) );

    MapBrushInit( out );

    shader = &bsp->shaders[in->shaderNum];
    Q_strncpyz( out->material, shader->shader, sizeof( out->material ) );

    for ( j = 0; j < in->numSides; j++ ) {
        brushSide = &bsp->brushSides[in->firstSide + j];

        outBrushSide = MapBrushSideInitFromBspBrushSide( bsp, brushSide );

        AddSideToMapBrush( out, outBrushSide );
    }

    return out;
}

static mapPrimitive_t *MapPatchInitFromBspPatch( const bspFile_t *bsp, dsurface_t *in ) {
    mapPrimitive_t *out;
    dshader_t *shader;
    drawVert_t *dv;
    surfVert_t *outVert;
    int c, j;

    assert( in->surfaceType == MST_PATCH );

    out = ( mapPrimitive_t * )ii.GetMemory( sizeof( *out ) );
    MapPatchInitWithSize( out, in->patchWidth, in->patchHeight );
    SurfacePatchSetSize( &out->patch, in->patchWidth, in->patchHeight );

    // load the full drawverts onto the stack
    c = in->patchWidth * in->patchHeight;
    if ( c > MAX_PATCH_VERTS ) {
        ii.Com_Error( ERR_DROP, "ParseMesh: MAX_PATCH_VERTS" );
    }

    dv = bsp->drawVerts + in->firstVert;
    outVert = out->patch.surf.verts;
    for ( j = 0 ; j < c ; j++, dv++, outVert++ ) {
        VectorCopy( dv->xyz, outVert->xyz );
        VectorCopy( dv->normal, outVert->normal );
    }

    shader = &bsp->shaders[in->shaderNum];
    Q_strncpyz( out->material, shader->shader, sizeof( out->material ) );

    out->vertSubdivisions = in->subdivisions;
    out->horzSubdivisions = in->subdivisions;
    out->explicitSubdivisions = qfalse;

    return out;
}

static mapPrimitive_t *MapMeshInitFromBspTerrain( const bspFile_t *bsp, dsurface_t *in ) {
    mapPrimitive_t *out;
    dshader_t *shader;
    drawVert_t *dv;
    surfVert_t *outVert;
    int *indexes;
    int j;

    assert( in->surfaceType == MST_TERRAIN );

    out = ( mapPrimitive_t * )ii.GetMemory( sizeof( *out ) );
    MapMeshInit( out );

    SurfaceResizeIndexes( &out->surf, in->numIndexes );
    dv = bsp->drawVerts + in->firstVert;
    outVert = out->surf.verts;
    for ( j = 0 ; j < in->numVerts ; j++, dv++, outVert++ ) {
        VectorCopy( dv->xyz, outVert->xyz );
        VectorCopy( dv->normal, outVert->normal );
    }

    indexes = bsp->drawIndexes + in->firstIndex;
    SurfaceResizeVerts( &out->surf, in->numVerts );
    for ( j = 0; j < in->numIndexes ; j++ ) {
        out->surf.indexes[j] = indexes[j];

        if ( indexes[j] < 0 || indexes[j] >= in->numVerts ) {
            Com_Error(ERR_DROP, "MapMeshInitFromBspTerrain: Bad index in trisoup surface");
        }
    }

    shader = &bsp->shaders[in->shaderNum];
    Q_strncpyz( out->material, shader->shader, sizeof( out->material ) );

    return out;
}

void TransferSubModelToMapEntity( const bspFile_t *bsp, dmodel_t *model, mapEntity_t *entity, int modelindex ) {
    int i, j;
    dbrush_t *brush;
    dbrushside_t *brushSide;
    dshader_t *shader;
    dplane_t *plane;
    mapPrimitive_t *outBrush;
    mapPrimitive_t *outPatch;
    mapBrushSide_t *outBrushSide;
    char            modelName[MAX_QPATH];

    if ( modelindex > 0 ) {
        sprintf( modelName, "*%d", modelindex );
        DictSet( &entity->epairs, "model", modelName );
    }

    for ( i = model->numBrushes - 1; i >= 0; i-- ) {
        brush = &bsp->brushes[model->firstBrush + i];
        shader = &bsp->shaders[brush->shaderNum];

        outBrush = MapBrushInitFromBspBrush( bsp, brush );

        AddPrimitiveToMapEntity( entity, outBrush );
    }

    for ( i = model->numSurfaces - 1; i >= model->firstSurface; i-- ) {
        dsurface_t *surf = &bsp->surfaces[i];
        if ( surf->surfaceType == MST_PATCH ) {
            outPatch = MapPatchInitFromBspPatch( bsp, surf );

            AddPrimitiveToMapEntity( entity, outPatch );
        } else if ( surf->surfaceType == MST_TERRAIN ) {
            outPatch = MapMeshInitFromBspTerrain( bsp, surf );

            AddPrimitiveToMapEntity( entity, outPatch );
        }

        // TODO: transfer triangle soups as render models
    }
}

/*
==================
CM_LoadBSP

Loads in the map and all submodels
==================
*/
void CM_LoadBSP( const char *filename, int *checksum ) {
    source_t *src;
    mapFile_t out;
    mapEntity_t *entity;
    dmodel_t *model;
    int entityNum;

	static unsigned	last_checksum;

	if ( !filename || !filename[0] ) {
		Com_Error( ERR_DROP, "CM_LoadBSP: NULL name" );
	}

	Com_DPrintf( "CM_LoadBSP( %s )\n", filename );

	if ( cm_bsp && !strcmp( cm_bsp->name, filename ) ) {
        // TODO: why is this client-only in spearmint?
		*checksum = cm_bsp->checksum;
		return;
	}

    CM_Clear();

    cm_bsp = BSP_Load( filename );
	if ( !cm_bsp ) {
		Com_Error( ERR_DROP, "Unable to load %s", filename );
	}

	*checksum = cm_bsp->checksum;

    MapFileInit( &out );
    Q_strncpyz( out.name, filename, sizeof( out.name ) );

    src = LoadSourceMemory( cm_bsp->entityString, cm_bsp->entityStringLength, "*bsp", NULL );

    if ( !ParseMapFileFromSource( &out, src ) ) {
        FreeSource( src );
        MapFileFree( &out );
        return;
    }

    FreeSource( src );
    src = NULL;

    entity = out.entities;
    entityNum = 0;
    while ( entity ) {
        if ( entityNum >= cm_bsp->numSubmodels ) {
            break;
        }
        model = &cm_bsp->submodels[entityNum];

        TransferSubModelToMapEntity( cm_bsp, model, entity, entityNum );

        entity = entity->next;
        entityNum++;
	
	}
	SetMapFileGeometryCRC( &out );

    MapFileWrite( &out, "out", ".map", qtrue );

	// load into heap
	CMod_LoadShaders();
	CMod_LoadLeafs();
	//CMod_LoadLeafBrushes();
	//CMod_LoadLeafSurfaces();
	CMod_LoadPlanes();
	//CMod_LoadBrushSides();
	//CMod_LoadBrushes();
	//CMod_LoadSubmodels();
	CMod_LoadNodes();
	CMod_LoadEntityString();
	CMod_LoadVisibility();

	CM_FloodAreaConnections ();

    CM_LoadMap2( &out, qfalse );

    // TODO: go over the entities one more time and load each surface as a render model...
    // this is for misc_model (compiled into the map)

	MapFileFree( &out );
}

/*
======================================================================

LEAF LISTING

======================================================================
*/

/*
==================
CM_PointLeafnum_r

==================
*/
int CM_PointLeafnum_r( const vec3_t p, int num ) {
	float		d;
	cNode_t		*node;
	cplane_t	*plane;

	while (num >= 0)
	{
		node = cmLocal.nodes + num;
		plane = node->plane;
		
		if (plane->type < 3)
			d = p[plane->type] - plane->dist;
		else
			d = DotProduct (plane->normal, p) - plane->dist;
		if (d < 0)
			num = node->children[1];
		else
			num = node->children[0];
	}

	c_pointcontents++;		// optimize counter

	return -1 - num;
}

int CM_PointLeafnum( const vec3_t p ) {
	if ( !cmLocal.numNodes ) {	// map not loaded
		return 0;
	}
	return CM_PointLeafnum_r (p, 0);
}


void CM_StoreLeafs( leafList_t *ll, int nodenum ) {
	int		leafNum;

	leafNum = -1 - nodenum;

	// store the lastLeaf even if the list is overflowed
	if ( cmLocal.leafs[ leafNum ].cluster != -1 ) {
		ll->lastLeaf = leafNum;
	}

	if ( ll->count >= ll->maxcount) {
		ll->overflowed = qtrue;
		return;
	}
	ll->list[ ll->count++ ] = leafNum;
}

/*
=============
CM_BoxLeafnums

Fills in a list of all the leafs touched
=============
*/
void CM_BoxLeafnums_r( leafList_t *ll, int nodenum ) {
	cplane_t	*plane;
	cNode_t		*node;
	int			s;

	while (1) {
		if (nodenum < 0) {
			ll->storeLeafs( ll, nodenum );
			return;
		}
	
		node = &cmLocal.nodes[nodenum];
		plane = node->plane;
        assert( node->plane >= cmLocal.planes && node->plane < cmLocal.planes + cmLocal.numPlanes );
		s = BoxOnPlaneSide( ll->bounds[0], ll->bounds[1], plane );
		if (s == 1) {
			nodenum = node->children[0];
		} else if (s == 2) {
			nodenum = node->children[1];
		} else {
			// go down both
			CM_BoxLeafnums_r( ll, node->children[0] );
			nodenum = node->children[1];
		}

	}
}

/*
==================
CM_BoxLeafnums
==================
*/
int	CM_BoxLeafnums( const vec3_t mins, const vec3_t maxs, int *list, int listsize, int *lastLeaf) {
	leafList_t	ll;

	//cmLocal.checkcount++;

	VectorCopy( mins, ll.bounds[0] );
	VectorCopy( maxs, ll.bounds[1] );
	ll.count = 0;
	ll.maxcount = listsize;
	ll.list = list;
	ll.storeLeafs = CM_StoreLeafs;
	ll.lastLeaf = 0;
	ll.overflowed = qfalse;

	CM_BoxLeafnums_r( &ll, 0 );

	*lastLeaf = ll.lastLeaf;
	return ll.count;
}

/*
===============================================================================

PVS

===============================================================================
*/

byte	*CM_ClusterPVS (int cluster) {
	if (cluster < 0 || cluster >= cmLocal.numClusters || !cmLocal.vised ) {
		return cmLocal.visibility;
	}

	return cmLocal.visibility + cluster * cmLocal.clusterBytes;
}



/*
===============================================================================

AREAPORTALS

===============================================================================
*/

void CM_FloodArea_r( int areaNum, int floodnum) {
	int		i;
	cmArea_t *area;
	int		*con;

	area = &cmLocal.areas[ areaNum ];

	if ( area->floodvalid == cmLocal.floodvalid ) {
		if (area->floodnum == floodnum)
			return;
		ii.Com_Error (ERR_DROP, "FloodArea_r: reflooded");
	}

	area->floodnum = floodnum;
	area->floodvalid = cmLocal.floodvalid;
	con = cmLocal.areaPortals + areaNum * cmLocal.numAreas;
	for ( i=0 ; i < cmLocal.numAreas  ; i++ ) {
		if ( con[i] > 0 ) {
			CM_FloodArea_r( i, floodnum );
		}
	}
}

/*
====================
CM_FloodAreaConnections

====================
*/
void	CM_FloodAreaConnections( void ) {
	int		i;
	cmArea_t	*area;
	int		floodnum;

	// all current floods are now invalid
	cmLocal.floodvalid++;
	floodnum = 0;

	for (i = 0 ; i < cmLocal.numAreas ; i++) {
		area = &cmLocal.areas[i];
		if (area->floodvalid == cmLocal.floodvalid) {
			continue;		// already flooded into
		}
		floodnum++;
		CM_FloodArea_r (i, floodnum);
	}

}

/*
====================
CM_AdjustAreaPortalState

====================
*/
void	CM_AdjustAreaPortalState( int area1, int area2, qboolean open ) {
	if ( area1 < 0 || area2 < 0 ) {
		return;
	}

	if ( area1 >= cmLocal.numAreas || area2 >= cmLocal.numAreas ) {
		ii.Com_Error (ERR_DROP, "CM_ChangeAreaPortalState: bad area number");
	}

	if ( open ) {
		cmLocal.areaPortals[ area1 * cmLocal.numAreas + area2 ]++;
		cmLocal.areaPortals[ area2 * cmLocal.numAreas + area1 ]++;
	} else {
		cmLocal.areaPortals[ area1 * cmLocal.numAreas + area2 ]--;
		cmLocal.areaPortals[ area2 * cmLocal.numAreas + area1 ]--;
		if ( cmLocal.areaPortals[ area2 * cmLocal.numAreas + area1 ] < 0 ) {
			ii.Com_Error (ERR_DROP, "CM_AdjustAreaPortalState: negative reference count");
		}
	}

	CM_FloodAreaConnections ();
}

/*
====================
CM_AreasConnected

====================
*/
qboolean	CM_AreasConnected( int area1, int area2 ) {
	if ( cm_noAreas->integer ) {
		return qtrue;
	}

	if ( area1 < 0 || area2 < 0 ) {
		return qfalse;
	}

	if (area1 >= cmLocal.numAreas || area2 >= cmLocal.numAreas) {
		ii.Com_Error (ERR_DROP, "area >= cm.numAreas");
	}

	if (cmLocal.areas[area1].floodnum == cmLocal.areas[area2].floodnum) {
		return qtrue;
	}
	return qfalse;
}


/*
=================
CM_WriteAreaBits

Writes a bit vector of all the areas
that are in the same flood as the area parameter
Returns the number of bytes needed to hold all the bits.

The bits are OR'd in, so you can CM_WriteAreaBits from multiple
viewpoints and get the union of all visible areas.

This is used to cull non-visible entities from snapshots
=================
*/
int CM_WriteAreaBits (byte *buffer, int area)
{
	int		i;
	int		floodnum;
	int		bytes;

	bytes = (cmLocal.numAreas+7)>>3;

	if (cm_noAreas->integer || area == -1)
	{	// for debugging, send everything
		Com_Memset (buffer, 255, bytes);
	}
	else
	{
		floodnum = cmLocal.areas[area].floodnum;
		for (i=0 ; i<cmLocal.numAreas ; i++)
		{
			if (cmLocal.areas[i].floodnum == floodnum || area == -1)
				buffer[i>>3] |= 1<<(i&7);
		}
	}

	return bytes;
}
