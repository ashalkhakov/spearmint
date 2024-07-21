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

#ifndef __MAPFILE_H__
#define __MAPFILE_H__

#include "../idlib/qfiles.h"
#include "../idlib/q_containers.h"
#include "../idlib/dict.h"
#include "../idlib/q_extramath.h"

/*
===============================================================================

	Reads or writes the contents of .map files into a standard internal
	format, which can then be moved into private formats for collision
	detection, map processing, or editor use.

	No validation (duplicate planes, null area brushes, etc) is performed.
	There are no limits to the number of any of the elements in maps.
	The order of entities, brushes, and sides is maintained.

===============================================================================
*/

#define OLD_MAP_VERSION					 1
#define CURRENT_MAP_VERSION				 2
#define DEFAULT_CURVE_SUBDIVISION		 4
#define DEFAULT_CURVE_MAX_ERROR			 4.0f
#define DEFAULT_CURVE_MAX_ERROR_CD		 24.0f
#define DEFAULT_CURVE_MAX_LENGTH		 -1.0f
#define DEFAULT_CURVE_MAX_LENGTH_CD		 -1.0;

typedef enum {
    PRIMTYPE_INVALID = -1,
    PRIMTYPE_BRUSH,
    PRIMTYPE_PATCH
} mapPrimitiveType_t;

typedef struct mapBrushSide_s {
	char					material[MAX_QPATH];
	plane_t				    plane;
	vec3_t					texMat[2];
	vec3_t					origin;
    struct mapBrushSide_s   *next;
} mapBrushSide_t;

typedef struct mapPrimitive_s {

	dict_t   				epairs;

	mapPrimitiveType_t		type;

    struct mapPrimitive_s * next;

    // brushes
	int						numSides;
	mapBrushSide_t *        sides;
    
    // patches
	char					material[MAX_QPATH];
	int						horzSubdivisions;
	int						vertSubdivisions;
	qboolean				explicitSubdivisions;
} mapPrimitive_t;

static ID_INLINE void MapBrushSideClear( mapBrushSide_t *bs ) {
    PlaneZero( bs->plane );
    VectorClear( bs->texMat[0] );
    VectorClear( bs->texMat[1] );
	VectorClear( bs->origin );
	bs->next = NULL;
    bs->material[0] = 0;
}

static ID_INLINE void MapBrushInit( mapPrimitive_t *p ) {
    DictInit( &p->epairs );
    p->type = PRIMTYPE_BRUSH;
    p->numSides = 0;
    p->sides = NULL;
    p->material[0] = 0;
    p->next = NULL;
}

static ID_INLINE void FreeMapBrush( mapPrimitive_t *p ) {
    mapBrushSide_t *side, *next;
    
    DictFree( &p->epairs );

    side = p->sides;
    while ( side ) {
        p->numSides--;
        next = side->next;
        ii.FreeMemory( side );
        side = next;
    }
    p->next = NULL;
}

static ID_INLINE void MapPrimitiveFree( mapPrimitive_t *p ) {
	switch ( p->type ) {
		case PRIMTYPE_BRUSH:
			FreeMapBrush( p );
			break;
		case PRIMTYPE_PATCH:
			// TODO: implement
			break;
		case PRIMTYPE_INVALID:
			break;
	}
}

typedef struct source_s source_t;

mapPrimitive_t *ParseMapBrush( source_t *src, const vec3_t origin, qboolean newFormat, float version );
mapPrimitive_t *ParseMapBrushQ3( source_t *src, const vec3_t origin );
qboolean MapBrushWrite( const mapPrimitive_t *p, fileHandle_t fp, int primitiveNum, const vec3_t origin );
unsigned int GetMapBrushGeometryCRC( const mapPrimitive_t *brush );

static ID_INLINE int GetNumSides( const mapPrimitive_t *p ) {
    assert( p->type == PRIMTYPE_BRUSH );
    return p->numSides;
}

static ID_INLINE int AddSideToMapBrush( mapPrimitive_t *p, mapBrushSide_t *side ) {
    assert( p->type == PRIMTYPE_BRUSH );
    side->next = p->sides;
    p->sides = side;
    p->numSides++;
    return p->numSides - 1;
}

static ID_INLINE mapBrushSide_t *GetMapBrushSide( const mapPrimitive_t *p, int i ) {
    mapBrushSide_t *side;

    assert( p->type == PRIMTYPE_BRUSH );
    side = p->sides;
    while ( i > 0 && side ) {
        side = side->next;
    }

    return side;
}

/*
class idMapPatch : public idMapPrimitive, public idSurface_Patch {
public:
							idMapPatch( void );
							idMapPatch( int maxPatchWidth, int maxPatchHeight );
							~idMapPatch( void ) { }
	static idMapPatch *		Parse( idLexer &src, const idVec3 &origin, bool patchDef3 = true, float version = CURRENT_MAP_VERSION );
	bool					Write( idFile *fp, int primitiveNum, const idVec3 &origin ) const;
	const char *			GetMaterial( void ) const { return material; }
	void					SetMaterial( const char *p ) { material = p; }
	int						GetHorzSubdivisions( void ) const { return horzSubdivisions; }
	int						GetVertSubdivisions( void ) const { return vertSubdivisions; }
	bool					GetExplicitlySubdivided( void ) const { return explicitSubdivisions; }
	void					SetHorzSubdivisions( int n ) { horzSubdivisions = n; }
	void					SetVertSubdivisions( int n ) { vertSubdivisions = n; }
	void					SetExplicitlySubdivided( bool b ) { explicitSubdivisions = b; }
	unsigned int			GetGeometryCRC( void ) const;

protected:
};

ID_INLINE idMapPatch::idMapPatch( void ) {
	type = TYPE_PATCH;
	horzSubdivisions = vertSubdivisions = 0;
	explicitSubdivisions = false;
	width = height = 0;
	maxWidth = maxHeight = 0;
	expanded = false;
}

ID_INLINE idMapPatch::idMapPatch( int maxPatchWidth, int maxPatchHeight ) {
	type = TYPE_PATCH;
	horzSubdivisions = vertSubdivisions = 0;
	explicitSubdivisions = false;
	width = height = 0;
	maxWidth = maxPatchWidth;
	maxHeight = maxPatchHeight;
	verts.SetNum( maxWidth * maxHeight );
	expanded = false;
}
*/

typedef struct mapEntity_s {
    dict_t              epairs;
    int                 numPrimitives;
	mapPrimitive_t *	primitives;
	struct mapEntity_s *next;
} mapEntity_t;

static ID_INLINE void MapEntityInit( mapEntity_t *e ) {
    DictInit( &e->epairs );
    e->primitives = NULL;
    e->numPrimitives = 0;
    e->next = NULL;
}

static ID_INLINE void MapEntityFree( mapEntity_t *e ) {
    mapPrimitive_t *p, *next;
	p = e->primitives;
	while ( p ) {
		next = p->next;

		MapPrimitiveFree( p );
		ii.FreeMemory( p );

		p = next;
	}
	e->primitives = NULL;
	e->numPrimitives = 0;
	DictFree( &e->epairs );
}

mapEntity_t *ParseMapEntity( source_t *src, qboolean worldSpawn, float version );
qboolean MapEntityWrite( const mapEntity_t *e, fileHandle_t fp, int entityNum );
void RemoveMapEntityPrimitiveData( mapEntity_t *e );
unsigned int GetMapEntityGeometryCRC( const mapEntity_t *e );

static ID_INLINE int GetMapEntityNumPrimitives( const mapEntity_t *e ) {
    return e->numPrimitives;
}

static ID_INLINE mapPrimitive_t *GetMapEntityPrimitive( const mapEntity_t *e, int i ) {
    mapPrimitive_t *prim;

    prim = e->primitives;
    while ( i > 0 && prim ) {
        prim = prim->next;
        i--;
    }

    return prim;
}

static ID_INLINE void AddPrimitiveToMapEntity( mapEntity_t *e, mapPrimitive_t *p ) {
    p->next = e->primitives;
    e->primitives = p;
    e->numPrimitives++;
}

void RemoveMapEntityPrimitiveData( mapEntity_t *e );

typedef struct mapFile_s {
	float					version;
	long					fileTime;
	unsigned int			geometryCRC;
	int						numEntities;
	mapEntity_t             *entities;
	char					name[MAX_QPATH];
	qboolean				hasPrimitiveData;
} mapFile_t;

qboolean ParseMapFileFromSource( mapFile_t *m, source_t *src );
qboolean ParseMapFile( mapFile_t *m, const char *filename, qboolean ignoreRegion, qboolean osPath );
qboolean MapFileWrite( const mapFile_t *m, const char *fileName, const char *ext, qboolean fromBasePath );
void SetMapFileGeometryCRC( mapFile_t *m );
int AddEntityToMapFile( mapFile_t *m, mapEntity_t *mapEnt );
mapEntity_t *FindEntityInMapFile( const mapFile_t *m, const char *name );
void RemoveEntityFromMapFile( mapFile_t *m, mapEntity_t *mapEnt );
void RemoveEntitiesFromMapFile( mapFile_t *m, const char *classname );
void RemoveAllEntitiesFromMapFile( mapFile_t *m );
void RemovePrimitiveDataFromMapFile( mapFile_t *m );
qboolean MapFileNeedsReload( mapFile_t *m );

static ID_INLINE void MapFileInit( mapFile_t *m ) {
	m->version = CURRENT_MAP_VERSION;
	m->fileTime = 0;
	m->geometryCRC = 0;
	m->entities = NULL;
	m->numEntities = 0;
	m->name[0] = 0;
	m->hasPrimitiveData = qfalse;
}

static ID_INLINE void MapFileFree( mapFile_t *m ) {
	RemoveAllEntitiesFromMapFile( m );
}

#endif /* !__MAPFILE_H__ */
