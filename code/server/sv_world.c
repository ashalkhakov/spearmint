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
// world.c -- world query functions

#include "server.h"
#include "../qcommon/debugvis.h"

/*
================
SV_ClipHandleForEntity

Returns a headnode that can be used for testing or clipping to a
given entity.  If the entity is a bsp model, the headnode will
be returned, otherwise a custom box tree or capsule will be constructed.
================
*/
clipHandle_t SV_ClipHandleForEntity( const sharedEntity_t *ent ) {
	char			modelName[MAX_QPATH];
	traceModel_t	trm;
	vec3_t			bounds[2];

	if ( ent->s.collisionType == CT_SUBMODEL ) {
		// explicit hulls in the BSP model
		snprintf( modelName, sizeof(modelName), "*%d", ent->s.modelindex );
		return cme.LoadModel( modelName, qfalse );
	}

	// create a temp tree or capsule from bounding box sizes
    trm.type = TRM_INVALID;
	VectorCopy( ent->s.mins, bounds[0] );
	VectorCopy( ent->s.maxs, bounds[1] );
	switch ( ent->s.collisionType ) {
		case CT_AABB:
			TraceModelSetupBox( &trm, bounds );
			break;
		case CT_CAPSULE:
			TraceModelSetupDodecahedron( &trm, bounds );
			break;
	}

	// FIXME: need a dummy material with the given contents
	return cme.SetupTrmModel( &trm, 0 );
}



/*
===============================================================================

ENTITY CHECKING

To avoid linearly searching through lists of entities during environment testing,
the world is carved up with an evenly spaced, axially aligned bsp tree.  Entities
are kept in chains either at the final leafs, or at the first node that splits
them, which prevents having to deal with multiple fragments of a single entity.

===============================================================================
*/

typedef struct worldSector_s {
	int		axis;		// -1 = leaf node
	float	dist;
	struct worldSector_s	*children[2];
	svEntity_t	*entities;
} worldSector_t;

#define	AREA_DEPTH	4
#define	AREA_NODES	64

worldSector_t	sv_worldSectors[AREA_NODES];
int			sv_numworldSectors;


/*
===============
SV_SectorList_f
===============
*/
void SV_SectorList_f( void ) {
#if USE_CLIP
    ClipPrintSectorList( &sv.clip );
#else
	int				i, c;
	worldSector_t	*sec;
	svEntity_t		*ent;

	for ( i = 0 ; i < AREA_NODES ; i++ ) {
		sec = &sv_worldSectors[i];

		c = 0;
		for ( ent = sec->entities ; ent ; ent = ent->nextEntityInWorldSector ) {
			c++;
		}
		Com_Printf( "sector %i: %i entities\n", i, c );
	}
#endif
}

#if !USE_CLIP
/*
===============
SV_CreateworldSector

Builds a uniformly subdivided tree for the given world size
===============
*/
static worldSector_t *SV_CreateworldSector( int depth, vec3_t mins, vec3_t maxs ) {
	worldSector_t	*anode;
	vec3_t		size;
	vec3_t		mins1, maxs1, mins2, maxs2;

	anode = &sv_worldSectors[sv_numworldSectors];
	sv_numworldSectors++;

	if (depth == AREA_DEPTH) {
		anode->axis = -1;
		anode->children[0] = anode->children[1] = NULL;
		return anode;
	}
	
	VectorSubtract (maxs, mins, size);
	if (size[0] > size[1]) {
		anode->axis = 0;
	} else {
		anode->axis = 1;
	}

	anode->dist = 0.5 * (maxs[anode->axis] + mins[anode->axis]);
	VectorCopy (mins, mins1);	
	VectorCopy (mins, mins2);	
	VectorCopy (maxs, maxs1);	
	VectorCopy (maxs, maxs2);	
	
	maxs1[anode->axis] = mins2[anode->axis] = anode->dist;
	
	anode->children[0] = SV_CreateworldSector (depth+1, mins2, maxs2);
	anode->children[1] = SV_CreateworldSector (depth+1, mins1, maxs1);

	return anode;
}
#endif

/*
===============
SV_ClearWorld

===============
*/
void SV_ClearWorld( void ) {
#if !USE_CLIP
	clipHandle_t	h;
	vec3_t			bounds[2];

	Com_Memset( sv_worldSectors, 0, sizeof(sv_worldSectors) );
	sv_numworldSectors = 0;

	// get world map bounds
	h = cme.GetModelBounds( 0, bounds );
	SV_CreateworldSector( 0, bounds[0], bounds[1] );
#endif
}


/*
===============
SV_UnlinkEntity

===============
*/
void SV_UnlinkEntity( sharedEntity_t *gEnt ) {
	svEntity_t		*ent;
	svEntity_t		*scan;
#if !USE_CLIP
	worldSector_t	*ws;
#endif
	sharedPlayerState_t	*ps;

	ent = SV_SvEntityForGentity( gEnt );

	gEnt->r.linked = qfalse;

	if (gEnt->s.number < MAX_CLIENTS) {
		ps = SV_GamePlayerNum(gEnt->s.number);
		ps->linked = qfalse;
	}

#if USE_CLIP
    ClipModelUnlink( &ent->model );
    ClipModelSetOwner( &ent->model, NULL );
#else
	ws = ent->worldSector;
	if ( !ws ) {
		return;		// not linked in anywhere
	}
	ent->worldSector = NULL;

	if ( ws->entities == ent ) {
		ws->entities = ent->nextEntityInWorldSector;
		return;
	}

	for ( scan = ws->entities ; scan ; scan = scan->nextEntityInWorldSector ) {
		if ( scan->nextEntityInWorldSector == ent ) {
			scan->nextEntityInWorldSector = ent->nextEntityInWorldSector;
			return;
		}
	}

	Com_Printf( "WARNING: SV_UnlinkEntity: not found in worldSector\n" );
#endif
}


/*
===============
SV_LinkEntity

===============
*/
#define MAX_TOTAL_ENT_LEAFS		128
void SV_LinkEntity( sharedEntity_t *gEnt ) {
	worldSector_t	*node;
	int			leafs[MAX_TOTAL_ENT_LEAFS];
	int			cluster;
	int			num_leafs;
	int			i;
	int			area;
	int			lastLeaf;
	float		*origin, *angles;
    vec3_t      axis[3];
	svEntity_t	*ent;
	sharedPlayerState_t	*ps;

	ent = SV_SvEntityForGentity( gEnt );

#if !USE_CLIP
	if ( ent->worldSector ) {
		SV_UnlinkEntity( gEnt );	// unlink from old position
	}
#endif

	// get the position
	origin = gEnt->r.currentOrigin;
	angles = gEnt->r.currentAngles;

	// set the abs box
	if ( gEnt->s.collisionType == CT_SUBMODEL && (angles[0] || angles[1] || angles[2]) ) {
		// expand for rotation
		float		max;

		max = RadiusFromBounds( gEnt->s.mins, gEnt->s.maxs );
		for (i=0 ; i<3 ; i++) {
			gEnt->r.absmin[i] = origin[i] - max;
			gEnt->r.absmax[i] = origin[i] + max;
		}
	} else {
		// normal
		VectorAdd (origin, gEnt->s.mins, gEnt->r.absmin);
		VectorAdd (origin, gEnt->s.maxs, gEnt->r.absmax);
	}

#if USE_CLIP
    // TODO: what should the lifetime of models be?
    ClipModelInit( &ent->model );

    // TODO: what should be put here?
    int newId = gEnt->s.number;
	if ( gEnt->s.collisionType == CT_SUBMODEL ) {
        char name[MAX_QPATH];

        AnglesToAxis( angles, axis );
        sprintf( name, "*%d", gEnt->s.modelindex );
        ClipModelInitFromModel( &ent->model, name );
    } else {
        traceModel_t trm;
        vec3_t bounds[2];

        // boxes don't rotate
        AxisClear( axis );

        VectorCopy( gEnt->s.mins, bounds[0] );
        VectorCopy( gEnt->s.maxs, bounds[1] );
        trm.type = TRM_INVALID;
        TraceModelSetupBox( &trm, bounds );
        ClipModelInitFromTraceModel( &ent->model, &trm );
    }
    ClipModelLink2( &ent->model, &sv.clip, gEnt, newId, origin, axis, -1 );
    if ( gEnt->r.ownerNum != -1 ) {
        sharedEntity_t *owner = SV_GentityNum( gEnt->r.ownerNum );
        ClipModelSetOwner( &ent->model, owner );
    } else {
        ClipModelSetOwner( &ent->model, NULL );
    }

#endif

	// because movement is clipped an epsilon away from an actual edge,
	// we must fully check even when bounding boxes don't quite touch
	gEnt->r.absmin[0] -= 1;
	gEnt->r.absmin[1] -= 1;
	gEnt->r.absmin[2] -= 1;
	gEnt->r.absmax[0] += 1;
	gEnt->r.absmax[1] += 1;
	gEnt->r.absmax[2] += 1;

	// link to PVS leafs
	ent->numClusters = 0;
	ent->lastCluster = 0;
	ent->areanum = -1;
	ent->areanum2 = -1;

	//get all leafs, including solids
	num_leafs = cme.BoxLeafnums( gEnt->r.absmin, gEnt->r.absmax,
		leafs, MAX_TOTAL_ENT_LEAFS, &lastLeaf );

	// if none of the leafs were inside the map, the
	// entity is outside the world and can be considered unlinked
	if ( !num_leafs ) {
		return;
	}

	// set areas, even from clusters that don't fit in the entity array
	for (i=0 ; i<num_leafs ; i++) {
		area = cme.LeafArea (leafs[i]);
		if (area != -1) {
			// doors may legally straggle two areas,
			// but nothing should evern need more than that
			if (ent->areanum != -1 && ent->areanum != area) {
				if (ent->areanum2 != -1 && ent->areanum2 != area && sv.state == SS_LOADING) {
					Com_DPrintf ("Object %i touching 3 areas at %f %f %f\n",
					gEnt->s.number,
					gEnt->r.absmin[0], gEnt->r.absmin[1], gEnt->r.absmin[2]);
				}
				ent->areanum2 = area;
			} else {
				ent->areanum = area;
			}
		}
	}

	// store as many explicit clusters as we can
	ent->numClusters = 0;
	for (i=0 ; i < num_leafs ; i++) {
		cluster = cme.LeafCluster( leafs[i] );
		if ( cluster != -1 ) {
			ent->clusternums[ent->numClusters++] = cluster;
			if ( ent->numClusters == MAX_ENT_CLUSTERS ) {
				break;
			}
		}
	}

	// store off a last cluster if we need to
	if ( i != num_leafs ) {
		ent->lastCluster = cme.LeafCluster( lastLeaf );
	}

#if !USE_CLIP
	gEnt->r.linkcount++;

	// find the first world sector node that the ent's box crosses
	node = sv_worldSectors;
	while (1)
	{
		if (node->axis == -1)
			break;
		if ( gEnt->r.absmin[node->axis] > node->dist)
			node = node->children[0];
		else if ( gEnt->r.absmax[node->axis] < node->dist)
			node = node->children[1];
		else
			break;		// crosses the node
	}
	
	// link it in
	ent->worldSector = node;
	ent->nextEntityInWorldSector = node->entities;
	node->entities = ent;
#endif
	gEnt->r.linked = qtrue;

	if (gEnt->s.number < MAX_CLIENTS) {
		ps = SV_GamePlayerNum(gEnt->s.number);
		ps->linked = qtrue;
	}
}

/*
============================================================================

AREA QUERY

Fills in a list of all entities who's absmin / absmax intersects the given
bounds.  This does NOT mean that they actually touch in the case of bmodels.
============================================================================
*/

typedef struct {
	const float	*mins;
	const float	*maxs;
	int			*list;
	int			count, maxcount;
} areaParms_t;

#if !USE_CLIP

/*
====================
SV_AreaEntities_r

====================
*/
static void SV_AreaEntities_r( worldSector_t *node, areaParms_t *ap ) {
	svEntity_t	*check, *next;
	sharedEntity_t *gcheck;

	for ( check = node->entities  ; check ; check = next ) {
		next = check->nextEntityInWorldSector;

		gcheck = SV_GEntityForSvEntity( check );

		if ( !gcheck->r.linked ) {
			continue;
		}

		if ( gcheck->r.absmin[0] > ap->maxs[0]
		|| gcheck->r.absmin[1] > ap->maxs[1]
		|| gcheck->r.absmin[2] > ap->maxs[2]
		|| gcheck->r.absmax[0] < ap->mins[0]
		|| gcheck->r.absmax[1] < ap->mins[1]
		|| gcheck->r.absmax[2] < ap->mins[2]) {
			continue;
		}

		if ( ap->count == ap->maxcount ) {
			Com_Printf ("SV_AreaEntities: MAXCOUNT\n");
			return;
		}

		ap->list[ap->count] = check - sv.svEntities;
		ap->count++;
	}
	
	if (node->axis == -1) {
		return;		// terminal node
	}

	// recurse down both sides
	if ( ap->maxs[node->axis] > node->dist ) {
		SV_AreaEntities_r ( node->children[0], ap );
	}
	if ( ap->mins[node->axis] < node->dist ) {
		SV_AreaEntities_r ( node->children[1], ap );
	}
}

/*
================
SV_AreaEntities
================
*/
int SV_AreaEntities( const vec3_t mins, const vec3_t maxs, int *entityList, int maxcount ) {
	areaParms_t		ap;
    
	ap.mins = mins;
	ap.maxs = maxs;
	ap.list = entityList;
	ap.count = 0;
	ap.maxcount = maxcount;

	SV_AreaEntities_r( sv_worldSectors, &ap );

	return ap.count;
}

#else

/*
================
SV_AreaEntities
================
*/
int SV_AreaEntities( const vec3_t mins, const vec3_t maxs, int *entityList, int maxcount ) {
	vec3_t bounds[2];
    sharedEntity_t **list;
    int i, num;

    VectorCopy( mins, bounds[0] );
    VectorCopy( maxs, bounds[1] );

    // TODO: get rid of memory allocation
    list = ii.GetMemory( maxcount * sizeof( *list ) );

    num = ClipEntitiesTouchingBounds( &sv.clip, bounds, -1, list, maxcount );
    for ( i = 0; i < num; i++ ) {
        entityList[i] = list[i]->s.number;
    }

    ii.FreeMemory( list );

	return num;
}


#endif

//===========================================================================

#if !USE_CLIP

typedef struct {
	vec3_t		boxmins, boxmaxs;// enclose the test object along entire move
	const float	*mins;
	const float *maxs;	// size of the moving object
	const float	*start;
	vec3_t		end;
	trace_t		trace;
	int			passEntityNum;
	int			contentmask;
	traceType_t	traceType;
} moveclip_t;


/*
====================
SV_ClipToEntity

====================
*/
void SV_ClipToEntity( trace_t *trace, const vec3_t start, const vec3_t mins, const vec3_t maxs, const vec3_t end, int entityNum, int contentmask, traceType_t type ) {
	sharedEntity_t	*touch;
	clipHandle_t	clipHandle;
	float			*origin, *angles;

	touch = SV_GentityNum( entityNum );

	Com_Memset(trace, 0, sizeof(trace_t));

	// if it doesn't have any brushes of a type we
	// are looking for, ignore it
	if ( ! ( contentmask & touch->s.contents ) ) {
		trace->fraction = 1.0;
		return;
	}

	// might intersect, so do an exact clip
	clipHandle = SV_ClipHandleForEntity (touch);

	origin = touch->r.currentOrigin;
	angles = touch->r.currentAngles;

	if ( touch->s.collisionType != CT_SUBMODEL ) {
		angles = vec3_origin;	// boxes don't rotate
	}

	CM_TransformedBoxTrace ( trace, start, end,
		mins, maxs, clipHandle,  contentmask,
		origin, angles, type);

	if ( trace->fraction < 1 ) {
		trace->entityNum = touch->s.number;
	}
}


/*
====================
SV_ClipMoveToEntities

====================
*/
static void SV_ClipMoveToEntities( moveclip_t *clip ) {
	int			i, num;
	int			touchlist[MAX_GENTITIES];
	sharedEntity_t *touch;
	int			passOwnerNum;
	trace_t		trace;
	clipHandle_t	clipHandle;
	float		*origin, *angles;

	num = SV_AreaEntities( clip->boxmins, clip->boxmaxs, touchlist, MAX_GENTITIES);

	if ( clip->passEntityNum != ENTITYNUM_NONE ) {
		passOwnerNum = ( SV_GentityNum( clip->passEntityNum ) )->r.ownerNum;
		if ( passOwnerNum == ENTITYNUM_NONE ) {
			passOwnerNum = -1;
		}
	} else {
		passOwnerNum = -1;
	}

	for ( i=0 ; i<num ; i++ ) {
		if ( clip->trace.allsolid ) {
			return;
		}
		touch = SV_GentityNum( touchlist[i] );

		// see if we should ignore this entity
		if ( clip->passEntityNum != ENTITYNUM_NONE ) {
			if ( touchlist[i] == clip->passEntityNum ) {
				continue;	// don't clip against the pass entity
			}
			if ( touch->r.ownerNum == clip->passEntityNum ) {
				continue;	// don't clip against own missiles
			}
			if ( touch->r.ownerNum == passOwnerNum ) {
				continue;	// don't clip against other missiles from our owner
			}
		}

		// if it doesn't have any brushes of a type we
		// are looking for, ignore it
		if ( ! ( clip->contentmask & touch->s.contents ) ) {
			continue;
		}

		// might intersect, so do an exact clip
		clipHandle = SV_ClipHandleForEntity (touch);

		// non-worldspawn entities must not use world as clip model!
		if ( clipHandle == 0 ) {
			continue;
		}

		origin = touch->r.currentOrigin;
		angles = touch->r.currentAngles;


		if ( touch->s.collisionType != CT_SUBMODEL ) {
			angles = vec3_origin;	// boxes don't rotate
		}

		CM_TransformedBoxTrace ( &trace, (float *)clip->start, (float *)clip->end,
			(float *)clip->mins, (float *)clip->maxs, clipHandle,  clip->contentmask,
			origin, angles, clip->traceType);

		if ( trace.allsolid ) {
			clip->trace.allsolid = qtrue;
			clip->trace.entityNum = touch->s.number;
		} else if ( trace.startsolid ) {
			clip->trace.startsolid = qtrue;
			clip->trace.entityNum = touch->s.number;
		}

		if ( trace.fraction < clip->trace.fraction ) {
			qboolean	oldStart;

			// make sure we keep a startsolid from a previous trace
			oldStart = clip->trace.startsolid;

			trace.entityNum = touch->s.number;
			clip->trace = trace;
			clip->trace.startsolid |= oldStart;
		}
	}
}


/*
==================
SV_Trace

Moves the given mins/maxs volume through the world from start to end.
passEntityNum and entities owned by passEntityNum are explicitly not checked.
==================
*/
void SV_Trace( trace_t *results, const vec3_t start, const vec3_t mins, const vec3_t maxs, const vec3_t end, int passEntityNum, int contentmask, traceType_t type ) {
	moveclip_t	clip;
	int			i;

	if ( !mins ) {
		mins = vec3_origin;
	}
	if ( !maxs ) {
		maxs = vec3_origin;
	}

	Com_Memset ( &clip, 0, sizeof ( moveclip_t ) );

	// clip to world
	CM_BoxTrace( &clip.trace, start, end, mins, maxs, 0, contentmask, type );
	clip.trace.entityNum = clip.trace.fraction != 1.0 ? ENTITYNUM_WORLD : ENTITYNUM_NONE;
	if ( clip.trace.fraction == 0 ) {
		*results = clip.trace;
		return;		// blocked immediately by the world
	}

	clip.contentmask = contentmask;
	clip.start = start;
//	VectorCopy( clip.trace.endpos, clip.end );
	VectorCopy( end, clip.end );
	clip.mins = mins;
	clip.maxs = maxs;
	clip.passEntityNum = passEntityNum;
	clip.traceType = type;

	// create the bounding box of the entire move
	// we can limit it to the part of the move not
	// already clipped off by the world, which can be
	// a significant savings for line of sight and shot traces
	for ( i=0 ; i<3 ; i++ ) {
		if ( end[i] > start[i] ) {
			clip.boxmins[i] = clip.start[i] + clip.mins[i] - 1;
			clip.boxmaxs[i] = clip.end[i] + clip.maxs[i] + 1;
		} else {
			clip.boxmins[i] = clip.end[i] + clip.mins[i] - 1;
			clip.boxmaxs[i] = clip.start[i] + clip.maxs[i] + 1;
		}
	}

	// clip to other solid entities
	SV_ClipMoveToEntities ( &clip );

	*results = clip.trace;
}


/*
=============
SV_ClipToEntities

SV_Trace that doesn't clip to world
=============
*/
void SV_ClipToEntities( trace_t *results, const vec3_t start, const vec3_t mins, const vec3_t maxs, const vec3_t end, int passEntityNum, int contentmask, traceType_t type ) {
	moveclip_t	clip;
	int			i;

	if ( !mins ) {
		mins = vec3_origin;
	}
	if ( !maxs ) {
		maxs = vec3_origin;
	}

	Com_Memset ( &clip, 0, sizeof ( moveclip_t ) );

	// Skip clipping to world
	clip.trace.fraction = 1; // assume it goes the entire distance until shown otherwise
	VectorCopy (end, clip.trace.endpos);
	clip.trace.entityNum = ENTITYNUM_NONE;

	clip.contentmask = contentmask;
	clip.start = start;
//	VectorCopy( clip.trace.endpos, clip.end );
	VectorCopy( end, clip.end );
	clip.mins = mins;
	clip.maxs = maxs;
	clip.passEntityNum = passEntityNum;
	clip.traceType = type;

	// create the bounding box of the entire move
	// we can limit it to the part of the move not
	// already clipped off by the world, which can be
	// a significant savings for line of sight and shot traces
	for ( i=0 ; i<3 ; i++ ) {
		if ( end[i] > start[i] ) {
			clip.boxmins[i] = clip.start[i] + clip.mins[i] - 1;
			clip.boxmaxs[i] = clip.end[i] + clip.maxs[i] + 1;
		} else {
			clip.boxmins[i] = clip.end[i] + clip.mins[i] - 1;
			clip.boxmaxs[i] = clip.start[i] + clip.maxs[i] + 1;
		}
	}

	// clip to other solid entities
	SV_ClipMoveToEntities ( &clip );

	*results = clip.trace;
}

/*
=============
SV_PointContents
=============
*/
int SV_PointContents( const vec3_t p, int passEntityNum ) {
	int			touch[MAX_GENTITIES];
	sharedEntity_t *hit;
	int			i, num;
	int			contents, c2;
	clipHandle_t	clipHandle;
	float		*angles;
	vec3_t		axis[3];

	// get base contents from world
	contents = cme.Contents( p, NULL, axisDefault, 0, 0, vec3_origin, axisDefault );

	// or in contents from all the other entities
	num = SV_AreaEntities( p, p, touch, MAX_GENTITIES );

	for ( i=0 ; i<num ; i++ ) {
		if ( touch[i] == passEntityNum ) {
			continue;
		}
		hit = SV_GentityNum( touch[i] );
		// might intersect, so do an exact clip
		clipHandle = SV_ClipHandleForEntity( hit );

		// non-worldspawn entities must not use world as clip model!
		if ( clipHandle == 0 ) {
			continue;
		}

		angles = hit->r.currentAngles;
		if ( hit->s.collisionType != CT_SUBMODEL ) {
			angles = vec3_origin;	// boxes don't rotate
		}

		AnglesToAxis( angles, axis );

		c2 = cme.Contents (p, NULL, axisDefault, 0, clipHandle, hit->r.currentOrigin, axis);

		contents |= c2;
	}

	return contents;
}

#else

static void CmTraceToTrace( cm_trace_t *t, const vec3_t start, const vec3_t end, const vec3_t mins, const vec3_t maxs, trace_t *trace, traceModel_t *trm, cmHandle_t model, vec3_t origin, vec3_t axis[3], int brushmask ) {
    int i;
    vec3_t dir;

    memset( trace, 0, sizeof( *trace ) );

	trace->fraction = t->fraction;
    VectorSubtract( end, start, dir );
    VectorMA( start, trace->fraction, dir, trace->endpos );
	//VectorCopy( t->endpos, trace->endpos ); // TODO: worldspace or not?

    if ( t->c.type != CONTACT_NONE ) {
        VectorCopy( t->c.point, trace->endpos );
        VectorCopy( t->c.normal, trace->plane.normal );
        trace->plane.dist = -t->c.dist;
		trace->plane.type = PlaneTypeForNormal( trace->plane.normal );
		SetPlaneSignbits( &trace->plane );
        trace->surfaceNum = t->c.modelFeature + 1;
        trace->surfaceFlags = cme.CM_GetMaterialSurfaceFlags( t->c.material );
        trace->contents = t->c.contents;
        trace->entityNum = t->c.entityNum;
    }
	trace->lateralFraction = 0.0f;
/*
    if ( cme.Contents( start, trm, axisDefault, -1, model, origin, axis ) & brushmask ) {
        trace->startsolid = qtrue;
    }
    if ( trace->startsolid ) {
        int c = cme.Contents( trace->endpos, trm, axisDefault, -1, model, origin, axis ) & brushmask;
        if ( c ) {
            trace->allsolid = qtrue;
            //trace->fraction = 0;
            //trace->contents = c;
        }
    }*/
}

/*
==================
SV_Trace

Moves the given mins/maxs volume through the world from start to end.
passEntityNum and entities owned by passEntityNum are explicitly not checked.
==================
*/
void SV_Trace( trace_t *results, const vec3_t start, const vec3_t mins, const vec3_t maxs, const vec3_t end, int passEntityNum, int contentmask, traceType_t type ) {
    cm_trace_t      trace;
    sharedEntity_t  *passEntity;
    vec3_t          bounds[2];

    passEntity = SV_GentityNum( passEntityNum );

	if ( !mins ) {
		mins = vec3_origin;
	}
	if ( !maxs ) {
		maxs = vec3_origin;
	}

    VectorCopy( mins, bounds[0] );
    VectorCopy( maxs, bounds[1] );

    ClipTraceBounds( &sv.clip, &trace, start, end, bounds, contentmask, passEntity );
    CmTraceToTrace( &trace, start, end, mins, maxs, results, NULL, 0, vec3_origin, axisDefault, contentmask );
}


/*
=============
SV_ClipToEntities

SV_Trace that doesn't clip to world
=============
*/
void SV_ClipToEntities( trace_t *results, const vec3_t start, const vec3_t mins, const vec3_t maxs, const vec3_t end, int passEntityNum, int contentmask, traceType_t type ) {
    cm_trace_t      trace;
    clipModel_t     *model;
    sharedEntity_t  *passEntity;
    vec3_t          bounds[2];
    traceModel_t    trm;

    passEntity = SV_GentityNum( passEntityNum );

	if ( !mins ) {
		mins = vec3_origin;
	}
	if ( !maxs ) {
		maxs = vec3_origin;
	}

    VectorCopy( mins, bounds[0] );
    VectorCopy( maxs, bounds[1] );
    trm.type = TRM_INVALID;
    TraceModelSetupBox( &trm, bounds );
    model = ClipDefaultClipModel( &sv.clip );
    ClipModelLoadTraceModel( model, &trm );

    ClipTranslationEntities( &sv.clip, &trace, start, end, model, axisDefault, contentmask, passEntity );
    CmTraceToTrace( &trace, start, end, mins, maxs, results, NULL, 0, vec3_origin, axisDefault, contentmask );
}

/*
=============
SV_PointContents
=============
*/
int SV_PointContents( const vec3_t p, int passEntityNum ) {
    clipModel_t *model;
    sharedEntity_t *passEntity;
    traceModel_t trm;

    passEntity = SV_GentityNum( passEntityNum );

    model = ClipDefaultClipModel( &sv.clip );

    trm.type = TRM_INVALID;
    TraceModelSetupBox2( &trm, 0 );
    ClipModelLoadTraceModel( model, &trm );

    return ClipContents( &sv.clip, p, model, axisDefault, -1, passEntity );
}

#endif
