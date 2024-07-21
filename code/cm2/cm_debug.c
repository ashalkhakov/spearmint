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
#include "../idlib/surfaceflags.h"
#include "../idlib/l_script.h"
#include "../idlib/l_precomp.h"

/*
===============================================================================

Visualisation code

===============================================================================
*/

const char *cm_contentsNameByIndex[] = {
	"none",							// 0
	"solid",						// 1
	//"opaque",						// 2
	"water",						// 3
	"playerclip",					// 4
	"monsterclip",					// 5
	//"moveableclip",					// 6
	//"ikclip",						// 7
	//"blood",						// 8
	"body",							// 9
	"corpse",						// 10
	"trigger",						// 11
	//"aas_solid",					// 12
	//"aas_obstacle",					// 13
	//"flashlight_trigger",			// 14
	NULL
};

int cm_contentsFlagByIndex[] = {
	-1,								// 0
	CONTENTS_SOLID,					// 1
	//CONTENTS_OPAQUE,				// 2
	CONTENTS_WATER,					// 3
	CONTENTS_PLAYERCLIP,			// 4
	CONTENTS_MONSTERCLIP,			// 5
	//CONTENTS_MOVEABLECLIP,			// 6
	//CONTENTS_IKCLIP,				// 7
	//CONTENTS_BLOOD,					// 8
	CONTENTS_BODY,					// 9
	CONTENTS_CORPSE,				// 10
	CONTENTS_TRIGGER,				// 11
	//CONTENTS_AAS_SOLID,				// 12
	//CONTENTS_AAS_OBSTACLE,			// 13
	//CONTENTS_FLASHLIGHT_TRIGGER,	// 14
	0
};

static vec4_t cm_color;

/*
================
CM_ContentsFromString
================
*/
int CM_ContentsFromString( const char *string ) {
	int i, contents = 0;
	source_t *src;
	token_t token;

	src = LoadSourceMemory( string, strlen( string ), "ContentsFromString", NULL );

	while( PC_ReadToken( src, &token ) ) {
		if ( !strcmp( token.string, "," ) ) {
			continue;
		}
		for ( i = 1; cm_contentsNameByIndex[i] != NULL; i++ ) {
			if ( Q_stricmp( token.string, cm_contentsNameByIndex[i] ) == 0 ) {
				contents |= cm_contentsFlagByIndex[i];
				break;
			}
		}
	}

	FreeSource( src );

	return contents;
}

/*
================
CM_StringFromContents
================
*/
const char *CM_StringFromContents( const int contents ) {
	int i, length = 0;
	static char contentsString[MAX_STRING_CHARS];

	contentsString[0] = '\0';

	for ( i = 1; cm_contentsFlagByIndex[i] != 0; i++ ) {
		if ( contents & cm_contentsFlagByIndex[i] ) {
			if ( length != 0 ) {
				length += snprintf( contentsString + length, sizeof( contentsString ) - length, "," );
			}
			length += snprintf( contentsString + length, sizeof( contentsString ) - length, cm_contentsNameByIndex[i] );
		}
	}

	return contentsString;
}

/*
================
CM_DrawEdge
================
*/
void CM_DrawEdge( cm_model_t *model, int edgeNum, const vec3_t origin, const vec3_t axis[3] ) {
	int side;
	cm_edge_t *edge;
	vec3_t start, end, mid, tmp;
	qboolean isRotated;

	isRotated = AxisIsRotated( axis );

	edge = model->edges + abs(edgeNum);
	side = edgeNum < 0;

	VectorCopy( model->vertices[edge->vertexNum[side]].p, start );
	VectorCopy( model->vertices[edge->vertexNum[!side]].p, end );
	if ( isRotated ) {
		VectorRotateSelf( start, axis );
		VectorRotateSelf( end, axis );
	}
	VectorAdd( start, origin, start );
	VectorAdd( end, origin, end );

	if ( edge->internal ) {
		if ( cm_drawInternal->integer ) {
			cmi.R_DebugArrow( colorGreen, start, end, 1, 0 );
		}
	} else {
		if ( edge->numUsers > 2 ) {
			cmi.R_DebugArrow( colorBlue, start, end, 1, 0 );
		} else {
			cmi.R_DebugArrow( cm_color, start, end, 1, 0 );
		}
	}

	if ( cm_drawNormals->integer ) {
		VectorAdd( start, end, mid );
		VectorScale( mid, 0.5f, mid );
		if ( isRotated ) {
			MatrixRotateVector( edge->normal, axis, tmp );
			VectorMA( mid, 5, tmp, end );
		} else {
			VectorMA( mid, 5, edge->normal, end );
		}
		cmi.R_DebugArrow( colorCyan, mid, end, 1, 0 );
	}
}

/*
================
CM_DrawPolygon
================
*/
void CM_DrawPolygon( cm_model_t *model, cm_polygon_t *p, const vec3_t origin, const vec3_t axis[3], const vec3_t viewOrigin ) {
	int i, edgeNum;
	cm_edge_t *edge;
	vec3_t center, end, dir, tmp;
	vec5_t tmp5;

	if ( cm_backFaceCull->integer ) {
		edgeNum = p->edges[0];
		edge = model->edges + abs(edgeNum);
		VectorSubtract( model->vertices[edge->vertexNum[0]].p, viewOrigin, dir );
		if ( DotProduct( dir, p->plane ) > 0.0f ) {
			return;
		}
	}

	if ( cm_drawNormals->integer ) {
		VectorClear( center );
		for ( i = 0; i < p->numEdges; i++ ) {
			edgeNum = p->edges[i];
			edge = model->edges + abs(edgeNum);
			VectorAdd( center, model->vertices[edge->vertexNum[edgeNum < 0]].p, center );
		}
		VectorScale( center, 1.0f / p->numEdges, center );
		if ( AxisIsRotated( axis ) ) {
			MatrixRotateVector( center, axis, tmp );
			VectorAdd( tmp, origin, center );
			MatrixRotateVector( p->plane, axis, tmp );
			VectorMA( center, 5, tmp, end );
		} else {
			VectorAdd( center, origin, center );
			VectorMA( center, 5, p->plane, end );
		}
		cmi.R_DebugArrow( colorMagenta, center, end, 1, 0 );
	}

	if ( cm_drawFilled->integer ) {
		fixedWinding_t winding;
		ClearFixedWinding( &winding );
		for ( i = p->numEdges - 1; i >= 0; i-- ) {
			edgeNum = p->edges[i];
			edge = model->edges + abs(edgeNum);
			MatrixRotateVector( model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]].p, axis, tmp );
			VectorAdd( origin, tmp, tmp );
			VectorCopy( tmp, tmp5 );
			tmp5[3] = tmp5[4] = 0.0f;
			AddPointToFixedWinding( &winding, tmp5 );
		}
        cmi.R_DebugPolygon( cm_color, &winding, 0, qfalse );
	} else {
		for ( i = 0; i < p->numEdges; i++ ) {
			edgeNum = p->edges[i];
			edge = model->edges + abs(edgeNum);
			if ( edge->checkcount == cmLocal.checkCount ) {
				continue;
			}
			edge->checkcount = cmLocal.checkCount;
			CM_DrawEdge( model, edgeNum, origin, axis );
		}
	}
}

/*
================
CM_DrawNodePolygons
================
*/
void CM_DrawNodePolygons( cm_model_t *model, cm_node_t *node,
										   const vec3_t origin, const vec3_t axis[3],
										   const vec3_t viewOrigin, const float radius ) {
	int i;
	cm_polygon_t *p;
	cm_polygonRef_t *pref;

	while (1) {
		for ( pref = node->polygons; pref; pref = pref->next ) {
			p = pref->p;
			if ( radius ) {
				// polygon bounds should overlap with trace bounds
				for ( i = 0; i < 3; i++ ) {
					if ( p->bounds[0][i] > viewOrigin[i] + radius ) {
						break;
					}
					if ( p->bounds[1][i] < viewOrigin[i] - radius ) {
						break;
					}
				}
				if ( i < 3 ) {
					continue;
				}
			}
			if ( p->checkcount == cmLocal.checkCount ) {
				continue;
			}
			if ( !( p->contents & cm_contentsFlagByIndex[cm_drawMask->integer] ) ) {
				continue;
			}

			CM_DrawPolygon( model, p, origin, axis, viewOrigin );
			p->checkcount = cmLocal.checkCount;
		}
		if ( node->planeType == -1 ) {
			break;
		}
		if ( radius && viewOrigin[node->planeType] > node->planeDist + radius  ) {
			node = node->children[0];
		} else if ( radius && viewOrigin[node->planeType] < node->planeDist - radius  ) {
			node = node->children[1];
		} else {
			CM_DrawNodePolygons( model, node->children[1], origin, axis, viewOrigin, radius );
			node = node->children[0];
		}
	}
}

/*
================
CM_DrawModel
================
*/
void CM_DrawModel( cmHandle_t handle, const vec3_t modelOrigin, const vec3_t modelAxis[3],
					const vec3_t viewOrigin, const float radius ) {

	cm_model_t *model;
	vec3_t viewPos, tmp;
	vec3_t invModelAxis[3];

	if ( handle < 0 && handle >= cmLocal.numModels ) {
		return;
	}

	if ( cm_drawColor->modified ) {
		sscanf( cm_drawColor->string, "%f %f %f %f", &cm_color[0], &cm_color[1], &cm_color[2], &cm_color[3] );
		cm_drawColor->modified = qfalse;
	}

	model = cmLocal.models[ handle ];
	TransposeAxis( modelAxis, invModelAxis );
	VectorSubtract( viewOrigin, modelOrigin, tmp );
	MatrixRotateVector( tmp, invModelAxis, viewPos );

	cmLocal.checkCount++;
	CM_DrawNodePolygons( model, model->node, modelOrigin, modelAxis, viewPos, radius );
}

/*
===============================================================================

Speed test code

===============================================================================
*/

extern cvar_t *cm_testCollision;
extern cvar_t *cm_testRotation;
extern cvar_t *cm_testModel;
extern cvar_t *cm_testTimes;
extern cvar_t *cm_testRandomMany;
extern cvar_t *cm_testOrigin;
extern cvar_t *cm_testReset;
extern cvar_t *cm_testBox;
extern cvar_t *cm_testBoxRotation;
extern cvar_t *cm_testWalk;
extern cvar_t *cm_testLength;
extern cvar_t *cm_testRadius;
extern cvar_t *cm_testAngle;

static int total_translation;
static int min_translation = 999999;
static int max_translation = -999999;
static int num_translation = 0;
static int total_rotation;
static int min_rotation = 999999;
static int max_rotation = -999999;
static int num_rotation = 0;
static vec3_t start;
static vec3_t *testend;

void CM_DebugOutput( const vec3_t origin ) {
	int i, k, t;
	char buf[128];
	vec3_t end;
	vec3_t boxAngles;
	vec3_t modelAxis[3], boxAxis[3];
	vec3_t bounds[2];
	cm_trace_t trace;

	if ( !cm_testCollision->integer ) {
		return;
	}

	testend = (vec3_t *) ii.GetMemory( cm_testTimes->integer * sizeof(vec3_t) );

	if ( cm_testReset->integer || ( cm_testWalk->integer && !VectorCompare( start, start ) ) ) {
		total_translation = total_rotation = 0;
		min_translation = min_rotation = 999999;
		max_translation = max_rotation = -999999;
		num_translation = num_rotation = 0;
		cmi.Cvar_Set( "cm_testReset", "0" );
	}

	if ( cm_testWalk->integer ) {
		VectorCopy( origin, start );
		cmi.Cvar_Set( "cm_testOrigin", va( "%1.2f %1.2f %1.2f", start[0], start[1], start[2] ) );
	} else {
		sscanf( cm_testOrigin->string, "%f %f %f", &start[0], &start[1], &start[2] );
	}

	sscanf( cm_testBox->string, "%f %f %f %f %f %f", &bounds[0][0], &bounds[0][1], &bounds[0][2],
										&bounds[1][0], &bounds[1][1], &bounds[1][2] );
	sscanf( cm_testBoxRotation->string, "%f %f %f", &boxAngles[0], &boxAngles[1], &boxAngles[2] );
	AnglesToAxis( boxAngles, boxAxis );
	AxisClear( modelAxis );

	traceModel_t itm;
	TraceModelSetupBox( &itm, bounds );
	int timer, timer1;

	if ( cm_testRandomMany->integer ) {
		// if many traces in one random direction
		for ( i = 0; i < 3; i++ ) {
			testend[0][i] = start[i] + crandom() * cm_testLength->value;
		}
		for ( k = 1; k < cm_testTimes->integer; k++ ) {
			VectorCopy( testend[0], testend[k] );
		}
	} else {
		// many traces each in a different random direction
		for ( k = 0; k < cm_testTimes->integer; k++ ) {
			for ( i = 0; i < 3; i++ ) {
				testend[k][i] = start[i] + crandom() * cm_testLength->validate;
			}
		}
	}

	// translational collision detection
	timer = ii.MilliSeconds();
	for ( i = 0; i < cm_testTimes->integer; i++ ) {
		CM_Translation( &trace, start, testend[i], &itm, boxAxis, CONTENTS_SOLID|CONTENTS_PLAYERCLIP, cm_testModel->integer, vec3_origin, modelAxis );
	}
	timer1 = ii.MilliSeconds();
	t = timer1 - timer;
	if ( t < min_translation ) min_translation = t;
	if ( t > max_translation ) max_translation = t;
	num_translation++;
	total_translation += t;
	if ( cm_testTimes->integer > 9999 ) {
		sprintf( buf, "%3dK", (int ) ( cm_testTimes->integer / 1000 ) );
	} else {
		sprintf( buf, "%4d", cm_testTimes->integer );
	}
	ii.Com_Printf("%s translations: %4d milliseconds, (min = %d, max = %d, av = %1.1f)\n", buf, t, min_translation, max_translation, (float) total_translation / num_translation );

	if ( cm_testRandomMany->integer ) {
		// if many traces in one random direction
		for ( i = 0; i < 3; i++ ) {
			testend[0][i] = start[i] + crandom() * cm_testRadius->value;
		}
		for ( k = 1; k < cm_testTimes->integer; k++ ) {
			VectorCopy( testend[0], testend[k] );
		}
	} else {
		// many traces each in a different random direction
		for ( k = 0; k < cm_testTimes->integer; k++ ) {
			for ( i = 0; i < 3; i++ ) {
				testend[k][i] = start[i] + crandom() * cm_testRadius->value;
			}
		}
	}

	if ( cm_testRotation->integer ) {
		// rotational collision detection
		vec3_t vec;
		VectorSet( vec, crandom(), crandom(), crandom() );
		VectorNormalize( vec );
		rotation_t rotation;
		
		RotationFromOriginAndVector( &rotation, vec3_origin, vec, cm_testAngle->value );

		timer = ii.MilliSeconds();
		for ( i = 0; i < cm_testTimes->integer; i++ ) {
			VectorCopy( testend[i], rotation.origin );
			CM_Rotation( &trace, start, &rotation, &itm, boxAxis, CONTENTS_SOLID|CONTENTS_PLAYERCLIP, cm_testModel->integer, vec3_origin, modelAxis );
		}
		timer1 = ii.MilliSeconds();
		t = timer1 - timer;
		if ( t < min_rotation ) min_rotation = t;
		if ( t > max_rotation ) max_rotation = t;
		num_rotation++;
		total_rotation += t;
		if ( cm_testTimes->integer > 9999 ) {
			sprintf( buf, "%3dK", (int ) ( cm_testTimes->integer / 1000 ) );
		} else {
			sprintf( buf, "%4d", cm_testTimes->integer );
		}
		ii.Com_Printf("%s rotation: %4d milliseconds, (min = %d, max = %d, av = %1.1f)\n", buf, t, min_rotation, max_rotation, (float) total_rotation / num_rotation );
	}

	ii.FreeMemory( testend );
	testend = NULL;
}
