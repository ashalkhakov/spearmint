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

#include "tr_local.h"
#include "../idlib/q_extramath.h"

#define MAX_DEBUG_LINES			16384

typedef struct debugLine_s {
	vec4_t		rgb;
	vec3_t		start;
	vec3_t		end;
	qboolean	depthTest;
	int			lifeTime;
} debugLine_t;

debugLine_t		rb_debugLines[ MAX_DEBUG_LINES ];
int				rb_numDebugLines = 0;
int				rb_debugLineTime = 0;

#define MAX_DEBUG_POLYGONS		8192

typedef struct debugPolygon_s {
	vec4_t  		rgb;
	fixedWinding_t	winding;
	qboolean    	depthTest;
	int		    	lifeTime;
} debugPolygon_t;

debugPolygon_t	rb_debugPolygons[ MAX_DEBUG_POLYGONS ];
int				rb_numDebugPolygons = 0;
int				rb_debugPolygonTime = 0;

/*
================
RB_SimpleWorldSetup
================
*/
void RB_SimpleWorldSetup( void ) {
/*
// FIXME: make this work
	backEnd.currentSpace = &backEnd.viewDef->worldSpace;
	qglLoadMatrixf( backEnd.viewDef->worldSpace.modelViewMatrix );

	backEnd.currentScissor = backEnd.viewDef->scissor;
	qglScissor( backEnd.viewDef->viewport.x1 + backEnd.currentScissor.x1, 
		backEnd.viewDef->viewport.y1 + backEnd.currentScissor.y1,
		backEnd.currentScissor.x2 + 1 - backEnd.currentScissor.x1,
		backEnd.currentScissor.y2 + 1 - backEnd.currentScissor.y1 );
*/
}

/*
================
RB_DrawBounds
================
*/
void RB_DrawBounds( const vec3_t bounds[2] ) {
	if ( bounds[0][0] > bounds[1][0] ) {
		return;
	}

	qglBegin( GL_LINE_LOOP );
	qglVertex3f( bounds[0][0], bounds[0][1], bounds[0][2] );
	qglVertex3f( bounds[0][0], bounds[1][1], bounds[0][2] );
	qglVertex3f( bounds[1][0], bounds[1][1], bounds[0][2] );
	qglVertex3f( bounds[1][0], bounds[0][1], bounds[0][2] );
	qglEnd();
	qglBegin( GL_LINE_LOOP );
	qglVertex3f( bounds[0][0], bounds[0][1], bounds[1][2] );
	qglVertex3f( bounds[0][0], bounds[1][1], bounds[1][2] );
	qglVertex3f( bounds[1][0], bounds[1][1], bounds[1][2] );
	qglVertex3f( bounds[1][0], bounds[0][1], bounds[1][2] );
	qglEnd();

	qglBegin( GL_LINES );
	qglVertex3f( bounds[0][0], bounds[0][1], bounds[0][2] );
	qglVertex3f( bounds[0][0], bounds[0][1], bounds[1][2] );

	qglVertex3f( bounds[0][0], bounds[1][1], bounds[0][2] );
	qglVertex3f( bounds[0][0], bounds[1][1], bounds[1][2] );

	qglVertex3f( bounds[1][0], bounds[0][1], bounds[0][2] );
	qglVertex3f( bounds[1][0], bounds[0][1], bounds[1][2] );

	qglVertex3f( bounds[1][0], bounds[1][1], bounds[0][2] );
	qglVertex3f( bounds[1][0], bounds[1][1], bounds[1][2] );
	qglEnd();
}

/*
================
RB_ClearDebugLines
================
*/
void RB_ClearDebugLines( int time ) {
	int			i;
	int			num;
	debugLine_t	*line;

	rb_debugLineTime = time;

	if ( !time ) {
		rb_numDebugLines = 0;
		return;
	}

	// copy any lines that still need to be drawn
	num	= 0;
	line = rb_debugLines;
	for ( i = 0 ; i < rb_numDebugLines; i++, line++ ) {
		if ( line->lifeTime > time ) {
			if ( num != i ) {
				rb_debugLines[ num ] = *line;
			}
			num++;
		}
	}
	rb_numDebugLines = num;
}

/*
================
RB_AddDebugLine
================
*/
void RB_AddDebugLine( const vec4_t color, const vec3_t start, const vec3_t end, const int lifeTime, const qboolean depthTest ) {
	debugLine_t *line;

	if ( rb_numDebugLines < MAX_DEBUG_LINES ) {
		line = &rb_debugLines[ rb_numDebugLines++ ];
        Vector4Copy( color, line->rgb );
        VectorCopy( start, line->start );
        VectorCopy( end, line->end );
		line->depthTest = depthTest;
		line->lifeTime	= rb_debugLineTime + lifeTime;
	}
}

/*
================
RB_ShowDebugLines
================
*/
void RB_ShowDebugLines( void ) {
	int			i;
	int			width;
	debugLine_t	*line;

	if ( !rb_numDebugLines ) {
		return;
	}

	// all lines are expressed in world coordinates
	RB_SimpleWorldSetup();

    // FIXME: Make this work
	//globalImages->BindNull();

	width = r_debugLineWidth->integer;
	if ( width < 1 ) {
		width = 1;
	} else if ( width > 10 ) {
		width = 10;
	}

	// draw lines
	GL_State( GLS_POLYMODE_LINE );//| GLS_DEPTHMASK ); //| GLS_SRCBLEND_ONE | GLS_DSTBLEND_ONE );
	qglLineWidth( width );

	if ( !r_debugLineDepthTest->integer ) {
		qglDisable( GL_DEPTH_TEST );
	}

	qglBegin( GL_LINES );

	line = rb_debugLines;
	for ( i = 0 ; i < rb_numDebugLines; i++, line++ ) {
		if ( !line->depthTest ) {
			qglColor3f( line->rgb[0], line->rgb[1], line->rgb[2] );
			qglVertex3fv( line->start );
			qglVertex3fv( line->end );
		}
	}
	qglEnd();

	if ( !r_debugLineDepthTest->integer ) {
		qglEnable( GL_DEPTH_TEST );
	}

	qglBegin( GL_LINES );

	line = rb_debugLines;
	for ( i = 0 ; i < rb_numDebugLines; i++, line++ ) {
		if ( line->depthTest ) {
			qglColor4f( line->rgb[0], line->rgb[1], line->rgb[2], line->rgb[3] );
			qglVertex3fv( line->start );
			qglVertex3fv( line->end );
		}
	}

	qglEnd();

	qglLineWidth( 1 );
	GL_State( GLS_DEFAULT );
}

/*
================
RB_ClearDebugPolygons
================
*/
void RB_ClearDebugPolygons( int time ) {
	int				i;
	int				num;
	debugPolygon_t	*poly;

	rb_debugPolygonTime = time;

	if ( !time ) {
		rb_numDebugPolygons = 0;
		return;
	}

	// copy any polygons that still need to be drawn
	num	= 0;

	poly = rb_debugPolygons;
	for ( i = 0 ; i < rb_numDebugPolygons; i++, poly++ ) {
		if ( poly->lifeTime > time ) {
			if ( num != i ) {
				rb_debugPolygons[ num ] = *poly;
			}
			num++;
		}
	}
	rb_numDebugPolygons = num;
}

/*
================
RB_AddDebugPolygon
================
*/
void RB_AddDebugPolygon( const vec4_t color, const fixedWinding_t *winding, const int lifeTime, const qboolean depthTest ) {
	debugPolygon_t *poly;

	if ( rb_numDebugPolygons < MAX_DEBUG_POLYGONS ) {
		poly = &rb_debugPolygons[ rb_numDebugPolygons++ ];
		Vector4Copy( color, poly->rgb );
		CopyFixedWinding( winding, &poly->winding );
		poly->depthTest = depthTest;
		poly->lifeTime	= rb_debugPolygonTime + lifeTime;
	}
}

/*
================
RB_ShowDebugPolygons
================
*/
void RB_ShowDebugPolygons( void ) {
	int				i, j;
	debugPolygon_t	*poly;

	if ( !rb_numDebugPolygons ) {
		return;
	}

	// all lines are expressed in world coordinates
	RB_SimpleWorldSetup();

    // FIXME: make this work
	//globalImages->BindNull();

	qglDisable( GL_TEXTURE_2D );
	qglDisable( GL_STENCIL_TEST );

	qglEnable( GL_DEPTH_TEST );

	if ( r_debugPolygonFilled->integer ) {
		GL_State( GLS_SRCBLEND_SRC_ALPHA | GLS_DSTBLEND_ONE_MINUS_SRC_ALPHA | GLS_DEPTHMASK_TRUE );
		qglPolygonOffset( -1, -2 );
		qglEnable( GL_POLYGON_OFFSET_FILL );
	} else {
		GL_State( GLS_POLYMODE_LINE );
		qglPolygonOffset( -1, -2 );
		qglEnable( GL_POLYGON_OFFSET_LINE );
	}

	poly = rb_debugPolygons;
	for ( i = 0 ; i < rb_numDebugPolygons; i++, poly++ ) {
//		if ( !poly->depthTest ) {

			qglColor4f( poly->rgb[0], poly->rgb[1], poly->rgb[2], poly->rgb[3] );

			qglBegin( GL_POLYGON );

			for ( j = 0; j < poly->winding.numPoints; j++) {
				qglVertex3fv( poly->winding.points[j] );
			}

			qglEnd();
//		}
	}

	GL_State( GLS_DEFAULT );

	if ( r_debugPolygonFilled->integer ) {
		qglDisable( GL_POLYGON_OFFSET_FILL );
	} else {
		qglDisable( GL_POLYGON_OFFSET_LINE );
	}

	qglDepthRange( 0, 1 );
	GL_State( GLS_DEFAULT );
}

/*
================
RB_DebugArrow
================
*/
void RB_DebugArrow( const vec4_t color, const vec3_t start, const vec3_t end, int size, const int lifetime ) {
	vec3_t forward, right, up, v1, v2;
	float a, s;
	int i;
	static float arrowCos[40];
	static float arrowSin[40];
	static int arrowStep;

	RB_AddDebugLine( color, start, end, lifetime, qfalse );

	if ( r_debugArrowStep->integer <= 10 ) {
		return;
	}
	// calculate sine and cosine when step size changes
	if ( arrowStep != r_debugArrowStep->integer ) {
		arrowStep = r_debugArrowStep->integer;
		for (i = 0, a = 0; a < 360.0f; a += arrowStep, i++) {
			arrowCos[i] = cos( DEG2RAD( a ) ); // TODO: Cos16
			arrowSin[i] = sin( DEG2RAD( a ) ); // TODO: Sin16
		}
		arrowCos[i] = arrowCos[0];
		arrowSin[i] = arrowSin[0];
	}
	// draw a nice arrow
	VectorSubtract( end, start, forward );
	VectorNormalize( forward );
	NormalVectors( forward, right, up );
	for (i = 0, a = 0; a < 360.0f; a += arrowStep, i++) {
		s = 0.5f * size * arrowCos[i];
        VectorMA( end, -size, forward, v1 );
        VectorMA( v1, s, right, v1 );
		s = 0.5f * size * arrowSin[i];
        VectorMA( v1, s, up, v1 );

		s = 0.5f * size * arrowCos[i+1];
        VectorMA( end, -size, forward, v2 );
        VectorMA( v2, s, right, v2 );
		s = 0.5f * size * arrowSin[i+1];
        VectorMA( v2, s, up, v2 );

		RB_AddDebugLine( color, v1, end, lifetime, qfalse );
		RB_AddDebugLine( color, v1, v2, lifetime, qfalse );
	}
}

/*
====================
RB_DebugWinding
====================
*/
void RB_DebugWinding( const vec4_t color, const fixedWinding_t *w, const vec3_t origin, const vec3_t axis[3], const int lifetime, const qboolean depthTest ) {
	int i;
	vec3_t point, lastPoint;

	if ( w->numPoints < 2 ) {
		return;
	}

    MatrixRotateVector( w->points[w->numPoints-1], axis, lastPoint );
    VectorAdd( lastPoint, origin, lastPoint );
	for ( i = 0; i < w->numPoints; i++ ) {
        MatrixRotateVector( w->points[i], axis, point );
        VectorAdd( point, origin, point );
		RB_AddDebugLine( color, lastPoint, point, lifetime, depthTest );
        VectorCopy( point, lastPoint );
	}
}

/*
====================
RB_DebugCircle
====================
*/
void RB_DebugCircle( const vec4_t color, const vec3_t origin, const vec3_t dir, const float radius, const int numSteps, const int lifetime, const qboolean depthTest ) {
	int i;
	float a;
	vec3_t left, up, point, lastPoint;

	VectorOrthogonalBasis( dir, left, up );
    VectorScale( left, radius, left );
    VectorScale( up, radius, up );
    VectorAdd( origin, up, lastPoint );
	for ( i = 1; i <= numSteps; i++ ) {
		a = 2.0 * M_PI * i / numSteps;
        // TODO: Sin16/Cos16
        VectorMA( origin, sin( a ), left, point );
        VectorMA( point, cos( a ), up, point );
		RB_AddDebugLine( color, lastPoint, point, lifetime, depthTest );
        VectorCopy( point, lastPoint );
	}
}

/*
============
RB_DebugSphere
============
*/
void RB_DebugSphere( const vec4_t color, const vec3_t origin, const float radius, const int lifetime, const qboolean depthTest ) {
	int i, j, n, num;
	float s, c;
	vec3_t p, lastp;
    vec3_t lastArray[24]; // 360 / 15

	num = 24; // 360 / 15;
    VectorCopy( origin, lastArray[0] );
    lastArray[0][2] += radius;
	for ( n = 1; n < num; n++ ) {
        VectorCopy( lastArray[0], lastArray[n] );
	}

	for ( i = 15; i <= 360; i += 15 ) {
        // TODO: Sin16/Cos16
		s = sin( DEG2RAD(i) );
		c = cos( DEG2RAD(i) );
		lastp[0] = origin[0];
		lastp[1] = origin[1] + radius * s;
		lastp[2] = origin[2] + radius * c;
		for ( n = 0, j = 15; j <= 360; j += 15, n++ ) {
            // TODO: Sin16/Cos16
			p[0] = origin[0] + sin( DEG2RAD(j) ) * radius * s;
			p[1] = origin[1] + cos( DEG2RAD(j) ) * radius * s;
			p[2] = lastp[2];

			RB_AddDebugLine( color, lastp, p, lifetime,depthTest );
			RB_AddDebugLine( color, lastp, lastArray[n], lifetime, depthTest );

            VectorCopy( lastp, lastArray[n] );
            VectorCopy( p, lastp );
		}
	}
}

/*
====================
RB_DebugBounds
====================
*/
void RB_DebugBounds( const vec4_t color, const vec3_t bounds[2], const vec3_t org, const int lifetime ) {
	int i;
	vec3_t v[8];

	if ( bounds[0][0] > bounds[1][0] ) {
		return;
	}

	for ( i = 0; i < 8; i++ ) {
		v[i][0] = org[0] + bounds[(i^(i>>1))&1][0];
		v[i][1] = org[1] + bounds[(i>>1)&1][1];
		v[i][2] = org[2] + bounds[(i>>2)&1][2];
	}
	for ( i = 0; i < 4; i++ ) {
		RB_AddDebugLine( color, v[i], v[(i+1)&3], lifetime, qfalse );
		RB_AddDebugLine( color, v[4+i], v[4+((i+1)&3)], lifetime, qfalse );
		RB_AddDebugLine( color, v[i], v[4+i], lifetime, qfalse );
	}
}

/*
================
RB_DebugAxis
================
*/
void RB_DebugAxis( const vec3_t origin, const vec3_t axis[3] ) {
	vec3_t start, end;
    
    VectorCopy( origin, start );
    VectorMA( start, 20.0f, axis[0], end );

	RB_DebugArrow( colorWhite, start, end, 2, 0 );
    VectorMA( start, -20.0f, axis[0], end );
	RB_DebugArrow( colorWhite, start, end, 2, 0 );
    VectorMA( start, 20.0f, axis[1], end );
	RB_DebugArrow( colorGreen, start, end, 2, 0 );
    VectorMA( start, -20.0f, axis[1], end );
	RB_DebugArrow( colorGreen, start, end, 2, 0 );
    VectorMA( start, 20.0f, axis[2], end );
	RB_DebugArrow( colorBlue, start, end, 2, 0 );
    VectorMA( start, -20.0f, axis[2], end );
	RB_DebugArrow( colorBlue, start, end, 2, 0 );
}

/*
=================
RB_RenderDebugTools
=================
*/
void RB_RenderDebugTools( void ) {
	// don't do anything if this was a 2D rendering
	if ( tr.refdef.rdflags & RDF_NOWORLDMODEL ) {
		return;
	}

	GL_State( GLS_DEFAULT );
    /*
	backEnd.currentScissor = backEnd.viewDef->scissor;
	qglScissor( backEnd.viewDef->viewport.x1 + backEnd.currentScissor.x1, 
		backEnd.viewDef->viewport.y1 + backEnd.currentScissor.y1,
		backEnd.currentScissor.x2 + 1 - backEnd.currentScissor.x1,
		backEnd.currentScissor.y2 + 1 - backEnd.currentScissor.y1 );
    */

	RB_ShowDebugLines();
	//RB_ShowDebugText();
	RB_ShowDebugPolygons();
}

/*
=================
RB_ShutdownDebugTools
=================
*/
void RB_ShutdownDebugTools( void ) {
	for ( int i = 0; i < MAX_DEBUG_POLYGONS; i++ ) {
		ClearFixedWinding( &rb_debugPolygons[i].winding );
	}
}
