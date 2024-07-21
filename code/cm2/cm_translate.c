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

/*
===============================================================================

Collision detection for translational motion

===============================================================================
*/

/*
================
CM_TranslateEdgeThroughEdge

  calculates fraction of the translation completed at which the edges collide
================
*/
ID_INLINE int CM_TranslateEdgeThroughEdge( vec3_t cross, vec6_t l1, vec6_t l2, float *fraction ) {

	float d, t;

	/*

	a = start of line
	b = end of line
	dir = movement direction
	l1 = pluecker coordinate for line
	l2 = pluecker coordinate for edge we might collide with
	a+dir = start of line after movement
	b+dir = end of line after movement
	t = scale factor
	solve pluecker inner product for t of line (a+t*dir : b+t*dir) and line l2

	v[0] = (a[0]+t*dir[0]) * (b[1]+t*dir[1]) - (b[0]+t*dir[0]) * (a[1]+t*dir[1]);
	v[1] = (a[0]+t*dir[0]) * (b[2]+t*dir[2]) - (b[0]+t*dir[0]) * (a[2]+t*dir[2]);
	v[2] = (a[0]+t*dir[0]) - (b[0]+t*dir[0]);
	v[3] = (a[1]+t*dir[1]) * (b[2]+t*dir[2]) - (b[1]+t*dir[1]) * (a[2]+t*dir[2]);
	v[4] = (a[2]+t*dir[2]) - (b[2]+t*dir[2]);
	v[5] = (b[1]+t*dir[1]) - (a[1]+t*dir[1]);

	l2[0] * v[4] + l2[1] * v[5] + l2[2] * v[3] + l2[4] * v[0] + l2[5] * v[1] + l2[3] * v[2] = 0;

	solve t

	v[0] = (a[0]+t*dir[0]) * (b[1]+t*dir[1]) - (b[0]+t*dir[0]) * (a[1]+t*dir[1]);
	v[0] = (a[0]*b[1]) + a[0]*t*dir[1] + b[1]*t*dir[0] + (t*t*dir[0]*dir[1]) -
			((b[0]*a[1]) + b[0]*t*dir[1] + a[1]*t*dir[0] + (t*t*dir[0]*dir[1]));
	v[0] = a[0]*b[1] + a[0]*t*dir[1] + b[1]*t*dir[0] - b[0]*a[1] - b[0]*t*dir[1] - a[1]*t*dir[0];

	v[1] = (a[0]+t*dir[0]) * (b[2]+t*dir[2]) - (b[0]+t*dir[0]) * (a[2]+t*dir[2]);
	v[1] = (a[0]*b[2]) + a[0]*t*dir[2] + b[2]*t*dir[0] + (t*t*dir[0]*dir[2]) -
			((b[0]*a[2]) + b[0]*t*dir[2] + a[2]*t*dir[0] + (t*t*dir[0]*dir[2]));
	v[1] = a[0]*b[2] + a[0]*t*dir[2] + b[2]*t*dir[0] - b[0]*a[2] - b[0]*t*dir[2] - a[2]*t*dir[0];

	v[2] = (a[0]+t*dir[0]) - (b[0]+t*dir[0]);
	v[2] = a[0] - b[0];

	v[3] = (a[1]+t*dir[1]) * (b[2]+t*dir[2]) - (b[1]+t*dir[1]) * (a[2]+t*dir[2]);
	v[3] = (a[1]*b[2]) + a[1]*t*dir[2] + b[2]*t*dir[1] + (t*t*dir[1]*dir[2]) -
			((b[1]*a[2]) + b[1]*t*dir[2] + a[2]*t*dir[1] + (t*t*dir[1]*dir[2]));
	v[3] = a[1]*b[2] + a[1]*t*dir[2] + b[2]*t*dir[1] - b[1]*a[2] - b[1]*t*dir[2] - a[2]*t*dir[1];

	v[4] = (a[2]+t*dir[2]) - (b[2]+t*dir[2]);
	v[4] = a[2] - b[2];

	v[5] = (b[1]+t*dir[1]) - (a[1]+t*dir[1]);
	v[5] = b[1] - a[1];


	v[0] = a[0]*b[1] + a[0]*t*dir[1] + b[1]*t*dir[0] - b[0]*a[1] - b[0]*t*dir[1] - a[1]*t*dir[0];
	v[1] = a[0]*b[2] + a[0]*t*dir[2] + b[2]*t*dir[0] - b[0]*a[2] - b[0]*t*dir[2] - a[2]*t*dir[0];
	v[2] = a[0] - b[0];
	v[3] = a[1]*b[2] + a[1]*t*dir[2] + b[2]*t*dir[1] - b[1]*a[2] - b[1]*t*dir[2] - a[2]*t*dir[1];
	v[4] = a[2] - b[2];
	v[5] = b[1] - a[1];

	v[0] = (a[0]*dir[1] + b[1]*dir[0] - b[0]*dir[1] - a[1]*dir[0]) * t + a[0]*b[1] - b[0]*a[1];
	v[1] = (a[0]*dir[2] + b[2]*dir[0] - b[0]*dir[2] - a[2]*dir[0]) * t + a[0]*b[2] - b[0]*a[2];
	v[2] = a[0] - b[0];
	v[3] = (a[1]*dir[2] + b[2]*dir[1] - b[1]*dir[2] - a[2]*dir[1]) * t + a[1]*b[2] - b[1]*a[2];
	v[4] = a[2] - b[2];
	v[5] = b[1] - a[1];

	l2[4] * (a[0]*dir[1] + b[1]*dir[0] - b[0]*dir[1] - a[1]*dir[0]) * t + l2[4] * (a[0]*b[1] - b[0]*a[1])
		+ l2[5] * (a[0]*dir[2] + b[2]*dir[0] - b[0]*dir[2] - a[2]*dir[0]) * t + l2[5] * (a[0]*b[2] - b[0]*a[2])
		+ l2[3] * (a[0] - b[0])
		+ l2[2] * (a[1]*dir[2] + b[2]*dir[1] - b[1]*dir[2] - a[2]*dir[1]) * t + l2[2] * (a[1]*b[2] - b[1]*a[2])
		+ l2[0] * (a[2] - b[2])
		+ l2[1] * (b[1] - a[1]) = 0

	t = (- l2[4] * (a[0]*b[1] - b[0]*a[1]) -
			l2[5] * (a[0]*b[2] - b[0]*a[2]) -
			l2[3] * (a[0] - b[0]) -
			l2[2] * (a[1]*b[2] - b[1]*a[2]) -
			l2[0] * (a[2] - b[2]) -
			l2[1] * (b[1] - a[1])) /
				(l2[4] * (a[0]*dir[1] + b[1]*dir[0] - b[0]*dir[1] - a[1]*dir[0]) +
				l2[5] * (a[0]*dir[2] + b[2]*dir[0] - b[0]*dir[2] - a[2]*dir[0]) +
				l2[2] * (a[1]*dir[2] + b[2]*dir[1] - b[1]*dir[2] - a[2]*dir[1]));

	d = l2[4] * (a[0]*dir[1] + b[1]*dir[0] - b[0]*dir[1] - a[1]*dir[0]) +
		l2[5] * (a[0]*dir[2] + b[2]*dir[0] - b[0]*dir[2] - a[2]*dir[0]) +
		l2[2] * (a[1]*dir[2] + b[2]*dir[1] - b[1]*dir[2] - a[2]*dir[1]);

	t = - ( l2[4] * (a[0]*b[1] - b[0]*a[1]) +
			l2[5] * (a[0]*b[2] - b[0]*a[2]) +
			l2[3] * (a[0] - b[0]) +
			l2[2] * (a[1]*b[2] - b[1]*a[2]) +
			l2[0] * (a[2] - b[2]) +
			l2[1] * (b[1] - a[1]));
	t /= d;

	MrE pats Pluecker on the head.. good monkey

	edgeDir = a - b;
	d = l2[4] * (edgeDir[0]*dir[1] - edgeDir[1]*dir[0]) +
		l2[5] * (edgeDir[0]*dir[2] - edgeDir[2]*dir[0]) +
		l2[2] * (edgeDir[1]*dir[2] - edgeDir[2]*dir[1]);
	*/

	d = l2[4] * cross[0] + l2[5] * cross[1] + l2[2] * cross[2];

	if ( d == 0.0f ) {
		*fraction = 1.0f;
		// no collision ever
		return qfalse;
	}

	t = -PlueckerPermutedInnerProduct( l1, l2 );
	// if the lines cross each other to begin with
	if ( t == 0.0f ) {
		*fraction = 0.0f;
		return qtrue;
	}
	// fraction of movement at the time the lines cross each other
	*fraction = t / d;
	return qtrue;
}

/*
================
CM_AddContact
================
*/
static ID_INLINE void CM_AddContact( cm_traceWork_t *tw ) {

	if ( tw->numContacts >= tw->maxContacts ) {
		return;
	}
	// copy contact information from cm_trace_t
	tw->contacts[tw->numContacts] = tw->trace.c;
	tw->numContacts++;
	// set fraction back to 1 to find all other contacts
	tw->trace.fraction = 1.0f;
}

/*
================
CM_SetVertexSidedness

  stores for the given model vertex at which side of one of the trm edges it passes
================
*/
static ID_INLINE void CM_SetVertexSidedness( cm_vertex_t *v, const vec6_t vpl, const vec6_t epl, const int bitNum ) {
	if ( !(v->sideSet & (1<<bitNum)) ) {
		float fl;
		fl = PlueckerPermutedInnerProduct( vpl, epl );
		v->side = (v->side & ~(1<<bitNum)) | (FLOATSIGNBITSET(fl) << bitNum);
		v->sideSet |= (1 << bitNum);
	}
}

/*
================
CM_SetEdgeSidedness

  stores for the given model edge at which side one of the trm vertices
================
*/
static ID_INLINE void CM_SetEdgeSidedness( cm_edge_t *edge, const vec6_t vpl, const vec6_t epl, const int bitNum ) {
	if ( !(edge->sideSet & (1<<bitNum)) ) {
		float fl;
		fl = PlueckerPermutedInnerProduct( vpl, epl );
		edge->side = (edge->side & ~(1<<bitNum)) | (FLOATSIGNBITSET(fl) << bitNum);
		edge->sideSet |= (1 << bitNum);
	}
}

/*
================
CM_TranslateTrmEdgeThroughPolygon
================
*/
void CM_TranslateTrmEdgeThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmEdge_t *trmEdge ) {
	int i, edgeNum;
	float f1, f2, dist, d1, d2;
	vec3_t start, end, normal, tmp1, tmp2;
	cm_edge_t *edge;
	cm_vertex_t *v1, *v2;
	vec6_t *pl, epsPl;

	// check edges for a collision
	for ( i = 0; i < poly->numEdges; i++) {
		edgeNum = poly->edges[i];
		edge = tw->model->edges + abs(edgeNum);
		// if this edge is already checked
		if ( edge->checkcount == cmLocal.checkCount ) {
			continue;
		}
		// can never collide with internal edges
		if ( edge->internal ) {
			continue;
		}
		pl = &tw->polygonEdgePlueckerCache[i];
		// get the sides at which the trm edge vertices pass the polygon edge
		CM_SetEdgeSidedness( edge, *pl, tw->vertices[trmEdge->vertexNum[0]].pl, trmEdge->vertexNum[0] );
		CM_SetEdgeSidedness( edge, *pl, tw->vertices[trmEdge->vertexNum[1]].pl, trmEdge->vertexNum[1] );
		// if the trm edge start and end vertex do not pass the polygon edge at different sides
		if ( !(((edge->side >> trmEdge->vertexNum[0]) ^ (edge->side >> trmEdge->vertexNum[1])) & 1) ) {
			continue;
		}
		// get the sides at which the polygon edge vertices pass the trm edge
		v1 = tw->model->vertices + edge->vertexNum[INTSIGNBITSET(edgeNum)];
		CM_SetVertexSidedness( v1, tw->polygonVertexPlueckerCache[i], trmEdge->pl, trmEdge->bitNum );
		v2 = tw->model->vertices + edge->vertexNum[INTSIGNBITNOTSET(edgeNum)];
		CM_SetVertexSidedness( v2, tw->polygonVertexPlueckerCache[i+1], trmEdge->pl, trmEdge->bitNum );
		// if the polygon edge start and end vertex do not pass the trm edge at different sides
		if ( !((v1->side ^ v2->side) & (1<<trmEdge->bitNum)) ) {
			continue;
		}
		// if there is no possible collision between the trm edge and the polygon edge
		if ( !CM_TranslateEdgeThroughEdge( trmEdge->cross, trmEdge->pl, *pl, &f1 ) ) {
			continue;
		}
		// if moving away from edge
		if ( f1 < 0.0f ) {
			continue;
		}

		// pluecker coordinate for epsilon expanded edge
        VectorMA( tw->model->vertices[edge->vertexNum[0]].p, CM_CLIP_EPSILON, edge->normal, start );
        VectorMA( tw->model->vertices[edge->vertexNum[1]].p, CM_CLIP_EPSILON, edge->normal, end );
        PlueckerFromLine( epsPl, start, end );
		// calculate collision fraction with epsilon expanded edge
		if ( !CM_TranslateEdgeThroughEdge( trmEdge->cross, trmEdge->pl, epsPl, &f2 ) ) {
			continue;
		}
		// if no collision with epsilon edge or moving away from edge
		if ( f2 > 1.0f || f1 < f2 ) {
			continue;
		}

		if ( f2 < 0.0f ) {
			f2 = 0.0f;
		}

		if ( f2 < tw->trace.fraction ) {
			tw->trace.fraction = f2;
			// create plane with normal vector orthogonal to both the polygon edge and the trm edge
            VectorCopy( tw->model->vertices[edge->vertexNum[0]].p, start );
			VectorCopy( tw->model->vertices[edge->vertexNum[1]].p, end );
            VectorSubtract( end, start, tmp1 );
            VectorSubtract( trmEdge->end, trmEdge->start, tmp2 );
            CrossProduct( tmp1, tmp2, tw->trace.c.normal );
			// FIXME: do this normalize when we know the first collision
			VectorNormalize( tw->trace.c.normal );
			tw->trace.c.dist = DotProduct( tw->trace.c.normal, start );
			// make sure the collision plane faces the trace model
			if ( DotProduct( tw->trace.c.normal, trmEdge->start ) - tw->trace.c.dist < 0.0f ) {
				VectorNegate( tw->trace.c.normal, tw->trace.c.normal );
				tw->trace.c.dist = -tw->trace.c.dist;
			}
			tw->trace.c.contents = poly->contents;
			tw->trace.c.material = poly->material;
			tw->trace.c.type = CONTACT_EDGE;
			tw->trace.c.modelFeature = edgeNum;
			tw->trace.c.trmFeature = trmEdge - tw->edges;
			// calculate collision point
			normal[0] = trmEdge->cross[2];
			normal[1] = -trmEdge->cross[1];
			normal[2] = trmEdge->cross[0];
			dist = DotProduct( normal, trmEdge->start );
			d1 = DotProduct( normal, start ) - dist;
			d2 = DotProduct( normal, end ) - dist;
			f1 = d1 / ( d1 - d2 );
			//assert( f1 >= 0.0f && f1 <= 1.0f );
            VectorSubtract( end, start, tmp1 );
            VectorMA( start, f1, tmp1, tw->trace.c.point );
			// if retrieving contacts
			if ( tw->getContacts ) {
				CM_AddContact( tw );
			}
		}
	}
}

/*
================
CM_TranslationPlaneFraction
================
*/

#if 0

float CM_TranslationPlaneFraction( idPlane &plane, idVec3 &start, idVec3 &end ) {
	float d1, d2;

	d2 = plane.Distance( end );
	// if the end point is closer to the plane than an epsilon we still take it for a collision
	if ( d2 >= CM_CLIP_EPSILON ) {
		return 1.0f;
	}
	d1 = plane.Distance( start );

	// if completely behind the polygon
	if ( d1 <= 0.0f ) {
		return 1.0f;
	}
	// leaves polygon
	if ( d1 <= d2 ) {
		return 1.0f;
	}
	return (d1-CM_CLIP_EPSILON) / (d1-d2);
}

#else

float CM_TranslationPlaneFraction( plane_t plane, vec3_t start, vec3_t end ) {
	float d1, d2, d2eps;

	d2 = PlaneDistance( plane, end );
	// if the end point is closer to the plane than an epsilon we still take it for a collision
	// if ( d2 >= CM_CLIP_EPSILON ) {
	d2eps = d2 - CM_CLIP_EPSILON;
	if ( FLOATSIGNBITNOTSET(d2eps) ) {
		return 1.0f;
	}
	d1 = PlaneDistance( plane, start );

	// if completely behind the polygon
	if ( FLOATSIGNBITSET(d1) ) {
		return 1.0f;
	}
	// if going towards the front of the plane and
	// the start and end point are not at equal distance from the plane
	// if ( d1 > d2 )
	d2 = d1 - d2;
	if ( d2 <= 0.0f ) {
		return 1.0f;
	}
	return (d1-CM_CLIP_EPSILON) / d2;
}

#endif

/*
================
CM_TranslateTrmVertexThroughPolygon
================
*/
void CM_TranslateTrmVertexThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmVertex_t *v, int bitNum ) {
	int i, edgeNum;
	float f;
	cm_edge_t *edge;
    vec3_t      dir;

	f = CM_TranslationPlaneFraction( poly->plane, v->p, v->endp );
	if ( f < tw->trace.fraction ) {

		for ( i = 0; i < poly->numEdges; i++ ) {
			edgeNum = poly->edges[i];
			edge = tw->model->edges + abs(edgeNum);
			CM_SetEdgeSidedness( edge, tw->polygonEdgePlueckerCache[i], v->pl, bitNum );
			if ( INTSIGNBITSET(edgeNum) ^ ((edge->side >> bitNum) & 1) ) {
				return;
			}
		}
		if ( f < 0.0f ) {
			f = 0.0f;
		}
		tw->trace.fraction = f;
		// collision plane is the polygon plane
		VectorCopy( poly->plane, tw->trace.c.normal );
		tw->trace.c.dist = PlaneGetDist( poly->plane );
		tw->trace.c.contents = poly->contents;
		tw->trace.c.material = poly->material;
		tw->trace.c.type = CONTACT_TRMVERTEX;
		tw->trace.c.modelFeature = *(int *)(&poly);
		tw->trace.c.trmFeature = v - tw->vertices;
        VectorSubtract( v->endp, v->p, dir );
        VectorMA( v->p, tw->trace.fraction, dir, tw->trace.c.point );
		// if retrieving contacts
		if ( tw->getContacts ) {
			CM_AddContact( tw );
			// no need to store the trm vertex more than once as a contact
			v->used = qfalse;
		}
	}
}

/*
================
CM_TranslatePointThroughPolygon
================
*/
void CM_TranslatePointThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmVertex_t *v ) {
	int i, edgeNum;
	float f;
	cm_edge_t *edge;
	vec6_t pl;
    vec3_t dir;

	f = CM_TranslationPlaneFraction( poly->plane, v->p, v->endp );
	if ( f < tw->trace.fraction ) {

		for ( i = 0; i < poly->numEdges; i++ ) {
			edgeNum = poly->edges[i];
			edge = tw->model->edges + abs(edgeNum);
			// if we didn't yet calculate the sidedness for this edge
			if ( edge->checkcount != cmLocal.checkCount ) {
				float fl;
				edge->checkcount = cmLocal.checkCount;
                PlueckerFromLine( pl, tw->model->vertices[edge->vertexNum[0]].p, tw->model->vertices[edge->vertexNum[1]].p );
				fl = PlueckerPermutedInnerProduct( v->pl, pl );
				edge->side = FLOATSIGNBITSET(fl);
			}
			// if the point passes the edge at the wrong side
			//if ( (edgeNum > 0) == edge->side ) {
			if ( INTSIGNBITSET(edgeNum) ^ edge->side ) {
				return;
			}
		}
		if ( f < 0.0f ) {
			f = 0.0f;
		}
		tw->trace.fraction = f;
		// collision plane is the polygon plane
		VectorCopy( poly->plane, tw->trace.c.normal );
		tw->trace.c.dist = PlaneGetDist( poly->plane );
		tw->trace.c.contents = poly->contents;
		tw->trace.c.material = poly->material;
		tw->trace.c.type = CONTACT_TRMVERTEX;
		tw->trace.c.modelFeature = *(int *)(&poly);
		tw->trace.c.trmFeature = v - tw->vertices;
        VectorSubtract( v->endp, v->p, dir );
        VectorMA( v->p, tw->trace.fraction, dir, tw->trace.c.point );
		// if retrieving contacts
		if ( tw->getContacts ) {
			CM_AddContact( tw );
			// no need to store the trm vertex more than once as a contact
			v->used = qfalse;
		}
	}
}

/*
================
CM_TranslateVertexThroughTrmPolygon
================
*/
void CM_TranslateVertexThroughTrmPolygon( cm_traceWork_t *tw, cm_trmPolygon_t *trmpoly, cm_polygon_t *poly, cm_vertex_t *v, vec3_t endp, vec6_t pl ) {
	int i, edgeNum;
	float f;
	cm_trmEdge_t *edge;
    vec3_t dir;

	f = CM_TranslationPlaneFraction( trmpoly->plane, v->p, endp );
	if ( f < tw->trace.fraction ) {

		for ( i = 0; i < trmpoly->numEdges; i++ ) {
			edgeNum = trmpoly->edges[i];
			edge = tw->edges + abs(edgeNum);

			CM_SetVertexSidedness( v, pl, edge->pl, edge->bitNum );
			if ( INTSIGNBITSET(edgeNum) ^ ((v->side >> edge->bitNum) & 1) ) {
				return;
			}
		}
		if ( f < 0.0f ) {
			f = 0.0f;
		}
		tw->trace.fraction = f;
		// collision plane is the inverse trm polygon plane
        VectorNegate( trmpoly->plane, tw->trace.c.normal );
		tw->trace.c.dist = -PlaneGetDist( trmpoly->plane );
		tw->trace.c.contents = poly->contents;
		tw->trace.c.material = poly->material;
		tw->trace.c.type = CONTACT_MODELVERTEX;
		tw->trace.c.modelFeature = v - tw->model->vertices;
		tw->trace.c.trmFeature = trmpoly - tw->polys;
        VectorSubtract( endp, v->p, dir );
        VectorMA( v->p, tw->trace.fraction, dir, tw->trace.c.point );
		// if retrieving contacts
		if ( tw->getContacts ) {
			CM_AddContact( tw );
		}
	}
}

/*
================
CM_TranslateTrmThroughPolygon

  returns true if the polygon blocks the complete translation
================
*/
qboolean CM_TranslateTrmThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *p ) {
	int i, j, k, edgeNum;
	float fraction, d;
	vec3_t endp;
	vec6_t *pl;
	cm_trmVertex_t *bv;
	cm_trmEdge_t *be;
	cm_trmPolygon_t *bp;
	cm_vertex_t *v;
	cm_edge_t *e;
    vec3_t  tmp;

	// if already checked this polygon
	if ( p->checkcount == cmLocal.checkCount ) {
		return qfalse;
	}
	p->checkcount = cmLocal.checkCount;

	// if this polygon does not have the right contents behind it
	if ( !(p->contents & tw->contents) ) {
		return qfalse;
	}

	// if the the trace bounds do not intersect the polygon bounds
	if ( !BoundsIntersect( tw->bounds[0], tw->bounds[1], p->bounds[0], p->bounds[1] ) ) {
		return qfalse;
	}

	// only collide with the polygon if approaching at the front
	if ( DotProduct( p->plane, tw->dir ) > 0.0f ) {
		return qfalse;
	}

	// if the polygon is too far from the first heart plane
	d = BoundsDistanceToPlane( p->bounds, tw->heartPlane1 );
	if ( fabs(d) > tw->maxDistFromHeartPlane1 ) {
		return qfalse;
	}

	// if the polygon is too far from the second heart plane
	d = BoundsDistanceToPlane( p->bounds, tw->heartPlane2 );
	if ( fabs(d) > tw->maxDistFromHeartPlane2 ) {
		return qfalse;
	}
	fraction = tw->trace.fraction;

	// fast point trace
	if ( tw->pointTrace ) {
		CM_TranslatePointThroughPolygon( tw, p, &tw->vertices[0] );
	}
	else {

		// trace bounds should cross polygon plane
		switch ( BoxOnPlaneSideSlow( tw->bounds[0], tw->bounds[1], p->plane, ON_EPSILON ) ) {
			case PLANESIDE_CROSS:
				break;
			case PLANESIDE_FRONT:
				if ( tw->model->isConvex ) {
					tw->quickExit = qtrue;
					return qtrue;
				}
			default:
				return qfalse;
		}

		// calculate pluecker coordinates for the polygon edges and polygon vertices
		for ( i = 0; i < p->numEdges; i++ ) {
			edgeNum = p->edges[i];
			e = tw->model->edges + abs(edgeNum);
			// reset sidedness cache if this is the first time we encounter this edge during this trace
			if ( e->checkcount != cmLocal.checkCount ) {
				e->sideSet = 0;
			}
			// pluecker coordinate for edge
			PlueckerFromLine( tw->polygonEdgePlueckerCache[i], tw->model->vertices[e->vertexNum[0]].p,
														tw->model->vertices[e->vertexNum[1]].p );

			v = &tw->model->vertices[e->vertexNum[INTSIGNBITSET(edgeNum)]];
			// reset sidedness cache if this is the first time we encounter this vertex during this trace
			if ( v->checkcount != cmLocal.checkCount ) {
				v->sideSet = 0;
			}
			// pluecker coordinate for vertex movement vector
            VectorNegate( tw->dir, tmp );
			PlueckerFromRay( tw->polygonVertexPlueckerCache[i], v->p, tmp );
		}
		// copy first to last so we can easily cycle through for the edges
        PlueckerCopy( tw->polygonVertexPlueckerCache[0], tw->polygonVertexPlueckerCache[p->numEdges] );

		// trace trm vertices through polygon
		for ( i = 0; i < tw->numVerts; i++ ) {
			bv = tw->vertices + i;
			if ( bv->used ) {
				CM_TranslateTrmVertexThroughPolygon( tw, p, bv, i );
			}
		}

		// trace trm edges through polygon
		for ( i = 1; i <= tw->numEdges; i++ ) {
			be = tw->edges + i;
			if ( be->used ) {
				CM_TranslateTrmEdgeThroughPolygon( tw, p, be);
			}
		}

		// trace all polygon vertices through the trm
		for ( i = 0; i < p->numEdges; i++ ) {
			edgeNum = p->edges[i];
			e = tw->model->edges + abs(edgeNum);

			if ( e->checkcount == cmLocal.checkCount ) {
				continue;
			}
			// set edge check count
			e->checkcount = cmLocal.checkCount;
			// can never collide with internal edges
			if ( e->internal ) {
				continue;
			}
			// got to check both vertices because we skip internal edges
			for ( k = 0; k < 2; k++ ) {

				v = tw->model->vertices + e->vertexNum[k ^ INTSIGNBITSET(edgeNum)];
				// if this vertex is already checked
				if ( v->checkcount == cmLocal.checkCount ) {
					continue;
				}
				// set vertex check count
				v->checkcount = cmLocal.checkCount;

				// if the vertex is outside the trace bounds
				if ( !BoundsIntersectPoint( tw->bounds[0], tw->bounds[1], v->p ) ) {
					continue;
				}

				// vertex end point after movement
                VectorSubtract( v->p, tw->dir, endp );
				// pluecker coordinate for vertex movement vector
				pl = &tw->polygonVertexPlueckerCache[i+k];

				for ( j = 0; j < tw->numPolys; j++ ) {
					bp = tw->polys + j;
					if ( bp->used ) {
						CM_TranslateVertexThroughTrmPolygon( tw, bp, p, v, endp, *pl );
					}
				}
			}
		}
	}

	// if there was a collision with this polygon and we are not retrieving contacts
	if ( tw->trace.fraction < fraction && !tw->getContacts ) {
		fraction = tw->trace.fraction;
		VectorMA( tw->start, fraction, tw->dir, endp );
		// decrease bounds
		for ( i = 0; i < 3; i++ ) {
			if ( tw->start[i] < endp[i] ) {
				tw->bounds[0][i] = tw->start[i] + tw->size[0][i] - CM_BOX_EPSILON;
				tw->bounds[1][i] = endp[i] + tw->size[1][i] + CM_BOX_EPSILON;
			}
			else {
				tw->bounds[0][i] = endp[i] + tw->size[0][i] - CM_BOX_EPSILON;
				tw->bounds[1][i] = tw->start[i] + tw->size[1][i] + CM_BOX_EPSILON;
			}
		}
	}

	return ( tw->trace.fraction == 0.0f );
}

/*
================
CM_SetupTrm
================
*/
void CM_SetupTrm( cm_traceWork_t *tw, const traceModel_t *trm ) {
	int i, j;

	// vertices
	tw->numVerts = trm->numVerts;
	for ( i = 0; i < trm->numVerts; i++ ) {
		VectorCopy( trm->verts[i], tw->vertices[i].p );
		tw->vertices[i].used = qfalse;
	}
	// edges
	tw->numEdges = trm->numEdges;
	for ( i = 1; i <= trm->numEdges; i++ ) {
		tw->edges[i].vertexNum[0] = trm->edges[i].v[0];
		tw->edges[i].vertexNum[1] = trm->edges[i].v[1];
		tw->edges[i].used = qfalse;
	}
	// polygons
	tw->numPolys = trm->numPolys;
	for ( i = 0; i < trm->numPolys; i++ ) {
		tw->polys[i].numEdges = trm->polys[i].numEdges;
		for ( j = 0; j < trm->polys[i].numEdges; j++ ) {
			tw->polys[i].edges[j] = trm->polys[i].edges[j];
		}
        VectorCopy( trm->polys[i].normal, tw->polys[i].plane );
		tw->polys[i].used = qfalse;
	}
	// is the trace model convex or not
	tw->isConvex = trm->isConvex;
}

/*
================
CM_SetupTranslationHeartPlanes
================
*/
void CM_SetupTranslationHeartPlanes( cm_traceWork_t *tw ) {
	vec3_t dir, normal1, normal2;

	// calculate trace heart planes
	VectorCopy( tw->dir, dir );
	VectorNormalize( dir );
    NormalVectors( dir, normal1, normal2 );
	VectorCopy( normal1, tw->heartPlane1 );
    PlaneFitThroughPoint( tw->heartPlane1, tw->start );
	VectorCopy( normal2, tw->heartPlane2 );
    PlaneFitThroughPoint( tw->heartPlane2, tw->start );
}

/*
================
CM_Translation
================
*/
#ifdef _DEBUG
static int entered = 0;
#endif

void CM_Translation( cm_trace_t *results, const vec3_t start, const vec3_t end,
										const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
										cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] ) {

	int i, j;
	float dist;
	qboolean model_rotated, trm_rotated;
	vec3_t dir, tmp;
	vec3_t invModelAxis[3], tmpAxis[3];
	cm_trmPolygon_t *poly;
	cm_trmEdge_t *edge;
	cm_trmVertex_t *vert;
	static cm_traceWork_t tw; // FIXME: ALIGN16

	assert( ((byte *)&start) < ((byte *)results) || ((byte *)&start) >= (((byte *)results) + sizeof( cm_trace_t )) );
	assert( ((byte *)&end) < ((byte *)results) || ((byte *)&end) >= (((byte *)results) + sizeof( cm_trace_t )) );
	assert( ((byte *)&trmAxis) < ((byte *)results) || ((byte *)&trmAxis) >= (((byte *)results) + sizeof( cm_trace_t )) );

	memset( results, 0, sizeof( *results ) );

	if ( model < 0 || model > MAX_SUBMODELS || model > cmLocal.maxModels ) {
		ii.Com_Printf("Translation: invalid model handle\n");
		return;
	}
	if ( !cmLocal.models[model] ) {
		ii.Com_Printf("Translation: invalid model\n");
		return;
	}

	// if case special position test
	if ( start[0] == end[0] && start[1] == end[1] && start[2] == end[2] ) {
		CM_ContentsTrm( results, start, trm, trmAxis, contentMask, model, modelOrigin, modelAxis );
		return;
	}

#ifdef _DEBUG
	qboolean startsolid = qfalse;
	// test whether or not stuck to begin with
	if ( cm_debugCollision->integer ) {
		if ( !entered && !cmLocal.getContacts ) {
			entered = 1;
			// if already messed up to begin with
			if ( CM_Contents( start, trm, trmAxis, -1, model, modelOrigin, modelAxis ) & contentMask ) {
				startsolid = qtrue;
			}
			entered = 0;
		}
	}
#endif

	cmLocal.checkCount++;

	tw.trace.fraction = 1.0f;
	tw.trace.c.contents = 0;
	tw.trace.c.type = CONTACT_NONE;
	tw.contents = contentMask;
	tw.isConvex = qtrue;
	tw.rotation = qfalse;
	tw.positionTest = qfalse;
	tw.quickExit = qfalse;
	tw.getContacts = cmLocal.getContacts;
	tw.contacts = cmLocal.contacts;
	tw.maxContacts = cmLocal.maxContacts;
	tw.numContacts = 0;
	tw.model = cmLocal.models[model];
    VectorSubtract( start, modelOrigin, tw.start );
    VectorSubtract( end, modelOrigin, tw.end );
    VectorSubtract( end, start, tw.dir );

	model_rotated = AxisIsRotated( modelAxis );
	if ( model_rotated ) {
        TransposeAxis( modelAxis, invModelAxis );
	}

	// if optimized point trace
	if ( !trm || ( trm->bounds[1][0] - trm->bounds[0][0] <= 0.0f &&
					trm->bounds[1][1] - trm->bounds[0][1] <= 0.0f &&
					trm->bounds[1][2] - trm->bounds[0][2] <= 0.0f ) ) {

		if ( model_rotated ) {
			// rotate trace instead of model
            VectorRotateSelf( tw.start, invModelAxis );
            VectorRotateSelf( tw.end, invModelAxis );
            VectorRotateSelf( tw.dir, invModelAxis );
		}

		// trace bounds
		for ( i = 0; i < 3; i++ ) {
			if ( tw.start[i] < tw.end[i] ) {
				tw.bounds[0][i] = tw.start[i] - CM_BOX_EPSILON;
				tw.bounds[1][i] = tw.end[i] + CM_BOX_EPSILON;
			}
			else {
				tw.bounds[0][i] = tw.end[i] - CM_BOX_EPSILON;
				tw.bounds[1][i] = tw.start[i] + CM_BOX_EPSILON;
			}
		}
		tw.extents[0] = tw.extents[1] = tw.extents[2] = CM_BOX_EPSILON;
        VectorClear( tw.size[0] );
        VectorClear( tw.size[1] );

		// setup trace heart planes
		CM_SetupTranslationHeartPlanes( &tw );
		tw.maxDistFromHeartPlane1 = CM_BOX_EPSILON;
		tw.maxDistFromHeartPlane2 = CM_BOX_EPSILON;
		// collision with single point
		tw.numVerts = 1;
		VectorCopy( tw.start, tw.vertices[0].p );
		VectorAdd( tw.vertices[0].p, tw.dir, tw.vertices[0].endp );
		PlueckerFromRay( tw.vertices[0].pl, tw.vertices[0].p, tw.dir );
		tw.numEdges = tw.numPolys = 0;
		tw.pointTrace = qtrue;
		// trace through the model
		CM_TraceThroughModel( &tw );
		// store results
		*results = tw.trace;
        VectorSubtract( end, start, tmp );
        VectorMA( start, results->fraction, tmp, results->endpos );
        AxisClear( results->endAxis );

		if ( results->fraction < 1.0f ) {
			// rotate trace plane normal if there was a collision with a rotated model
			if ( model_rotated ) {
                VectorRotateSelf( results->c.normal, modelAxis );
                VectorRotateSelf( results->c.point, modelAxis );
			}
			VectorAdd( results->c.point, modelOrigin, results->c.point );
			results->c.dist += DotProduct( modelOrigin, results->c.normal );
		}
		cmLocal.numContacts = tw.numContacts;
		return;
	}

	// the trace fraction is too inaccurate to describe translations over huge distances
	if ( VectorLengthSquared( tw.dir ) > Square( CM_MAX_TRACE_DIST ) ) {
		results->fraction = 0.0f;
		VectorCopy( start, results->endpos );
        AxisCopy( trmAxis, results->endAxis );
		VectorClear( results->c.normal );
		results->c.material = 0;
		VectorCopy( start, results->c.point );
		cmi.R_DebugArrow( colorRed, start, end, 1, 0 );
		ii.Com_Printf( "Translation: huge translation\n" );
		return;
	}

	tw.pointTrace = qfalse;
    ClearBounds( tw.size[0], tw.size[1] );

	// setup trm structure
	CM_SetupTrm( &tw, trm );

	trm_rotated = AxisIsRotated( trmAxis );

	// calculate vertex positions
	if ( trm_rotated ) {
		for ( i = 0; i < tw.numVerts; i++ ) {
			// rotate trm around the start position
            VectorRotateSelf( tw.vertices[i].p, trmAxis );
		}
	}
	for ( i = 0; i < tw.numVerts; i++ ) {
		// set trm at start position
        VectorAdd( tw.vertices[i].p, tw.start, tw.vertices[i].p );
	}
	if ( model_rotated ) {
		for ( i = 0; i < tw.numVerts; i++ ) {
			// rotate trm around model instead of rotating the model
            VectorRotateSelf( tw.vertices[i].p, invModelAxis );
		}
	}

	// add offset to start point
	if ( trm_rotated ) {
        MatrixRotateVector( trm->offset, trmAxis, dir );
        VectorAdd( tw.start, dir, tw.start );
        VectorAdd( tw.end, dir, tw.end );
	} else {
        VectorAdd( tw.start, trm->offset, tw.start );
        VectorAdd( tw.end, trm->offset, tw.end );
	}
	if ( model_rotated ) {
		// rotate trace instead of model
        VectorRotateSelf( tw.start, invModelAxis );
        VectorRotateSelf( tw.end, invModelAxis );
        VectorRotateSelf( tw.dir, invModelAxis );
	}

	// rotate trm polygon planes
	if ( trm_rotated & model_rotated ) {
        MatrixMultiply( trmAxis, invModelAxis, tmpAxis );
		for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
            PlaneRotateSelf( poly->plane, tmpAxis );
		}
	} else if ( trm_rotated ) {
		for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
            PlaneRotateSelf( poly->plane, trmAxis );
		}
	} else if ( model_rotated ) {
		for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
            PlaneRotateSelf( poly->plane, invModelAxis );
		}
	}

	// setup trm polygons
	for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
		// if the trm poly plane is facing in the movement direction
		dist = DotProduct( poly->plane, tw.dir );
		if ( dist > 0.0f || ( !trm->isConvex && dist == 0.0f ) ) {
			// this trm poly and it's edges and vertices need to be used for collision
			poly->used = qtrue;
			for ( j = 0; j < poly->numEdges; j++ ) {
				edge = &tw.edges[abs( poly->edges[j] )];
				edge->used = qtrue;
				tw.vertices[edge->vertexNum[0]].used = qtrue;
				tw.vertices[edge->vertexNum[1]].used = qtrue;
			}
		}
	}

	// setup trm vertices
	for ( vert = tw.vertices, i = 0; i < tw.numVerts; i++, vert++ ) {
		if ( !vert->used ) {
			continue;
		}
		// get axial trm size after rotations
		VectorSubtract( vert->p, tw.start, tmp );
        AddPointToBounds( tmp, tw.size[0], tw.size[1] );
		// calculate the end position of each vertex for a full trace
        VectorAdd( vert->p, tw.dir, vert->endp );
		// pluecker coordinate for vertex movement line
		PlueckerFromRay( vert->pl, vert->p, tw.dir );
	}

	// setup trm edges
	for ( edge = tw.edges + 1, i = 1; i <= tw.numEdges; i++, edge++ ) {
		if ( !edge->used ) {
			continue;
		}
		// edge start, end and pluecker coordinate
		VectorCopy( tw.vertices[edge->vertexNum[0]].p, edge->start );
		VectorCopy( tw.vertices[edge->vertexNum[1]].p, edge->end );
		PlueckerFromLine( edge->pl, edge->start, edge->end );
		// calculate normal of plane through movement plane created by the edge
        VectorSubtract( edge->start, edge->end, dir );
		edge->cross[0] = dir[0] * tw.dir[1] - dir[1] * tw.dir[0];
		edge->cross[1] = dir[0] * tw.dir[2] - dir[2] * tw.dir[0];
		edge->cross[2] = dir[1] * tw.dir[2] - dir[2] * tw.dir[1];
		// bit for vertex sidedness bit cache
		edge->bitNum = i;
	}

	// set trm plane distances
	for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
		if ( poly->used ) {
            PlaneFitThroughPoint( poly->plane, tw.edges[abs(poly->edges[0])].start );
		}
	}

	// bounds for full trace, a little bit larger for epsilons
	for ( i = 0; i < 3; i++ ) {
		if ( tw.start[i] < tw.end[i] ) {
			tw.bounds[0][i] = tw.start[i] + tw.size[0][i] - CM_BOX_EPSILON;
			tw.bounds[1][i] = tw.end[i] + tw.size[1][i] + CM_BOX_EPSILON;
		} else {
			tw.bounds[0][i] = tw.end[i] + tw.size[0][i] - CM_BOX_EPSILON;
			tw.bounds[1][i] = tw.start[i] + tw.size[1][i] + CM_BOX_EPSILON;
		}
		if ( fabs( tw.size[0][i] ) > fabs( tw.size[1][i] ) ) {
			tw.extents[i] = fabs( tw.size[0][i] ) + CM_BOX_EPSILON;
		} else {
			tw.extents[i] = fabs( tw.size[1][i] ) + CM_BOX_EPSILON;
		}
	}

	// setup trace heart planes
	CM_SetupTranslationHeartPlanes( &tw );
	tw.maxDistFromHeartPlane1 = 0;
	tw.maxDistFromHeartPlane2 = 0;
	// calculate maximum trm vertex distance from both heart planes
	for ( vert = tw.vertices, i = 0; i < tw.numVerts; i++, vert++ ) {
		if ( !vert->used ) {
			continue;
		}
		dist = fabs( PlaneDistance( tw.heartPlane1, vert->p ) );
		if ( dist > tw.maxDistFromHeartPlane1 ) {
			tw.maxDistFromHeartPlane1 = dist;
		}
		dist = fabs( PlaneDistance( tw.heartPlane2, vert->p ) );
		if ( dist > tw.maxDistFromHeartPlane2 ) {
			tw.maxDistFromHeartPlane2 = dist;
		}
	}
	// for epsilons
	tw.maxDistFromHeartPlane1 += CM_BOX_EPSILON;
	tw.maxDistFromHeartPlane2 += CM_BOX_EPSILON;

	// trace through the model
	CM_TraceThroughModel( &tw );

	// if we're getting contacts
	if ( tw.getContacts ) {
		// move all contacts to world space
		if ( model_rotated ) {
			for ( i = 0; i < tw.numContacts; i++ ) {
                VectorRotateSelf( tw.contacts[i].normal, modelAxis );
                VectorRotateSelf( tw.contacts[i].point, modelAxis );
			}
		}
		if ( modelOrigin != vec3_origin ) {
			for ( i = 0; i < tw.numContacts; i++ ) {
				VectorAdd( tw.contacts[i].point, modelOrigin, tw.contacts[i].point );
				tw.contacts[i].dist += DotProduct( modelOrigin, tw.contacts[i].normal );
			}
		}
		cmLocal.numContacts = tw.numContacts;
	} else {
		// store results
		*results = tw.trace;
        VectorSubtract( end, start, tmp );
        VectorMA( start, results->fraction, tmp, results->endpos );
        AxisCopy( trmAxis, results->endAxis );

		if ( results->fraction < 1.0f ) {
			// if the fraction is tiny the actual movement could end up zero
			if ( results->fraction > 0.0f && VectorCompare( results->endpos, start ) ) {
				results->fraction = 0.0f;
			}
			// rotate trace plane normal if there was a collision with a rotated model
			if ( model_rotated ) {
                VectorRotateSelf( results->c.normal, modelAxis );
                VectorRotateSelf( results->c.point, modelAxis );
			}
            VectorAdd( results->c.point, modelOrigin, results->c.point );
			results->c.dist += DotProduct( modelOrigin, results->c.normal );
		}
	}

#ifdef _DEBUG
	// test for missed collisions
	if ( cm_debugCollision.integerValue ) {
		if ( !entered && !cmLocal.getContacts ) {
			entered = 1;
			// if the trm is stuck in the model
			if ( CM_Contents( results->endpos, trm, trmAxis, -1, model, modelOrigin, modelAxis ) & contentMask ) {
				cm_trace_t tr;

				// test where the trm is stuck in the model
				CM_Contents( results->endpos, trm, trmAxis, -1, model, modelOrigin, modelAxis );
				// re-run collision detection to find out where it failed
				CM_Translation( &tr, start, end, trm, trmAxis, contentMask, model, modelOrigin, modelAxis );
			}
			entered = 0;
		}
	}
#endif
}
