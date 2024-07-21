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

Collision detection for rotational motion

===============================================================================
*/

// epsilon for round-off errors in epsilon calculations
#define CM_PL_RANGE_EPSILON			1e-4f
// if the collision point is this close to the rotation axis it is not considered a collision
#define ROTATION_AXIS_EPSILON		(CM_CLIP_EPSILON*0.25f)


/*
================
CM_RotatePoint

  rotates a point about an arbitrary axis using the tangent of half the rotation angle
================
*/
void CM_RotatePoint( vec3_t point, const vec3_t origin, const vec3_t axis, const float tanHalfAngle ) {
	double d, t, s, c;
	vec3_t proj, v1, v2;
    int i;

    VectorSubtract( point, origin, point );
    VectorScale( axis, DotProduct( point, axis ), proj );

    VectorSubtract( point, proj, v1 );
    CrossProduct( axis, v1, v2 );

	// r = tan( a / 2 );
	// sin(a) = 2*r/(1+r*r);
	// cos(a) = (1-r*r)/(1+r*r);
	t = tanHalfAngle * tanHalfAngle;
	d = 1.0f / ( 1.0f + t );
	s = 2.0f * tanHalfAngle * d;
	c = ( 1.0f - t ) * d;

    for ( i = 0; i < 3; i++ ) {
	    point[i] = v1[i] * c - v2[i] * s + proj[i] + origin[i];
    }
}

/*
================
CM_RotateEdge

  rotates an edge about an arbitrary axis using the tangent of half the rotation angle
================
*/
void CM_RotateEdge( vec3_t start, vec3_t end, const vec3_t origin, const vec3_t axis, const float tanHalfAngle ) {
	double d, t, s, c;
	vec3_t proj, v1, v2;
    int i;

	// r = tan( a / 2 );
	// sin(a) = 2*r/(1+r*r);
	// cos(a) = (1-r*r)/(1+r*r);
	t = tanHalfAngle * tanHalfAngle;
	d = 1.0f / ( 1.0f + t );
	s = 2.0f * tanHalfAngle * d;
	c = ( 1.0f - t ) * d;

    VectorSubtract( start, origin, start );
    VectorScale( axis, DotProduct( start, axis ), proj );
    VectorSubtract( start, proj, v1 );
    CrossProduct( axis, v1, v2 );
    for ( i = 0; i < 3; i++ ) {
	    start[i] = v1[i] * c - v2[i] * s + proj[i] + origin[i];
    }

    VectorSubtract( end, origin, end );
    VectorScale( axis, DotProduct( end, axis ), proj );
	VectorSubtract( end, proj, v1 );
    CrossProduct( axis, v1, v2 );
    for ( i = 0; i < 3; i++ ) {
	    end[i] = v1[i] * c - v2[i] * s + proj[i] + origin[i];
    }
}

/*
================
CM_CollisionBetweenEdgeBounds

  verifies if the collision of two edges occurs between the edge bounds
  also calculates the collision point and collision plane normal if the collision occurs between the bounds
================
*/
int CM_CollisionBetweenEdgeBounds( cm_traceWork_t *tw, const vec3_t va, const vec3_t vb,
												   const vec3_t vc, const vec3_t vd, float tanHalfAngle,
												   vec3_t collisionPoint, vec3_t collisionNormal ) {
	float d1, d2, d;
	vec3_t at, bt, dir, dir1, dir2, tmp;
	vec6_t pl1, pl2;
    int i;

    VectorCopy( va, at );
    VectorCopy( vb, bt );
	if ( tanHalfAngle != 0.0f ) {
		CM_RotateEdge( at, bt, tw->origin, tw->axis, tanHalfAngle );
	}

    VectorSubtract( at, tw->origin, tmp );
    CrossProduct( tmp, tw->axis, dir1 );
    VectorSubtract( bt, tw->origin, tmp );
    CrossProduct( tmp, tw->axis, dir2 );
	if ( DotProduct( dir1, dir1 ) > DotProduct( dir2, dir2 ) ) {
        VectorCopy( dir1, dir );
	}
	else {
        VectorCopy( dir2, dir );
	}
	if ( tw->angle < 0.0f ) {
        VectorNegate( dir, dir );
	}

    PlueckerFromLine( pl1, at, bt );
    PlueckerFromRay( pl2, vc, dir );
	d1 = PlueckerPermutedInnerProduct( pl1, pl2 );
	PlueckerFromRay( pl2, vd, dir );
	d2 = PlueckerPermutedInnerProduct( pl1, pl2 );
	if ( ( d1 > 0.0f && d2 > 0.0f ) || ( d1 < 0.0f && d2 < 0.0f ) ) {
		return qfalse;
	}

	PlueckerFromLine( pl1, vc, vd );
	PlueckerFromRay( pl2, at, dir );
	d1 = PlueckerPermutedInnerProduct( pl1, pl2 );
	PlueckerFromRay( pl2, bt, dir );
	d2 = PlueckerPermutedInnerProduct( pl1, pl2 );
	if ( ( d1 > 0.0f && d2 > 0.0f ) || ( d1 < 0.0f && d2 < 0.0f ) ) {
		return qfalse;
	}

	// collision point on the edge at-bt
    VectorSubtract( vd, vc, tmp );
    CrossProduct( tmp, dir, dir1 );
	d = DotProduct( dir1, vc );
	d1 = DotProduct( dir1, at ) - d;
	d2 = DotProduct( dir1, bt ) - d;
	if ( d1 == d2 ) {
		return qfalse;
	}
    for ( i = 0; i < 3; i++ ) {
	    collisionPoint[i] = at[i] + ( d1 / (d1 - d2) ) * ( bt[i] - at[i] );
    }

	// normal is cross product of the rotated edge va-vb and the edge vc-vd
    VectorSubtract( bt, at, dir1 );
    VectorSubtract( vd, vc, dir2 );
    CrossProduct( dir1, dir2, collisionNormal );

	return qtrue;
}

/*
================
CM_RotateEdgeThroughEdge

  calculates the tangent of half the rotation angle at which the edges collide
================
*/
int CM_RotateEdgeThroughEdge( cm_traceWork_t *tw, const vec6_t pl1,
												const vec3_t vc, const vec3_t vd,
												const float minTan, float *tanHalfAngle ) {
	double v0, v1, v2, a, b, c, d, sqrtd, q, frac1, frac2;
	vec3_t ct, dt;
	vec6_t pl2;
    vec3_t tmp;

	/*

	a = start of line being rotated
	b = end of line being rotated
	pl1 = pluecker coordinate for line (a - b)
	pl2 = pluecker coordinate for edge we might collide with (c - d)
	t = rotation angle around the z-axis
	solve pluecker inner product for t of rotating line a-b and line l2

	// start point of rotated line during rotation
	an[0] = a[0] * cos(t) + a[1] * sin(t)
	an[1] = a[0] * -sin(t) + a[1] * cos(t)
	an[2] = a[2];
	// end point of rotated line during rotation
	bn[0] = b[0] * cos(t) + b[1] * sin(t)
	bn[1] = b[0] * -sin(t) + b[1] * cos(t)
	bn[2] = b[2];

	pl1[0] = a[0] * b[1] - b[0] * a[1];
	pl1[1] = a[0] * b[2] - b[0] * a[2];
	pl1[2] = a[0] - b[0];
	pl1[3] = a[1] * b[2] - b[1] * a[2];
	pl1[4] = a[2] - b[2];
	pl1[5] = b[1] - a[1];

	v[0] = (a[0] * cos(t) + a[1] * sin(t)) * (b[0] * -sin(t) + b[1] * cos(t)) - (b[0] * cos(t) + b[1] * sin(t)) * (a[0] * -sin(t) + a[1] * cos(t));
	v[1] = (a[0] * cos(t) + a[1] * sin(t)) * b[2] - (b[0] * cos(t) + b[1] * sin(t)) * a[2];
	v[2] = (a[0] * cos(t) + a[1] * sin(t)) - (b[0] * cos(t) + b[1] * sin(t));
	v[3] = (a[0] * -sin(t) + a[1] * cos(t)) * b[2] - (b[0] * -sin(t) + b[1] * cos(t)) * a[2];
	v[4] = a[2] - b[2];
	v[5] = (b[0] * -sin(t) + b[1] * cos(t)) - (a[0] * -sin(t) + a[1] * cos(t));

	pl2[0] * v[4] + pl2[1] * v[5] + pl2[2] * v[3] + pl2[4] * v[0] + pl2[5] * v[1] + pl2[3] * v[2] = 0;

	v[0] = (a[0] * cos(t) + a[1] * sin(t)) * (b[0] * -sin(t) + b[1] * cos(t)) - (b[0] * cos(t) + b[1] * sin(t)) * (a[0] * -sin(t) + a[1] * cos(t));
	v[0] = (a[1] * b[1] - a[0] * b[0]) * cos(t) * sin(t) + (a[0] * b[1] + a[1] * b[0] * cos(t)^2) - (a[1] * b[0]) - ((b[1] * a[1] - b[0] * a[0]) * cos(t) * sin(t) + (b[0] * a[1] + b[1] * a[0]) * cos(t)^2 - (b[1] * a[0]))
	v[0] = - (a[1] * b[0]) - ( - (b[1] * a[0]))
	v[0] = (b[1] * a[0]) - (a[1] * b[0])

	v[0] = (a[0]*b[1]) - (a[1]*b[0]);
	v[1] = (a[0]*b[2] - b[0]*a[2]) * cos(t) + (a[1]*b[2] - b[1]*a[2]) * sin(t);
	v[2] = (a[0]-b[0]) * cos(t) + (a[1]-b[1]) * sin(t);
	v[3] = (b[0]*a[2] - a[0]*b[2]) * sin(t) + (a[1]*b[2] - b[1]*a[2]) * cos(t);
	v[4] = a[2] - b[2];
	v[5] = (a[0]-b[0]) * sin(t) + (b[1]-a[1]) * cos(t);

	v[0] = (a[0]*b[1]) - (a[1]*b[0]);
	v[1] = (a[0]*b[2] - b[0]*a[2]) * cos(t) + (a[1]*b[2] - b[1]*a[2]) * sin(t);
	v[2] = (a[0]-b[0]) * cos(t) - (b[1]-a[1]) * sin(t);
	v[3] = (a[0]*b[2] - b[0]*a[2]) * -sin(t) + (a[1]*b[2] - b[1]*a[2]) * cos(t);
	v[4] = a[2] - b[2];
	v[5] = (a[0]-b[0]) * sin(t) + (b[1]-a[1]) * cos(t);

	v[0] = pl1[0];
	v[1] = pl1[1] * cos(t) + pl1[3] * sin(t);
	v[2] = pl1[2] * cos(t) - pl1[5] * sin(t);
	v[3] = pl1[3] * cos(t) - pl1[1] * sin(t);
	v[4] = pl1[4];
	v[5] = pl1[5] * cos(t) + pl1[2] * sin(t);

	pl2[0] * v[4] + pl2[1] * v[5] + pl2[2] * v[3] + pl2[4] * v[0] + pl2[5] * v[1] + pl2[3] * v[2] = 0;

	0 =	pl2[0] * pl1[4] +
		pl2[1] * (pl1[5] * cos(t) + pl1[2] * sin(t)) +
		pl2[2] * (pl1[3] * cos(t) - pl1[1] * sin(t)) +
		pl2[4] * pl1[0] +
		pl2[5] * (pl1[1] * cos(t) + pl1[3] * sin(t)) +
		pl2[3] * (pl1[2] * cos(t) - pl1[5] * sin(t));

	v2 * cos(t) + v1 * sin(t) + v0 = 0;

	// rotation about the z-axis
	v0 = pl2[0] * pl1[4] + pl2[4] * pl1[0];
	v1 = pl2[1] * pl1[2] - pl2[2] * pl1[1] + pl2[5] * pl1[3] - pl2[3] * pl1[5];
	v2 = pl2[1] * pl1[5] + pl2[2] * pl1[3] + pl2[5] * pl1[1] + pl2[3] * pl1[2];

	// rotation about the x-axis
	//v0 = pl2[3] * pl1[2] + pl2[2] * pl1[3];
	//v1 = -pl2[5] * pl1[0] + pl2[4] * pl1[1] - pl2[1] * pl1[4] + pl2[0] * pl1[5];
	//v2 = pl2[4] * pl1[0] + pl2[5] * pl1[1] + pl2[0] * pl1[4] + pl2[1] * pl1[5];

	r = tan(t / 2);
	sin(t) = 2*r/(1+r*r);
	cos(t) = (1-r*r)/(1+r*r);

	v1 * 2 * r / (1 + r*r) + v2 * (1 - r*r) / (1 + r*r) + v0 = 0
	(v1 * 2 * r + v2 * (1 - r*r)) / (1 + r*r) = -v0
	(v1 * 2 * r + v2 - v2 * r*r) / (1 + r*r) = -v0
	v1 * 2 * r + v2 - v2 * r*r = -v0 * (1 + r*r)
	v1 * 2 * r + v2 - v2 * r*r = -v0 + -v0 * r*r
	(v0 - v2) * r * r + (2 * v1) * r + (v0 + v2) = 0;

	MrE gives Pluecker a banana.. good monkey

	*/

	*tanHalfAngle = tw->maxTan;

	// transform rotation axis to z-axis
    VectorSubtract( vc, tw->origin, tmp );
    MatrixRotateVector( tmp, tw->matrix, ct );
    VectorSubtract( vd, tw->origin, tmp );
    MatrixRotateVector( tmp, tw->matrix, dt );

	PlueckerFromLine( pl2, ct, dt );

	v0 = pl2[0] * pl1[4] + pl2[4] * pl1[0];
	v1 = pl2[1] * pl1[2] - pl2[2] * pl1[1] + pl2[5] * pl1[3] - pl2[3] * pl1[5];
	v2 = pl2[1] * pl1[5] + pl2[2] * pl1[3] + pl2[5] * pl1[1] + pl2[3] * pl1[2];

	a = v0 - v2;
	b = v1;
	c = v0 + v2;
	if ( a == 0.0f ) {
		if ( b == 0.0f ) {
			return qfalse;
		}
		frac1 = -c / ( 2.0f * b );
		frac2 = 1e10;	// = tan( idMath::HALF_PI )
	}
	else {
		d = b * b - c * a;
		if ( d <= 0.0f ) {
			return qfalse;
		}
		sqrtd = sqrt( d );
		if ( b > 0.0f ) {
			q = - b + sqrtd;
		}
		else {
			q = - b - sqrtd;
		}
		frac1 = q / a;
		frac2 = c / q;
	}

	if ( tw->angle < 0.0f ) {
		frac1 = -frac1;
		frac2 = -frac2;
	}

	// get smallest tangent for which a collision occurs
	if ( frac1 >= minTan && frac1 < *tanHalfAngle ) {
		*tanHalfAngle = frac1;
	}
	if ( frac2 >= minTan && frac2 < *tanHalfAngle ) {
		*tanHalfAngle = frac2;
	}

	if ( tw->angle < 0.0f ) {
		*tanHalfAngle = -*tanHalfAngle;
	}

	return qtrue;
}

/*
================
CM_EdgeFurthestFromEdge

  calculates the direction of motion at the initial position, where dir < 0 means the edges move towards each other
  if the edges move away from each other the tangent of half the rotation angle at which
  the edges are furthest apart is also calculated
================
*/
int CM_EdgeFurthestFromEdge( cm_traceWork_t *tw, const vec6_t pl1,
												const vec3_t vc, const vec3_t vd,
												float *tanHalfAngle, float *dir ) {
	double v0, v1, v2, a, b, c, d, sqrtd, q, frac1, frac2;
	vec3_t ct, dt;
	vec6_t pl2;
    vec3_t tmp;

	/*

	v2 * cos(t) + v1 * sin(t) + v0 = 0;

	// rotation about the z-axis
	v0 = pl2[0] * pl1[4] + pl2[4] * pl1[0];
	v1 = pl2[1] * pl1[2] - pl2[2] * pl1[1] + pl2[5] * pl1[3] - pl2[3] * pl1[5];
	v2 = pl2[1] * pl1[5] + pl2[2] * pl1[3] + pl2[5] * pl1[1] + pl2[3] * pl1[2];

	derivative:
	v1 * cos(t) - v2 * sin(t) = 0;

	r = tan(t / 2);
	sin(t) = 2*r/(1+r*r);
	cos(t) = (1-r*r)/(1+r*r);

	-v2 * 2 * r / (1 + r*r) + v1 * (1 - r*r)/(1+r*r);
	-v2 * 2 * r + v1 * (1 - r*r) / (1 + r*r) = 0;
	-v2 * 2 * r + v1 * (1 - r*r) = 0;
	(-v1) * r * r + (-2 * v2) * r + (v1) = 0;

	*/

	*tanHalfAngle = 0.0f;

	// transform rotation axis to z-axis
    VectorSubtract( vc, tw->origin, tmp );
    MatrixRotateVector( tmp, tw->matrix, ct );
    VectorSubtract( vd, tw->origin, tmp );
    MatrixRotateVector( tmp, tw->matrix, dt );

	PlueckerFromLine( pl2, ct, dt );

	v0 = pl2[0] * pl1[4] + pl2[4] * pl1[0];
	v1 = pl2[1] * pl1[2] - pl2[2] * pl1[1] + pl2[5] * pl1[3] - pl2[3] * pl1[5];
	v2 = pl2[1] * pl1[5] + pl2[2] * pl1[3] + pl2[5] * pl1[1] + pl2[3] * pl1[2];

	// get the direction of motion at the initial position
	c = v0 + v2;
	if ( tw->angle > 0.0f ) {
		if ( c > 0.0f ) {
			*dir = v1;
		}
		else {
			*dir = -v1;
		}
	}
	else {
		if ( c > 0.0f ) {
			*dir = -v1;
		}
		else {
			*dir = v1;
		}
	}
	// negative direction means the edges move towards each other at the initial position
	if ( *dir <= 0.0f ) {
		return qtrue;
	}

	a = -v1;
	b = -v2;
	c = v1;
	if ( a == 0.0f ) {
		if ( b == 0.0f ) {
			return qfalse;
		}
		frac1 = -c / ( 2.0f * b );
		frac2 = 1e10;	// = tan( idMath::HALF_PI )
	}
	else {
		d = b * b - c * a;
		if ( d <= 0.0f ) {
			return qfalse;
		}
		sqrtd = sqrt( d );
		if ( b > 0.0f ) {
			q = - b + sqrtd;
		}
		else {
			q = - b - sqrtd;
		}
		frac1 = q / a;
		frac2 = c / q;
	}

	if ( tw->angle < 0.0f ) {
		frac1 = -frac1;
		frac2 = -frac2;
	}

	if ( frac1 < 0.0f && frac2 < 0.0f ) {
		return qfalse;
	}

	if ( frac1 > frac2 ) {
		*tanHalfAngle = frac1;
	}
	else {
		*tanHalfAngle = frac2;
	}

	if ( tw->angle < 0.0f ) {
		*tanHalfAngle = -*tanHalfAngle;
	}

	return qtrue;
}

/*
================
CM_RotateTrmEdgeThroughPolygon
================
*/
void CM_RotateTrmEdgeThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmEdge_t *trmEdge ) {
	int i, j, edgeNum;
	float f1, f2, startTan, dir, tanHalfAngle;
	cm_edge_t *edge;
	cm_vertex_t *v1, *v2;
	vec3_t collisionPoint, collisionNormal, origin, epsDir;
	vec6_t epsPl;
	vec3_t bounds[2];
    vec3_t tmp1, tmp2;

	// if the trm is convex and the rotation axis intersects the trm
	if ( tw->isConvex && tw->axisIntersectsTrm ) {
		// if both points are behind the polygon the edge cannot collide within a 180 degrees rotation
		if ( tw->vertices[trmEdge->vertexNum[0]].polygonSide & tw->vertices[trmEdge->vertexNum[1]].polygonSide ) {
			return;
		}
	}

	// if the trace model edge rotation bounds do not intersect the polygon bounds
	if ( !BoundsIntersect( trmEdge->rotationBounds[0], trmEdge->rotationBounds[1], poly->bounds[0], poly->bounds[1] ) ) {
		return;
	}

	// edge rotation bounds should cross polygon plane
	if ( BoxOnPlaneSideSlow( trmEdge->rotationBounds[0], trmEdge->rotationBounds[1], poly->plane, ON_EPSILON ) != PLANESIDE_CROSS ) {
		return;
	}

	// check edges for a collision
	for ( i = 0; i < poly->numEdges; i++ ) {
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

		v1 = tw->model->vertices + edge->vertexNum[INTSIGNBITSET(edgeNum)];
		v2 = tw->model->vertices + edge->vertexNum[INTSIGNBITNOTSET(edgeNum)];

		// edge bounds
		for ( j = 0; j < 3; j++ ) {
			if ( v1->p[j] > v2->p[j] ) {
				bounds[0][j] = v2->p[j];
				bounds[1][j] = v1->p[j];
			}
			else {
				bounds[0][j] = v1->p[j];
				bounds[1][j] = v2->p[j];
			}
		}

		// if the trace model edge rotation bounds do not intersect the polygon edge bounds
		if ( !BoundsIntersect( trmEdge->rotationBounds[0], trmEdge->rotationBounds[1], bounds[0], bounds[1] ) ) {
			continue;
		}

		f1 = PlueckerPermutedInnerProduct( trmEdge->pl, tw->polygonEdgePlueckerCache[i] );

		// pluecker coordinate for epsilon expanded edge
        VectorScale( edge->normal, (CM_CLIP_EPSILON+CM_PL_RANGE_EPSILON), epsDir );
        VectorAdd( tw->model->vertices[edge->vertexNum[0]].p, epsDir, tmp1 );
        VectorAdd( tw->model->vertices[edge->vertexNum[1]].p, epsDir, tmp2 );
		PlueckerFromLine( epsPl, tmp1, tmp2 );

		f2 = PlueckerPermutedInnerProduct( trmEdge->pl, epsPl );

		// if the rotating edge is inbetween the polygon edge and the epsilon expanded edge
		if ( ( f1 < 0.0f && f2 > 0.0f ) || ( f1 > 0.0f && f2 < 0.0f ) ) {

			if ( !CM_EdgeFurthestFromEdge( tw, trmEdge->plzaxis, v1->p, v2->p, &startTan, &dir ) ) {
				continue;
			}

			if ( dir <= 0.0f ) {
				// moving towards the polygon edge so stop immediately
				tanHalfAngle = 0.0f;
			}
			else if ( fabs( startTan ) >= tw->maxTan ) {
				// never going to get beyond the start tangent during the current rotation
				continue;
			}
			else {
				// collide with the epsilon expanded edge
                VectorAdd( v1->p, epsDir, tmp1 );
                VectorAdd( v2->p, epsDir, tmp2 );
				if ( !CM_RotateEdgeThroughEdge(tw, trmEdge->plzaxis, tmp1, tmp2, fabs( startTan ), &tanHalfAngle ) ) {
					tanHalfAngle = startTan;
				}
			}
		}
		else {
			// collide with the epsilon expanded edge
            VectorScale( edge->normal, CM_CLIP_EPSILON, epsDir );
            VectorAdd( v1->p, epsDir, tmp1 );
            VectorAdd( v2->p, epsDir, tmp2 );
			if ( !CM_RotateEdgeThroughEdge(tw, trmEdge->plzaxis, tmp1, tmp2, 0.0f, &tanHalfAngle ) ) {
				continue;
			}
		}

		if ( fabs( tanHalfAngle ) >= tw->maxTan ) {
			continue;
		}

		// check if the collision is between the edge bounds
		if ( !CM_CollisionBetweenEdgeBounds( tw, trmEdge->start, trmEdge->end, v1->p, v2->p,
												tanHalfAngle, collisionPoint, collisionNormal ) ) {
			continue;
		}

		// allow rotation if the rotation axis goes through the collisionPoint
        VectorSubtract( collisionPoint, tw->origin, tmp1 );
        VectorMA( tw->origin, DotProduct( tw->axis, tmp1 ), tw->axis, origin );
        VectorSubtract( collisionPoint, origin, tmp1 );
		if ( VectorLengthSquared( tmp1 ) < ROTATION_AXIS_EPSILON * ROTATION_AXIS_EPSILON ) {
			continue;
		}

		// fill in trace structure
		tw->maxTan = fabs( tanHalfAngle );
		VectorCopy( collisionNormal, tw->trace.c.normal );
		VectorNormalize( tw->trace.c.normal );
		tw->trace.c.dist = DotProduct( tw->trace.c.normal, v1->p );
		// make sure the collision plane faces the trace model
		if ( DotProduct( tw->trace.c.normal, trmEdge->start ) - tw->trace.c.dist < 0 ) {
			VectorNegate( tw->trace.c.normal, tw->trace.c.normal );
			tw->trace.c.dist = -tw->trace.c.dist;
		}
		tw->trace.c.contents = poly->contents;
		tw->trace.c.material = poly->material;
		tw->trace.c.type = CONTACT_EDGE;
		tw->trace.c.modelFeature = edgeNum;
		tw->trace.c.trmFeature = trmEdge - tw->edges;
        VectorCopy( collisionPoint, tw->trace.c.point );
		// if no collision can be closer
		if ( tw->maxTan == 0.0f ) {
			break;
		}
	}
}

/*
================
CM_RotatePointThroughPlane

  calculates the tangent of half the rotation angle at which the point collides with the plane
================
*/
int CM_RotatePointThroughPlane( const cm_traceWork_t *tw, const vec3_t point, const plane_t plane,
													const float angle, const float minTan, float *tanHalfAngle ) {
	double v0, v1, v2, a, b, c, d, sqrtd, q, frac1, frac2;
	vec3_t p, normal, tmp;

	/*

	p[0] = point[0] * cos(t) + point[1] * sin(t)
	p[1] = point[0] * -sin(t) + point[1] * cos(t)
	p[2] = point[2];

	normal[0] * (p[0] * cos(t) + p[1] * sin(t)) +
		normal[1] * (p[0] * -sin(t) + p[1] * cos(t)) +
			normal[2] * p[2] + dist = 0

	normal[0] * p[0] * cos(t) + normal[0] * p[1] * sin(t) +
		-normal[1] * p[0] * sin(t) + normal[1] * p[1] * cos(t) +
			normal[2] * p[2] + dist = 0

	v2 * cos(t) + v1 * sin(t) + v0

	// rotation about the z-axis
	v0 = normal[2] * p[2] + dist
	v1 = normal[0] * p[1] - normal[1] * p[0]
	v2 = normal[0] * p[0] + normal[1] * p[1]

	r = tan(t / 2);
	sin(t) = 2*r/(1+r*r);
	cos(t) = (1-r*r)/(1+r*r);

	v1 * 2 * r / (1 + r*r) + v2 * (1 - r*r) / (1 + r*r) + v0 = 0
	(v1 * 2 * r + v2 * (1 - r*r)) / (1 + r*r) = -v0
	(v1 * 2 * r + v2 - v2 * r*r) / (1 + r*r) = -v0
	v1 * 2 * r + v2 - v2 * r*r = -v0 * (1 + r*r)
	v1 * 2 * r + v2 - v2 * r*r = -v0 + -v0 * r*r
	(v0 - v2) * r * r + (2 * v1) * r + (v0 + v2) = 0;

	*/

	*tanHalfAngle = tw->maxTan;

	// transform rotation axis to z-axis
    VectorSubtract( point, tw->origin, tmp );
    MatrixRotateVector( tmp, tw->matrix, p );
	d = plane[3] + DotProduct( plane, tw->origin );
    MatrixRotateVector( plane, tw->matrix, normal );

	v0 = normal[2] * p[2] + d;
	v1 = normal[0] * p[1] - normal[1] * p[0];
	v2 = normal[0] * p[0] + normal[1] * p[1];

	a = v0 - v2;
	b = v1;
	c = v0 + v2;
	if ( a == 0.0f ) {
		if ( b == 0.0f ) {
			return qfalse;
		}
		frac1 = -c / ( 2.0f * b );
		frac2 = 1e10;	// = tan( idMath::HALF_PI )
	}
	else {
		d = b * b - c * a;
		if ( d <= 0.0f ) {
			return qfalse;
		}
		sqrtd = sqrt( d );
		if ( b > 0.0f ) {
			q = - b + sqrtd;
		}
		else {
			q = - b - sqrtd;
		}
		frac1 = q / a;
		frac2 = c / q;
	}

	if ( angle < 0.0f ) {
		frac1 = -frac1;
		frac2 = -frac2;
	}

	// get smallest tangent for which a collision occurs
	if ( frac1 >= minTan && frac1 < *tanHalfAngle ) {
		*tanHalfAngle = frac1;
	}
	if ( frac2 >= minTan && frac2 < *tanHalfAngle ) {
		*tanHalfAngle = frac2;
	}

	if ( angle < 0.0f ) {
		*tanHalfAngle = -*tanHalfAngle;
	}

	return qtrue;
}

/*
================
CM_PointFurthestFromPlane

  calculates the direction of motion at the initial position, where dir < 0 means the point moves towards the plane
  if the point moves away from the plane the tangent of half the rotation angle at which
  the point is furthest away from the plane is also calculated
================
*/
int CM_PointFurthestFromPlane( const cm_traceWork_t *tw, const vec3_t point, const plane_t plane,
													const float angle, float *tanHalfAngle, float *dir ) {

	double v1, v2, a, b, c, d, sqrtd, q, frac1, frac2;
	vec3_t p, normal, tmp;

	/*

	v2 * cos(t) + v1 * sin(t) + v0 = 0;

	// rotation about the z-axis
	v0 = normal[2] * p[2] + dist
	v1 = normal[0] * p[1] - normal[1] * p[0]
	v2 = normal[0] * p[0] + normal[1] * p[1]

	derivative:
	v1 * cos(t) - v2 * sin(t) = 0;

	r = tan(t / 2);
	sin(t) = 2*r/(1+r*r);
	cos(t) = (1-r*r)/(1+r*r);

	-v2 * 2 * r / (1 + r*r) + v1 * (1 - r*r)/(1+r*r);
	-v2 * 2 * r + v1 * (1 - r*r) / (1 + r*r) = 0;
	-v2 * 2 * r + v1 * (1 - r*r) = 0;
	(-v1) * r * r + (-2 * v2) * r + (v1) = 0;

	*/

	*tanHalfAngle = 0.0f;

	// transform rotation axis to z-axis
    VectorSubtract( point, tw->origin, tmp );
    MatrixRotateVector( tmp, tw->matrix, p );
    MatrixRotateVector( plane, tw->matrix, normal );

	v1 = normal[0] * p[1] - normal[1] * p[0];
	v2 = normal[0] * p[0] + normal[1] * p[1];

	// the point will always start at the front of the plane, therefore v0 + v2 > 0 is always true
	if ( angle < 0.0f ) {
		*dir = -v1;
	}
	else {
		*dir = v1;
	}
	// negative direction means the point moves towards the plane at the initial position
	if ( *dir <= 0.0f ) {
		return qtrue;
	}

	a = -v1;
	b = -v2;
	c = v1;
	if ( a == 0.0f ) {
		if ( b == 0.0f ) {
			return qfalse;
		}
		frac1 = -c / ( 2.0f * b );
		frac2 = 1e10;	// = tan( idMath::HALF_PI )
	}
	else {
		d = b * b - c * a;
		if ( d <= 0.0f ) {
			return qfalse;
		}
		sqrtd = sqrt( d );
		if ( b > 0.0f ) {
			q = - b + sqrtd;
		}
		else {
			q = - b - sqrtd;
		}
		frac1 = q / a;
		frac2 = c / q;
	}

	if ( angle < 0.0f ) {
		frac1 = -frac1;
		frac2 = -frac2;
	}

	if ( frac1 < 0.0f && frac2 < 0.0f ) {
		return qfalse;
	}

	if ( frac1 > frac2 ) {
		*tanHalfAngle = frac1;
	}
	else {
		*tanHalfAngle = frac2;
	}

	if ( angle < 0.0f ) {
		*tanHalfAngle = -*tanHalfAngle;
	}

	return qtrue;
}

/*
================
CM_RotatePointThroughEpsilonPlane
================
*/
int CM_RotatePointThroughEpsilonPlane( const cm_traceWork_t *tw, const vec3_t point, const vec3_t endPoint,
							const plane_t plane, const float angle, const vec3_t origin,
							float *tanHalfAngle, vec3_t collisionPoint, vec3_t endDir ) {
	float d, dir, startTan;
	vec3_t vec, startDir, tmp;
	plane_t epsPlane;

	// epsilon expanded plane
    PlaneCopy( plane, epsPlane );
	epsPlane[3] = ( epsPlane[3] + CM_CLIP_EPSILON );

	// if the rotation sphere at the rotation origin is too far away from the polygon plane
	d = PlaneDistance( epsPlane, origin );
	VectorSubtract( point, origin, vec );
	if ( d * d > DotProduct( vec, vec ) ) {
		return qfalse;
	}

	// calculate direction of motion at vertex start position
	VectorSubtract( point, origin, tmp );
	CrossProduct( tmp, tw->axis, startDir );
	if ( angle < 0.0f ) {
		VectorNegate( startDir, startDir );
	}
	// if moving away from plane at start position
	if ( DotProduct( startDir, epsPlane ) >= 0.0f ) {
		// if end position is outside epsilon range
		d = PlaneDistance( epsPlane, endPoint );
		if ( d >= 0.0f ) {
			return qfalse;	// no collision
		}
		// calculate direction of motion at vertex end position
		VectorSubtract( endPoint, origin, tmp );
		CrossProduct( tmp, tw->axis, endDir );
		if ( angle < 0.0f ) {
			VectorNegate( endDir, endDir );
		}
		// if also moving away from plane at end position
		if ( DotProduct( endDir, epsPlane ) > 0.0f ) {
			return qfalse; // no collision
		}
	}

	// if the start position is in the epsilon range
	d = PlaneDistance( epsPlane, point );
	if ( d <= CM_PL_RANGE_EPSILON ) {

		// calculate tangent of half the rotation for which the vertex is furthest away from the plane
		if ( !CM_PointFurthestFromPlane( tw, point, plane, angle, &startTan, &dir ) ) {
			return qfalse;
		}

		if ( dir <= 0.0f ) {
			// moving towards the polygon plane so stop immediately
			*tanHalfAngle = 0.0f;
		}
		else if ( fabs( startTan ) >= tw->maxTan ) {
			// never going to get beyond the start tangent during the current rotation
			return qfalse;
		}
		else {
			// calculate collision with epsilon expanded plane
			if ( !CM_RotatePointThroughPlane( tw, point, epsPlane, angle, fabs( startTan ), tanHalfAngle ) ) {
				*tanHalfAngle = startTan;
			}
		}
	}
	else {
		// calculate collision with epsilon expanded plane
		if ( !CM_RotatePointThroughPlane( tw, point, epsPlane, angle, 0.0f, tanHalfAngle ) ) {
			return qfalse;
		}
	}

	// calculate collision point
	VectorCopy( point, collisionPoint );
	if ( *tanHalfAngle != 0.0f ) {
		CM_RotatePoint( collisionPoint, tw->origin, tw->axis, *tanHalfAngle );
	}
	// calculate direction of motion at collision point
	VectorSubtract( collisionPoint, origin, tmp );
	CrossProduct( tmp, tw->axis, endDir );
	if ( angle < 0.0f ) {
		VectorNegate( endDir, endDir );
	}
	return qtrue;
}

/*
================
CM_RotateTrmVertexThroughPolygon
================
*/
void CM_RotateTrmVertexThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *poly, cm_trmVertex_t *v, int vertexNum ) {
	int i;
	float tanHalfAngle;
	vec3_t endDir, collisionPoint;
	vec6_t pl;

	// if the trm vertex is behind the polygon plane it cannot collide with the polygon within a 180 degrees rotation
	if ( tw->isConvex && tw->axisIntersectsTrm && v->polygonSide ) {
		return;
	}

	// if the trace model vertex rotation bounds do not intersect the polygon bounds
	if ( !BoundsIntersect( v->rotationBounds[0], v->rotationBounds[1], poly->bounds[0], poly->bounds[1] ) ) {
		return;
	}

	// vertex rotation bounds should cross polygon plane
	if ( BoxOnPlaneSideSlow( v->rotationBounds[0], v->rotationBounds[1], poly->plane, ON_EPSILON ) != PLANESIDE_CROSS ) {
		return;
	}

	// rotate the vertex through the epsilon plane
	if ( !CM_RotatePointThroughEpsilonPlane( tw, v->p, v->endp, poly->plane, tw->angle, v->rotationOrigin,
											&tanHalfAngle, collisionPoint, endDir ) ) {
		return;
	}

	if ( fabs( tanHalfAngle ) < tw->maxTan ) {
		// verify if 'collisionPoint' moving along 'endDir' moves between polygon edges
		PlueckerFromRay( pl, collisionPoint, endDir );
		for ( i = 0; i < poly->numEdges; i++ ) {
			if ( poly->edges[i] < 0 ) {
				if ( PlueckerPermutedInnerProduct( pl, tw->polygonEdgePlueckerCache[i] ) > 0.0f ) {
					return;
				}
			}
			else {
				if ( PlueckerPermutedInnerProduct( pl, tw->polygonEdgePlueckerCache[i] ) < 0.0f ) {
					return;
				}
			}
		}
		tw->maxTan = fabs( tanHalfAngle );
		// collision plane is the polygon plane
		VectorCopy( poly->plane, tw->trace.c.normal );
		tw->trace.c.dist = PlaneGetDist( poly->plane );
		tw->trace.c.contents = poly->contents;
		tw->trace.c.material = poly->material;
		tw->trace.c.type = CONTACT_TRMVERTEX;
		tw->trace.c.modelFeature = *(int *)(&poly);
		tw->trace.c.trmFeature = v - tw->vertices;
		VectorCopy( collisionPoint, tw->trace.c.point );
	}
}

/*
================
CM_RotateVertexThroughTrmPolygon
================
*/
void CM_RotateVertexThroughTrmPolygon( cm_traceWork_t *tw, cm_trmPolygon_t *trmpoly, cm_polygon_t *poly, cm_vertex_t *v, vec3_t rotationOrigin ) {
	int i, edgeNum;
	float tanHalfAngle;
	vec3_t dir, endp, endDir, collisionPoint;
	vec6_t pl;
	cm_trmEdge_t *edge;

	// if the polygon vertex is behind the trm plane it cannot collide with the trm polygon within a 180 degrees rotation
	if ( tw->isConvex && tw->axisIntersectsTrm
			&& PlaneDistance( trmpoly->plane, v->p ) < 0.0f ) {
		return;
	}

	// if the model vertex is outside the trm polygon rotation bounds
	if ( !BoundsIntersectPoint( trmpoly->rotationBounds[0], trmpoly->rotationBounds[1], v->p ) ) {
		return;
	}

	// if the rotation axis goes through the polygon vertex
	VectorSubtract( v->p, rotationOrigin, dir );
	if ( DotProduct( dir, dir ) < ROTATION_AXIS_EPSILON * ROTATION_AXIS_EPSILON ) {
		return;
	}

	// calculate vertex end position
	VectorCopy( v->p, endp );
	RotationRotatePoint( &tw->modelVertexRotation, endp );

	// rotate the vertex through the epsilon plane
	if ( !CM_RotatePointThroughEpsilonPlane( tw, v->p, endp, trmpoly->plane, -tw->angle, rotationOrigin,
											&tanHalfAngle, collisionPoint, endDir ) ) {
		return;
	}

	if ( fabs( tanHalfAngle ) < tw->maxTan ) {
		// verify if 'collisionPoint' moving along 'endDir' moves between polygon edges
		PlueckerFromRay( pl, collisionPoint, endDir );
		for ( i = 0; i < trmpoly->numEdges; i++ ) {
			edgeNum = trmpoly->edges[i];
			edge = tw->edges + abs(edgeNum);
			if ( edgeNum < 0 ) {
				if ( PlueckerPermutedInnerProduct( pl, edge->pl ) > 0.0f ) {
					return;
				}
			}
			else {
				if ( PlueckerPermutedInnerProduct( pl, edge->pl ) < 0.0f ) {
					return;
				}
			}
		}
		tw->maxTan = fabs( tanHalfAngle );
		// collision plane is the flipped trm polygon plane
		VectorNegate( trmpoly->plane, tw->trace.c.normal );
		tw->trace.c.dist = DotProduct( tw->trace.c.normal, v->p );
		tw->trace.c.contents = poly->contents;
		tw->trace.c.material = poly->material;
		tw->trace.c.type = CONTACT_MODELVERTEX;
		tw->trace.c.modelFeature = v - tw->model->vertices;
		tw->trace.c.trmFeature = trmpoly - tw->polys;
		VectorCopy( v->p, tw->trace.c.point );
	}
}

/*
================
CM_RotateTrmThroughPolygon

  returns true if the polygon blocks the complete rotation
================
*/
qboolean CM_RotateTrmThroughPolygon( cm_traceWork_t *tw, cm_polygon_t *p ) {
	int i, j, k, edgeNum;
	float d;
	cm_trmVertex_t *bv;
	cm_trmEdge_t *be;
	cm_trmPolygon_t *bp;
	cm_vertex_t *v;
	cm_edge_t *e;
	vec3_t *rotationOrigin;
	vec3_t tmp1;

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

	// back face culling
	if ( tw->isConvex ) {
		// if the center of the convex trm is behind the polygon plane
		if ( PlaneDistance( p->plane, tw->start ) < 0.0f ) {
			// if the rotation axis intersects the trace model
			if ( tw->axisIntersectsTrm ) {
				return qfalse;
			}
			else {
				// if the direction of motion at the start and end position of the
				// center of the trm both go towards or away from the polygon plane
				// or if the intersections of the rotation axis with the expanded heart planes
				// are both in front of the polygon plane
			}
		}
	}

	// if the polygon is too far from the first heart plane
	d = BoundsDistanceToPlane( p->bounds, tw->heartPlane1 );
	if ( fabs(d) > tw->maxDistFromHeartPlane1 ) {
		return qfalse;
	}

	// rotation bounds should cross polygon plane
	switch( BoxOnPlaneSideSlow( tw->bounds[0], tw->bounds[1], p->plane, ON_EPSILON ) ) {
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

	for ( i = 0; i < tw->numVerts; i++ ) {
		bv = tw->vertices + i;
		// calculate polygon side this vertex is on
		d = PlaneDistance( p->plane, bv->p );
		bv->polygonSide = FLOATSIGNBITSET( d );
	}

	for ( i = 0; i < p->numEdges; i++ ) {
		edgeNum = p->edges[i];
		e = tw->model->edges + abs(edgeNum);
		v = tw->model->vertices + e->vertexNum[INTSIGNBITSET(edgeNum)];

		// pluecker coordinate for edge
		PlueckerFromLine( tw->polygonEdgePlueckerCache[i],
						tw->model->vertices[e->vertexNum[0]].p,
						tw->model->vertices[e->vertexNum[1]].p );

		// calculate rotation origin projected into rotation plane through the vertex
		VectorSubtract( v->p, tw->origin, tmp1 );
		VectorMA( tw->origin, DotProduct( tw->axis, tmp1 ), tw->axis, tw->polygonRotationOriginCache[i] );
	}
	// copy first to last so we can easily cycle through
	VectorCopy( tw->polygonRotationOriginCache[0], tw->polygonRotationOriginCache[p->numEdges] );

	// fast point rotation
	if ( tw->pointTrace ) {
		CM_RotateTrmVertexThroughPolygon( tw, p, &tw->vertices[0], 0 );
	}
	else {
		// rotate trm vertices through polygon
		for ( i = 0; i < tw->numVerts; i++ ) {
			bv = tw->vertices + i;
			if ( bv->used ) {
				CM_RotateTrmVertexThroughPolygon( tw, p, bv, i );
			}
		}

		// rotate trm edges through polygon
		for ( i = 1; i <= tw->numEdges; i++ ) {
			be = tw->edges + i;
			if ( be->used ) {
				CM_RotateTrmEdgeThroughPolygon( tw, p, be );
			}
		}

		// rotate all polygon vertices through the trm
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

				// if the vertex is outside the trm rotation bounds
				if ( !BoundsIntersectPoint( tw->bounds[0], tw->bounds[1], v->p ) ) {
					continue;
				}

				rotationOrigin = &tw->polygonRotationOriginCache[i+k];

				for ( j = 0; j < tw->numPolys; j++ ) {
					bp = tw->polys + j;
					if ( bp->used ) {
						CM_RotateVertexThroughTrmPolygon( tw, bp, p, v, *rotationOrigin );
					}
				}
			}
		}
	}

	return ( tw->maxTan == 0.0f );
}

/*
================
CM_BoundsForRotation

  only for rotations < 180 degrees
================
*/
void CM_BoundsForRotation( const vec3_t origin, const vec3_t axis, const vec3_t start, const vec3_t end, vec3_t bounds[2] ) {
	int i;
	float radiusSqr;
	vec3_t v1, v2, tmp;

	VectorSubtract( start, origin, tmp );
	radiusSqr = VectorLengthSquared( tmp );
	CrossProduct( tmp, axis, v1 );
	VectorSubtract( end, origin, tmp );
	CrossProduct( tmp, axis, v2 );

	for ( i = 0; i < 3; i++ ) {
		// if the derivative changes sign along this axis during the rotation from start to end
		if ( ( v1[i] > 0.0f && v2[i] < 0.0f ) || ( v1[i] < 0.0f && v2[i] > 0.0f ) ) {
			if ( ( 0.5f * (start[i] + end[i]) - origin[i] ) > 0.0f ) {
				bounds[0][i] = ( start[i] < end[i] ? start[i] : end[i] );
				bounds[1][i] = origin[i] + sqrt( radiusSqr * ( 1.0f - axis[i] * axis[i] ) );
			}
			else {
				bounds[0][i] = origin[i] - sqrt( radiusSqr * ( 1.0f - axis[i] * axis[i] ) );
				bounds[1][i] = ( start[i] > end[i] ? start[i] : end[i] );
			}
		}
		else if ( start[i] > end[i] ) {
			bounds[0][i] = end[i];
			bounds[1][i] = start[i];
		}
		else {
			bounds[0][i] = start[i];
			bounds[1][i] = end[i];
		}
		// expand for epsilons
		bounds[0][i] -= CM_BOX_EPSILON;
		bounds[1][i] += CM_BOX_EPSILON;
	}
}

/*
================
CM_Rotation180
================
*/
void CM_Rotation180( cm_trace_t *results, const vec3_t rorg, const vec3_t axis,
										const float startAngle, const float endAngle, const vec3_t start,
										const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
										cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] ) {
	int i, j, edgeNum;
	float d, maxErr, initialTan;
	qboolean model_rotated, trm_rotated;
	vec3_t tmp, vr, vup, at, bt;
	vec3_t invModelAxis[3], tmpAxis[3];
	rotation_t startRotation, endRotation;
	vec6_t plaxis;
	cm_trmPolygon_t *poly;
	cm_trmEdge_t *edge;
	cm_trmVertex_t *vert;
	static cm_traceWork_t tw; // TODO: ALIGN16

	if ( model < 0 || model > MAX_SUBMODELS || model > cmLocal.maxModels ) {
		ii.Com_Printf("CM_Rotation180: invalid model handle\n");
		return;
	}
	if ( !cmLocal.models[model] ) {
		ii.Com_Printf("CM_Rotation180: invalid model\n");
		return;
	}

	cmLocal.checkCount++;

	tw.trace.fraction = 1.0f;
	tw.trace.c.contents = 0;
	tw.trace.c.type = CONTACT_NONE;
	tw.contents = contentMask;
	tw.isConvex = qtrue;
	tw.rotation = qtrue;
	tw.positionTest = qfalse;
	tw.axisIntersectsTrm = qfalse;
	tw.quickExit = qfalse;
	tw.angle = endAngle - startAngle;
	assert( tw.angle > -180.0f && tw.angle < 180.0f );
	tw.maxTan = initialTan = fabs( tan( ( M_PI / 360.0f ) * tw.angle ) );
	tw.model = cmLocal.models[model];
	VectorSubtract( start, modelOrigin, tw.start );
	// rotation axis, axis is assumed to be normalized
	VectorCopy( axis, tw.axis );

	assert( tw.axis[0] * tw.axis[0] + tw.axis[1] * tw.axis[1] + tw.axis[2] * tw.axis[2] > 0.99f );
	// rotation origin projected into rotation plane through tw.start
	VectorSubtract( rorg, modelOrigin, tw.origin );
	d = DotProduct( tw.axis, tw.origin ) - DotProduct( tw.axis, tw.start );
	VectorScale( tw.axis, d, tmp );
	VectorMA( tw.origin, -d, tw.axis, tw.origin );
	// radius of rotation
	VectorSubtract( tw.start, tw.origin, tmp );
	tw.radius = VectorLength( tmp );
	// maximum error of the circle approximation traced through the axial BSP tree
	d = tw.radius * tw.radius - (CIRCLE_APPROXIMATION_LENGTH*CIRCLE_APPROXIMATION_LENGTH*0.25f);
	if ( d > 0.0f ) {
		maxErr = tw.radius - sqrt( d );
	} else {
		maxErr = tw.radius;
	}

	model_rotated = AxisIsRotated( modelAxis );
	if ( model_rotated ) {
		TransposeAxis( modelAxis, invModelAxis );
		VectorRotateSelf( tw.axis, invModelAxis );
		VectorRotateSelf( tw.origin, invModelAxis );
	}

	RotationFromOriginAndVector( &startRotation, tw.origin, tw.axis, startAngle );
	RotationFromOriginAndVector( &endRotation, tw.origin, tw.axis, endAngle );

	// create matrix which rotates the rotation axis to the z-axis
	NormalVectors( tw.axis, vr, vup );
	
	tw.matrix[0][0] = vr[0];
	tw.matrix[1][0] = vr[1];
	tw.matrix[2][0] = vr[2];
	tw.matrix[0][1] = -vup[0];
	tw.matrix[1][1] = -vup[1];
	tw.matrix[2][1] = -vup[2];
	tw.matrix[0][2] = tw.axis[0];
	tw.matrix[1][2] = tw.axis[1];
	tw.matrix[2][2] = tw.axis[2];

	// if optimized point trace
	if ( !trm || ( trm->bounds[1][0] - trm->bounds[0][0] <= 0.0f &&
					trm->bounds[1][1] - trm->bounds[0][1] <= 0.0f &&
					trm->bounds[1][2] - trm->bounds[0][2] <= 0.0f ) ) {

		if ( model_rotated ) {
			// rotate trace instead of model
			VectorRotateSelf( tw.start, invModelAxis );
		}
		VectorCopy( tw.start, tw.end );
		// if we start at a specific angle
		if ( startAngle != 0.0f ) {
			RotationRotatePoint( &startRotation, tw.start );
		}
		// calculate end position of rotation
		RotationRotatePoint( &endRotation, tw.end );

		// calculate rotation origin projected into rotation plane through the vertex
		tw.numVerts = 1;
		VectorCopy( tw.start, tw.vertices[0].p );
		VectorCopy( tw.end, tw.vertices[0].endp );
		tw.vertices[0].used = qtrue;
		VectorSubtract( tw.vertices[0].p, tw.origin, tmp );
		VectorMA( tw.origin, DotProduct( tw.axis, tmp ), tw.axis, tw.vertices[0].rotationOrigin );
		CM_BoundsForRotation( tw.vertices[0].rotationOrigin, tw.axis, tw.start, tw.end, tw.vertices[0].rotationBounds );
		// rotation bounds
		VectorCopy( tw.vertices[0].rotationBounds[0], tw.bounds[0] );
		VectorCopy( tw.vertices[0].rotationBounds[1], tw.bounds[1] );
		tw.numEdges = tw.numPolys = 0;

		// collision with single point
		tw.pointTrace = qtrue;

		// extents is set to maximum error of the circle approximation traced through the axial BSP tree
		tw.extents[0] = tw.extents[1] = tw.extents[2] = maxErr + CM_BOX_EPSILON;

		// setup rotation heart plane
		VectorCopy( tw.axis, tw.heartPlane1 );
		PlaneFitThroughPoint( tw.heartPlane1, tw.start );
		tw.maxDistFromHeartPlane1 = CM_BOX_EPSILON;

		// trace through the model
		CM_TraceThroughModel( &tw );

		// store results
		*results = tw.trace;
		VectorCopy( start, results->endpos );
		if ( tw.maxTan == initialTan ) {
			results->fraction = 1.0f;
		} else {
			results->fraction = fabs( atan( tw.maxTan ) * ( 2.0f * 180.0f / M_PI ) / tw.angle );
		}
		assert( results->fraction <= 1.0f );
		RotationFromOriginAndVector( &endRotation, rorg, axis, startAngle + (endAngle-startAngle) * results->fraction );
		RotationRotatePoint( &endRotation, results->endpos );
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
		return;
	}

	tw.pointTrace = qfalse;

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
			VectorRotateSelf( tw.vertices[i].p, invModelAxis );
		}
	}
	for ( i = 0; i < tw.numVerts; i++ ) {
		VectorCopy( tw.vertices[i].p, tw.vertices[i].endp );
	}
	// if we start at a specific angle
	if ( startAngle != 0.0f ) {
		for ( i = 0; i < tw.numVerts; i++ ) {
			RotationRotatePoint( &startRotation, tw.vertices[i].p );
		}
	}
	for ( i = 0; i < tw.numVerts; i++ ) {
		// end position of vertex
		RotationRotatePoint( &endRotation, tw.vertices[i].endp );
	}

	// add offset to start point
	if ( trm_rotated ) {
        MatrixRotateVector( trm->offset, trmAxis, tmp );
		VectorAdd( tw.start, tmp, tw.start );
	} else {
		VectorAdd( tw.start, trm->offset, tw.start );
	}
	// if the model is rotated
	if ( model_rotated ) {
		// rotate trace instead of model
		VectorRotateSelf( tw.start, invModelAxis );
	}
	VectorCopy( tw.start, tw.end );
	// if we start at a specific angle
	if ( startAngle != 0.0f ) {
		RotationRotatePoint( &startRotation, tw.start );
	}
	// calculate end position of rotation
	RotationRotatePoint( &endRotation, tw.end );

	// setup trm vertices
	for ( vert = tw.vertices, i = 0; i < tw.numVerts; i++, vert++ ) {
		// calculate rotation origin projected into rotation plane through the vertex
		VectorSubtract( vert->p, tw.origin, tmp );
		VectorMA( tw.origin, DotProduct( tw.axis, tmp ), tw.axis, vert->rotationOrigin );
		// calculate rotation bounds for this vertex
		CM_BoundsForRotation( vert->rotationOrigin, tw.axis, vert->p, vert->endp, vert->rotationBounds );
		// if the rotation axis goes through the vertex then the vertex is not used
		VectorSubtract( vert->p, vert->rotationOrigin, tmp );
		d = VectorLengthSquared( tmp );
		if ( d > ROTATION_AXIS_EPSILON * ROTATION_AXIS_EPSILON ) {
			vert->used = qtrue;
		}
	}

	// setup trm edges
	for ( edge = tw.edges + 1, i = 1; i <= tw.numEdges; i++, edge++ ) {
		// if the rotation axis goes through both the edge vertices then the edge is not used
		if ( tw.vertices[edge->vertexNum[0]].used | tw.vertices[edge->vertexNum[1]].used ) {
			edge->used = qtrue;
		}
		// edge start, end and pluecker coordinate
		VectorCopy( tw.vertices[edge->vertexNum[0]].p, edge->start );
		VectorCopy( tw.vertices[edge->vertexNum[1]].p, edge->end );
		PlueckerFromLine( edge->pl, edge->start, edge->end );
		// pluecker coordinate for edge being rotated about the z-axis
		VectorSubtract( edge->start, tw.origin, tmp );
		MatrixRotateVector( tmp, tw.matrix, at );
		VectorSubtract( edge->end, tw.origin, tmp );
		MatrixRotateVector( tmp, tw.matrix, bt );
		PlueckerFromLine( edge->plzaxis, at, bt );
		// get edge rotation bounds from the rotation bounds of both vertices
		VectorCopy( tw.vertices[edge->vertexNum[0]].rotationBounds[0], edge->rotationBounds[0] );
		VectorCopy( tw.vertices[edge->vertexNum[0]].rotationBounds[1], edge->rotationBounds[1] );
		AddBoundsToBounds( tw.vertices[edge->vertexNum[1]].rotationBounds, edge->rotationBounds );
		// used to calculate if the rotation axis intersects the trm
		edge->bitNum = 0;
	}

	ClearBounds( tw.bounds[0], tw.bounds[1] );

	// rotate trm polygon planes
	if ( trm_rotated & model_rotated ) {
		MatrixMultiply( trmAxis, invModelAxis, tmpAxis );
		for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
			VectorRotateSelf( poly->plane, tmpAxis );
		}
	} else if ( trm_rotated ) {
		for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
			VectorRotateSelf( poly->plane, trmAxis );
		}
	} else if ( model_rotated ) {
		for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
			VectorRotateSelf( poly->plane, invModelAxis );
		}
	}

	// setup trm polygons
	for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
		poly->used = qtrue;
		// set trm polygon plane distance
		PlaneFitThroughPoint( poly->plane, tw.edges[abs(poly->edges[0])].start );
		// get polygon bounds from edge bounds
		ClearBounds( poly->rotationBounds[0], poly->rotationBounds[1] );
		for ( j = 0; j < poly->numEdges; j++ ) {
			// add edge rotation bounds to polygon rotation bounds
			edge = &tw.edges[abs( poly->edges[j] )];
			AddBoundsToBounds( edge->rotationBounds, poly->rotationBounds );
		}
		// get trace bounds from polygon bounds
		AddBoundsToBounds( poly->rotationBounds, tw.bounds );
	}

	// extents including the maximum error of the circle approximation traced through the axial BSP tree
	for ( i = 0; i < 3; i++ ) {
		tw.size[0][i] = tw.bounds[0][i] - tw.start[i];
		tw.size[1][i] = tw.bounds[1][i] - tw.start[i];
		if ( fabs( tw.size[0][i] ) > fabs( tw.size[1][i] ) ) {
			tw.extents[i] = fabs( tw.size[0][i] ) + maxErr + CM_BOX_EPSILON;
		} else {
			tw.extents[i] = fabs( tw.size[1][i] ) + maxErr + CM_BOX_EPSILON;
		}
	}

	// for back-face culling
	if ( tw.isConvex ) {
		if ( tw.start == tw.origin ) {
			tw.axisIntersectsTrm = qtrue;
		} else {
			// determine if the rotation axis intersects the trm
			PlueckerFromRay( plaxis, tw.origin, tw.axis );
			for ( poly = tw.polys, i = 0; i < tw.numPolys; i++, poly++ ) {
				// back face cull polygons
				if ( DotProduct( poly->plane, tw.axis ) > 0.0f ) {
					continue;
				}
				// test if the axis goes between the polygon edges
				for ( j = 0; j < poly->numEdges; j++ ) {
					edgeNum = poly->edges[j];
					edge = tw.edges + abs(edgeNum);
					if ( !(edge->bitNum & 2) ) {
						d = PlueckerPermutedInnerProduct( plaxis, edge->pl );
						edge->bitNum = FLOATSIGNBITSET( d ) | 2;
					}
					if ( ( edge->bitNum ^ INTSIGNBITSET( edgeNum ) ) & 1 ) {
						break;
					}
				}
				if ( j >= poly->numEdges ) {
					tw.axisIntersectsTrm = qtrue;
					break;
				}
			}
		}
	}

	// setup rotation heart plane
	VectorCopy( tw.axis, tw.heartPlane1 );
	PlaneFitThroughPoint( tw.heartPlane1, tw.start );
	tw.maxDistFromHeartPlane1 = 0.0f;
	for ( i = 0; i < tw.numVerts; i++ ) {
		d = fabs( PlaneDistance( tw.heartPlane1, tw.vertices[i].p ) );
		if ( d > tw.maxDistFromHeartPlane1 ) {
			tw.maxDistFromHeartPlane1 = d;
		}
	}
	tw.maxDistFromHeartPlane1 += CM_BOX_EPSILON;

	// inverse rotation to rotate model vertices towards trace model
	RotationFromOriginAndVector( &tw.modelVertexRotation, tw.origin, tw.axis, -tw.angle );

	// trace through the model
	CM_TraceThroughModel( &tw );

	// store results
	*results = tw.trace;
	VectorCompare( start, results->endpos );
	if ( tw.maxTan == initialTan ) {
		results->fraction = 1.0f;
	} else {
		results->fraction = fabs( atan( tw.maxTan ) * ( 2.0f * 180.0f / M_PI ) / tw.angle );
	}
	assert( results->fraction <= 1.0f );
	RotationFromOriginAndVector( &endRotation, rorg, axis, startAngle + (endAngle-startAngle) * results->fraction );
	RotationReCalculateMatrix( &endRotation );
	RotationRotatePoint( &endRotation, results->endpos );
	MatrixMultiply( trmAxis, endRotation.axis, results->endAxis );

	if ( results->fraction < 1.0f ) {
		// rotate trace plane normal if there was a collision with a rotated model
		if ( model_rotated ) {
			VectorRotateSelf( results->c.normal, modelAxis );
			VectorRotateSelf( results->c.point, modelAxis );
		}
		VectorAdd( results->c.point, modelOrigin, results->c.point );
		results->c.dist += DotProduct( modelOrigin, results->c.normal );
	}
}

/*
================
CM_Rotation
================
*/
#ifdef _DEBUG
static int entered = 0;
#endif

void CM_Rotation( cm_trace_t *results, const vec3_t start, const rotation_t *rotation,
										const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
										cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] ) {
	vec3_t tmp;
	float maxa, stepa, a, lasta;

	assert( ((byte *)&start) < ((byte *)results) || ((byte *)&start) > (((byte *)results) + sizeof( cm_trace_t )) );
	assert( ((byte *)&trmAxis) < ((byte *)results) || ((byte *)&trmAxis) > (((byte *)results) + sizeof( cm_trace_t )) );

	memset( results, 0, sizeof( *results ) );

	// if special position test
	if ( rotation->angle == 0.0f ) {
		CM_ContentsTrm( results, start, trm, trmAxis, contentMask, model, modelOrigin, modelAxis );
		return;
	}

#ifdef _DEBUG
	qboolean startsolid = qfalse;
	// test whether or not stuck to begin with
	if ( cm_debugCollision.integerValue ) {
		if ( !entered ) {
			entered = 1;
			// if already messed up to begin with
			if ( CM_Contents( start, trm, trmAxis, -1, model, modelOrigin, modelAxis ) & contentMask ) {
				startsolid = qtrue;
			}
			entered = 0;
		}
	}
#endif

	if ( rotation->angle >= 180.0f || rotation->angle <= -180.0f) {
		if ( rotation->angle >= 360.0f ) {
			maxa = 360.0f;
			stepa = 120.0f;			// three steps strictly < 180 degrees
		} else if ( rotation->angle <= -360.0f ) {
			maxa = -360.0f;
			stepa = -120.0f;		// three steps strictly < 180 degrees
		} else {
			maxa = rotation->angle;
			stepa = rotation->angle * 0.5f;	// two steps strictly < 180 degrees
		}
		for ( lasta = 0.0f, a = stepa; fabs( a ) < fabs( maxa ) + 1.0f; lasta = a, a += stepa ) {
			// partial rotation
			CM_Rotation180( results, rotation->origin, rotation->vec, lasta, a, start, trm, trmAxis, contentMask, model, modelOrigin, modelAxis );
			// if there is a collision
			if ( results->fraction < 1.0f ) {
				// fraction of total rotation
				results->fraction = (lasta + stepa * results->fraction) / rotation->angle;
				return;
			}
		}
		results->fraction = 1.0f;
		return;
	}

	CM_Rotation180( results, rotation->origin, rotation->vec, 0.0f, rotation->angle, start, trm, trmAxis, contentMask, model, modelOrigin, modelAxis );

#ifdef _DEBUG
	// test for missed collisions
	if ( cm_debugCollision.integerValue) ) {
		if ( !entered ) {
			entered = 1;
			// if the trm is stuck in the model
			if ( CM_Contents( results->endpos, trm, results->endAxis, -1, model, modelOrigin, modelAxis ) & contentMask ) {
				cm_trace_t tr;

				// test where the trm is stuck in the model
				CM_Contents( results->endpos, trm, results->endAxis, -1, model, modelOrigin, modelAxis );
				// re-run collision detection to find out where it failed
				CM_Rotation( &tr, start, rotation, trm, trmAxis, contentMask, model, modelOrigin, modelAxis );
			}
			entered = 0;
		}
	}
#endif
}
