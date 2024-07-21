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

Contents test

===============================================================================
*/

/*
================
CM_TestTrmVertsInBrush

  returns true if any of the trm vertices is inside the brush
================
*/
qboolean CM_TestTrmVertsInBrush( cm_traceWork_t *tw, cm_brush_t *b ) {
	int i, j, numVerts, bestPlane;
	float d, bestd;
	vec3_t *p;

	if ( b->checkcount == cmLocal.checkCount ) {
		return qfalse;
	}
	b->checkcount = cmLocal.checkCount;

	if ( !(b->contents & tw->contents) ) {
		return qfalse;
	}

	// if the brush bounds don't intersect the trace bounds
	if ( !BoundsIntersect( b->bounds[0], b->bounds[1], tw->bounds[0], tw->bounds[1] ) ) {
		return qfalse;
	}

	if ( tw->pointTrace ) {
		numVerts = 1;
	}
	else {
		numVerts = tw->numVerts;
	}

	for ( j = 0; j < numVerts; j++ ) {
		p = &tw->vertices[j].p;

		// see if the point is inside the brush
		bestPlane = 0;
		bestd = -Q_INFINITY;
		for ( i = 0; i < b->numPlanes; i++ ) {
			d = PlaneDistance( b->planes[i], *p );
			if ( d >= 0.0f ) {
				break;
			}
			if ( d > bestd ) {
				bestd = d;
				bestPlane = i;
			}
		}
		if ( i >= b->numPlanes ) {
			tw->trace.fraction = 0.0f;
			tw->trace.c.type = CONTACT_TRMVERTEX;
            VectorCopy( b->planes[bestPlane], tw->trace.c.normal );
			tw->trace.c.dist = PlaneGetDist( b->planes[bestPlane] );
			tw->trace.c.contents = b->contents;
			tw->trace.c.material = b->material;
            VectorCopy( *p, tw->trace.c.point );
			tw->trace.c.modelFeature = 0;
			tw->trace.c.trmFeature = j;
			return qtrue;
		}
	}
	return qfalse;
}

/*
================
CM_SetTrmEdgeSidedness
================
*/
#define CM_SetTrmEdgeSidedness( edge, bpl, epl, bitNum ) {							\
	if ( !(edge->sideSet & (1<<bitNum)) ) {											\
		float fl;																	\
		fl = PlueckerPermutedInnerProduct( bpl, epl );							    \
		edge->side = (edge->side & ~(1<<bitNum)) | (FLOATSIGNBITSET(fl) << bitNum);	\
		edge->sideSet |= (1 << bitNum);												\
	}																				\
}

/*
================
CM_SetTrmPolygonSidedness
================
*/
#define CM_SetTrmPolygonSidedness( v, plane, bitNum ) {								\
	if ( !((v)->sideSet & (1<<bitNum)) ) {											\
		float fl;																	\
		fl = PlaneDistance( (plane), (v)->p );					                    \
		/* cannot use float sign bit because it is undetermined when fl == 0.0f */	\
		if ( fl < 0.0f ) {															\
			(v)->side |= (1 << bitNum);												\
		}																			\
		else {																		\
			(v)->side &= ~(1 << bitNum);											\
		}																			\
		(v)->sideSet |= (1 << bitNum);												\
	}																				\
}

/*
================
CM_TestTrmInPolygon

  returns true if the trm intersects the polygon
================
*/
qboolean CM_TestTrmInPolygon( cm_traceWork_t *tw, cm_polygon_t *p ) {
	int i, j, k, edgeNum, flip, trmEdgeNum, bitNum, bestPlane;
	int sides[MAX_TRACEMODEL_VERTS];
	float d, bestd;
	cm_trmEdge_t *trmEdge;
	cm_edge_t *edge;
	cm_vertex_t *v, *v1, *v2;

	// if already checked this polygon
	if ( p->checkcount == cmLocal.checkCount ) {
		return qfalse;
	}
	p->checkcount = cmLocal.checkCount;

	// if this polygon does not have the right contents behind it
	if ( !(p->contents & tw->contents) ) {
		return qfalse;
	}

	// if the polygon bounds don't intersect the trace bounds
	if ( !BoundsIntersect( p->bounds[0], p->bounds[1], tw->bounds[0], tw->bounds[1] ) ) {
		return qfalse;
	}

	// bounds should cross polygon plane
	switch( BoxOnPlaneSideSlow( tw->bounds[0], tw->bounds[1], p->plane, ON_EPSILON ) ) {
		case PLANESIDE_CROSS: // cross
			break;
		case PLANESIDE_FRONT: // front
			if ( tw->model->isConvex ) {
				tw->quickExit = qtrue;
				return qtrue;
			}
		default:
			return qfalse;
	}

	// if the trace model is convex
	if ( tw->isConvex ) {
		// test if any polygon vertices are inside the trm
		for ( i = 0; i < p->numEdges; i++ ) {
			edgeNum = p->edges[i];
			edge = tw->model->edges + abs(edgeNum);
			// if this edge is already tested
			if ( edge->checkcount == cmLocal.checkCount ) {
				continue;
			}

			for ( j = 0; j < 2; j++ ) {
				v = &tw->model->vertices[edge->vertexNum[j]];
				// if this vertex is already tested
				if ( v->checkcount == cmLocal.checkCount ) {
					continue;
				}

				bestPlane = 0;
				bestd = -Q_INFINITY;
				for ( k = 0; k < tw->numPolys; k++ ) {
					d = PlaneDistance( tw->polys[k].plane, v->p );
					if ( d >= 0.0f ) {
						break;
					}
					if ( d > bestd ) {
						bestd = d;
						bestPlane = k;
					}
				}
				if ( k >= tw->numPolys ) {
					tw->trace.fraction = 0.0f;
					tw->trace.c.type = CONTACT_MODELVERTEX;
                    VectorNegate( tw->polys[bestPlane].plane, tw->trace.c.normal );
					tw->trace.c.dist = -PlaneGetDist( tw->polys[bestPlane].plane );
					tw->trace.c.contents = p->contents;
					tw->trace.c.material = p->material;
                    VectorCopy( v->p, tw->trace.c.point );
					tw->trace.c.modelFeature = edge->vertexNum[j];
					tw->trace.c.trmFeature = 0;
					return qtrue;
				}
			}
		}
	}

	for ( i = 0; i < p->numEdges; i++ ) {
		edgeNum = p->edges[i];
		edge = tw->model->edges + abs(edgeNum);
		// reset sidedness cache if this is the first time we encounter this edge
		if ( edge->checkcount != cmLocal.checkCount ) {
			edge->sideSet = 0;
		}
		// pluecker coordinate for edge
		PlueckerFromLine( tw->polygonEdgePlueckerCache[i],
                         tw->model->vertices[edge->vertexNum[0]].p,
						 tw->model->vertices[edge->vertexNum[1]].p );
		v = &tw->model->vertices[edge->vertexNum[INTSIGNBITSET(edgeNum)]];
		// reset sidedness cache if this is the first time we encounter this vertex
		if ( v->checkcount != cmLocal.checkCount ) {
			v->sideSet = 0;
		}
		v->checkcount = cmLocal.checkCount;
	}

	// get side of polygon for each trm vertex
	for ( i = 0; i < tw->numVerts; i++ ) {
		d = PlaneDistance( p->plane, tw->vertices[i].p );
		sides[i] = d < 0.0f ? -1 : 1;
	}

	// test if any trm edges go through the polygon
	for ( i = 1; i <= tw->numEdges; i++ ) {
		// if the trm edge does not cross the polygon plane
		if ( sides[tw->edges[i].vertexNum[0]] == sides[tw->edges[i].vertexNum[1]] ) {
			continue;
		}
		// check from which side to which side the trm edge goes
		flip = INTSIGNBITSET( sides[tw->edges[i].vertexNum[0]] );
		// test if trm edge goes through the polygon between the polygon edges
		for ( j = 0; j < p->numEdges; j++ ) {
			edgeNum = p->edges[j];
			edge = tw->model->edges + abs(edgeNum);
#if 1
			CM_SetTrmEdgeSidedness( edge, tw->edges[i].pl, tw->polygonEdgePlueckerCache[j], i );
			if ( INTSIGNBITSET(edgeNum) ^ ((edge->side >> i) & 1) ^ flip ) {
				break;
			}
#else
			d = PlueckerPermutedInnerProduct( tw->edges[i].pl, tw->polygonEdgePlueckerCache[j] );
			if ( flip ) {
				d = -d;
			}
			if ( edgeNum > 0 ) {
				if ( d <= 0.0f ) {
					break;
				}
			}
			else {
				if ( d >= 0.0f ) {
					break;
				}
			}
#endif
		}
		if ( j >= p->numEdges ) {
			tw->trace.fraction = 0.0f;
			tw->trace.c.type = CONTACT_EDGE;
            VectorCopy( p->plane, tw->trace.c.normal );
			tw->trace.c.dist = PlaneGetDist( p->plane );
			tw->trace.c.contents = p->contents;
			tw->trace.c.material = p->material;
            VectorCopy( tw->vertices[tw->edges[i].vertexNum[ !flip ]].p, tw->trace.c.point );
			tw->trace.c.modelFeature = *(int *)(&p);
			tw->trace.c.trmFeature = i;
			return qtrue;
		}
	}

	// test if any polygon edges go through the trm polygons
	for ( i = 0; i < p->numEdges; i++ ) {
		edgeNum = p->edges[i];
		edge = tw->model->edges + abs(edgeNum);
		if ( edge->checkcount == cmLocal.checkCount ) {
			continue;
		}
		edge->checkcount = cmLocal.checkCount;

		for ( j = 0; j < tw->numPolys; j++ ) {
#if 1
			v1 = tw->model->vertices + edge->vertexNum[0];
			CM_SetTrmPolygonSidedness( v1, tw->polys[j].plane, j );
			v2 = tw->model->vertices + edge->vertexNum[1];
			CM_SetTrmPolygonSidedness( v2, tw->polys[j].plane, j );
			// if the polygon edge does not cross the trm polygon plane
			if ( !(((v1->side ^ v2->side) >> j) & 1) ) {
				continue;
			}
			flip = (v1->side >> j) & 1;
#else
			float d1, d2;

			v1 = tw->model->vertices + edge->vertexNum[0];
			d1 = tw->polys[j].plane.Distance( v1->p );
			v2 = tw->model->vertices + edge->vertexNum[1];
			d2 = tw->polys[j].plane.Distance( v2->p );
			// if the polygon edge does not cross the trm polygon plane
			if ( (d1 >= 0.0f && d2 >= 0.0f) || (d1 <= 0.0f && d2 <= 0.0f) ) {
				continue;
			}
			flip = false;
			if ( d1 < 0.0f ) {
				flip = true;
			}
#endif
			// test if polygon edge goes through the trm polygon between the trm polygon edges
			for ( k = 0; k < tw->polys[j].numEdges; k++ ) {
				trmEdgeNum = tw->polys[j].edges[k];
				trmEdge = tw->edges + abs(trmEdgeNum);
#if 1
				bitNum = abs(trmEdgeNum);
				CM_SetTrmEdgeSidedness( edge, trmEdge->pl, tw->polygonEdgePlueckerCache[i], bitNum );
				if ( INTSIGNBITSET(trmEdgeNum) ^ ((edge->side >> bitNum) & 1) ^ flip ) {
					break;
				}
#else
				d = trmEdge->pl.PermutedInnerProduct( tw->polygonEdgePlueckerCache[i] );
				if ( flip ) {
					d = -d;
				}
				if ( trmEdgeNum > 0 ) {
					if ( d <= 0.0f ) {
						break;
					}
				}
				else {
					if ( d >= 0.0f ) {
						break;
					}
				}
#endif
			}
			if ( k >= tw->polys[j].numEdges ) {
				tw->trace.fraction = 0.0f;
				tw->trace.c.type = CONTACT_EDGE;
                VectorNegate( tw->polys[j].plane, tw->trace.c.normal );
				tw->trace.c.dist = -PlaneGetDist( tw->polys[j].plane );
				tw->trace.c.contents = p->contents;
				tw->trace.c.material = p->material;
                VectorCopy( tw->model->vertices[edge->vertexNum[ !flip ]].p, tw->trace.c.point );
				tw->trace.c.modelFeature = edgeNum;
				tw->trace.c.trmFeature = j;
				return qtrue;
			}
		}
	}
	return qfalse;
}

/*
================
CM_PointNode
================
*/
cm_node_t *CM_PointNode( const vec3_t p, cm_model_t *model ) {
	cm_node_t *node;

	node = model->node;
	while ( node->planeType != -1 ) {
		if (p[node->planeType] > node->planeDist) {
			node = node->children[0];
		}
		else {
			node = node->children[1];
		}

		assert( node != NULL );
	}
	return node;
}

/*
================
CM_PointContents
================
*/
int CM_PointContents( const vec3_t p, cmHandle_t model ) {
	int i;
	float d;
	cm_node_t *node;
	cm_brushRef_t *bref;
	cm_brush_t *b;
	plane_t *plane;

	node = CM_PointNode( p, cmLocal.models[model] );
	for ( bref = node->brushes; bref; bref = bref->next ) {
		b = bref->b;
		// test if the point is within the brush bounds
		for ( i = 0; i < 3; i++ ) {
			if ( p[i] < b->bounds[0][i] ) {
				break;
			}
			if ( p[i] > b->bounds[1][i] ) {
				break;
			}
		}
		if ( i < 3 ) {
			continue;
		}
		// test if the point is inside the brush
		plane = &b->planes[0];
		for ( i = 0; i < b->numPlanes; i++, plane++ ) {
			d = PlaneDistance( *plane, p );
			if ( d >= 0 ) {
				break;
			}
		}
		if ( i >= b->numPlanes ) {
			return b->contents;
		}
	}
	return 0;
}

/*
==================
CM_TransformedPointContents
==================
*/
int	CM_TransformedPointContents( const vec3_t p, cmHandle_t model, const vec3_t origin, const vec3_t modelAxis[3] ) {
	vec3_t p_l, tmp;

	// subtract origin offset
    VectorSubtract( p, origin, p_l );
	if ( AxisIsRotated( modelAxis ) ) {
        VectorRotateSelf( p_l, modelAxis );
	}
	return CM_PointContents( p_l, model );
}


/*
==================
CM_ContentsTrm
==================
*/
int CM_ContentsTrm( cm_trace_t *results, const vec3_t start,
									const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
									cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] ) {
	int i;
	qboolean model_rotated, trm_rotated;
	vec3_t invModelAxis[3], tmpAxis[3];
	vec3_t dir, tmp;
	cm_traceWork_t tw; // TODO: ALIGN16

	// fast point case
	if ( !trm || ( trm->bounds[1][0] - trm->bounds[0][0] <= 0.0f &&
					trm->bounds[1][1] - trm->bounds[0][1] <= 0.0f &&
					trm->bounds[1][2] - trm->bounds[0][2] <= 0.0f ) ) {

		results->c.contents = CM_TransformedPointContents( start, model, modelOrigin, modelAxis );
		results->fraction = ( results->c.contents == 0 );
        VectorCopy( start, results->endpos );
        AxisCopy( trmAxis, results->endAxis );

		return results->c.contents;
	}

	cmLocal.checkCount++;

	tw.trace.fraction = 1.0f;
	tw.trace.c.contents = 0;
	tw.trace.c.type = CONTACT_NONE;
	tw.contents = contentMask;
	tw.isConvex = qtrue;
	tw.rotation = qfalse;
	tw.positionTest = qtrue;
	tw.pointTrace = qfalse;
	tw.quickExit = qfalse;
	tw.numContacts = 0;
	tw.model = cmLocal.models[model];
    VectorSubtract( start, modelOrigin, tw.start );
    VectorCopy( tw.start, tw.end );

	model_rotated = AxisIsRotated( modelAxis );
	if ( model_rotated ) {
        TransposeAxis( modelAxis, invModelAxis );
	}

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
        MatrixRotateVector(  trm->offset, trmAxis, dir );
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
	}


	// setup trm vertices
	ClearBounds( tw.size[0], tw.size[1] );
	for ( i = 0; i < tw.numVerts; i++ ) {
		// get axial trm size after rotations
        VectorSubtract( tw.vertices[i].p, tw.start, tmp );
        AddPointToBounds( tmp, tw.size[0], tw.size[1] );
	}

	// setup trm edges
	for ( i = 1; i <= tw.numEdges; i++ ) {
		// edge start, end and pluecker coordinate
        VectorCopy( tw.vertices[tw.edges[i].vertexNum[0]].p, tw.edges[i].start );
        VectorCopy( tw.vertices[tw.edges[i].vertexNum[1]].p, tw.edges[i].end );
		PlueckerFromLine( tw.edges[i].pl, tw.edges[i].start, tw.edges[i].end );
	}

	// setup trm polygons
	if ( trm_rotated & model_rotated ) {
        MatrixMultiply( trmAxis, invModelAxis, tmpAxis );
		for ( i = 0; i < tw.numPolys; i++ ) {
            PlaneRotateSelf( tw.polys[i].plane, tmpAxis );
		}
	} else if ( trm_rotated ) {
		for ( i = 0; i < tw.numPolys; i++ ) {
            PlaneRotateSelf( tw.polys[i].plane, trmAxis );
		}
	} else if ( model_rotated ) {
		for ( i = 0; i < tw.numPolys; i++ ) {
            PlaneRotateSelf( tw.polys[i].plane, invModelAxis );
		}
	}
	for ( i = 0; i < tw.numPolys; i++ ) {
		PlaneFitThroughPoint( tw.polys[i].plane, tw.edges[abs(tw.polys[i].edges[0])].start );
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
		if ( fabs(tw.size[0][i]) > fabs(tw.size[1][i]) ) {
			tw.extents[i] = fabs( tw.size[0][i] ) + CM_BOX_EPSILON;
		} else {
			tw.extents[i] = fabs( tw.size[1][i] ) + CM_BOX_EPSILON;
		}
	}

	// trace through the model
	CM_TraceThroughModel( &tw );

	*results = tw.trace;
	results->fraction = ( results->c.contents == 0 );
    VectorCopy( start, results->endpos );
    AxisCopy( trmAxis, results->endAxis );

	return results->c.contents;
}

/*
==================
CM_Contents
==================
*/
int CM_Contents( const vec3_t start,
									const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
									cmHandle_t model, const vec3_t modelOrigin, const vec3_t modelAxis[3] ) {
	cm_trace_t results;

	if ( model < 0 || model > cmLocal.maxModels || model > MAX_SUBMODELS ) {
		ii.Com_Printf( PRINT_ALL, "CM_Contents: invalid model handle\n");
		return 0;
	}
	if ( !cmLocal.models || !cmLocal.models[model] ) {
		ii.Com_Printf( PRINT_ALL, "CM_Contents: invalid model\n");
		return 0;
	}

	return CM_ContentsTrm( &results, start, trm, trmAxis, contentMask, model, modelOrigin, modelAxis );
}
