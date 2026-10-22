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
// tr_portals.c -- Doom 3's portal visibility on a Quake 3 BSP: the tree is
// cut into portals between its leaves when the map loads (as qbsp's
// MakeTreePortals), and each view flows from the eye's cluster through the
// portals between clusters (Doom 3's areas; leafs are too many and small),
// narrowed to each portal it goes through (as Doom 3's
// FlowViewThroughPortals, RenderWorld_portals.cpp), marking the clusters'
// leafs it reaches. A view through a portal is kept as a screen rectangle:
// a cluster is flowed into again only when that grows what of it is seen,
// so it ends without limits, and never sees less than it should. The PVS
// and area portals (doors) still prune it.

#include "tr_local.h"


#define MAX_PORTAL_POINTS	64
#define PORTAL_EYE_EPSILON	4.0f	// the eye this near a portal's plane looks through it
#define RECT_GRID			16.0f	// a view's rectangles snap out to 1/32 of the screen
#define MAX_FLOW_STEPS		65536	// portals gone through in a view at most (a safety)

typedef struct {
	int		numPoints;
	vec3_t	points[MAX_PORTAL_POINTS];
} pwinding_t;

typedef struct {
	float	x0, y0, x1, y1;		// on the screen, -1 to 1
} flowRect_t;

typedef struct {
	int			cluster;
	flowRect_t	rect;
} flowStep_t;

static world_t		*s_world;		// being cut
static rportal_t	*s_portals;		// growing while the tree is cut
static int			s_numPortals, s_maxPortals;
static vec3_t		*s_points;
static int			s_numPoints, s_maxPoints;

/*
================
ClipWinding

the part of the winding on the plane's front (dot - dist >= 0) or back
================
*/
static void ClipWinding( const pwinding_t *in, const vec3_t normal, float dist, qboolean front, pwinding_t *out ) {
	float	dists[MAX_PORTAL_POINTS + 1];
	int		i;

	for ( i = 0; i < in->numPoints; i++ ) {
		dists[i] = DotProduct( in->points[i], normal ) - dist;
		if ( !front ) {
			dists[i] = -dists[i];
		}
	}
	dists[in->numPoints] = dists[0];

	out->numPoints = 0;
	for ( i = 0; i < in->numPoints && out->numPoints < MAX_PORTAL_POINTS - 1; i++ ) {
		const float *p1 = in->points[i], *p2 = in->points[( i + 1 ) % in->numPoints];

		if ( dists[i] >= 0 ) {
			VectorCopy( p1, out->points[out->numPoints] );
			out->numPoints++;
		}
		if ( ( dists[i] < 0 && dists[i + 1] > 0 ) || ( dists[i] > 0 && dists[i + 1] < 0 ) ) {
			float	frac = dists[i] / ( dists[i] - dists[i + 1] );
			float	*o = out->points[out->numPoints];

			o[0] = p1[0] + frac * ( p2[0] - p1[0] );
			o[1] = p1[1] + frac * ( p2[1] - p1[1] );
			o[2] = p1[2] + frac * ( p2[2] - p1[2] );
			out->numPoints++;
		}
	}
	if ( out->numPoints < 3 ) {
		out->numPoints = 0;
	}
}

/*
================
BaseWinding

a square on the plane, as big as the world
================
*/
static void BaseWinding( const cplane_t *plane, const vec3_t mins, const vec3_t maxs, pwinding_t *w ) {
	vec3_t	right, up, center, mid;
	float	size;

	// centred on the world, as big as it: a smaller square clips with less
	// rounding than one the size of MAX_WORLD_COORD
	VectorAdd( mins, maxs, mid );
	VectorScale( mid, 0.5f, mid );
	size = Distance( mins, maxs ) + 64;
	VectorMA( mid, plane->dist - DotProduct( mid, plane->normal ), plane->normal, center );
	PerpendicularVector( right, plane->normal );
	CrossProduct( plane->normal, right, up );
	VectorScale( right, size, right );
	VectorScale( up, size, up );
	VectorSubtract( center, right, w->points[0] );
	VectorAdd( w->points[0], up, w->points[0] );
	VectorAdd( center, right, w->points[1] );
	VectorAdd( w->points[1], up, w->points[1] );
	VectorAdd( center, right, w->points[2] );
	VectorSubtract( w->points[2], up, w->points[2] );
	VectorSubtract( center, right, w->points[3] );
	VectorSubtract( w->points[3], up, w->points[3] );
	w->numPoints = 4;
}

/*
================
AddPortal
================
*/
static void AddPortal( const pwinding_t *w, const cplane_t *plane, int front, int back ) {
	rportal_t *p;

	if ( s_numPortals == s_maxPortals ) {
		rportal_t *grown;

		s_maxPortals = MAX( 1024, s_maxPortals * 2 );
		grown = ri.Malloc( s_maxPortals * sizeof( *grown ) );
		if ( s_numPortals ) {
			Com_Memcpy( grown, s_portals, s_numPortals * sizeof( *grown ) );
			ri.Free( s_portals );
		}
		s_portals = grown;
	}
	while ( s_numPoints + w->numPoints > s_maxPoints ) {
		vec3_t *grown;

		s_maxPoints = MAX( 8192, s_maxPoints * 2 );
		grown = ri.Malloc( s_maxPoints * sizeof( *grown ) );
		if ( s_numPoints ) {
			Com_Memcpy( grown, s_points, s_numPoints * sizeof( *grown ) );
			ri.Free( s_points );
		}
		s_points = grown;
	}

	p = &s_portals[s_numPortals++];
	p->firstPoint = s_numPoints;
	p->numPoints = w->numPoints;
	Com_Memcpy( s_points + s_numPoints, w->points, w->numPoints * sizeof( vec3_t ) );
	s_numPoints += w->numPoints;
	VectorCopy( plane->normal, p->plane );
	p->plane[3] = plane->dist;
	p->leafs[0] = front;
	p->leafs[1] = back;
}

/*
================
FilterBack

a node's winding, already down to a leaf of its front side, down its back
side: a portal between the two leaves for each piece that ends in one
================
*/
static void FilterBack( mnode_t *node, const pwinding_t *w, const cplane_t *plane, int front ) {
	pwinding_t part;

	if ( node->isLeaf ) {
		if ( node->cluster >= 0 ) {
			AddPortal( w, plane, front, node - s_world->nodes );
		}
		return;
	}
	ClipWinding( w, node->plane->normal, node->plane->dist, qtrue, &part );
	if ( part.numPoints ) {
		FilterBack( node->children[0], &part, plane, front );
	}
	ClipWinding( w, node->plane->normal, node->plane->dist, qfalse, &part );
	if ( part.numPoints ) {
		FilterBack( node->children[1], &part, plane, front );
	}
}

/*
================
FilterFront

a node's winding down its front side, each piece then down its back side
================
*/
static void FilterFront( mnode_t *node, const pwinding_t *w, mnode_t *split ) {
	pwinding_t part;

	if ( node->isLeaf ) {
		if ( node->cluster >= 0 ) {
			FilterBack( split->children[1], w, split->plane, node - s_world->nodes );
		}
		return;
	}
	ClipWinding( w, node->plane->normal, node->plane->dist, qtrue, &part );
	if ( part.numPoints ) {
		FilterFront( node->children[0], &part, split );
	}
	ClipWinding( w, node->plane->normal, node->plane->dist, qfalse, &part );
	if ( part.numPoints ) {
		FilterFront( node->children[1], &part, split );
	}
}

/*
================
R_PortalizeWorld

the portals between the world's leaves: each node's plane, cut to its cell
(by its parents' planes and the world's bounds), split by both its sides'
subtrees down to their leaves
================
*/
void R_PortalizeWorld( world_t *world ) {
	mnode_t		*node;
	pwinding_t	w, clipped;
	int			i, j, *counts;

	world->numPortals = 0;
	if ( !world->numDecisionNodes ) {
		return;
	}
	s_world = world;
	s_portals = NULL;
	s_numPortals = s_maxPortals = 0;
	s_points = NULL;
	s_numPoints = s_maxPoints = 0;

	for ( i = 0; i < world->numDecisionNodes; i++ ) {
		mnode_t *child, *parent;

		node = &world->nodes[i];
		BaseWinding( node->plane, world->nodes[0].mins, world->nodes[0].maxs, &w );

		// within its cell: its parents' sides it's on
		for ( child = node, parent = node->parent; parent && w.numPoints; child = parent, parent = parent->parent ) {
			ClipWinding( &w, parent->plane->normal, parent->plane->dist, parent->children[0] == child, &clipped );
			w = clipped;
		}
		// and the world's bounds
		for ( j = 0; j < 3 && w.numPoints; j++ ) {
			vec3_t normal = { 0, 0, 0 };

			normal[j] = 1;
			ClipWinding( &w, normal, world->nodes[0].mins[j] - 8, qtrue, &clipped );
			w = clipped;
			if ( w.numPoints ) {
				ClipWinding( &w, normal, world->nodes[0].maxs[j] + 8, qfalse, &clipped );
				w = clipped;
			}
		}
		if ( w.numPoints ) {
			FilterFront( node->children[0], &w, node );
		}
	}

	// into the hunk, with each leaf's list of its portals
	world->numPortals = s_numPortals;
	world->portals = ri.Hunk_Alloc( MAX( 1, s_numPortals ) * sizeof( rportal_t ), h_low );
	world->portalPoints = ri.Hunk_Alloc( MAX( 1, s_numPoints ) * sizeof( vec3_t ), h_low );
	if ( s_numPortals ) {
		Com_Memcpy( world->portals, s_portals, s_numPortals * sizeof( rportal_t ) );
		Com_Memcpy( world->portalPoints, s_points, s_numPoints * sizeof( vec3_t ) );
	}
	world->leafFirstPortal = ri.Hunk_Alloc( world->numnodes * sizeof( int ), h_low );
	world->leafNumPortals = ri.Hunk_Alloc( world->numnodes * sizeof( int ), h_low );
	world->leafPortals = ri.Hunk_Alloc( MAX( 1, s_numPortals * 2 ) * sizeof( int ), h_low );
	counts = ri.Hunk_AllocateTempMemory( world->numnodes * sizeof( int ) );
	Com_Memset( counts, 0, world->numnodes * sizeof( int ) );
	for ( i = 0; i < s_numPortals; i++ ) {
		world->leafNumPortals[s_portals[i].leafs[0]]++;
		world->leafNumPortals[s_portals[i].leafs[1]]++;
	}
	for ( i = 0, j = 0; i < world->numnodes; i++ ) {
		world->leafFirstPortal[i] = j;
		j += world->leafNumPortals[i];
	}
	for ( i = 0; i < s_numPortals; i++ ) {
		for ( j = 0; j < 2; j++ ) {
			int leaf = s_portals[i].leafs[j];

			world->leafPortals[world->leafFirstPortal[leaf] + counts[leaf]++] = i;
		}
	}
	ri.Hunk_FreeTempMemory( counts );

	// the portals between clusters, each cluster's, and each cluster's leafs
	world->clusterFirstPortal = ri.Hunk_Alloc( MAX( 1, world->numClusters ) * sizeof( int ), h_low );
	world->clusterNumPortals = ri.Hunk_Alloc( MAX( 1, world->numClusters ) * sizeof( int ), h_low );
	world->clusterPortals = ri.Hunk_Alloc( MAX( 1, s_numPortals * 2 ) * sizeof( int ), h_low );
	world->clusterFirstLeaf = ri.Hunk_Alloc( MAX( 1, world->numClusters ) * sizeof( int ), h_low );
	world->clusterNumLeafs = ri.Hunk_Alloc( MAX( 1, world->numClusters ) * sizeof( int ), h_low );
	world->clusterLeafs = ri.Hunk_Alloc( MAX( 1, world->numnodes ) * sizeof( int ), h_low );
	counts = ri.Hunk_AllocateTempMemory( MAX( 1, world->numClusters ) * sizeof( int ) );
	for ( i = 0; i < s_numPortals; i++ ) {
		int c0 = world->nodes[s_portals[i].leafs[0]].cluster, c1 = world->nodes[s_portals[i].leafs[1]].cluster;

		if ( c0 != c1 && c0 < world->numClusters && c1 < world->numClusters ) {
			world->clusterNumPortals[c0]++;
			world->clusterNumPortals[c1]++;
		}
	}
	for ( i = world->numDecisionNodes; i < world->numnodes; i++ ) {
		if ( world->nodes[i].cluster >= 0 && world->nodes[i].cluster < world->numClusters ) {
			world->clusterNumLeafs[world->nodes[i].cluster]++;
		}
	}
	{
		int p = 0, l = 0;

		for ( i = 0; i < world->numClusters; i++ ) {
			world->clusterFirstPortal[i] = p;
			p += world->clusterNumPortals[i];
			world->clusterFirstLeaf[i] = l;
			l += world->clusterNumLeafs[i];
		}
	}
	Com_Memset( counts, 0, MAX( 1, world->numClusters ) * sizeof( int ) );
	for ( i = 0; i < s_numPortals; i++ ) {
		int c0 = world->nodes[s_portals[i].leafs[0]].cluster, c1 = world->nodes[s_portals[i].leafs[1]].cluster;

		if ( c0 != c1 && c0 < world->numClusters && c1 < world->numClusters ) {
			world->clusterPortals[world->clusterFirstPortal[c0] + counts[c0]++] = i;
			world->clusterPortals[world->clusterFirstPortal[c1] + counts[c1]++] = i;
		}
	}
	Com_Memset( counts, 0, MAX( 1, world->numClusters ) * sizeof( int ) );
	for ( i = world->numDecisionNodes; i < world->numnodes; i++ ) {
		int c = world->nodes[i].cluster;

		if ( c >= 0 && c < world->numClusters ) {
			world->clusterLeafs[world->clusterFirstLeaf[c] + counts[c]++] = i;
		}
	}
	ri.Hunk_FreeTempMemory( counts );

	// the flow's scratch
	world->maxFlowSteps = world->numClusters * 4 + 64;
	world->flowRects = ri.Hunk_Alloc( MAX( 1, world->numClusters ) * sizeof( flowRect_t ), h_low );
	world->flowSeen = ri.Hunk_Alloc( MAX( 1, world->numClusters ), h_low );
	world->flowSteps = ri.Hunk_Alloc( world->maxFlowSteps * sizeof( flowStep_t ), h_low );

	if ( s_portals ) {
		ri.Free( s_portals );
	}
	if ( s_points ) {
		ri.Free( s_points );
	}
	ri.Printf( PRINT_DEVELOPER, "%d portals between %d leafs\n", world->numPortals, world->numnodes - world->numDecisionNodes );
}

/*
============================================================

FLOW

============================================================
*/

static const byte	*s_flowPVS;
static int			s_flowLeafs;		// marked this view, for r_portalVis 2

/*
================
MarkLeaf

as R_MarkLeaves marks a leaf in view: it and its parents
================
*/
static void MarkLeaf( mnode_t *leaf ) {
	mnode_t *parent = leaf;

	if ( leaf->visCounts[tr.visIndex] != tr.visCounts[tr.visIndex] ) {
		s_flowLeafs++;
	}
	do {
		if ( parent->visCounts[tr.visIndex] == tr.visCounts[tr.visIndex] ) {
			break;
		}
		parent->visCounts[tr.visIndex] = tr.visCounts[tr.visIndex];
		parent = parent->parent;
	} while ( parent );
}


static flowRect_t	*s_clusterRects;	// what of each cluster has been seen, this view
static byte			*s_clusterSeen;
static flowStep_t	*s_steps;			// the list, a ring
static int			s_maxSteps;
static float		s_tanX, s_tanY;

/*
================
ProjectWinding

the screen rectangle of a winding in front of the eye
================
*/
static void ProjectWinding( const pwinding_t *w, flowRect_t *r ) {
	int i;

	r->x0 = r->y0 = 1e9f;
	r->x1 = r->y1 = -1e9f;
	for ( i = 0; i < w->numPoints; i++ ) {
		vec3_t	v;
		float	z, x, y;

		VectorSubtract( w->points[i], tr.viewParms.or.origin, v );
		z = MAX( DotProduct( v, tr.viewParms.or.axis[0] ), 0.5f );
		x = -DotProduct( v, tr.viewParms.or.axis[1] ) / ( z * s_tanX );
		y = DotProduct( v, tr.viewParms.or.axis[2] ) / ( z * s_tanY );
		r->x0 = MIN( r->x0, x );
		r->x1 = MAX( r->x1, x );
		r->y0 = MIN( r->y0, y );
		r->y1 = MAX( r->y1, y );
	}
}

/*
================
RectPlanes

the four sides of the frustum through a screen rectangle, facing in
================
*/
static void RectPlanes( const flowRect_t *r, vec4_t planes[4] ) {
	const float	*eye = tr.viewParms.or.origin;
	float		xs[2] = { r->x0, r->x1 }, ys[2] = { r->y0, r->y1 };
	vec3_t		corner[4], center;
	int			i;

	// the corners' directions: screen x right (against axis[1], left), y up
	for ( i = 0; i < 4; i++ ) {
		float x = xs[( i == 1 || i == 2 ) ? 1 : 0], y = ys[( i >= 2 ) ? 1 : 0];

		VectorCopy( tr.viewParms.or.axis[0], corner[i] );
		VectorMA( corner[i], -x * s_tanX, tr.viewParms.or.axis[1], corner[i] );
		VectorMA( corner[i], y * s_tanY, tr.viewParms.or.axis[2], corner[i] );
	}
	VectorCopy( tr.viewParms.or.axis[0], center );
	VectorMA( center, -0.5f * ( r->x0 + r->x1 ) * s_tanX, tr.viewParms.or.axis[1], center );
	VectorMA( center, 0.5f * ( r->y0 + r->y1 ) * s_tanY, tr.viewParms.or.axis[2], center );
	for ( i = 0; i < 4; i++ ) {
		vec3_t normal;

		CrossProduct( corner[i], corner[( i + 1 ) & 3], normal );
		VectorNormalize( normal );
		if ( DotProduct( normal, center ) < 0 ) {
			VectorNegate( normal, normal );
		}
		VectorCopy( normal, planes[i] );
		planes[i][3] = DotProduct( normal, eye );
	}
}

/*
================
MarkCluster

all a cluster's leafs, but those behind closed doors
================
*/
static void MarkCluster( int cluster ) {
	world_t	*w = tr.world;
	int		i;

	for ( i = 0; i < w->clusterNumLeafs[cluster]; i++ ) {
		mnode_t *leaf = &w->nodes[w->clusterLeafs[w->clusterFirstLeaf[cluster] + i]];

		if ( !( tr.refdef.areamask[leaf->area >> 3] & ( 1 << ( leaf->area & 7 ) ) ) ) {
			MarkLeaf( leaf );
		}
	}
}

/*
================
FlowClusters

from the eye's cluster, seen through the whole screen: each cluster taken
off the list flows through its portals to other clusters, the portal clipped
to the part of the screen it's seen through and in front of the eye; the
other cluster goes on the list when what of it is seen grows
================
*/
static void FlowClusters( int eyeCluster, const flowRect_t *screen ) {
	world_t	*w = tr.world;
	int		numSteps = 0, head = 0, total = 0;

	s_steps[numSteps].cluster = eyeCluster;
	s_steps[numSteps].rect = *screen;
	numSteps++;
	s_clusterRects[eyeCluster] = *screen;
	s_clusterSeen[eyeCluster] = 1;
	MarkCluster( eyeCluster );

	while ( head != numSteps && total < MAX_FLOW_STEPS ) {
		flowStep_t	step = s_steps[head];
		vec4_t		planes[4];
		int			i, j;

		head = ( head + 1 ) % s_maxSteps;
		RectPlanes( &step.rect, planes );

		for ( i = 0; i < w->clusterNumPortals[step.cluster]; i++ ) {
			rportal_t	*p = &w->portals[w->clusterPortals[w->clusterFirstPortal[step.cluster] + i]];
			qboolean	fromFront = w->nodes[p->leafs[0]].cluster == step.cluster;
			mnode_t		*otherLeaf = &w->nodes[fromFront ? p->leafs[1] : p->leafs[0]];
			int			other = otherLeaf->cluster;
			float		eyeDist = DotProduct( tr.viewParms.or.origin, p->plane ) - p->plane[3];
			pwinding_t	win, clipped;
			flowRect_t	r, *seen;

			total++;
			// a sightline crosses a portal only from the eye's side of it
			if ( fromFront ? eyeDist < -PORTAL_EYE_EPSILON : eyeDist > PORTAL_EYE_EPSILON ) {
				continue;
			}
			// the PVS and doors still say what can't be seen
			if ( s_flowPVS && !( s_flowPVS[other >> 3] & ( 1 << ( other & 7 ) ) ) ) {
				continue;
			}
			if ( tr.refdef.areamask[otherLeaf->area >> 3] & ( 1 << ( otherLeaf->area & 7 ) ) ) {
				continue;
			}

			// the eye (nearly) in the portal's plane sees it edge on, as a line
			// on the screen, but looks right through it: all it sees goes on
			if ( fabs( eyeDist ) < PORTAL_EYE_EPSILON ) {
				r = step.rect;
				goto through;
			}

			// the part of the portal in front of the eye, in what's seen
			win.numPoints = p->numPoints;
			Com_Memcpy( win.points, w->portalPoints + p->firstPoint, p->numPoints * sizeof( vec3_t ) );
			ClipWinding( &win, tr.viewParms.or.axis[0], DotProduct( tr.viewParms.or.axis[0], tr.viewParms.or.origin ) + 1.0f, qtrue, &clipped );
			win = clipped;
			for ( j = 0; j < 4 && win.numPoints; j++ ) {
				ClipWinding( &win, planes[j], planes[j][3], qtrue, &clipped );
				win = clipped;
			}
			if ( !win.numPoints ) {
				continue;
			}
			ProjectWinding( &win, &r );
			r.x0 = MAX( r.x0, step.rect.x0 );
			r.y0 = MAX( r.y0, step.rect.y0 );
			r.x1 = MIN( r.x1, step.rect.x1 );
			r.y1 = MIN( r.y1, step.rect.y1 );
			if ( r.x0 >= r.x1 || r.y0 >= r.y1 ) {
				continue;
			}
			// outward to a coarse grid: slivers of growth don't send a cluster
			// round again
			r.x0 = MAX( -1.0f, floor( r.x0 * RECT_GRID ) / RECT_GRID );
			r.y0 = MAX( -1.0f, floor( r.y0 * RECT_GRID ) / RECT_GRID );
			r.x1 = MIN( 1.0f, ceil( r.x1 * RECT_GRID ) / RECT_GRID );
			r.y1 = MIN( 1.0f, ceil( r.y1 * RECT_GRID ) / RECT_GRID );

through:
			// on only if it shows more of the other cluster
			seen = &s_clusterRects[other];
			if ( s_clusterSeen[other] ) {
				if ( r.x0 >= seen->x0 && r.y0 >= seen->y0 && r.x1 <= seen->x1 && r.y1 <= seen->y1 ) {
					continue;
				}
				seen->x0 = MIN( seen->x0, r.x0 );
				seen->y0 = MIN( seen->y0, r.y0 );
				seen->x1 = MAX( seen->x1, r.x1 );
				seen->y1 = MAX( seen->y1, r.y1 );
				r = *seen;
			} else {
				*seen = r;
				s_clusterSeen[other] = 1;
				MarkCluster( other );
			}
			if ( ( numSteps + 1 ) % s_maxSteps == head ) {
				continue;		// (the list is full: it has been seen anyway)
			}
			s_steps[numSteps].cluster = other;
			s_steps[numSteps].rect = r;
			numSteps = ( numSteps + 1 ) % s_maxSteps;
		}
	}
}

/*
================
R_FlowView

the leaves the view sees through the portals from the eye's leaf, marked
as R_MarkLeaves marks the PVS's; qfalse if it can't (no portals, the eye in
something solid), and R_MarkLeaves does it its way
================
*/
qboolean R_FlowView( mnode_t *eyeLeaf, const byte *pvs ) {
	flowRect_t	screen = { -1, -1, 1, 1 };
	int			i;

	if ( !r_portalVis->integer || !tr.world->numPortals || !eyeLeaf || eyeLeaf->cluster < 0
		|| eyeLeaf->cluster >= tr.world->numClusters ) {
		return qfalse;
	}

	s_tanX = tan( DEG2RAD( tr.viewParms.fovX * 0.5f ) );
	s_tanY = tan( DEG2RAD( tr.viewParms.fovY * 0.5f ) );
	s_clusterRects = tr.world->flowRects;
	s_clusterSeen = tr.world->flowSeen;
	s_steps = tr.world->flowSteps;
	s_maxSteps = tr.world->maxFlowSteps;

	// a new marking every view: it depends on where the eye looks
	tr.visIndex = ( tr.visIndex + 1 ) % MAX_VISCOUNTS;
	tr.visCounts[tr.visIndex]++;
	tr.visClusters[tr.visIndex] = -2;

	Com_Memset( s_clusterSeen, 0, tr.world->numClusters );
	s_flowPVS = pvs;
	s_flowLeafs = 0;
	FlowClusters( eyeLeaf->cluster, &screen );

	if ( r_portalVis->integer == 2 ) {
		int pvsLeafs = 0;

		for ( i = tr.world->numDecisionNodes; i < tr.world->numnodes; i++ ) {
			int cluster = tr.world->nodes[i].cluster;

			if ( cluster >= 0 && ( !pvs || ( pvs[cluster >> 3] & ( 1 << ( cluster & 7 ) ) ) ) ) {
				pvsLeafs++;
			}
		}
		ri.Printf( PRINT_ALL, "portal flow: %d leafs (the PVS has %d)\n", s_flowLeafs, pvsLeafs );
	}
	return qtrue;
}
