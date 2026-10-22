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
// tr_fogvolume.c -- fog volumes drawn as Doom 3's fog lights: a fog brush
// whose shader says "fogLight" is fogged through its hull per pixel, by
// how much of each view ray is inside it (fogvolume_fp.glsl), instead of
// on the surfaces in its bounds as a Quake 3 fog volume. Daikatana's fog
// volumes are these (bsp-sekai).

#include "tr_local.h"

#define	MAX_HULL_POINTS		64

/*
================
ClipHullPolygon

the part of the polygon behind the plane (dot - dist <= 0)
================
*/
static int ClipHullPolygon( vec3_t *in, int numIn, const vec4_t plane, vec3_t *out ) {
	float	dists[MAX_HULL_POINTS + 1];
	int		i, numOut;

	for ( i = 0; i < numIn; i++ ) {
		dists[i] = DotProduct( in[i], plane ) - plane[3];
	}
	dists[numIn] = dists[0];

	numOut = 0;
	for ( i = 0; i < numIn && numOut < MAX_HULL_POINTS - 1; i++ ) {
		float *p1 = in[i], *p2 = in[( i + 1 ) % numIn];

		if ( dists[i] <= 0 ) {
			VectorCopy( p1, out[numOut] );
			numOut++;
		}
		if ( ( dists[i] < 0 && dists[i + 1] > 0 ) || ( dists[i] > 0 && dists[i + 1] < 0 ) ) {
			float frac = dists[i] / ( dists[i] - dists[i + 1] );

			out[numOut][0] = p1[0] + frac * ( p2[0] - p1[0] );
			out[numOut][1] = p1[1] + frac * ( p2[1] - p1[1] );
			out[numOut][2] = p1[2] + frac * ( p2[2] - p1[2] );
			numOut++;
		}
	}
	return numOut;
}

/*
================
R_FogLightHull

a fog light's hull from its brush's planes (each pointing out): a face on
each plane, clipped by the others, as triangles; its bounds from them
================
*/
void R_FogLightHull( fog_t *fog, const vec4_t *planes, int numPlanes ) {
	vec3_t	poly[2][MAX_HULL_POINTS];
	vec3_t	*verts;
	vec4_t	*faces;
	int		i, j, k, n, cur, numVerts;

	numPlanes = MIN( numPlanes, MAX_FOGLIGHT_PLANES );
	fog->fogLight = qtrue;
	fog->numPlanes = numPlanes;
	fog->planes = ri.Hunk_Alloc( numPlanes * sizeof( vec4_t ), h_low );
	Com_Memcpy( fog->planes, planes, numPlanes * sizeof( vec4_t ) );

	// at most a fan of MAX_HULL_POINTS - 2 triangles a face
	verts = ri.Hunk_AllocateTempMemory( numPlanes * ( MAX_HULL_POINTS - 2 ) * 3 * sizeof( vec3_t ) );
	faces = ri.Hunk_AllocateTempMemory( numPlanes * ( MAX_HULL_POINTS - 2 ) * sizeof( vec4_t ) );
	numVerts = 0;
	ClearBounds( fog->bounds[0], fog->bounds[1] );

	for ( i = 0; i < numPlanes; i++ ) {
		vec3_t	right, up, center;

		// a big square on the plane
		PerpendicularVector( right, planes[i] );
		CrossProduct( planes[i], right, up );
		VectorScale( planes[i], planes[i][3], center );
		VectorScale( right, MAX_WORLD_COORD * 2, right );
		VectorScale( up, MAX_WORLD_COORD * 2, up );
		VectorSubtract( center, right, poly[0][0] );
		VectorAdd( poly[0][0], up, poly[0][0] );
		VectorAdd( center, right, poly[0][1] );
		VectorAdd( poly[0][1], up, poly[0][1] );
		VectorAdd( center, right, poly[0][2] );
		VectorSubtract( poly[0][2], up, poly[0][2] );
		VectorSubtract( center, right, poly[0][3] );
		VectorSubtract( poly[0][3], up, poly[0][3] );
		n = 4;
		cur = 0;

		for ( j = 0; j < numPlanes && n >= 3; j++ ) {
			if ( j == i ) {
				continue;
			}
			n = ClipHullPolygon( poly[cur], n, planes[j], poly[cur ^ 1] );
			cur ^= 1;
		}
		if ( n < 3 ) {
			continue;		// a bevel: it touches the hull, it has no face
		}

		for ( k = 2; k < n; k++ ) {
			VectorCopy( poly[cur][0], verts[numVerts + 0] );
			VectorCopy( poly[cur][k - 1], verts[numVerts + 1] );
			VectorCopy( poly[cur][k], verts[numVerts + 2] );
			Vector4Copy( planes[i], faces[numVerts / 3] );
			numVerts += 3;
		}
		for ( k = 0; k < n; k++ ) {
			AddPointToBounds( poly[cur][k], fog->bounds[0], fog->bounds[1] );
		}
	}

	fog->numHullVerts = numVerts;
	fog->hullVerts = ri.Hunk_Alloc( MAX( numVerts, 1 ) * sizeof( vec3_t ), h_low );
	fog->hullPlanes = ri.Hunk_Alloc( MAX( numVerts / 3, 1 ) * sizeof( vec4_t ), h_low );
	Com_Memcpy( fog->hullVerts, verts, numVerts * sizeof( vec3_t ) );
	Com_Memcpy( fog->hullPlanes, faces, ( numVerts / 3 ) * sizeof( vec4_t ) );

	ri.Hunk_FreeTempMemory( faces );
	ri.Hunk_FreeTempMemory( verts );

	ri.Printf( PRINT_DEVELOPER, "fog light %s: %d planes, %d triangles, bounds (%.0f %.0f %.0f) (%.0f %.0f %.0f)\n",
		fog->shader ? fog->shader->name : "?", numPlanes, numVerts / 3,
		fog->bounds[0][0], fog->bounds[0][1], fog->bounds[0][2], fog->bounds[1][0], fog->bounds[1][1], fog->bounds[1][2] );
}

/*
================
RB_DrawFogVolumes

after a view's surfaces: each fog light in view, its faces that look back
at the eye (all of them from inside it), fogging what's behind them
================
*/
void RB_DrawFogVolumes( void ) {
	shaderProgram_t	*sp = &tr.fogVolumeShader;
	vec4_t			viewInfo;
	vec3_t			eye;
	int				i, j, numPlanes;
	qboolean		bound = qfalse;

	if ( !tr.world || !r_fogLights->integer || !tr.hdrDepthFbo || !r_depthPrepass->integer
		|| ( backEnd.refdef.rdflags & RDF_NOWORLDMODEL ) ) {
		return;		// (it needs the depth prepass's depth, as a float)
	}

	VectorCopy( backEnd.viewParms.or.origin, eye );
	numPlanes = ( backEnd.viewParms.flags & VPF_FARPLANEFRUSTUM ) ? 5 : 4;

	for ( i = 1; i < tr.world->numfogs; i++ ) {
		fog_t	*fog = &tr.world->fogs[i];

		if ( !fog->fogLight || !fog->numHullVerts ) {
			continue;
		}
		for ( j = 0; j < numPlanes; j++ ) {
			if ( BoxOnPlaneSide( fog->bounds[0], fog->bounds[1], &backEnd.viewParms.frustum[j] ) == 2 ) {
				break;
			}
		}
		if ( j < numPlanes ) {
			continue;		// out of view
		}

		// its faces seen from behind
		tess.numVertexes = 0;
		tess.numIndexes = 0;
		tess.firstIndex = 0;
		for ( j = 0; j < fog->numHullVerts && tess.numVertexes + 3 <= SHADER_MAX_VERTEXES; j += 3 ) {
			const float *face = fog->hullPlanes[j / 3];
			int k;

			if ( DotProduct( eye, face ) - face[3] >= 0 ) {
				continue;
			}
			for ( k = 0; k < 3; k++ ) {
				VectorCopy( fog->hullVerts[j + k], tess.xyz[tess.numVertexes] );
				tess.xyz[tess.numVertexes][3] = 1.0f;
				tess.indexes[tess.numIndexes++] = tess.numVertexes++;
			}
		}
		if ( !tess.numIndexes ) {
			continue;
		}

		if ( !bound ) {
			bound = qtrue;
			GL_SetModelviewMatrix( backEnd.viewParms.world.modelMatrix );
			GL_State( GLS_SRCBLEND_SRC_ALPHA | GLS_DSTBLEND_ONE_MINUS_SRC_ALPHA | GLS_DEPTHTEST_DISABLE );
			GL_Cull( CT_TWO_SIDED );

			GLSL_BindProgram( sp );
			GL_BindToTMU( tr.hdrDepthImage, TB_COLORMAP );

			VectorSet4( viewInfo, backEnd.viewParms.zFar / r_znear->value, backEnd.viewParms.zFar,
						1.0f / tr.hdrDepthImage->width, 1.0f / tr.hdrDepthImage->height );
			GLSL_SetUniformVec4( sp, UNIFORM_VIEWINFO, viewInfo );
			GLSL_SetUniformVec3( sp, UNIFORM_VIEWORIGIN, eye );
			GLSL_SetUniformVec3( sp, UNIFORM_VIEWFORWARD, backEnd.viewParms.or.axis[0] );
			GLSL_SetUniformMat4( sp, UNIFORM_MODELVIEWPROJECTIONMATRIX, glState.modelviewProjection );
		}

		GLSL_SetUniformVec4FogPlanes( sp, UNIFORM_FOGPLANES, (const vec4_t *)fog->planes, fog->numPlanes );
		GLSL_SetUniformInt( sp, UNIFORM_FOGNUMPLANES, fog->numPlanes );
		{
			vec4_t info, color;

			VectorSet4( info, fog->bounds[0][2], fog->bounds[1][2],
						MAX( 1.0f, fog->shader->fogParms.depthForOpaque ) / MAX( 0.01f, r_fogLightDensity->value ), 0 );
			GLSL_SetUniformVec4( sp, UNIFORM_FOGVOLUMEINFO, info );
			VectorCopy( fog->shader->fogParms.color, color );
			color[3] = 1.0f;
			GLSL_SetUniformVec4( sp, UNIFORM_COLOR, color );
		}

		RB_UpdateTessVao( ATTR_POSITION );
		R_DrawElements( tess.numIndexes, tess.firstIndex );

		tess.numIndexes = 0;
		tess.numVertexes = 0;
		tess.firstIndex = 0;
	}
}
