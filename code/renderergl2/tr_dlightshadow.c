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
// tr_dlightshadow.c -- dynamic lights' shadow maps, after The Dark Mod's
// atlas (ShadowMapStage.cpp, tr_light.cpp's atlas pages): the lights that
// matter most to a view get places in one depth texture, six faces for a
// point light and one for a spot light, sized by how much they matter, each
// face rendered as a depth-only view from the light (so the portal flow
// culls its casters as it does the eye's view, tr_portals.c). The light
// passes (ForwardDlight, lightall_fp.glsl) pick the face and soften it as
// The Dark Mod does (PCSS: wider the farther from the caster); the fog
// lights (fogvolume_fp.glsl) read it too, for shafts in the fog.

#include "tr_local.h"

// a point light's faces, in the order the shader picks them: -X +X -Y +Y -Z +Z
static const float cubeAxes[6][3][3] = {		// forward, left (up x forward), up
	{ { -1,  0,  0 }, {  0, -1,  0 }, {  0,  0,  1 } },
	{ {  1,  0,  0 }, {  0,  1,  0 }, {  0,  0,  1 } },
	{ {  0, -1,  0 }, {  1,  0,  0 }, {  0,  0,  1 } },
	{ {  0,  1,  0 }, { -1,  0,  0 }, {  0,  0,  1 } },
	{ {  0,  0, -1 }, {  0,  1,  0 }, {  1,  0,  0 } },
	{ {  0,  0,  1 }, {  0, -1,  0 }, {  1,  0,  0 } },
};

/*
================
R_DlightShadowScore

how much a light's shadow matters to the view: bright, big and near
================
*/
static float R_DlightShadowScore( const dlight_t *dl, const vec3_t eye ) {
	float dist = MAX( Distance( dl->origin, eye ) - dl->radius, 1.0f );

	return dl->intensity * MAX( dl->color[0], MAX( dl->color[1], dl->color[2] ) ) * dl->radius / dist;
}

/*
================
R_RenderDlightShadows

the shadow maps of the r_dlightShadows lights that matter most to the view,
into the atlas: its first one a quarter of it a face, the next ones smaller
================
*/
void R_RenderDlightShadows( const refdef_t *fd ) {
	int			order[MAX_DLIGHTS];
	float		score[MAX_DLIGHTS];
	int			i, j, n, numLights, atlas;
	int			cursorX = 0, cursorY = 0, rowHeight = 0;

	tr.refdef.numDlightShadows = 0;
	for ( i = 0; i < tr.refdef.num_dlights; i++ ) {
		tr.refdef.dlights[i].shadow = -1;
	}
	if ( r_dlightShadows->integer <= 0 || !tr.dlightShadowAtlasFbo || tr.dlightShadowAtlasImage->width < 512
		|| !tr.world || ( fd->rdflags & RDF_NOWORLDMODEL ) ) {
		return;
	}
	atlas = tr.dlightShadowAtlasImage->width;

	// the lights that can have one, best first
	for ( i = 0, numLights = 0; i < tr.refdef.num_dlights; i++ ) {
		dlight_t *dl = &tr.refdef.dlights[i];

		if ( ( dl->flags & ( REF_DIRECTED_DLIGHT | REF_VERTEX_DLIGHT ) ) || !( dl->flags & REF_SURFACE_DLIGHT ) || dl->radius < 32 ) {
			continue;
		}
		score[numLights] = R_DlightShadowScore( dl, fd->vieworg );
		order[numLights++] = i;
	}
	for ( i = 0; i < numLights; i++ ) {
		for ( j = i + 1; j < numLights; j++ ) {
			if ( score[j] > score[i] ) {
				float	s = score[i];
				int		o = order[i];

				score[i] = score[j]; order[i] = order[j];
				score[j] = s; order[j] = o;
			}
		}
	}
	numLights = MIN( numLights, MIN( r_dlightShadows->integer, MAX_DLIGHT_SHADOWS ) );

	for ( n = 0; n < numLights; n++ ) {
		dlight_t	*dl = &tr.refdef.dlights[order[n]];
		int			numFaces = ( dl->flags & REF_SPOT_DLIGHT ) ? 1 : 6;
		int			size = atlas / ( n == 0 ? 4 : n < 3 ? 8 : 16 );
		int			shadowNum = tr.refdef.numDlightShadows;

		// room for its faces, row by row
		for ( j = 0; j < numFaces; j++ ) {
			if ( cursorX + size > atlas ) {
				cursorX = 0;
				cursorY += rowHeight;
				rowHeight = 0;
			}
			if ( cursorY + size > atlas ) {
				break;
			}
			cursorX += size;
			rowHeight = MAX( rowHeight, size );
		}
		if ( j < numFaces ) {
			break;		// the atlas is full
		}

		tr.refdef.dlightShadows[shadowNum].numFaces = numFaces;
		tr.refdef.numDlightShadows++;
		dl->shadow = shadowNum;
	}

	// render them: the same placing again, now face by face
	cursorX = cursorY = rowHeight = 0;
	for ( n = 0; n < tr.refdef.numDlightShadows; n++ ) {
		dlight_t	*dl = NULL;
		int			size = atlas / ( n == 0 ? 4 : n < 3 ? 8 : 16 );

		for ( i = 0; i < tr.refdef.num_dlights; i++ ) {
			if ( tr.refdef.dlights[i].shadow == n ) {
				dl = &tr.refdef.dlights[i];
				break;
			}
		}
		if ( !dl ) {
			continue;
		}

		for ( j = 0; j < tr.refdef.dlightShadows[n].numFaces; j++ ) {
			viewParms_t	parms;
			float		*tile = tr.refdef.dlightShadows[n].tiles[j];

			if ( cursorX + size > atlas ) {
				cursorX = 0;
				cursorY += rowHeight;
				rowHeight = 0;
			}

			Com_Memset( &parms, 0, sizeof( parms ) );
			parms.viewportX = cursorX;
			parms.viewportY = cursorY;
			parms.viewportWidth = size;
			parms.viewportHeight = size;
			parms.targetFbo = tr.dlightShadowAtlasFbo;
			parms.flags = VPF_DEPTHSHADOW | VPF_NOVIEWMODEL | VPF_FIXEDFAR;
			parms.zFar = dl->radius;
			VectorCopy( dl->origin, parms.or.origin );
			VectorCopy( dl->origin, parms.pvsOrigin );

			if ( dl->flags & REF_SPOT_DLIGHT ) {
				float angle = RAD2DEG( acos( Com_Clamp( -1.0f, 1.0f, dl->spotCos ) ) );

				parms.fovX = parms.fovY = MIN( 170.0f, angle * 2.0f + 4.0f );
				VectorCopy( dl->dir, parms.or.axis[0] );
				PerpendicularVector( parms.or.axis[2], parms.or.axis[0] );
				CrossProduct( parms.or.axis[2], parms.or.axis[0], parms.or.axis[1] );
			} else {
				parms.fovX = parms.fovY = 90;
				VectorCopy( cubeAxes[j][0], parms.or.axis[0] );
				VectorCopy( cubeAxes[j][1], parms.or.axis[1] );
				VectorCopy( cubeAxes[j][2], parms.or.axis[2] );
			}

			R_RenderView( &parms );
			tr.refdef.dlightShadows[n].zNear = r_znear->value;
			tr.refdef.dlightShadows[n].zFar = dl->radius;
			tr.refdef.dlightShadows[n].tanHalfFov = tan( DEG2RAD( parms.fovX * 0.5f ) );

			// world to this face's clip space, and its place in the atlas
			Mat4Multiply( tr.viewParms.projectionMatrix, tr.viewParms.world.modelMatrix, tr.refdef.dlightShadows[n].mvp[j] );
			tile[0] = (float)cursorX / atlas;
			tile[1] = (float)cursorY / atlas;
			tile[2] = (float)size / atlas;
			tile[3] = (float)size / atlas;

			cursorX += size;
			rowHeight = MAX( rowHeight, size );
		}
	}
}
