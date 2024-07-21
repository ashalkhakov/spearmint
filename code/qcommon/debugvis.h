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

#ifndef __DEBUGVIS_H__
#define __DEBUGVIS_H__

#include "../idlib/q_shared.h"
#include "../idlib/q_extramath.h"

void	Com_DebugClearLines( int time );
void    Com_DebugLine( const vec4_t color, const vec3_t start, const vec3_t end, const int lifetime, const qboolean depthTest );
void	Com_DebugArrow( const vec4_t color, const vec3_t start, const vec3_t end, int size, const int lifetime );
void	Com_DebugWinding( const vec4_t color, const fixedWinding_t *w, const vec3_t origin, const vec3_t axis[3], const int lifetime, const qboolean depthTest );
void	Com_DebugCircle( const vec4_t color, const vec3_t origin, const vec3_t dir, const float radius, const int numSteps, const int lifetime, const qboolean depthTest );
void	Com_DebugBounds( const vec4_t color, const vec3_t bounds[2], const vec3_t org, const int lifetime );
void	Com_DebugAxis( const vec3_t origin, const vec3_t axis[3] );

void	Com_DebugClearPolygons( int time );
void	Com_DebugPolygon( const vec4_t color, const fixedWinding_t *winding, const int lifeTime, const qboolean depthTest );

#endif /* !__DEBUGVIS_H__ */
