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

#include "debugvis.h"
#include "../renderercommon/tr_public.h"

// Structure containing functions exported from refresh DLL
extern refexport_t	re;

/*
============
Com_DebugClearLines
============
*/
void Com_DebugClearLines( int time ) {
    if ( !re.DebugClearLines ) {
        return;
    }

    re.DebugClearLines( time );
}

/*
============
Com_DebugLine
============
*/
void Com_DebugLine( const vec4_t color, const vec3_t start, const vec3_t end, const int lifetime, const qboolean depthTest ) {
    if ( !re.DebugLine ) {
        return;
    }

    re.DebugLine( color, start, end, lifetime, depthTest );
}

/*
============
Com_DebugArrow
============
*/
void Com_DebugArrow( const vec4_t color, const vec3_t start, const vec3_t end, int size, const int lifetime ) {
    if ( !re.DebugArrow ) {
        return;
    }
    
    re.DebugArrow( color, start, end, size, lifetime );
}

/*
============
Com_DebugWinding
============
*/
void Com_DebugWinding( const vec4_t color, const fixedWinding_t *w, const vec3_t origin, const vec3_t axis[3], const int lifetime, const qboolean depthTest ) {
    if ( !re.DebugWinding ) {
        return;
    }

    re.DebugWinding( color, w, origin, axis, lifetime, depthTest );
}

/*
============
Com_DebugCircle
============
*/
void Com_DebugCircle( const vec4_t color, const vec3_t origin, const vec3_t dir, const float radius, const int numSteps, const int lifetime, const qboolean depthTest ) {
    if ( !re.DebugCircle ) {
        return;
    }

    re.DebugCircle( color, origin, dir, radius, numSteps, lifetime, depthTest );
}

/*
============
Com_DebugBounds
============
*/
void Com_DebugBounds( const vec4_t color, const vec3_t bounds[2], const vec3_t org, const int lifetime ) {
    if ( !re.DebugBounds ) {
        return;
    }

    re.DebugBounds( color, bounds, org, lifetime );
}

/*
============
Com_DebugAxis
============
*/
void Com_DebugAxis( const vec3_t origin, const vec3_t axis[3] ) {
    if ( !re.DebugAxis ) {
        return;
    }

    re.DebugAxis( origin, axis );
}

/*
============
Com_DebugClearPolygons
============
*/
void Com_DebugClearPolygons( int time ) {
    if ( !re.DebugClearPolygons ) {
        return;
    }

    re.DebugClearPolygons( time );
}

/*
============
Com_DebugPolygon
============
*/
void Com_DebugPolygon( const vec4_t color, const fixedWinding_t *winding, const int lifeTime, const qboolean depthTest ) {
    if ( !re.DebugPolygon ) {
        return;
    }

    re.DebugPolygon( color, winding, lifeTime, depthTest );
}
