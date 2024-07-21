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

Retrieving contacts

===============================================================================
*/

/*
==================
CM_Contacts
==================
*/
int CM_Contacts( contactInfo_t *contacts, const int maxContacts, const vec3_t start, const vec6_t dir, const float depth,
								const traceModel_t *trm, const vec3_t trmAxis[3], int contentMask,
								cmHandle_t model, const vec3_t origin, const vec3_t modelAxis[3] ) {
	cm_trace_t results;
	vec3_t end;

	// same as Translation but instead of storing the first collision we store all collisions as contacts
	cmLocal.getContacts = qtrue;
	cmLocal.contacts = contacts;
	cmLocal.maxContacts = maxContacts;
	cmLocal.numContacts = 0;
    VectorMA( start, depth, dir, end );
	CM_Translation( &results, start, end, trm, trmAxis, contentMask, model, origin, modelAxis );
	if ( VectorLengthSquared( dir + 3 ) != 0.0f ) {
		// FIXME: rotational contacts
	}
	cmLocal.getContacts = qfalse;
	cmLocal.maxContacts = 0;

	return cmLocal.numContacts;
}
