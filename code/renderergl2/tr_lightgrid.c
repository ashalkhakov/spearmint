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
// tr_lightgrid.c -- light grids a room, after RBDOOM-3-BFG's (one an area):
// a lattice over the world, every LATTICE_SPACING units, each open point
// labelled with its room, so a model lit from it takes no light from the
// other side of a wall.
//
// The rooms are a watershed over voxels of the open space: each voxel's
// depth is its distance from the nearest wall, and the voxels flood in from
// the deepest down. Each room fills from its middle; where two floods meet
// at a narrow place (a doorway, shallow beside the rooms' middles) is the
// border between them. Floods that meet nearly as deep as the shallower one
// (r_roomSaddle of it) are one room: a wide opening, or a flat floor.
//
// The light at each point is baked (bakeLightGrid) from a cube of the
// lightmaps seen from there, untextured: what falls on the room around it.
// The ambient is their average; the light comes from where it lands from, so
// a lit floor below is light from above, and lit walls all round cancel out.
// Saved as maps/<map>.lightgrid, with the rooms.

#include "tr_local.h"

#define LATTICE_SPACING		64.0f
#define VOXEL_SIZE			32.0f		// LATTICE_SPACING over a whole number
#define LATTICE_OFFSET		31.0f		// off the planes, which are mostly on even units

#define LIGHTGRID_IDENT		( ( 'G' << 24 ) + ( 'L' << 16 ) + ( 'K' << 8 ) + 'D' )
#define LIGHTGRID_VERSION	1

/*
================
LatticeOpen

whether a point is in the open (a leaf with a cluster)
================
*/
static qboolean LatticeOpen( world_t *w, const vec3_t p ) {
	mnode_t *node = w->nodes;

	while ( !node->isLeaf ) {
		node = node->children[DotProduct( p, node->plane->normal ) - node->plane->dist > 0 ? 0 : 1];
	}
	return node->cluster >= 0;
}

/*
================
LatticeClear

whether a line is in the open all along (every 8 units)
================
*/
static qboolean LatticeClear( world_t *w, const vec3_t a, const vec3_t b ) {
	vec3_t	d, p;
	float	len;
	int		i, steps;

	VectorSubtract( b, a, d );
	len = VectorLength( d );
	steps = MAX( 1, (int)( len / 8.0f ) );
	for ( i = 1; i <= steps; i++ ) {
		VectorMA( a, (float)i / steps, d, p );
		if ( !LatticeOpen( w, p ) ) {
			return qfalse;
		}
	}
	return qtrue;
}

static void LatticePointOf( const world_t *w, int index, vec3_t p ) {
	int x = index % w->latticeBounds[0];
	int y = ( index / w->latticeBounds[0] ) % w->latticeBounds[1];
	int z = index / ( w->latticeBounds[0] * w->latticeBounds[1] );

	p[0] = w->latticeOrigin[0] + x * LATTICE_SPACING;
	p[1] = w->latticeOrigin[1] + y * LATTICE_SPACING;
	p[2] = w->latticeOrigin[2] + z * LATTICE_SPACING;
}

/*
================
BakeAlloc

the bake's scratch memory, which outgrows the zone on big maps (millions of
voxels): the system's, a drop if there's none
================
*/
static void *BakeAlloc( size_t size ) {
	void *p = malloc( size );

	if ( !p ) {
		ri.Error( ERR_DROP, "bakeLightGrid: out of memory for %lu bytes", (unsigned long)size );
	}
	return p;
}

static int RoomRoot( int *parent, int i ) {
	while ( parent[i] != i ) {
		parent[i] = parent[parent[i]];
		i = parent[i];
	}
	return i;
}

/*
================
R_BuildRoomLattice

the lattice and its rooms (into the lattice already there, if any: the
world's doesn't change)
================
*/
static void R_BuildRoomLattice( world_t *w ) {
	int		vb[3], step[3], n, i, j, k, x, y, z, head, tail, maxDepth, numRooms;
	byte	*edges;		// bit k: open to the neighbour +k (k 0..2)
	byte	*depth;		// 255 solid
	int		*queue, *parent, *peak, *roomOf;
	float	saddle = r_roomSaddle->value;
	int		start = ri.Milliseconds();

	for ( i = 0; i < 3; i++ ) {
		w->latticeOrigin[i] = LATTICE_SPACING * floor( w->nodes[0].mins[i] / LATTICE_SPACING ) + LATTICE_OFFSET - LATTICE_SPACING;
		w->latticeBounds[i] = (int)ceil( ( w->nodes[0].maxs[i] - w->latticeOrigin[i] ) / LATTICE_SPACING ) + 1;
		vb[i] = ( w->latticeBounds[i] - 1 ) * (int)( LATTICE_SPACING / VOXEL_SIZE ) + 1;
	}
	step[0] = 1;
	step[1] = vb[0];
	step[2] = vb[0] * vb[1];
	n = vb[0] * vb[1] * vb[2];
	edges = BakeAlloc( n );
	depth = BakeAlloc( n );
	queue = BakeAlloc( n * sizeof( int ) );
	parent = BakeAlloc( n * sizeof( int ) );

	// the voxels in the open, and which neighbours each reaches
	for ( i = 0, z = 0; z < vb[2]; z++ ) {
		for ( y = 0; y < vb[1]; y++ ) {
			for ( x = 0; x < vb[0]; x++, i++ ) {
				vec3_t p;

				p[0] = w->latticeOrigin[0] + x * VOXEL_SIZE;
				p[1] = w->latticeOrigin[1] + y * VOXEL_SIZE;
				p[2] = w->latticeOrigin[2] + z * VOXEL_SIZE;
				depth[i] = LatticeOpen( w, p ) ? 254 : 255;
				edges[i] = 0;
			}
		}
	}
	for ( i = 0, z = 0; z < vb[2]; z++ ) {
		for ( y = 0; y < vb[1]; y++ ) {
			for ( x = 0; x < vb[0]; x++, i++ ) {
				int		c[3] = { x, y, z };
				vec3_t	p, q;

				if ( depth[i] == 255 ) {
					continue;
				}
				p[0] = w->latticeOrigin[0] + x * VOXEL_SIZE;
				p[1] = w->latticeOrigin[1] + y * VOXEL_SIZE;
				p[2] = w->latticeOrigin[2] + z * VOXEL_SIZE;
				for ( k = 0; k < 3; k++ ) {
					if ( c[k] + 1 >= vb[k] || depth[i + step[k]] == 255 ) {
						continue;
					}
					VectorCopy( p, q );
					q[k] += VOXEL_SIZE;
					if ( LatticeClear( w, p, q ) ) {
						edges[i] |= 1 << k;
					}
				}
			}
		}
	}

	// depth: steps from the nearest voxel beside a wall, outward from all of them
#define	VOXEL_LINKED(i, k, ni) ( (k) & 1 ? ( (ni) = (i) - step[(k) >> 1], ( (i) / step[(k) >> 1] ) % vb[(k) >> 1] > 0 && ( edges[ni] & ( 1 << ( (k) >> 1 ) ) ) ) \
									 : ( (ni) = (i) + step[(k) >> 1], ( edges[i] & ( 1 << ( (k) >> 1 ) ) ) ) )
	head = tail = 0;
	for ( i = 0; i < n; i++ ) {
		int ni, links = 0;

		if ( depth[i] == 255 ) {
			continue;
		}
		for ( k = 0; k < 6; k++ ) {
			if ( VOXEL_LINKED( i, k, ni ) ) {
				links++;
			}
		}
		if ( links < 6 ) {
			depth[i] = 0;
			queue[tail++] = i;
		}
	}
	maxDepth = 0;
	while ( head < tail ) {
		int c = queue[head++], ni;

		for ( k = 0; k < 6; k++ ) {
			if ( VOXEL_LINKED( c, k, ni ) && depth[ni] == 254 ) {
				depth[ni] = depth[c] + 1;
				maxDepth = MAX( maxDepth, depth[ni] );
				queue[tail++] = ni;
			}
		}
	}

	// the flood, deepest first (the queue holds the voxels shallowest first)
	for ( i = 0; i < n; i++ ) {
		parent[i] = -1;
	}
	peak = BakeAlloc( n * sizeof( int ) );		// a basin's deepest, at its root
	for ( j = tail - 1; j >= 0; j-- ) {
		int c = queue[j], ni, d = depth[c];

		parent[c] = c;
		peak[c] = d;
		for ( k = 0; k < 6; k++ ) {
			int a, b;

			if ( !VOXEL_LINKED( c, k, ni ) || parent[ni] < 0 ) {
				continue;
			}
			a = RoomRoot( parent, c );
			b = RoomRoot( parent, ni );
			if ( a == b ) {
				continue;
			}
			if ( a == c && peak[a] == d ) {
				parent[a] = b;		// a voxel new to the flood joins it
				continue;
			}
			// two floods meet here: one room if this is nearly as deep as the shallower
			if ( d + 1 >= saddle * ( MIN( peak[a], peak[b] ) + 1 ) ) {
				if ( peak[a] > peak[b] ) {
					parent[b] = a;
				} else {
					parent[a] = b;
				}
			}
		}
	}
#undef VOXEL_LINKED

	// the rooms, numbered; the lattice's points are voxels too
	roomOf = peak;
	for ( i = 0; i < n; i++ ) {
		roomOf[i] = -1;
	}
	numRooms = 0;
	for ( i = 0; i < n; i++ ) {
		if ( parent[i] >= 0 && RoomRoot( parent, i ) == i ) {
			roomOf[i] = numRooms++;
		}
	}
	w->numRooms = numRooms;
	if ( !w->latticeRoom ) {
		w->latticeRoom = ri.Hunk_Alloc( w->latticeBounds[0] * w->latticeBounds[1] * w->latticeBounds[2] * sizeof( short ), h_low );
	}
	j = (int)( LATTICE_SPACING / VOXEL_SIZE );
	for ( i = 0, z = 0; z < w->latticeBounds[2]; z++ ) {
		for ( y = 0; y < w->latticeBounds[1]; y++ ) {
			for ( x = 0; x < w->latticeBounds[0]; x++, i++ ) {
				int v = x * j * step[0] + y * j * step[1] + z * j * step[2];

				w->latticeRoom[i] = parent[v] < 0 ? -1 : roomOf[RoomRoot( parent, v )];
			}
		}
	}

	// developer: the voxels' rooms, for looking at (maps/<map>.rooms: the
	// size, three ints, then a short a voxel, -1 solid)
	if ( ri.Cvar_VariableIntegerValue( "developer" ) ) {
		int		size = 12 + n * 2;
		byte	*buf = BakeAlloc( size );
		short	*out = (short *)( buf + 12 );
		int		*sizes = BakeAlloc( MAX( 1, numRooms ) * sizeof( int ) );
		int		big = 0, best = 0;
		char	name[MAX_QPATH];

		Com_Memcpy( buf, vb, 12 );
		Com_Memset( sizes, 0, MAX( 1, numRooms ) * sizeof( int ) );
		for ( i = 0; i < n; i++ ) {
			out[i] = parent[i] < 0 ? -1 : roomOf[RoomRoot( parent, i )];
			if ( out[i] >= 0 ) {
				sizes[out[i]]++;
			}
		}
		COM_StripExtension( w->name, name, sizeof( name ) );
		Q_strcat( name, sizeof( name ), ".rooms" );
		ri.FS_WriteFile( name, buf, size );
		for ( i = 0; i < numRooms; i++ ) {
			big += sizes[i] >= 64;
			best = MAX( best, sizes[i] );
		}
		ri.Printf( PRINT_ALL, "%d rooms in %dx%dx%d voxels, %d of 64 voxels or more, the biggest %d, the deepest %d, %d ms\n",
			numRooms, vb[0], vb[1], vb[2], big, best, maxDepth, ri.Milliseconds() - start );
		free( sizes );
		free( buf );
	}

	free( peak );
	free( parent );
	free( queue );
	free( depth );
	free( edges );
}

/*
================
R_LoadRoomLattice

the rooms and their light, baked, if there are (maps/<map>.lightgrid: the
ident and version, the lattice's origin, bounds and rooms, a short a point
for its room, -1 solid, then the light, 8 bytes an open point)
================
*/
void R_LoadRoomLattice( world_t *w ) {
	char		name[MAX_QPATH];
	byte		*buf, *p;
	int			size, n, i, numOpen;
	const int	*header;

	w->numRooms = 0;
	w->latticeRoom = NULL;
	w->latticeLight = NULL;

	Com_sprintf( name, sizeof( name ), "maps/%s.lightgrid", w->baseName );
	size = ri.FS_ReadFile( name, (void **)&buf );
	if ( !buf ) {
		return;
	}
	header = (const int *)buf;
	if ( size < 32 || LittleLong( header[0] ) != LIGHTGRID_IDENT || LittleLong( header[1] ) != LIGHTGRID_VERSION ) {
		ri.Printf( PRINT_WARNING, "WARNING: %s is not a light grid of this version\n", name );
		ri.FS_FreeFile( buf );
		return;
	}
	for ( i = 0; i < 3; i++ ) {
		w->latticeOrigin[i] = LittleFloat( ( (const float *)header )[2 + i] );
		w->latticeBounds[i] = LittleLong( header[5 + i] );
	}
	w->numRooms = LittleLong( header[8] );
	n = w->latticeBounds[0] * w->latticeBounds[1] * w->latticeBounds[2];
	if ( n <= 0 || size < 36 + n * 2 ) {
		ri.Printf( PRINT_WARNING, "WARNING: %s is cut short\n", name );
		w->numRooms = 0;
		ri.FS_FreeFile( buf );
		return;
	}
	w->latticeRoom = ri.Hunk_Alloc( n * sizeof( short ), h_low );
	w->latticeLight = ri.Hunk_Alloc( n * 8, h_low );
	p = buf + 36;
	for ( i = 0, numOpen = 0; i < n; i++, p += 2 ) {
		w->latticeRoom[i] = LittleShort( *(short *)p );
		numOpen += w->latticeRoom[i] >= 0;
	}
	if ( size < 36 + n * 2 + numOpen * 8 ) {
		ri.Printf( PRINT_WARNING, "WARNING: %s is cut short\n", name );
		w->numRooms = 0;
		w->latticeRoom = NULL;
		w->latticeLight = NULL;
		ri.FS_FreeFile( buf );
		return;
	}
	for ( i = 0; i < n; i++ ) {
		if ( w->latticeRoom[i] >= 0 ) {
			Com_Memcpy( w->latticeLight + i * 8, p, 8 );
			p += 8;
		}
	}
	ri.FS_FreeFile( buf );
	ri.Printf( PRINT_DEVELOPER, "%s: %d rooms, %d points lit\n", name, w->numRooms, numOpen );
}

/*
================
R_SetupEntityLightingRooms

an entity's light from the points of its room around it: the room of the
nearest it sees. qfalse when it sees none (in a wall, or off the lattice)
================
*/
qboolean R_SetupEntityLightingRooms( trRefEntity_t *ent, world_t *w ) {
	vec3_t	org, v, direction;
	int		pos[3], corner[8], i, j, room = -1;
	float	frac[3], dist[8], totalFactor = 0;
	byte	tried = 0;

	if ( ent->e.renderfx & RF_LIGHTING_ORIGIN ) {
		VectorCopy( ent->e.lightingOrigin, org );
	} else {
		VectorCopy( ent->e.origin, org );
	}
	for ( i = 0; i < 3; i++ ) {
		float f = ( org[i] - w->latticeOrigin[i] ) / LATTICE_SPACING;

		pos[i] = floor( f );
		if ( pos[i] < 0 || pos[i] >= w->latticeBounds[i] - 1 ) {
			return qfalse;
		}
		frac[i] = f - pos[i];
	}
	for ( i = 0; i < 8; i++ ) {
		corner[i] = ( pos[0] + ( i & 1 ) ) + ( pos[1] + ( ( i >> 1 ) & 1 ) ) * w->latticeBounds[0]
			+ ( pos[2] + ( i >> 2 ) ) * w->latticeBounds[0] * w->latticeBounds[1];
		LatticePointOf( w, corner[i], v );
		dist[i] = DistanceSquared( v, org );
	}

	// the room: the nearest open corner seen from here
	while ( room < 0 ) {
		int best = -1;

		for ( i = 0; i < 8; i++ ) {
			if ( !( tried & ( 1 << i ) ) && w->latticeRoom[corner[i]] >= 0 && ( best < 0 || dist[i] < dist[best] ) ) {
				best = i;
			}
		}
		if ( best < 0 ) {
			return qfalse;
		}
		tried |= 1 << best;
		LatticePointOf( w, corner[best], v );
		if ( LatticeClear( w, org, v ) ) {
			room = w->latticeRoom[corner[best]];
		}
	}

	VectorClear( ent->ambientLight );
	VectorClear( ent->directedLight );
	VectorClear( direction );
	for ( i = 0; i < 8; i++ ) {
		const byte	*data = w->latticeLight + corner[i] * 8;
		float		factor = 1.0f, lat, lng;

		if ( w->latticeRoom[corner[i]] != room ) {
			continue;
		}
		for ( j = 0; j < 3; j++ ) {
			factor *= ( i >> j ) & 1 ? frac[j] : 1.0f - frac[j];
		}
		totalFactor += factor;
		for ( j = 0; j < 3; j++ ) {
			// 128 a unit of light, 255 a unit for the shaders
			ent->ambientLight[j] += factor * data[j] * ( 255.0f / 128.0f );
			ent->directedLight[j] += factor * data[3 + j] * ( 255.0f / 128.0f );
		}
		lat = data[7] * ( 2.0f * M_PI / 256.0f );
		lng = data[6] * ( 2.0f * M_PI / 256.0f );
		v[0] = cos( lat ) * sin( lng );
		v[1] = sin( lat ) * sin( lng );
		v[2] = cos( lng );
		VectorMA( direction, factor, v, direction );
	}
	if ( totalFactor <= 0 ) {
		return qfalse;
	}
	VectorScale( ent->ambientLight, 1.0f / totalFactor, ent->ambientLight );
	VectorScale( ent->directedLight, 1.0f / totalFactor, ent->directedLight );
	if ( VectorNormalize2( direction, ent->lightDir ) == 0 ) {
		VectorSet( ent->lightDir, 0, 0, 1 );
	}
	return qtrue;
}

/*
================
BakePoint

a point's light from its cube: rgba floats, the faces side by side
================
*/
static void BakePoint( const float *pixels, int row, const vec3_t axes[6][3], byte *out ) {
	const int	S = LIGHTGRID_BAKE_FACE, W = LIGHTGRID_BAKE_FACE * 6;
	vec3_t		mean = { 0, 0, 0 }, flow[3], d;
	float		total = 0, lat, lng;
	int			f, x, y, c;

	Com_Memset( flow, 0, sizeof( flow ) );
	for ( f = 0; f < 6; f++ ) {
		for ( y = 0; y < S; y++ ) {
			const float *px = pixels + ( ( row * S + y ) * W + f * S ) * 4;
			float sy = ( y + 0.5f ) / S * 2.0f - 1.0f;

			for ( x = 0; x < S; x++, px += 4 ) {
				float	sx = ( x + 0.5f ) / S * 2.0f - 1.0f;
				float	r2 = 1.0f + sx * sx + sy * sy;
				float	weight = 1.0f / ( r2 * sqrt( r2 ) );	// the pixel's solid angle
				vec3_t	dir;

				// the way it looks (screen right is the face's -left), then
				// the way the light came: the other way
				VectorScale( axes[f][0], 1.0f, dir );
				VectorMA( dir, -sx, axes[f][1], dir );
				VectorMA( dir, sy, axes[f][2], dir );
				VectorNormalize( dir );
				for ( c = 0; c < 3; c++ ) {
					mean[c] += weight * px[c];
					VectorMA( flow[c], -weight * px[c], dir, flow[c] );
				}
				total += weight;
			}
		}
	}
	VectorScale( mean, 1.0f / total, mean );
	for ( c = 0; c < 3; c++ ) {
		VectorScale( flow[c], 1.0f / total, flow[c] );
	}

	// the light's way, by brightness; then ambient + directed * max(0, n.d)
	// with the same average, and as much of it coming that way
	VectorScale( flow[0], 0.299f, d );
	VectorMA( d, 0.587f, flow[1], d );
	VectorMA( d, 0.114f, flow[2], d );
	if ( VectorNormalize( d ) == 0 ) {
		VectorSet( d, 0, 0, 1 );
	}
	for ( c = 0; c < 3; c++ ) {
		float directed = 2.0f * MAX( 0.0f, DotProduct( flow[c], d ) );
		float ambient = MAX( 0.0f, mean[c] - directed * 0.5f );

		out[c] = (byte)Com_Clamp( 0, 255, ambient * 128.0f + 0.5f );
		out[3 + c] = (byte)Com_Clamp( 0, 255, directed * 128.0f + 0.5f );
	}
	lat = atan2( d[1], d[0] );
	if ( lat < 0 ) {
		lat += 2.0f * M_PI;
	}
	lng = acos( Com_Clamp( -1, 1, d[2] ) );
	out[6] = (byte)( (int)( lng * 256.0f / ( 2.0f * M_PI ) + 0.5f ) & 255 );
	out[7] = (byte)( (int)( lat * 256.0f / ( 2.0f * M_PI ) + 0.5f ) & 255 );
}

/*
================
R_BakeLightGrid_f

bakeLightGrid: the rooms and their light, for the map that's loaded
================
*/
void R_BakeLightGrid_f( void ) {
	// faces: forward, left, up (left = up x forward)
	static const vec3_t axes[6][3] = {
		{ {  1, 0, 0 }, { 0,  1, 0 }, { 0, 0, 1 } },
		{ { -1, 0, 0 }, { 0, -1, 0 }, { 0, 0, 1 } },
		{ { 0,  1, 0 }, { -1, 0, 0 }, { 0, 0, 1 } },
		{ { 0, -1, 0 }, {  1, 0, 0 }, { 0, 0, 1 } },
		{ { 0, 0,  1 }, { 0,  1, 0 }, { -1, 0, 0 } },
		{ { 0, 0, -1 }, { 0,  1, 0 }, {  1, 0, 0 } },
	};
	world_t	*w = tr.world;
	int		n, i, j, f, numOpen, *open, size, start;
	int		savedLightmap, savedPortalVis, savedFogLights;
	float	*pixels;
	byte	*buf, *p;
	char	name[MAX_QPATH];

	if ( !w || !w->numDecisionNodes ) {
		ri.Printf( PRINT_ALL, "bakeLightGrid: no map loaded\n" );
		return;
	}
	if ( !tr.lightGridBakeFbo ) {
		ri.Printf( PRINT_ALL, "bakeLightGrid: needs framebuffer objects\n" );
		return;
	}
	start = ri.Milliseconds();

	R_BuildRoomLattice( w );
	n = w->latticeBounds[0] * w->latticeBounds[1] * w->latticeBounds[2];
	if ( !w->latticeLight ) {
		w->latticeLight = ri.Hunk_Alloc( n * 8, h_low );
	}
	Com_Memset( w->latticeLight, 0, n * 8 );
	open = BakeAlloc( n * sizeof( int ) );
	for ( i = 0, numOpen = 0; i < n; i++ ) {
		if ( w->latticeRoom[i] >= 0 ) {
			open[numOpen++] = i;
		}
	}
	ri.Printf( PRINT_ALL, "bakeLightGrid: %d rooms, %d points\n", w->numRooms, numOpen );

	// the lightmaps alone, from everywhere: no portal flow, no lit fog
	savedLightmap = r_lightmap->integer;
	savedPortalVis = r_portalVis->integer;
	savedFogLights = r_fogLights->integer;
	r_lightmap->integer = 1;
	r_portalVis->integer = 0;
	r_fogLights->integer = 0;

	pixels = BakeAlloc( tr.lightGridBakeFbo->width * tr.lightGridBakeFbo->height * 4 * sizeof( float ) );
	for ( i = 0; i < numOpen; i += LIGHTGRID_BAKE_BATCH ) {
		refdef_t	fd;
		int			batch = MIN( LIGHTGRID_BAKE_BATCH, numOpen - i );

		Com_Memset( &fd, 0, sizeof( fd ) );
		fd.width = tr.lightGridBakeFbo->width;
		fd.height = tr.lightGridBakeFbo->height;
		fd.fov_x = fd.fov_y = 90;
		LatticePointOf( w, open[i], fd.vieworg );
		AxisCopy( axisDefault, fd.viewaxis );

		RE_ClearScene();
		RE_BeginScene( &fd );
		for ( j = 0; j < batch; j++ ) {
			vec3_t org;

			LatticePointOf( w, open[i + j], org );
			for ( f = 0; f < 6; f++ ) {
				viewParms_t parms;

				Com_Memset( &parms, 0, sizeof( parms ) );
				parms.viewportX = f * LIGHTGRID_BAKE_FACE;
				parms.viewportY = j * LIGHTGRID_BAKE_FACE;
				parms.viewportWidth = parms.viewportHeight = LIGHTGRID_BAKE_FACE;
				parms.targetFbo = tr.lightGridBakeFbo;
				parms.flags = VPF_NOVIEWMODEL | VPF_NOCUBEMAPS;
				parms.fovX = parms.fovY = 90;
				VectorCopy( org, parms.or.origin );
				VectorCopy( org, parms.pvsOrigin );
				VectorCopy( axes[f][0], parms.or.axis[0] );
				VectorCopy( axes[f][1], parms.or.axis[1] );
				VectorCopy( axes[f][2], parms.or.axis[2] );
				R_RenderView( &parms );
			}
		}
		RE_EndScene();
		R_IssuePendingRenderCommands();
		R_InitNextFrame();

		FBO_Bind( tr.lightGridBakeFbo );
		qglReadPixels( 0, 0, tr.lightGridBakeFbo->width, tr.lightGridBakeFbo->height, GL_RGBA, GL_FLOAT, pixels );
		FBO_Bind( NULL );
		for ( j = 0; j < batch; j++ ) {
			BakePoint( pixels, j, axes, w->latticeLight + open[i + j] * 8 );
		}
		if ( ( i / LIGHTGRID_BAKE_BATCH ) % 2000 == 1999 ) {
			ri.Printf( PRINT_ALL, "bakeLightGrid: %d of %d\n", i + batch, numOpen );
		}
	}
	free( pixels );
	r_lightmap->integer = savedLightmap;
	r_portalVis->integer = savedPortalVis;
	r_fogLights->integer = savedFogLights;

	// saved
	size = 36 + n * 2 + numOpen * 8;
	buf = BakeAlloc( size );
	( (int *)buf )[0] = LittleLong( LIGHTGRID_IDENT );
	( (int *)buf )[1] = LittleLong( LIGHTGRID_VERSION );
	for ( i = 0; i < 3; i++ ) {
		( (float *)buf )[2 + i] = LittleFloat( w->latticeOrigin[i] );
		( (int *)buf )[5 + i] = LittleLong( w->latticeBounds[i] );
	}
	( (int *)buf )[8] = LittleLong( w->numRooms );
	p = buf + 36;
	for ( i = 0; i < n; i++, p += 2 ) {
		*(short *)p = LittleShort( w->latticeRoom[i] );
	}
	for ( i = 0; i < numOpen; i++, p += 8 ) {
		Com_Memcpy( p, w->latticeLight + open[i] * 8, 8 );
	}
	Com_sprintf( name, sizeof( name ), "maps/%s.lightgrid", w->baseName );
	ri.FS_WriteFile( name, buf, size );
	free( buf );
	free( open );
	ri.Printf( PRINT_ALL, "bakeLightGrid: %s, %d points in %.1f s\n", name, numOpen, ( ri.Milliseconds() - start ) / 1000.0f );
}
