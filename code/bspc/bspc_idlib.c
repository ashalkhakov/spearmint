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

// bspc_idlib.c -- what idlib, the BSP loaders and the old collision code bspc
// links expect from the engine: idlib's import table (ii) over bspc's
// allocator and stdio, and the zone

#include <stdarg.h>
#include <stdio.h>
#include "../idlib/q_shared.h"
#include "l_cmd.h"
#include "l_mem.h"
#include "../idlib/idlib_public.h"
#undef GetMemory
#undef FreeMemory

static void *BspcIdGetMemory( int bytes, char *label, char *file, int line ) { return GetMemory( bytes ); }
static void BspcIdFreeMemory( void *buf, char *label, char *file, int line ) { FreeMemory( buf ); }
static void *BspcIdHunkAlloc( int size ) { return GetClearedMemory( size ); }
static void BspcIdFreeFile( void *buf ) { FreeMemory( buf ); }

static void QDECL BspcIdPrintf( const char *fmt, ... ) {
	va_list ap;
	va_start( ap, fmt );
	vprintf( fmt, ap );
	va_end( ap );
}

static void QDECL __attribute__((noreturn)) BspcIdError( int level, const char *fmt, ... ) {
	char text[1024];
	va_list ap;
	va_start( ap, fmt );
	vsnprintf( text, sizeof( text ), fmt, ap );
	va_end( ap );
	Error( "%s", text );
}

// files by stdio, handles 1..
#define BSPC_MAX_FILES 16
static FILE *bspcFiles[BSPC_MAX_FILES];

static int BspcIdFOpenFile( const char *qpath, fileHandle_t *f, fsMode_t mode ) {
	int i;
	long len;

	for ( i = 1; i < BSPC_MAX_FILES && bspcFiles[i]; i++ ) {
	}
	*f = 0;
	if ( i == BSPC_MAX_FILES || !( bspcFiles[i] = fopen( qpath, mode == FS_READ ? "rb" : "wb" ) ) ) {
		return -1;
	}
	*f = i;
	if ( mode != FS_READ ) {
		return 0;
	}
	fseek( bspcFiles[i], 0, SEEK_END );
	len = ftell( bspcFiles[i] );
	fseek( bspcFiles[i], 0, SEEK_SET );
	return len;
}

static int BspcIdRead( void *buffer, int len, fileHandle_t f ) { return fread( buffer, 1, len, bspcFiles[f] ); }
static int BspcIdWrite( const void *buffer, int len, fileHandle_t f ) { return fwrite( buffer, 1, len, bspcFiles[f] ); }
static void BspcIdFCloseFile( fileHandle_t f ) { if ( bspcFiles[f] ) fclose( bspcFiles[f] ); bspcFiles[f] = NULL; }

idlib_import_t ii = {
	.Com_Error = BspcIdError,
	.Com_Printf = BspcIdPrintf,
	.Com_DPrintf = BspcIdPrintf,
	.GetMemoryDebug = BspcIdGetMemory,
	.FreeMemoryDebug = BspcIdFreeMemory,
	.HunkAlloc = BspcIdHunkAlloc,
	.FS_FOpenFile = BspcIdFOpenFile,
	.FS_Read = BspcIdRead,
	.FS_Write = BspcIdWrite,
	.FS_FCloseFile = BspcIdFCloseFile,
	.FS_FreeFile = BspcIdFreeFile,
};

// what the engine's zone gives the BSP loaders and the old collision code
void *Z_MallocDebug( int size, char *label, char *file, int line ) { return GetClearedMemory( size ); }
void Z_FreeDebug( void *ptr, char *label, char *file, int line ) { FreeMemory( ptr ); }

// l_struct.c writes structures through it; bspc never does
int FS_WriteFloatString( fileHandle_t fs, const char *fmt, ... ) { return 0; }
