#include "cm_local.h"

cmLocal_t cmLocal;
const bspFile_t	*cm_bsp = NULL;
cmimport_t cmi;

cvar_t *cm_drawMask;
cvar_t *cm_drawColor;
cvar_t *cm_drawFilled;
cvar_t *cm_drawInternal;
cvar_t *cm_drawNormals;
cvar_t *cm_backFaceCull;
cvar_t *cm_debugCollision;

cvar_t *cm_testCollision;
cvar_t *cm_testRotation;
cvar_t *cm_testModel;
cvar_t *cm_testTimes;
cvar_t *cm_testRandomMany;
cvar_t *cm_testOrigin;
cvar_t *cm_testReset;
cvar_t *cm_testBox;
cvar_t *cm_testBoxRotation;
cvar_t *cm_testWalk;
cvar_t *cm_testLength;
cvar_t *cm_testRadius;
cvar_t *cm_testAngle;

cvar_t *cm_noAreas;
cvar_t *cm_noCurves;
cvar_t *cm_playerCurveClip;


#ifdef USE_COLLISIONMODEL_DLOPEN
Q_EXPORT cmexport_t* QDECL GetCollisionModelAPI ( int apiVersion, cmimport_t *cmimp ) {
#else
cmexport_t *GetCollisionModelAPI ( int apiVersion, cmimport_t *cmimp ) {
#endif

	static cmexport_t	cme;

	cmi = *cmimp;

	Com_Memset( &cme, 0, sizeof( cme ) );

	if ( apiVersion != CM_API_VERSION ) {
		ii.Com_Printf( "Mismatched CM_API_VERSION: expected %i, got %i\n", 
			CM_API_VERSION, apiVersion );
		return NULL;
	}

	cme.Init = CM_Init;
	cme.Shutdown = CM_Shutdown;

	cme.LoadMap = CM_LoadMap;
	cme.LoadBSP = CM_LoadBSP;
	cme.FreeMap = CM_FreeMap;

	cme.LoadModel = CM_LoadModel;
	cme.SetupTrmModel = CM_SetupTrmModel;
	cme.TrmFromModel = CM_TrmFromModel;

	cme.GetModelName = CM_GetModelName;
	cme.GetModelBounds = CM_GetModelBounds;
	cme.GetModelContents = CM_GetModelContents;
	cme.GetModelVertex = CM_GetModelVertex;
	cme.GetModelEdge = CM_GetModelEdge;
	cme.GetModelPolygon = CM_GetModelPolygon;

	cme.Translation = CM_Translation;
	cme.Rotation = CM_Rotation;
	cme.Contents = CM_Contents;
	cme.Contacts = CM_Contacts;

	cme.DebugOutput = CM_DebugOutput;
	cme.DrawModel = CM_DrawModel;
	cme.ModelInfo = CM_ModelInfo;
	cme.ListModels = CM_ListModels;
	cme.WriteCollisionModelForMapEntity = CM_WriteCollisionModelForMapEntity;

    cme.CM_RegisterMaterial = CM_RegisterMaterial;
    cme.CM_GetMaterialName = CM_GetMaterialName;
    cme.CM_GetMaterialContentFlags = CM_GetMaterialContentFlags;
    cme.CM_GetMaterialSurfaceFlags = CM_GetMaterialSurfaceFlags;

	cme.NumClusters = CM_NumClusters;
	cme.NumInlineModels = CM_NumInlineModels;
	cme.EntityString = CM_EntityString;
	cme.GetEntityToken = CM_GetEntityToken;

	cme.ClusterPVS = CM_ClusterPVS;

	cme.PointLeafnum = CM_PointLeafnum;

	cme.BoxLeafnums = CM_BoxLeafnums;

	cme.LeafCluster = CM_LeafCluster;
	cme.LeafArea = CM_LeafArea;

	cme.AdjustAreaPortalState = CM_AdjustAreaPortalState;
	cme.AreasConnected = CM_AreasConnected;

	cme.WriteAreaBits = CM_WriteAreaBits;

	return &cme;
}

void CM_Init( void ) {

    cm_drawMask = cmi.Cvar_Get( "cm_drawMask",			"none",		CVAR_ARCHIVE );
    cm_drawColor = cmi.Cvar_Get( "cm_drawColor",			"1 0 0 .5",	CVAR_ARCHIVE );
    cm_drawFilled = cmi.Cvar_Get( "cm_drawFilled",		"0",		CVAR_ARCHIVE );
    cm_drawInternal = cmi.Cvar_Get( "cm_drawInternal",		"1",		CVAR_ARCHIVE );
    cm_drawNormals = cmi.Cvar_Get(	"cm_drawNormals",		"0",		CVAR_ARCHIVE );
    cm_backFaceCull = cmi.Cvar_Get( "cm_backFaceCull",		"0",		CVAR_ARCHIVE );
    cm_debugCollision = cmi.Cvar_Get( "cm_debugCollision",	"0",		CVAR_ARCHIVE );

	cm_testCollision = cmi.Cvar_Get( "cm_testCollision",		"0",					0 );
	cm_testRotation = cmi.Cvar_Get( "cm_testRotation",		"1",					0 );
	cm_testModel = cmi.Cvar_Get( "cm_testModel",			"0",					0 );
	cm_testTimes = cmi.Cvar_Get( "cm_testTimes",			"1000",					0 );
	cm_testRandomMany = cmi.Cvar_Get( "cm_testRandomMany",	"0",					0 );
	cm_testOrigin = cmi.Cvar_Get( "cm_testOrigin",		"0 0 0",				0 );
	cm_testReset = cmi.Cvar_Get( "cm_testReset",			"0",					0 );
	cm_testBox = cmi.Cvar_Get( "cm_testBox",			"-16 -16 0 16 16 64",	0 );
	cm_testBoxRotation = cmi.Cvar_Get( "cm_testBoxRotation",	"0 0 0",				0 );
	cm_testWalk = cmi.Cvar_Get( "cm_testWalk",			"1",					0 );
	cm_testLength = cmi.Cvar_Get( "cm_testLength",		"1024",					0 );
	cm_testRadius = cmi.Cvar_Get( "cm_testRadius",		"64",					0 );
	cm_testAngle = cmi.Cvar_Get( "cm_testAngle",			"60",					0 );

	cm_noAreas = cmi.Cvar_Get ("cm_noAreas", "0", CVAR_CHEAT);
	cm_noCurves = cmi.Cvar_Get ("cm_noCurves", "0", CVAR_CHEAT);
	cm_playerCurveClip = cmi.Cvar_Get ("cm_playerCurveClip", "1", CVAR_ARCHIVE|CVAR_CHEAT );

    // HACK: find a better place to put this
    HashIndexInit( &cmLocal.materialsHash, 1024, 1024 );
}

void CM_Shutdown() {

}

