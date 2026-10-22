// A fog volume (a shader with fogLight), drawn on its hull's faces that look
// back at the eye: each pixel fogged by the length of its view ray inside
// the hull, up to what's drawn there, as Doom 3's fog lights. The fog is
// Daikatana's (gl_fogsurf.cpp): that length over the volume's distance,
// thicker the lower the ray goes in it.

uniform sampler2D u_ScreenDepthMap;	// the depth prepass's, as a float

uniform vec4   u_ViewInfo;			// zfar / znear, zfar, 1 / width, 1 / height
uniform vec3   u_ViewOrigin;
uniform vec3   u_ViewForward;		// unit

uniform vec4   u_FogPlanes[MAX_FOGLIGHT_PLANES];	// pointing out of the hull
uniform int    u_FogNumPlanes;
uniform vec4   u_FogVolumeInfo;		// lowest z, highest z, distance to opaque
uniform vec4   u_Color;

varying vec3   var_Position;

// depth is GL_DEPTH_COMPONENT24
// so the maximum error is 1.0 / 2^24
#define DEPTH_MAX_ERROR 0.000000059604644775390625

void main()
{
	vec3 toFace = var_Position - u_ViewOrigin;
	float tExit = length(toFace);
	vec3 dir = toFace / tExit;

	// how far along the ray what's drawn is
	float sampleZDivW = texture2D(u_ScreenDepthMap, gl_FragCoord.xy * u_ViewInfo.zw).r - DEPTH_MAX_ERROR;
	float viewZ = u_ViewInfo.y / mix(u_ViewInfo.x, 1.0, sampleZDivW);
	float tScene = viewZ / max(dot(dir, u_ViewForward), 0.0001);

	// where the ray goes into the hull: the last plane it crosses inward;
	// and whether the eye is in it
	float tEnter = 0.0;
	float eyeOut = -1.0;
	for (int i = 0; i < MAX_FOGLIGHT_PLANES; i++)
	{
		if (i >= u_FogNumPlanes)
			break;

		float eyeDist = dot(u_ViewOrigin, u_FogPlanes[i].xyz) - u_FogPlanes[i].w;
		float facing = dot(dir, u_FogPlanes[i].xyz);
		eyeOut = max(eyeOut, eyeDist);
		if (facing < 0.0)
			tEnter = max(tEnter, -eyeDist / facing);
	}

	// from outside, only what's drawn in it is fogged, as Daikatana fogs the
	// surfaces in its fog's leaves: what's seen through it (a pond's floor
	// through its murky top) isn't
	if (eyeOut > 0.0 && tScene > tExit)
		discard;

	float tEnd = min(tExit, tScene);
	float len = tEnd - tEnter;
	if (len <= 0.0)
		discard;

	// Daikatana's R_CalcDensityAdjustment, between where the ray is in it
	float z1 = u_ViewOrigin.z + dir.z * tEnter;
	float z2 = u_ViewOrigin.z + dir.z * tEnd;
	float density = (u_FogVolumeInfo.y - 0.75 * min(z1, z2) - 0.25 * max(z1, z2)) / max(u_FogVolumeInfo.y - u_FogVolumeInfo.x, 1.0);

	gl_FragColor = vec4(u_Color.rgb, clamp(len * max(density, 0.0) / u_FogVolumeInfo.z, 0.0, 1.0));
}
