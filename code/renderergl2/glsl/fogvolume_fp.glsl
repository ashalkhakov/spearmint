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

// the dynamic lights reaching it: scattered toward the eye, they brighten it
uniform vec4   u_FogLights[MAX_FOG_DLIGHTS];		// origin, radius
uniform vec4   u_FogLightColors[MAX_FOG_DLIGHTS];	// colour times intensity and scattering
uniform vec4   u_FogLightDirs[MAX_FOG_DLIGHTS];	// a spot's direction and cone cosine, w -2 a point
uniform int    u_FogNumLights;

// the shadows of some of them (tr_dlightshadow.c): dark shafts in the fog
uniform sampler2D u_DlightShadowMap;
uniform vec4   u_FogShadowMvps[MAX_FOG_SHADOWS * 24];	// a face's world to clip matrix, 4 columns each
uniform vec4   u_FogShadowTiles[MAX_FOG_SHADOWS * 6];
uniform vec4   u_FogShadowLights[MAX_FOG_SHADOWS];		// which of the lights (-1 none), faces, near, far

// whether the light reaches a point in the fog, by its shadow map (one tap:
// the fog's eight samples soften it)
float FogShadow(int s, vec3 p)
{
	vec3 v = p - u_FogLights[int(u_FogShadowLights[s].x)].xyz;
	int face = 0;

	if (u_FogShadowLights[s].y > 1.5)
	{
		vec3 a = abs(v);

		if (a.x >= a.y && a.x >= a.z)
			face = v.x > 0.0 ? 1 : 0;
		else if (a.y >= a.z)
			face = v.y > 0.0 ? 3 : 2;
		else
			face = v.z > 0.0 ? 5 : 4;
	}

	int m0 = s * 24 + face * 4;
	mat4 m = mat4(u_FogShadowMvps[m0], u_FogShadowMvps[m0 + 1], u_FogShadowMvps[m0 + 2], u_FogShadowMvps[m0 + 3]);
	vec4 tile = u_FogShadowTiles[s * 6 + face];
	vec4 clip = m * vec4(p, 1.0);
	vec2 uv = tile.xy + (clip.xy / clip.w * 0.5 + 0.5) * tile.zw;
	float z = texture2D(u_DlightShadowMap, uv).r * 2.0 - 1.0;
	float n = u_FogShadowLights[s].z, f = u_FogShadowLights[s].w;
	float d = 2.0 * n * f / (f + n - z * (f - n));

	return float(d >= clip.w - 2.0);
}

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

	// where the ray goes into the hull: the last plane it crosses inward.
	// What's seen through it is fogged by the whole run, as Daikatana draws
	// its fog volume's far faces (a mist slab's bottom) as fog over what's
	// behind them
	float tEnter = 0.0;
	for (int i = 0; i < MAX_FOGLIGHT_PLANES; i++)
	{
		if (i >= u_FogNumPlanes)
			break;

		float eyeDist = dot(u_ViewOrigin, u_FogPlanes[i].xyz) - u_FogPlanes[i].w;
		float facing = dot(dir, u_FogPlanes[i].xyz);
		if (facing < 0.0)
			tEnter = max(tEnter, -eyeDist / facing);
	}


	float tEnd = min(tExit, tScene);
	float len = tEnd - tEnter;
	if (len <= 0.0)
		discard;

	// Daikatana's R_CalcDensityAdjustment, between where the ray is in it
	float z1 = u_ViewOrigin.z + dir.z * tEnter;
	float z2 = u_ViewOrigin.z + dir.z * tEnd;
	float density = (u_FogVolumeInfo.y - 0.75 * min(z1, z2) - 0.25 * max(z1, z2)) / max(u_FogVolumeInfo.y - u_FogVolumeInfo.x, 1.0);

	// the lights' light along the run, in the fog's own hue
	vec3 fogColor = u_Color.rgb;
	if (u_FogNumLights > 0)
	{
		vec3 lit = vec3(0.0);

		for (int k = 0; k < 8; k++)
		{
			vec3 p = u_ViewOrigin + dir * mix(tEnter, tEnd, (float(k) + 0.5) / 8.0);

			for (int l = 0; l < MAX_FOG_DLIGHTS; l++)
			{
				if (l >= u_FogNumLights)
					break;

				vec3 toLight = u_FogLights[l].xyz - p;
				float dist = length(toLight);
				float a = clamp(1.0 - dist / u_FogLights[l].w, 0.0, 1.0);

				a *= a;
				if (u_FogLightDirs[l].w > -1.5)
					a *= smoothstep(u_FogLightDirs[l].w, mix(u_FogLightDirs[l].w, 1.0, 0.3), dot(-toLight / max(dist, 0.001), u_FogLightDirs[l].xyz));
				for (int sh = 0; sh < MAX_FOG_SHADOWS; sh++)
				{
					if (int(u_FogShadowLights[sh].x) == l && a > 0.0)
						a *= FogShadow(sh, p);
				}
				lit += u_FogLightColors[l].rgb * a;
			}
		}
		fogColor += u_Color.rgb / max(max(u_Color.r, u_Color.g), max(u_Color.b, 0.001)) * lit * 0.125;
	}

	gl_FragColor = vec4(fogColor, clamp(len * max(density, 0.0) / u_FogVolumeInfo.z, 0.0, 1.0));
}
