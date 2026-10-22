uniform sampler2D u_DiffuseMap;

uniform int       u_AlphaTest;
uniform float     u_AlphaTestRef;

#if defined(USE_LIGHTMAP)
uniform sampler2D u_LightMap;

uniform int       u_Texture1Env;
#endif

varying vec2      var_DiffuseTex;

#if defined(USE_LIGHTMAP)
varying vec2      var_LightTex;
#endif

varying vec4      var_Color;

#if defined(USE_FOG)
uniform int       u_FogType;
uniform float     u_FogEyeT;
uniform vec4      u_FogColorMask;

varying vec2      var_FogST;

float CalcFog(vec2 st)
{
	float s = st.x;
	float t = st.y;

	float eyeOutside = float(u_FogEyeT < 0.0);

	if ( u_FogType == FT_LINEAR ) {
		if ( eyeOutside == 0.0 )
			t += u_FogEyeT;

		return s * t;
	} else {
		// FT_EXP
		s *= 8.0;

		float fogged = float(t >= eyeOutside);

		t += 1e-6;
		t *= fogged / (t - u_FogEyeT * eyeOutside);

		return s * t;
	}
}
#endif


void main()
{
	vec4 color  = texture2D(u_DiffuseMap, var_DiffuseTex);
#if defined(USE_LIGHTMAP)
	vec4 color2 = texture2D(u_LightMap, var_LightTex);
  #if defined(RGBM_LIGHTMAP)
	color2.rgb *= color2.a;
	color2.a = 1.0;
  #endif

	if (u_Texture1Env == TEXENV_MODULATE)
	{
		color *= color2;
	}
	else if (u_Texture1Env == TEXENV_ADD)
	{
		color += color2;
	}
	else if (u_Texture1Env == TEXENV_REPLACE)
	{
		color = color2;
	}
	
	//color = color * (u_Texture1Env.xxxx + color2 * u_Texture1Env.z) + color2 * u_Texture1Env.y;
#endif

	vec4 vertColor = var_Color;
#if defined(USE_FOG)
	vertColor *= vec4(1.0) - u_FogColorMask * sqrt(clamp(CalcFog(var_FogST), 0.0, 1.0));
#endif

	float alpha = color.a * vertColor.a;
	if (u_AlphaTest == U_ATEST_EQUAL)
	{
		if (alpha != u_AlphaTestRef)
			discard;
	}
	else if (u_AlphaTest == U_ATEST_GREATEREQUAL)
	{
		if (alpha < u_AlphaTestRef)
			discard;
	}
	else if (u_AlphaTest == U_ATEST_LESS)
	{
		if (alpha >= u_AlphaTestRef)
			discard;
	}
	else if (u_AlphaTest == U_ATEST_LESSEQUAL)
	{
		if (alpha > u_AlphaTestRef)
			discard;
	}
	else if (u_AlphaTest == U_ATEST_NOTEQUAL)
	{
		if (alpha == u_AlphaTestRef)
			discard;
	}
	else if (u_AlphaTest == U_ATEST_GREATER)
	{
		if (alpha <= u_AlphaTestRef)
			discard;
	}
	
	gl_FragColor.rgb = color.rgb * vertColor.rgb;
	gl_FragColor.a = alpha;
}
