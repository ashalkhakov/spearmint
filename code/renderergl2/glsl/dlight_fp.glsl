uniform sampler2D u_DiffuseMap;

uniform int       u_AlphaTest;
uniform float     u_AlphaTestRef;

uniform vec4      u_DlightInfo;
uniform float     u_Intensity;
uniform vec4      u_SpotInfo;		// a spot light's direction and cone cosine; w -2 for a point light

varying vec2      var_Tex1;
varying vec4      var_Color;
varying vec3      var_DlightDist;
varying vec3      var_Normal;
varying float     var_PerPixel;


void main()
{
	vec4 color = texture2D(u_DiffuseMap, var_Tex1);
	vec4 lightColor = var_Color;

	// a Q3 dlight: its depth falloff and the surface's facing here, not
	// interpolated from the vertexes (which breaks up on big polygons)
	if (var_PerPixel > 0.5)
	{
		float dlightmod = step(0.0, dot(var_DlightDist, var_Normal));
		dlightmod *= clamp(u_Intensity * 2.0 * (1.0 - abs(var_DlightDist.z) * u_DlightInfo.a), 0.0, 1.0);
		if (u_SpotInfo.w > -1.5)
			dlightmod *= smoothstep(u_SpotInfo.w, mix(u_SpotInfo.w, 1.0, 0.3), dot(-normalize(var_DlightDist), u_SpotInfo.xyz));
		lightColor *= dlightmod;
	}

	float alpha = color.a * lightColor.a;
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
	
	gl_FragColor.rgb = color.rgb * lightColor.rgb;
	gl_FragColor.a = alpha;
}
