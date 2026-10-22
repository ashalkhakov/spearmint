uniform vec4  u_Color;

uniform int   u_FogType;
uniform float u_FogEyeT;

varying vec2  var_FogST;

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

void main()
{
	gl_FragColor = u_Color;
	gl_FragColor.a = sqrt(clamp(CalcFog(var_FogST) * u_Color.a * u_Color.a, 0.0, 1.0));
}
