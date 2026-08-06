#version 410 core

in vec3 v_color;
in vec3 v_normal_view;
in vec3 v_position_view;
in float v_alpha;

uniform vec3 u_light_dir;
uniform vec3 u_fog_color;
uniform float u_fog_near;
uniform float u_fog_far;
uniform float u_fog_strength;
uniform float u_ambient;
uniform float u_shininess;
uniform int u_two_sided;

out vec4 f_color;

void main()
{
    vec3 normal = normalize(v_normal_view);
    vec3 view_dir = normalize(-v_position_view);

    // Cartoon ribbons and membrane sheets are open surfaces, so back faces
    // must be lit as if their normal were flipped rather than left black.
    if (u_two_sided != 0 && dot(normal, view_dir) < 0.0) normal = -normal;

    vec3 L = normalize(u_light_dir);
    float ndl = dot(normal, L);
    float wrap = clamp((ndl + 0.35) / 1.35, 0.0, 1.0);
    vec3 half_vec = normalize(L + view_dir);
    float spec = pow(max(dot(normal, half_vec), 0.0), u_shininess);

    vec3 color = v_color * (u_ambient + (1.0 - u_ambient) * wrap)
               + vec3(0.9, 0.93, 1.0) * spec * 0.30;

    float rim = pow(1.0 - max(dot(normal, view_dir), 0.0), 3.0);
    color += v_color * rim * 0.20;

    float depth = -v_position_view.z;
    float fog = clamp((depth - u_fog_near) / max(u_fog_far - u_fog_near, 1e-3),
                      0.0, 1.0) * u_fog_strength;
    color = mix(color, u_fog_color, fog);

    f_color = vec4(color, v_alpha);
}
