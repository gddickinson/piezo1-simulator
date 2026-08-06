#version 410 core
// Cylinder impostor: bonds, backbone traces and mode-displacement arrows.
// One instance per cylinder; a six-vertex triangle strip forms a screen-facing
// box that bounds the cylinder, and the fragment stage ray-casts it.

layout(location = 0) in vec3 in_start;
layout(location = 1) in vec3 in_end;
layout(location = 2) in float in_radius;
layout(location = 3) in vec3 in_color_a;
layout(location = 4) in vec3 in_color_b;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_radius_scale;

out vec3 v_start_view;
out vec3 v_end_view;
out float v_radius;
out vec3 v_color_a;
out vec3 v_color_b;
out vec3 v_quad_view;

void main()
{
    vec3 a = (u_view * vec4(in_start, 1.0)).xyz;
    vec3 b = (u_view * vec4(in_end, 1.0)).xyz;
    float radius = in_radius * u_radius_scale;

    vec3 axis = b - a;
    float len = length(axis);
    vec3 dir = (len > 1e-6) ? axis / len : vec3(0.0, 0.0, 1.0);

    // Build a frame whose x lies across the screen-projected axis.
    vec3 to_eye = normalize(-0.5 * (a + b));
    vec3 side = normalize(cross(dir, to_eye));
    if (any(isnan(side))) side = normalize(cross(dir, vec3(0.0, 0.0, 1.0)));

    // Four corners of a screen-facing rectangle, padded by the radius.
    int id = gl_VertexID;
    float along = ((id & 1) == 0) ? -1.0 : 1.0;
    float across = ((id & 2) == 0) ? -1.0 : 1.0;
    vec3 base = (along < 0.0) ? a : b;
    vec3 pos = base + dir * (along * radius * 1.05) + side * (across * radius * 1.05);

    v_start_view = a;
    v_end_view = b;
    v_radius = radius;
    v_color_a = in_color_a;
    v_color_b = in_color_b;
    v_quad_view = pos;

    gl_Position = u_proj * vec4(pos, 1.0);
}
