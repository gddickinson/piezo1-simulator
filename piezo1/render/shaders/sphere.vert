#version 410 core
// Sphere impostor, vertex stage.
//
// Each atom is one instance. Rather than tessellating a sphere, we emit a
// single view-aligned quad that is guaranteed to cover the sphere's
// silhouette; the fragment stage then ray-casts the true sphere inside it.
// Four vertices per atom instead of a few hundred triangles, and the result is
// pixel-exact at any zoom.

layout(location = 0) in vec3 in_center;
layout(location = 1) in float in_radius;
layout(location = 2) in vec3 in_color;
layout(location = 3) in float in_flags;   // bit 0: selected, bit 1: dimmed

uniform mat4 u_view;
uniform mat4 u_proj;
uniform float u_radius_scale;

out vec3 v_color;
out vec3 v_center_view;
out float v_radius;
out vec2 v_offset;
out float v_flags;

void main()
{
    // Corners of a unit quad from gl_VertexID, drawn as a triangle strip.
    vec2 corner = vec2(((gl_VertexID & 1) == 0) ? -1.0 : 1.0,
                       ((gl_VertexID & 2) == 0) ? -1.0 : 1.0);

    float radius = in_radius * u_radius_scale;
    vec4 center_view = u_view * vec4(in_center, 1.0);

    // Oversize slightly: under perspective the silhouette of a sphere is
    // marginally larger than its equatorial disc, and clipping it looks awful.
    vec2 offset = corner * radius * 1.12;

    v_color = in_color;
    v_center_view = center_view.xyz;
    v_radius = radius;
    v_offset = offset;
    v_flags = in_flags;

    gl_Position = u_proj * vec4(center_view.xyz + vec3(offset, 0.0), 1.0);
}
