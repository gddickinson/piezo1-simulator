#version 410 core
// Generic lit triangle mesh: cartoons, tubes, molecular surfaces, the membrane.

layout(location = 0) in vec3 in_position;
layout(location = 1) in vec3 in_normal;
layout(location = 2) in vec3 in_color;
layout(location = 3) in float in_alpha;

uniform mat4 u_view;
uniform mat4 u_proj;
uniform mat3 u_normal_matrix;

out vec3 v_color;
out vec3 v_normal_view;
out vec3 v_position_view;
out float v_alpha;

void main()
{
    vec4 pos_view = u_view * vec4(in_position, 1.0);
    v_position_view = pos_view.xyz;
    v_normal_view = normalize(u_normal_matrix * in_normal);
    v_color = in_color;
    v_alpha = in_alpha;
    gl_Position = u_proj * pos_view;
}
