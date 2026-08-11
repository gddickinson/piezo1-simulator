#version 410 core
// Analytic ray-cylinder intersection, capped at both ends.

in vec3 v_start_view;
in vec3 v_end_view;
in float v_radius;
in vec3 v_color_a;
in vec3 v_color_b;
in vec3 v_quad_view;

uniform mat4 u_proj;
uniform vec3 u_light_dir;
uniform vec3 u_fog_color;
uniform float u_fog_near;
uniform float u_fog_far;
uniform float u_fog_strength;
uniform float u_ambient;
uniform float u_shininess;

out vec4 f_color;

void main()
{
    vec3 ray = normalize(v_quad_view);
    vec3 axis = v_end_view - v_start_view;
    float height = length(axis);
    if (height < 1e-6) discard;
    vec3 dir = axis / height;

    // The eye is the origin in view space, so (eye - base) is just -base.
    vec3 oc = -v_start_view;
    // Components perpendicular to the cylinder axis.
    vec3 rp = ray - dir * dot(ray, dir);
    // This used `-oc`, which negates B and therefore negates both roots of the
    // quadratic — so the near hit came out behind the eye and the `t < 0`
    // guard below discarded every fragment. Cylinders have been invisible
    // since the renderer was written: ball-and-stick drew balls only, and the
    // HaloTag seam drew nothing.
    vec3 op = oc - dir * dot(oc, dir);

    float A = dot(rp, rp);
    if (A < 1e-9) discard;
    float B = 2.0 * dot(rp, op);
    float C = dot(op, op) - v_radius * v_radius;
    float disc = B * B - 4.0 * A * C;
    if (disc < 0.0) discard;

    float t = (-B - sqrt(disc)) / (2.0 * A);
    if (t < 0.0) discard;

    vec3 hit = ray * t;
    float along = dot(hit - v_start_view, dir);
    vec3 normal;
    if (along < 0.0 || along > height) {
        // Missed the barrel — try the flat end cap facing us.
        float plane_d = (along < 0.0) ? 0.0 : height;
        vec3 cap_point = v_start_view + dir * plane_d;
        float denom = dot(ray, dir);
        if (abs(denom) < 1e-6) discard;
        t = dot(cap_point, dir) / denom;
        if (t < 0.0) discard;
        hit = ray * t;
        if (length(hit - cap_point) > v_radius) discard;
        normal = (along < 0.0) ? -dir : dir;
        along = clamp(plane_d, 0.0, height);
    } else {
        normal = normalize(hit - (v_start_view + dir * along));
    }

    vec4 clip = u_proj * vec4(hit, 1.0);
    gl_FragDepth = 0.5 * (clip.z / clip.w) + 0.5;

    vec3 base = mix(v_color_a, v_color_b, clamp(along / height, 0.0, 1.0));

    vec3 L = normalize(u_light_dir);
    vec3 view_dir = -ray;
    float wrap = clamp((dot(normal, L) + 0.35) / 1.35, 0.0, 1.0);
    vec3 half_vec = normalize(L + view_dir);
    float spec = pow(max(dot(normal, half_vec), 0.0), u_shininess);
    vec3 color = base * (u_ambient + (1.0 - u_ambient) * wrap)
               + vec3(0.9, 0.93, 1.0) * spec * 0.30;

    float depth = -hit.z;
    float fog = clamp((depth - u_fog_near) / max(u_fog_far - u_fog_near, 1e-3),
                      0.0, 1.0) * u_fog_strength;
    f_color = vec4(mix(color, u_fog_color, fog), 1.0);
}
