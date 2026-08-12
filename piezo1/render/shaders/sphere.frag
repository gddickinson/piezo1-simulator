#version 410 core
// Sphere impostor, fragment stage: analytic ray-sphere intersection.

in vec3 v_color;
in vec3 v_center_view;
in float v_radius;
in vec2 v_offset;
in float v_flags;

uniform mat4 u_proj;
uniform vec3 u_light_dir;      // in view space, pointing towards the light
uniform vec3 u_fog_color;
uniform float u_fog_near;
uniform float u_fog_far;
uniform float u_fog_strength;
uniform float u_ambient;
uniform float u_shininess;
uniform int u_outline;
// Per-batch opacity. 1.0 for every molecular representation; the pore
// probe spheres are drawn translucent so the lining stays visible
// through them, which is the whole point of drawing them at all.
uniform float u_alpha;

out vec4 f_color;

void main()
{
    // The eye sits at the origin in view space looking down -z.
    vec3 quad_point = v_center_view + vec3(v_offset, 0.0);
    vec3 ray = normalize(quad_point);

    float b = dot(ray, v_center_view);
    float c = dot(v_center_view, v_center_view) - v_radius * v_radius;
    float disc = b * b - c;
    if (disc < 0.0) discard;               // the quad corner missed the sphere

    float t = b - sqrt(disc);
    vec3 hit = ray * t;
    vec3 normal = (hit - v_center_view) / v_radius;

    // Write the true sphere depth so impostors interpenetrate correctly with
    // each other and with ordinary triangle geometry.
    vec4 clip = u_proj * vec4(hit, 1.0);
    gl_FragDepth = 0.5 * (clip.z / clip.w) + 0.5;

    vec3 base = v_color;
    bool selected = (int(v_flags) & 1) != 0;
    bool dimmed   = (int(v_flags) & 2) != 0;
    if (dimmed) base = mix(base, u_fog_color, 0.65);

    // Blinn-Phong with a soft wrap term, which keeps the underside of a large
    // assembly readable instead of pure black.
    vec3 L = normalize(u_light_dir);
    float ndl = dot(normal, L);
    float wrap = clamp((ndl + 0.35) / 1.35, 0.0, 1.0);
    vec3 view_dir = -ray;
    vec3 half_vec = normalize(L + view_dir);
    float spec = pow(max(dot(normal, half_vec), 0.0), u_shininess);

    vec3 color = base * (u_ambient + (1.0 - u_ambient) * wrap)
               + vec3(0.9, 0.93, 1.0) * spec * 0.35;

    // Rim light picks out silhouettes against the dark background.
    float rim = pow(1.0 - max(dot(normal, view_dir), 0.0), 3.0);
    color += base * rim * 0.25;

    if (selected) {
        float edge = smoothstep(0.35, 0.0, max(dot(normal, view_dir), 0.0));
        color = mix(color, vec3(1.0, 0.85, 0.2), edge * 0.9);
    }

    // Depth cue: fade distant atoms into the background colour.
    float depth = -hit.z;
    float fog = clamp((depth - u_fog_near) / max(u_fog_far - u_fog_near, 1e-3),
                      0.0, 1.0) * u_fog_strength;
    color = mix(color, u_fog_color, fog);

    f_color = vec4(color, u_alpha);
}
