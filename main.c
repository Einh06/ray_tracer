#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <time.h>

#include "common.h"
#include "vec3.h"
#include "material.h"

#include "vec3.c"
#include "material.c"

v3 random_in_unit_sphere(void) {
    v3 p;
    do {
        p = v3_Sub(
            v3_Mul((v3){drand48(), drand48(), drand48()}, 2.0),
            (v3){1.0, 1.0, 1.0}
            );
    } while (v3_SquaredLength(p) >= 1.0);
    return p;
}

v3 random_in_unit_disk(void) {
    v3 p;
    do {
        p = v3_Sub(
            v3_Mul((v3){drand48(), drand48(), 0.f}, 2.f), 
            (v3){1.f, 1.f, 0.f});
    } while( v3_Dot(p, p) >= 1.f);
    return p;
}

typedef struct ray3 {
    v3 origin;
    v3 direction;
    float time;
} ray3;

inline v3 ray3_PointAt(ray3 r, float t) {
    return v3_Add(r.origin, v3_Mul(r.direction, t));
}

typedef struct Camera {
    v3 origin;
    v3 lower_left;
    v3 horizontal;
    v3 vertical;
    v3 u, v, w;
    float lens_radius;
    float time0, time1;
} Camera;

Camera Camera_Make(v3 lookFrom, v3 lookAt, v3 vup, float vfov, float aspect, float aperture, float focus_dist, float time0, float time1) {
    
    float lens_radius = aperture / 2.f;
    
    float theta = vfov*M_PI/180.f;
    float half_height = tan(theta/2.f);
    float half_width = aspect * half_height;
    
    v3 w = v3_Normalize(v3_Sub(lookFrom, lookAt));
    v3 u = v3_Normalize(v3_Cross(vup, w));
    v3 v = v3_Cross(w, u);
    
    v3 horizontal = v3_Mul(u, 2*half_width*focus_dist);
    v3 vertical = v3_Mul(v, 2*half_height*focus_dist);
    //((origin - u*half_width) - v * half_height) - w
    v3 lower_left =
        v3_Sub(v3_Sub(v3_Sub(
        lookFrom, 
        v3_Mul(u, half_width * focus_dist)), 
                      v3_Mul(v, half_height * focus_dist)), 
               v3_Mul(w, focus_dist));
    
    return (Camera){
        .origin     = lookFrom,
        .horizontal = horizontal,
        .vertical   = vertical,
        .lower_left = lower_left,
        .u = u, .v = v, .w = w,
        .lens_radius = lens_radius,
        .time0 = time0, .time1 = time1,
    };
}

ray3 Camera_GetRay(Camera *c, float s, float t) {
    v3 rd = v3_Mul(random_in_unit_disk(), c->lens_radius);
    v3 offset = v3_Add(
        v3_Mul(c->u, rd.x),
        v3_Mul(c->v, rd.y));
    
    v3 h = v3_Mul(c->horizontal, s);
    v3 v = v3_Mul(c->vertical, t);
    v3 o = v3_Add(c->origin, offset);
    
    float time = c->time0 + (drand48() * (c->time1 - c->time0));

    v3 direction = v3_Sub(v3_Add(c->lower_left, v3_Add(h, v)), o);
    return  (ray3) {
        .origin = o,
        .direction = direction,
        .time = time,
    };
}

typedef enum ObjectKind {
    OBJECTKIND_SPHERE,
}ObjectKind;

typedef struct Object {
    
    ObjectKind kind;
    v3 pos0, pos1; // managing motion blur
    float time0, time1;
    union {
        struct {
            float radius;
        } sphere;
    };
    Material mat;
} Object;

v3 Object_CenterAtTime(Object *o, float time) {
    float t = ((time - o->time0) / (o->time1 - o->time0));
    return v3_Add(o->pos0, v3_Mul(v3_Sub(o->pos1, o->pos0), t));
}

typedef struct HitInformation {
    float t;
    v3 pos;
    v3 norm;
    Material mat;
} HitInformation;

v3 Reflected(v3 v, v3 n) {
    return v3_Sub(v,
                  v3_Mul(n, 2.0f * v3_Dot(v, n))
                  );
}

bool Refracted(v3 v, v3 n, float ni_over_nt, v3* refracted) 
{
    v3 uv = v3_Normalize(v);
    
    float dt = v3_Dot(uv, n);
    float discriminant = 1.0 - ni_over_nt * ni_over_nt * (1 - dt * dt);
    if (discriminant > 0) {
        //ni_over_nt*(uv - n*dt)
        v3 lhs = v3_Mul(
            v3_Sub(uv, v3_Mul(n, dt)), 
            ni_over_nt);
        //n*sqrt(discriminant)
        v3 rhs = v3_Mul(n, sqrt(discriminant));
        //(ni_over_nt*(uv-n*dt) - (n*sqrt(discriminant))
        *refracted = v3_Sub(lhs, rhs);
        return true;
    } else {
        return false;
    }
}

float Schilk(float cosine, float ref_idx) {
    float r0 = (1.f-ref_idx) / (1.f+ref_idx);
    r0 = r0 * r0;
    return r0 + (1.f-r0) * powf((1.f-cosine), 5);
}

bool Scatter(ray3* r, HitInformation* h, v3* attenuation, ray3* scattered) {
    switch(h->mat.kind) {
        case MATKIND_LAMBERTIAN:
        {
            v3 target = v3_Add(
                random_in_unit_sphere(),
                v3_Add(h->pos, h->norm)
                );
            
            *scattered = (ray3){.origin = h->pos, .direction = v3_Sub(target, h->pos), .time = r->time};
            *attenuation = h->mat.lambertian.albedo;
            return true;
        }
        
        case MATKIND_METAL: 
        {
            v3 reflected = Reflected(v3_Normalize(r->direction), h->norm);
            v3 random_in_sphere = v3_Mul(random_in_unit_sphere(), h->mat.metal.fuzz);
            *scattered = (ray3){
                .origin = h->pos, 
                .direction = v3_Add(reflected, random_in_sphere),
                .time = r->time,
            };
            *attenuation = h->mat.metal.albedo;
            return v3_Dot(scattered->direction, h->norm) > 0.f;
        }
        
        case MATKIND_DIELETRIC:
        {
            v3 outward_normal;
            v3 reflected = Reflected(r->direction, h->norm);
            float ni_over_nl;
            *attenuation = (v3){1.0, 1.0, 1.0};
            v3 refracted;
            float cosine;
            float reflect_prob;
            
            float ref_idx = h->mat.dieletric.ref_idx;
            
            float dot = v3_Dot(r->direction, h->norm);
            if (dot > 0) {
                outward_normal = v3_Neg(h->norm);
                ni_over_nl = ref_idx;
                cosine = dot / v3_Length(r->direction);
                cosine = sqrt( 1 - ref_idx*ref_idx * (1 - cosine*cosine));
            } else {
                outward_normal = h->norm;
                ni_over_nl = 1.0 / ref_idx;
                cosine = -dot / v3_Length(r->direction);
            }
            
            if (Refracted(r->direction, outward_normal, ni_over_nl, &refracted)) {
                reflect_prob = Schilk(cosine, ref_idx);
            } else {
                reflect_prob = 1.0;
            }
            
            if (drand48() < reflect_prob) {
                *scattered = (ray3){h->pos, reflected, r->time};
            } else {
                *scattered = (ray3){h->pos, refracted, r->time};
            }
            return true;
        }
        default:
        exit(1);
    }
    return false;
}

bool hit_sphere(ray3 r, Object *o, float tmin, float tmax, HitInformation *h) {
    assert(o->kind == OBJECTKIND_SPHERE);
    
    v3 sphere_center = Object_CenterAtTime(o, r.time);
    float radius = o->sphere.radius;
    
    v3 oc = v3_Sub(r.origin, sphere_center);
    float a = v3_Dot(r.direction, r.direction);
    float b = v3_Dot(oc, r.direction);
    float c = v3_Dot(oc, oc) - radius * radius;
    float discriminant = b*b - a*c;
    
    if (discriminant >= 0.f) {
        float temp = (-b - sqrt(discriminant)) / a;
        if (temp >= tmin && temp <= tmax) {
            h->t = temp;
            h->pos = ray3_PointAt(r, temp);
            h->norm = v3_Normalize(v3_Mul(v3_Sub(h->pos, sphere_center), 1.f/radius));
            return true;
        }
        
        temp = (-b + sqrt(discriminant)) / a;
        if (temp >= tmin && temp <= tmax) {
            h->t = temp;
            h->pos = ray3_PointAt(r, temp);
            h->norm = v3_Normalize(v3_Mul(v3_Sub(h->pos, sphere_center), 1.f/radius));
            return true;
        }
    }
    
    return false;
}

v3 Color(ray3 r, Object *objects, int object_count, int depth) {
    HitInformation h;
    
    float closest_so_far = FLT_MAX;
    bool hit_something = false;
    for (Object *o = objects; o < objects + object_count; o++){
        
        HitInformation temp_hit;
        
        switch (o->kind) {
            case OBJECTKIND_SPHERE: {
                if (!hit_sphere(r, o, 0.001, closest_so_far,  &temp_hit)) { continue; }
                closest_so_far = temp_hit.t;
                hit_something = true;
                h = temp_hit;
                h.mat = o->mat;
                break;
            }
            
            default: {
                perror("Unrecognized object type");
                exit(1);
            }
        }
    }
    
    if (hit_something) {
        ray3 new_ray;
        v3 attenutation;
        if (depth < 50 && Scatter(&r, &h,&attenutation, &new_ray)) {
            v3 color = Color(new_ray, objects, object_count, depth + 1);
            return (v3){color.r * attenutation.r, color.g * attenutation.g, color.b * attenutation.b};
        } else {
            return (v3){ 0.f, 0.f, 0.f};
        }
    } else {
        v3 unit_direction = v3_Normalize(r.direction);
        float t = (unit_direction.y + 1.f) * 0.5f;
        return v3_Add(
            v3_Mul(v3_Make(1.f, 1.f, 1.f), (1.f-t)),
            v3_Mul(v3_Make(0.5f, 0.7f, 1.f), t)
            );
    }
}

#define OUT_PRINT(str, ...) (fprintf(output, str, __VA_ARGS__))
int main(int argc, char** argv) {
    
    srand(time(NULL));
    FILE* output = fopen("output/result.ppm", "w");
    if (output == NULL) {
        return 1;
    }
    
    int width = 1200;
    int height = 800;
    
    float inv_w = 1.0f / ((float)width);
    float inv_h = 1.0f / ((float)height);
    
    OUT_PRINT("P3\n%d %d\n255\n", width, height);
    
#if 0
    v3 lookfrom = v3_Make(3.f, 3.f, 2.f);
    v3 lookat = v3_Make(0.f, 0.f, -1.f);
    
    float distance_to_focus = v3_Length(v3_Sub(lookfrom, lookat));
    float aperture = 2.0f;
    
    Camera c = Camera_Make(lookfrom, lookat, v3_Make(0.f, 1.f, 0.f), 20.f, ((float)width) / ((float) height), aperture, distance_to_focus);
    
    Object world[] = {
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(0.f, 0.f, -1.f), 
            .sphere = {.radius = 0.5f},
            .mat = Material_MakeLambertian(v3_Make(0.1f, 0.2f, 0.5f)),
        },
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(1.f, 0.f, -1.f), 
            .sphere = {.radius = 0.5f},
            .mat = Material_MakeMetal(v3_Make(0.8f, 0.6f, 0.2f), 0.0f),
        },
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(-1.f, 0.f, -1.f), 
            .sphere = {.radius = 0.5f},
            .mat = Material_MakeDieletric(1.5f),
        },
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(-1.f, 0.f, -1.f), 
            .sphere = {.radius = -0.45f},
            .mat = Material_MakeDieletric(1.5f),s
        },
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(0.f, -100.5f, -1.f), 
            .sphere = {.radius = 100.0f},
            .mat = Material_MakeLambertian(v3_Make(0.8f, 0.8f, 0.0f)),
        },
    };
    
    int world_size = ArrayCount(world);
#else
    v3 lookfrom = v3_Make(13.f, 2.f, 3.f);
    v3 lookat = v3_Make(0.f, 0.f, 0.f);
    
    float distance_to_focus = 10.f;
    float aperture = 0.1f;
    
    Camera c = Camera_Make(lookfrom, lookat, v3_Make(0.f, 1.f, 0.f), 20.f, ((float)width) / ((float) height), aperture, distance_to_focus,
        0.0f, 1.0f);
    
    Object world[501];
    world[0] = (Object){
        .kind = OBJECTKIND_SPHERE,
        .pos0 = (v3){0, -1000, 0}, .pos1 = (v3){0, -1000, 0},
        .time0 = 0.0f, .time1 = 1.0f,
        .sphere = {.radius = 1000},
        .mat = Material_MakeLambertian((v3){0.5, 0.5, 0.5}),
    };
    
    int o = 1;
    for (int a = -11; a < 11; a++) {
        for (int b = -11; b < 11; b++) {
            float choose_mat = drand48();
            v3 center0 = (v3){a + 0.9 * drand48(), 0.2, b + 0.9 * drand48()};
            v3 center1 = center0;

            if ( v3_Length(v3_Sub(center0, (v3){4, 0.2, 0})) > 0.9 ) {
                Material m;
                if (choose_mat < 0.8) { // diffuse
                    center1 = v3_Add(center1, (v3){0, 0.3 * (1 + drand48()), 0});
                    m = Material_MakeLambertian((v3){ drand48() * drand48(), drand48() * drand48(), drand48() * drand48()});
                } else if (choose_mat < 0.95) { // Metal
                    m = Material_MakeMetal((v3){0.5 * ( 1 * drand48() ),  0.5 * ( 1 * drand48() ), 0.5 * ( 1 * drand48() )}, 0.5 * drand48());
               q } else {
                    m = Material_MakeDieletric(1.5);
                }
                world[o++] = (Object){
                    .kind = OBJECTKIND_SPHERE,
                    .pos0 = center0, .pos1 = center1,
                    .time0 = 0.0f, .time1 = 1.0f,
                    .sphere = {.radius = 0.2},
                    .mat = m,
                };
            }
        }
    }
    
    world[o++] = (Object){
        .kind = OBJECTKIND_SPHERE,
        .pos0 = (v3){0, 1, 0}, .pos1 = (v3){0, 1, 0},
        .time0 = 0.0f, .time1 = 1.0f,
        .sphere = {.radius = 1.0},
        .mat = Material_MakeDieletric(1.5),
    };
    
    world[o++] = (Object){
        .kind = OBJECTKIND_SPHERE,
        .pos0 = (v3){-4, 1, 0}, .pos1 = (v3){-4, 1, 0},
        .time0 = 0.0f, .time1 = 1.0f,
        .sphere = {.radius = 1.0},
        .mat = Material_MakeLambertian((v3){0.4, 0.2, 0.1}),
    };
    
    world[o++] = (Object){
        .kind = OBJECTKIND_SPHERE,
        .pos0 = (v3){4, 1, 0}, .pos1 = (v3){4, 1, 0},
        .time0 = 0.0f, .time1 = 1.0f,
        .sphere = {.radius = 1.0},
        .mat = Material_MakeMetal((v3){0.7, 0.6, 0.5}, 0.0),
    };
    
    int world_size = o;
#endif
    
    int rand_step = 100;
    for (int j = height - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            OUT_PRINT(" ");
            v3 color = {0};
            for (int s = 0; s < rand_step; s++) {
                float u = ((float) i + drand48()) * inv_w;
                float v = ((float) j + drand48()) * inv_h;
                
                ray3 r = Camera_GetRay(&c, u, v);
                
                v3_Added(&color, Color(r, world, world_size, 0)); // get color mapped from 0.0 - 1.0
            }
            v3_Scaled(&color, (1.f/((float)rand_step)));
            color = (v3){sqrt(color.r), sqrt(color.g), sqrt(color.b)};
            color = v3_Mul(color, 255.99); // map colors from 0.0 - 1.0 to 0 - 255
            OUT_PRINT("%d %d %d", (int)color.r, (int)color.g, (int)color.b);
        }
        OUT_PRINT("\n");
    }
    return 0;
}

