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
#include "ray.h"
#include "camera.h"

#include "vec3.c"
#include "math_util.c"
#include "material.c"
#include "ray.c"
#include "camera.c"

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

int Object_XCmp(const void *a, const void *b) {
    Object *o1 = (Object *)a, *o2 = (Object *)b;
    if (o1->pos0.x < o2->pos0.x) return -1;
    if (o1->pos0.x > o2->pos0.x) return 1;
    return 0;
}

int Object_YCmp(const void *a, const void *b) {
    Object *o1 = (Object *)a, *o2 = (Object *)b;
    if (o1->pos0.y < o2->pos0.y) return -1;
    if (o1->pos0.y > o2->pos0.y) return 1;
    return 0;
}

int Object_ZCmp(const void *a, const void *b) {
    Object *o1 = (Object *)a, *o2 = (Object *)b;
    if (o1->pos0.z < o2->pos0.z) return -1;
    if (o1->pos0.z > o2->pos0.z) return 1;
    return 0;
}

typedef struct Aabb {
    v3 min, max;
} Aabb;

typedef enum BvhNodeKind {
    BVHNODEKIND_OBJECTS,
    BVHNODEKIND_NODES,
} BvhNodeKind;

typedef struct BvhNode {
    BvhNodeKind kind;
    union {
        struct {
            Object *left, *right;
        } objects;
        struct {
            struct BvhNode *left, *right;
        } nodes;
    };
    Aabb box;
} BvhNode;

bool Aabb_MakeFromObject(Object *o, float t0, float t1, Aabb *b) {

    v3 posMin = v3_Min(o->pos0, o->pos1);
    v3 posMax = v3_Max(o->pos0, o->pos1);
    switch (o->kind) {
        case OBJECTKIND_SPHERE: {
            *b = (Aabb){
                .min = v3_Sub(posMin, (v3){o->sphere.radius, o->sphere.radius, o->sphere.radius}),
                .max = v3_Add(posMax, (v3){o->sphere.radius, o->sphere.radius, o->sphere.radius}),
            };
            return true;
        } break;
    }
    return false;
}

bool Aabb_Compound(Aabb b1, Aabb b2, Aabb *r) {
    *r = (Aabb) {
        .min = (v3){fmin(b1.min.x, b2.min.x), fmin(b1.min.y, b2.min.y), fmin(b1.min.z, b2.min.z)},
        .max = (v3){fmax(b1.max.x, b2.max.x), fmax(b1.max.y, b2.max.y), fmax(b1.max.z, b2.max.z)}
    };
    return true;
}

bool hit_aabb(Aabb b, ray3 r) {
    v3 invRayDir = v3_Inv(r.direction);

    v3 t0 = v3_Sub(b.min, r.origin);
    t0 = (v3){t0.x * invRayDir.x, t0.y * invRayDir.y, t0.z * invRayDir.z};

    v3 t1 = v3_Sub(b.max, r.origin);
    t1 = (v3){t1.x * invRayDir.x, t1.y * invRayDir.y, t1.z * invRayDir.z}; 

    v3 tmin = v3_Min(t0, t1), tmax = v3_Max(t0, t1);
    return v3_MaxCompo(tmin) <= v3_MinCompo(tmax);
}

enum {
    SortOrderX = 0,
    SortOrderY = 1,
    SortOrderZ = 2,
    SortOrder_COUNT,
};

BvhNode *BvhNode_Alloc(void) {
    return malloc(sizeof(BvhNode));
}

BvhNode *BvhNode_MakeWithObjects(Object *left, Object *right) {
    BvhNode *result = BvhNode_Alloc();
    result->kind = BVHNODEKIND_OBJECTS;
    result->objects.left = left;
    result->objects.right = right;

    Aabb b1, b2;
    Aabb_MakeFromObject(left, 0, 1, &b1);
    Aabb_MakeFromObject(right, 0, 1, &b2);
    Aabb_Compound(b1, b2, &result->box);

    return result;
}

bool BvhTree_MakeFromObjects(Object *o, int o_count, float t0, float t1, BvhNode *n) {

    int sort_order = ((int)(drand48() * (double)SortOrder_COUNT)) % SortOrder_COUNT;
    switch (sort_order) {
        case SortOrderX: {
            qsort(o, o_count, sizeof(Object), Object_XCmp);
        } break;
        case SortOrderY: {
            qsort(o, o_count, sizeof(Object), Object_YCmp);
        } break;
        case SortOrderZ: {
            qsort(o, o_count, sizeof(Object), Object_ZCmp);
        } break;
    }

    Aabb bl, br;
    if (o_count == 1) { 
        n->kind = BVHNODEKIND_OBJECTS;
        n->objects.left = n->objects.right = o;

        Aabb_MakeFromObject(o, t0, t1, &bl);
        br = bl;
    } else if (o_count == 2) { 
        n->kind = BVHNODEKIND_OBJECTS;
        n->objects.left = o;
        n->objects.right = o + 1;

        if (!Aabb_MakeFromObject(o, t0, t1, &bl)) return false;
        if (!Aabb_MakeFromObject(o + 1, t0, t1, &br)) return false;
    } else {
        n->kind = BVHNODEKIND_NODES;
        n->nodes.left = BvhNode_Alloc();
        n->nodes.right = BvhNode_Alloc();

        int middle = o_count / 2;
        if (!BvhTree_MakeFromObjects(o, middle, t0, t1, n->nodes.left)) return false;
        if (!BvhTree_MakeFromObjects(o + middle, o_count - middle, t0, t1, n->nodes.right)) return false;

        bl = n->nodes.left->box;
        br = n->nodes.right->box;
    }
    Aabb_Compound(bl, br, &n->box);

    return true;
}

v3 Object_PosAtTime(Object *o, float time) {
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
                RandomInUnitSphere(),
                v3_Add(h->pos, h->norm)
                );
            
            *scattered = (ray3){.origin = h->pos, .direction = v3_Sub(target, h->pos), .time = r->time};
            *attenuation = h->mat.lambertian.albedo;
            return true;
        }
        
        case MATKIND_METAL: 
        {
            v3 reflected = Reflected(v3_Normalize(r->direction), h->norm);
            v3 random_in_sphere = v3_Mul(RandomInUnitSphere(), h->mat.metal.fuzz);
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
    
    v3 sphere_center = Object_PosAtTime(o, r.time);
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
            h->mat = o->mat;
            return true;
        }
        
        temp = (-b + sqrt(discriminant)) / a;
        if (temp >= tmin && temp <= tmax) {
            h->t = temp;
            h->pos = ray3_PointAt(r, temp);
            h->norm = v3_Normalize(v3_Mul(v3_Sub(h->pos, sphere_center), 1.f/radius));
            h->mat = o->mat;
            return true;
        }
    }
    
    return false;
}

bool hit_bvh(ray3 r, BvhNode *n, float tmin, float tmax, HitInformation *h) {
#if 1
    if (hit_aabb(n->box, r)) {
#else
    if (hit_aabb1(n->box, r, tmin, tmax)) {
#endif
        bool did_hit_left = false, did_hit_right = false;
        HitInformation left_info, right_info;
        switch (n->kind) {
            case BVHNODEKIND_NODES: {
                did_hit_left = hit_bvh(r, n->nodes.left, tmin, tmax, &left_info);
                did_hit_right = hit_bvh(r, n->nodes.right, tmin, tmax, &right_info);
            } break;
            case BVHNODEKIND_OBJECTS: {
                did_hit_left = hit_sphere(r, n->objects.left, tmin, tmax, &left_info);
                did_hit_right = hit_sphere(r, n->objects.right, tmin, tmax, &right_info);
            } break;
        }

        if (did_hit_left && did_hit_right) {
            *h = left_info.t < right_info.t ? left_info : right_info;
            return true;
        } else if (did_hit_left) {
            *h = left_info;
            return true;
        } else if (did_hit_right) { 
            *h = right_info;
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

v3 ColorAABB(ray3 r, BvhNode *tree, int depth) {

    HitInformation h;
    if (hit_bvh(r, tree, 0.001, FLT_MAX, &h)) {
        ray3 new_ray;
        v3 attenutation;
        if (depth < 50 && Scatter(&r, &h,&attenutation, &new_ray)) {
            v3 color = ColorAABB(new_ray, tree, depth + 1);
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
    
    v3 lookfrom = v3_Make(13.f, 2.f, 3.f);
    v3 lookat = v3_Make(0.f, 0.f, 0.f);
    
    float distance_to_focus = 10.f;
    float aperture = 0.1f;
    
    Camera c = Camera_Make(lookfrom, lookat, v3_Make(0.f, 1.f, 0.f), 20.f, ((float)width) / ((float) height), aperture, distance_to_focus,
        0.0f, 1.0f);
    
    Object world[501] = {0};
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
//                    center1 = v3_Add(center1, (v3){0, 0.3 * (1 + drand48()), 0});
                    m = Material_MakeLambertian((v3){ drand48() * drand48(), drand48() * drand48(), drand48() * drand48()});
                } else if (choose_mat < 0.95) { // Metal
                    m = Material_MakeMetal((v3){0.5 * ( 1 * drand48() ),  0.5 * ( 1 * drand48() ), 0.5 * ( 1 * drand48() )}, 0.5 * drand48());
                } else {
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

    BvhNode tree;
    if (!BvhTree_MakeFromObjects(world, world_size, 0.0f, 1.0f, &tree)) {
        return 1;
    }

    int rand_step = 100;
    for (int j = height - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            OUT_PRINT(" ");
            v3 color = {0};
            for (int s = 0; s < rand_step; s++) {
                float u = ((float) i + drand48()) * inv_w;
                float v = ((float) j + drand48()) * inv_h;
                
                ray3 r = Camera_GetRay(&c, u, v);
                
                v3_Added(&color, ColorAABB(r, &tree, 0)); // get color mapped from 0.0 - 1.0
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

