#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <time.h>

#define ArrayCount(a) (sizeof(a) / sizeof(a[0]))
#define MIN(a,b) (a) < (b) ? (a) : (b)
double drand48() {
    return ((double)rand()/(RAND_MAX+1));
}

typedef struct v3 {
    union {
        float v[3];
        struct {
            float x;
            float y;
            float z;
        };
        struct {
            float r;
            float g;
            float b;
        };
    };
} v3;

inline v3 v3_Make(float x, float y, float z) {
    return (v3){.x = x, .y = y, .z = z};
}

inline float v3_Length( v3 v) {
    return sqrt( v.x * v.x + v.y * v.y + v.z * v.z);
}

inline float v3_SquaredLength(v3 v) {
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

inline v3 v3_Normalize(v3 v) {
    float inv_l = 1.0f / v3_Length(v);
    return (v3) {.x = v.x * inv_l, .y = v.y * inv_l, .z = v.z * inv_l};
}

inline void v3_Normalized(v3* v) {
    float inv_l = 1.0f / v3_Length(*v);
    v->x *= inv_l;
    v->y *= inv_l;
    v->z *= inv_l;
}

inline v3 v3_Add(v3 a, v3 b) {
    return (v3){.x = a.x + b.x, .y = a.y + b.y, .z = a.z + b.z };
}

inline void v3_Added(v3* v, v3 a) {
    v->x += a.x;
    v->y += a.y;
    v->z += a.z;
}

inline v3 v3_Sub(v3 a, v3 b) {
    return (v3){.x = a.x - b.x, .y = a.y - b.y, .z = a.z - b.z };
}

inline v3 v3_Subed(v3* v, v3 a) {
    v->x -= a.x;
    v->y -= a.y;
    v->z -= a.z;
}

inline v3 v3_Scale(v3 v, float s) {
    return (v3){.x = v.x * s, .y = v.y * s, .z = v.z * s };
}

inline v3_Scaled(v3* v, float s) {
    v->x *= s;
    v->y *= s;
    v->z *= s;
}

inline float v3_Dot(v3 a, v3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline v3 v3_Cross(v3 a, v3 b) {
    return (v3){ .x = a.y * b.z - a.z * b.y, .y = a.z * b.x - a.x * b.z, .z = a.x * b.y - a.y * b.x};
}

typedef struct ray3 {
    v3 origin;
    v3 direction;
} ray3;

inline v3 ray3_PointAt(ray3 r, float t) {
    return v3_Add(r.origin, v3_Scale(r.direction, t));
}

typedef struct Camera {
    v3 origin;
    v3 lower_left;
    v3 horizontal;
    v3 vertical;
} Camera;

ray3 Camera_GetRay(Camera *c, float u, float v) {
    return  (ray3) {
        .origin = c->origin,
        .direction = v3_Add(c->lower_left,
                            v3_Add(
            v3_Scale(c->horizontal, u),
            v3_Scale(c->vertical, v)
            )
                            )
    };
}

v3 random_in_unit_sphere(void) {
    v3 p;
    do {
        p = v3_Sub(
            v3_Scale((v3){drand48(), drand48(), drand48()}, 2.0),
            (v3){1.0, 1.0, 1.0}
            );
    } while (v3_SquaredLength(p) >= 1.0);
    return p;
}

typedef enum MaterialKind {
    MATKIND_LAMBERTIAN,
    MATKIND_METAL,
    MATKIND_COUNT,
} MaterialKind;

typedef struct Material {
    MaterialKind kind;
    union {
        struct {
            v3 albedo;
        } lambertian;
        
        struct {
            v3 albedo;
            float fuzz;
        } metal;
    };
} Material;

Material Material_MakeLambertian(v3 albedo) {
    return (Material){.kind = MATKIND_LAMBERTIAN, .lambertian = {albedo}};
}

Material Material_MakeMetal(v3 albedo, float fuzz) {
    fuzz = MIN(fuzz, 1.0);
    return (Material){.kind = MATKIND_METAL, .metal = {albedo, fuzz}};
}

typedef enum ObjectKind {
    OBJECTKIND_SPHERE,
}ObjectKind;

typedef struct Object {
    
    ObjectKind kind;
    v3 pos;
    union {
        struct {
            float radius;
        } sphere;
    };
    Material mat;
} Object;

typedef struct HitInformation {
    float t;
    v3 pos;
    v3 norm;
    Material mat;
} HitInformation;

bool Scatter(ray3* r, HitInformation* h, v3* attenuation, ray3* scattered) {
    switch(h->mat.kind) {
        case MATKIND_LAMBERTIAN:
        {
            v3 target = v3_Add(
                random_in_unit_sphere(),
                v3_Add(h->pos, h->norm)
                );
            
            *scattered = (ray3){.origin = h->pos, .direction = v3_Sub(target, h->pos)};
            *attenuation = h->mat.lambertian.albedo;
            return true;
        }
        case MATKIND_METAL: 
        {
            v3 reflected = v3_Sub(r->direction,
                                  v3_Scale(h->norm, 2.0f * v3_Dot(r->direction, h->norm))
                                  );
            v3 random_in_sphere = v3_Scale(random_in_unit_sphere(), h->mat.metal.fuzz);
            *scattered = (ray3){.origin = h->pos, .direction = v3_Add(reflected, random_in_sphere)};
            *attenuation = h->mat.metal.albedo;
            return v3_Dot(scattered->direction, h->norm) > 0.f;
        }
        default:
        exit(1);
    }
    return false;
}


bool hit_sphere(ray3 r, Object *o, float tmin, float tmax, HitInformation *h) {
    assert(o->kind == OBJECTKIND_SPHERE);
    
    v3 sphere_center = o->pos;
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
            h->norm = v3_Normalize( v3_Sub(h->pos, sphere_center) );
            return true;
        }
        
        temp = (-b + sqrt(discriminant)) / a;
        if (temp >= tmin && temp <= tmax) {
            h->t = temp;
            h->pos = ray3_PointAt(r, temp);
            h->norm = v3_Normalize( v3_Sub(h->pos, sphere_center) );
            return true;
        }
    }
    
    return false;
}

v3 Color(ray3 r, Object *objects, int object_count, int depth) {
    
    HitInformation h;
    
    for (Object *o = objects; o < objects + object_count; o++){
        switch (o->kind) {
            case OBJECTKIND_SPHERE: {
                
                if (!hit_sphere(r, o, 0.001, FLT_MAX,  &h)) { continue; }
                h.mat = o->mat;
                ray3 new_ray;
                v3 attenutation;
                if (depth < 50 && Scatter(&r, &h,&attenutation, &new_ray)) {
                    v3 color = Color(new_ray, objects, object_count, depth + 1);
                    return (v3){color.r * attenutation.r, color.g * attenutation.g, color.b * attenutation.b};
                } else {
                    return (v3){ 0.f, 0.f, 0.f};
                }
            }
            
            default: {
                perror("Unrecognized object type");
                exit(1);
            }
        }
    }
    
    v3 unit_direction = v3_Normalize(r.direction);
    float t = (unit_direction.y + 1.f) * 0.5f;
    return v3_Add(
        v3_Scale(v3_Make(1.f, 1.f, 1.f), (1.f-t)),
        v3_Scale(v3_Make(0.5f, 0.7f, 1.f), t)
        );
}

#define OUT_PRINT(str, ...) (fprintf(output, str, __VA_ARGS__))
int main(int argc, char** argv) {
    
    srand(time(NULL));
    FILE* output = fopen("output/result.ppm", "w");
    if (output == NULL) {
        goto end;
    }
    
    int width = 400;
    int height = 200;
    int rand_step = 100;
    
    OUT_PRINT("P3\n%d %d\n255\n", width, height);
    
    Camera c = {
        .origin = (v3){.x = 0.0f, .y = 0.0f, .z = 0.0f},
        .lower_left = (v3){.x = -4.0f, .y = -2.0f, .z = -1.0f},
        .horizontal = (v3){.x = 8.0f, .y = 0.0f, .z = 0.0f},
        .vertical = (v3){.x = 0.0f, .y = 4.0f, .z = 0.0f},
    };
    
    float inv_w = 1.0f / ((float)width);
    float inv_h = 1.0f / ((float)height);
    
    Object world[] = {
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(0.f, 0.f, -1.f), 
            .sphere = {.radius = 0.5f},
            .mat = Material_MakeLambertian(v3_Make(0.8f, 0.3f, 0.3f)),
        },
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(1.f, 0.f, -1.f), 
            .sphere = {.radius = 0.5f},
            .mat = Material_MakeMetal(v3_Make(0.8f, 0.6f, 0.2f), 1.0f),
        },
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(-1.f, 0.f, -1.f), 
            .sphere = {.radius = 0.5f},
            .mat = Material_MakeMetal(v3_Make(0.8f, 0.8f, 0.8f), 0.3f),
        },
        {
            .kind = OBJECTKIND_SPHERE, 
            .pos = v3_Make(0.f, -100.5f, -1.f), 
            .sphere = {.radius = 100.0f},
            .mat = Material_MakeLambertian(v3_Make(0.8f, 0.8f, 0.0f)),
        },
    };
    
    for (int j = height - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            OUT_PRINT(" ");
            v3 color = {0};
            for (int s = 0; s < rand_step; s++) {
                float u = ((float) i + ((double)rand()/(RAND_MAX+1))) * inv_w;
                float v = ((float) j + ((double)rand()/(RAND_MAX+1))) * inv_h;
                
                ray3 r = Camera_GetRay(&c, u, v);
                
                v3_Added(&color, Color(r, world, ArrayCount(world), 0)); // get color mapped from 0.0 - 1.0
            }
            v3_Scaled(&color, (1.f/((float)rand_step)));
            color = (v3){sqrt(color.r), sqrt(color.g), sqrt(color.b)};
            color = v3_Scale(color, 255.99); // map colors from 0.0 - 1.0 to 0 - 255
            OUT_PRINT("%d %d %d", (int)color.r, (int)color.g, (int)color.b);
        }
        OUT_PRINT("\n");
    }
    end:
    return 0;
}

