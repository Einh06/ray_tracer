#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

#define OUT_PRINT(str, ...) (fprintf(output, str, __VA_ARGS__))
int main(int argc, char** argv) {
    FILE* output = fopen("result.ppm", "w");
    if (output == NULL) {
        goto end;
    }
    
    int width = 400;
    int height = 200;
    
    OUT_PRINT("P3\n%d %d\n255\n", width, height);
    
    float inv_w = 1.0f / ((float) width);
    float inv_h = 1.0f / ((float) height);
    
    for (int j = height - 1; j >= 0; j--) {
        for (int i = 0; i < width; i++) {
            
            OUT_PRINT(" ");
            v3 color = {.r = ((float) i) * inv_w, .g = ((float) j) * inv_h, .b = 0.2};
            v3_Scaled(&color, 255.99);
            
            OUT_PRINT("%d %d %d", (int)color.r, (int)color.g, (int)color.b);
        }
        OUT_PRINT("\n");
    }
    end:
    return 0;
}

