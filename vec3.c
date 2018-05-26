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

inline v3 v3_Neg(v3 v) {
    return (v3){-v.x, -v.y, -v.z};
}

inline v3 v3_Mul(v3 v, float s) {
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


