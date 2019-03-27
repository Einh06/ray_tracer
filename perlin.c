#define PERLIN_PERM_COUNT 256

static float ranfloat[PERLIN_PERM_COUNT];
static int perm_x[PERLIN_PERM_COUNT];
static int perm_y[PERLIN_PERM_COUNT];
static int perm_z[PERLIN_PERM_COUNT];

inline float trilinear_interp(float c[2][2][2], float u, float v, float w) {
    float result = 0;
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                result += ((i*u + (1-i)*(1-u))*
                           (j*v + (1-j)*(1-v))*
                           (k*w + (1-k)*(1-w)) ) * c[i][j][k];
            }
        }
    }
    return result;
}

void Perlin_Init(void) {
    // generate rand float
    for (int i = 0; i < PERLIN_PERM_COUNT; i++) {
        ranfloat[i] = drand48();
    }

    // gernerate perm for x, y z
    int *perms[3] = {perm_x, perm_y, perm_z};
    for (int p = 0; p < 3; p++) {
        int *perm = perms[p];
        for (int i = 0; i < PERLIN_PERM_COUNT; i++) {
            perm[i] = i;
        }
        for (int i = 0; i < PERLIN_PERM_COUNT; i++) {
            int target = (int)(drand48() * (i + 1));
            int tmp = perm[i];
            perm[i] = perm[target];
            perm[target] = tmp;
        }
    }
}

float Perlin3(v3 p) {
    float u = p.x - floor(p.x);
    float v = p.y - floor(p.y);
    float w = p.z - floor(p.z);
    int i = ((int)(4.f * p.x)) & 255;
    int j = ((int)(4.f * p.y)) & 255;
    int k = ((int)(4.f * p.z)) & 255;
    return ranfloat[perm_x[i] ^ perm_y[j] ^ perm_z[k]];
}

float Perlin3_Interp(v3 p) {
    float u = p.x - floor(p.x);
    float v = p.y - floor(p.y);
    float w = p.z - floor(p.z);
    u = u * u * (3 - 2*u);
    v = v * v * (3 - 2*v);
    w = w * w * (3 - 2*w);
    int i = floor(p.x);
    int j = floor(p.y);
    int k = floor(p.z);
    float c[2][2][2];
    for (int di = 0; di < 2; di++) {
        for (int dj = 0; dj < 2; dj++) {
            for (int dk = 0; dk < 2; dk++) {
                c[di][dj][dk] = ranfloat[perm_x[(i+di) & 255] ^ perm_y[(j+dj) & 255] ^ perm_z[(k+dk) & 255]];
            }
        }
    }
    return trilinear_interp(c, u, v, w);
}

float Turbulence(v3 p, int depth) {
    float acc = 0;
    v3 temp_p = p;
    float weight = 1.0;
    for (int i = 0; i < depth; ++i) {
        acc += weight * Perlin3_Interp(temp_p); 
        weight *= 0.5;
        v3_Scaled(&temp_p, 2);
    }
    return fabs(acc);
}
