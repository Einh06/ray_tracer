typedef struct Camera {
    v3 origin;
    v3 lower_left;
    v3 horizontal;
    v3 vertical;
    v3 u, v, w;
    float lens_radius;
    float time0, time1;
} Camera;
