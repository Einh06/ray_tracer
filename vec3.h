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
