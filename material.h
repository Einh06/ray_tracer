typedef enum MaterialKind {
    MATKIND_LAMBERTIAN,
    MATKIND_METAL,
    MATKIND_DIELETRIC,
    MATKIND_COUNT,
} MaterialKind;

typedef struct Material {
    MaterialKind kind;
    union {
        struct {
            Texture albedo;
        } lambertian;
        
        struct {
            v3 albedo;
            float fuzz;
        } metal;
        
        struct { 
            float ref_idx;
        } dieletric;
    };
} Material;

