typedef enum TextureKind {
    TextureKind_Color = 0,
    TextureKind_Checker,
} TextureKind;

typedef struct Texture {
    TextureKind kind;
    union {
        v3 color;
        struct {
            struct Texture *even, *odd;
        } checker;
    };
} Texture;
