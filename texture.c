Texture Texture_MakeColor(v3 c) {
    return (Texture){ .kind = TextureKind_Color, .color = c};
}

Texture Texture_MakeChecker(Texture *even, Texture *odd) {
    return (Texture){ .kind = TextureKind_Checker, .checker = {even, odd}};
}

v3 Texture_ColorAt(Texture *t, float u, float v, v3 *p) {
    switch (t->kind) {
        case TextureKind_Color: {
            return t->color;
        } break;
        case TextureKind_Checker: {
            float sines = sin(10 * p->x) * sin(10 * p->y) * sin(10 * p->z); 
            if (sines < 0) {
                return Texture_ColorAt(t->checker.odd, u, v, p);
            } else {
                return Texture_ColorAt(t->checker.even, u, v, p);
            }
        } break;
        default: {
            assert(false);
        } break;
    }
    return (v3){0};
}
