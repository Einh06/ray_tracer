Material Material_MakeLambertian(Texture texture) {
    return (Material){.kind = MATKIND_LAMBERTIAN, .lambertian = {texture}};
}

Material Material_MakeMetal(v3 albedo, float fuzz) {
    fuzz = CLAMP_MAX(fuzz, 1.0);
    return (Material){.kind = MATKIND_METAL, .metal = {albedo, fuzz}};
}

Material Material_MakeDieletric(float ref_idx) {
    return (Material){.kind = MATKIND_DIELETRIC, .dieletric = {ref_idx}};
}

