Material Material_MakeLambertian(v3 albedo) {
    return (Material){.kind = MATKIND_LAMBERTIAN, .lambertian = {albedo}};
}

Material Material_MakeMetal(v3 albedo, float fuzz) {
    fuzz = CLAMP_MAX(fuzz, 1.0);
    return (Material){.kind = MATKIND_METAL, .metal = {albedo, fuzz}};
}

Material Material_MakeDieletric(float ref_idx) {
    return (Material){.kind = MATKIND_DIELETRIC, .dieletric = {ref_idx}};
}


