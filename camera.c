Camera Camera_Make(v3 lookFrom, v3 lookAt, v3 vup, float vfov, float aspect, float aperture, float focus_dist, float time0, float time1) {
    float lens_radius = aperture / 2.f;
    
    float theta = vfov*M_PI/180.f;
    float half_height = tan(theta/2.f);
    float half_width = aspect * half_height;
    
    v3 w = v3_Normalize(v3_Sub(lookFrom, lookAt));
    v3 u = v3_Normalize(v3_Cross(vup, w));
    v3 v = v3_Cross(w, u);
    
    v3 horizontal = v3_Mul(u, 2*half_width*focus_dist);
    v3 vertical = v3_Mul(v, 2*half_height*focus_dist);
    //((origin - u*half_width) - v * half_height) - w
    v3 lower_left =
        v3_Sub(v3_Sub(v3_Sub(
        lookFrom, 
        v3_Mul(u, half_width * focus_dist)), 
                      v3_Mul(v, half_height * focus_dist)), 
               v3_Mul(w, focus_dist));
    
    return (Camera){
        .origin     = lookFrom,
        .horizontal = horizontal,
        .vertical   = vertical,
        .lower_left = lower_left,
        .u = u, .v = v, .w = w,
        .lens_radius = lens_radius,
        .time0 = time0, .time1 = time1,
    };
}

ray3 Camera_GetRay(Camera *c, float s, float t) {
    v3 rd = v3_Mul(RandomInUnitDisk(), c->lens_radius);
    v3 offset = v3_Add(
        v3_Mul(c->u, rd.x),
        v3_Mul(c->v, rd.y));
    
    v3 h = v3_Mul(c->horizontal, s);
    v3 v = v3_Mul(c->vertical, t);
    v3 o = v3_Add(c->origin, offset);
    
    float time = c->time0 + (drand48() * (c->time1 - c->time0));

    v3 direction = v3_Sub(v3_Add(c->lower_left, v3_Add(h, v)), o);
    return  (ray3) {
        .origin = o,
        .direction = direction,
        .time = time,
    };
}
