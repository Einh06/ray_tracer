inline v3 ray3_PointAt(ray3 r, float t) {
    return v3_Add(r.origin, v3_Mul(r.direction, t));
}
