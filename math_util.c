v3 RandomInUnitSphere(void) {
    v3 p;
    do {
        p = v3_Sub(
            v3_Mul((v3){drand48(), drand48(), drand48()}, 2.0),
            (v3){1.0, 1.0, 1.0}
            );
    } while (v3_SquaredLength(p) >= 1.0);
    return p;
}

v3 RandomInUnitDisk(void) {
    v3 p;
    do {
        p = v3_Sub(
            v3_Mul((v3){drand48(), drand48(), 0.f}, 2.f), 
            (v3){1.f, 1.f, 0.f});
    } while( v3_Dot(p, p) >= 1.f);
    return p;
}
