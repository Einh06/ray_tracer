#define ArrayCount(a) (sizeof(a) / sizeof(a[0]))
#define MIN(a,b) (a) < (b) ? (a) : (b)
#define MAX(a,b) (a) > (b) ? (a) : (b)
#define CLAMP_MAX(a,b) MIN(a,b)
#define CLAMP_MIN(a,b) MAX(a,b)
double drand48() {
    return ((double)rand()/(RAND_MAX+1.0f));
}
