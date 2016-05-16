/* An implementation of the Histogram Sketch data structure described in
 * Ben-Haim and Tom-Tov's "A Streaming Parallel Decision Tree Algorithm" in
 * Journal of Machine Learning Research 11
 * (http://www.jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf).
 *
 * This implementation is mostly a port of github.com/aaw/histosketch to C as
 * a Redis module.
 *
 * The histogram sketch consists of a fixed number of (value, count) centroids
 * sorted by increasing value. Each time a value V is added to the sketch, it's
 * either merged into an existing centroid with the same value or added as a
 * (V, 1) centroid and then the two centroids in the sketch that have values
 * closest together are merged. Quantiles and counts are estimated by finding
 * the two centroids (V1, C1) and (V2, C2) in the sketch bordering the desired
 * value and using the area under the trapezoid defined by (V1, 0), (V1, C1),
 * (V2, 0), (V2, C2) to guide the estimate.
 */

#include "float.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

#include "redismodule.h"

#define HISTK_MODULE_VERSION 1
#define HISTK_ENCODING_VERSION 0
#define HISTK_DEFAULT_NUM_CENTROIDS 64
#define HISTK_MAX_NUM_CENTROIDS 2048
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
#define HISTK_STR_MAX_CENTROIDS STR(HISTK_MAX_NUM_CENTROIDS)
#define HISTK_ERRORMSG_COUNTNOTINT    "ERR count is not an integer."
#define HISTK_ERRORMSG_VALUENOTDOUBLE "ERR value is not a double."
#define HISTK_ERRORMSG_BADEMPTYKEYSET "ERR problem setting value for empty key."
#define HISTK_ERRORMSG_BADQUANTILE    "ERR argument must be in the range " \
                                      "[0.0, 1.0]."
#define HISTK_ERRORMSG_EMPTYSKETCH    "ERR empty histogram."
#define HISTK_ERRORMSG_CENTROIDLIMIT  "ERR invalid size: number of centroids " \
                                      "must be at most " \
                                      HISTK_STR_MAX_CENTROIDS "."
#define UNUSED(x) (void)(x)

static RedisModuleType *HistKType;

struct Centroid {
    double value;
    long long count;
};

struct HistK {
    // Array of centroids, sorted by value.
    struct Centroid *cs;
    // Total number of values observed by the sketch.
    unsigned long long totalCount;
    // Minimum value ever added to the sketch..
    double min;
    // Maximum value ever added to the sketch.
    double max;
    // Current number of centroids in this sketch..
    unsigned short int numCentroids;
    // Maximum number of centroids allowed in this sketch.
    unsigned short int maxCentroids;
};

struct HistK *createHistK(unsigned short int maxCentroids) {
    struct HistK *h;
    h = RedisModule_Alloc(sizeof(*h));
    h->totalCount = 0;
    h->numCentroids = 0;
    h->min = DBL_MAX;
    h->max = DBL_MIN;
    // Allocate one more centroid than we need as a workspace for adding new
    // values to the sketch. We'll add the centroid as a singleton then merge
    // the two closest centroids.
    h->cs = RedisModule_Alloc((maxCentroids + 1) * sizeof(struct Centroid));
    h->maxCentroids = maxCentroids;
    return h;
}

void freeHistK(struct HistK *o) {
    RedisModule_Free(o->cs);
    RedisModule_Free(o);
}

// Merge the centroid at h->cs[i+1] into the centroid at h->cs[i].
inline void mergeCentroidWithNext(struct HistK* h, unsigned int i) {
    long long s = h->cs[i].count + h->cs[i+1].count;
    h->cs[i].value = ((h->cs[i].value * h->cs[i].count) +
                      (h->cs[i+1].value * h->cs[i+1].count)) / s;
    h->cs[i].count = s;
}

// Add <count> <value>s to the sketch.
void add(struct HistK *h, double value, unsigned long long count) {
    if (value < h->min) { h->min = value; }
    if (value > h->max) { h->max = value; }
    // Find the index k in the sorted list of centroids where (value, count)
    // belongs.
    int i = h->numCentroids - 1;
    unsigned char exact_match = 0;
    for (; i >= 0; i--) {
        if (h->cs[i].value == value) {
            exact_match = 1;
	    break;
	} else if (h->cs[i].value < value) {
            break;
        }
        h->cs[i+1] = h->cs[i];
    }
    h->cs[i+1].value = value;
    h->cs[i+1].count = count;
    h->numCentroids++;
    h->totalCount += count;

    if (h->numCentroids <= h->maxCentroids && !exact_match) {
        return;
    }

    // Find index where |h->cs[i].value - h->cs[i+1].value| is minimized.
    unsigned int mi = h->numCentroids - 1;
    double md = DBL_MAX;
    unsigned int n = 1;
    for(int i = 0; i < h->numCentroids - 1; i++) {
      double d = fabs(h->cs[i+1].value - h->cs[i].value);
      if (d < md || (d == md && rand() * n++ < 1.0)) {
	mi = i;
	md = d;
      }
    }

    // Merge centroids in h->cs[i] and h->cs[i+1].
    mergeCentroidWithNext(h, mi);
    h->numCentroids--;
    for(int i = mi + 1; i < h->numCentroids; i++) {
      h->cs[i] = h->cs[i+1];
    }
}

// Populate ci and cj with the two centroids in h before and after index i. If
// i is 0, we populate ci with a dummy centroid holding the min value observed
// by the sketch. If i is h->numCentroids, we populate cj with a dummy centroid
// holding the max value observed by the sketch. Using dummy centroids like
// this makes the quantiles and counts at extreme values more accurate.
void getBorderingCentroids(const struct HistK *h, unsigned int i,
                           struct Centroid *ci, struct Centroid *cj) {
    if (i == 0) {
        ci->value = h->min;
        ci->count = 0;
        *cj = h->cs[0];
    } else if (i == h->numCentroids) {
        *ci = h->cs[h->numCentroids-1];
        cj->value = h->max;
        cj->count = 0;
    } else {
        *ci = h->cs[i-1];
        *cj = h->cs[i];
    }
}

// Return an estimate of the smallest value V observed by the sketch such that
// q * h->totalCount elements observed were less than or equal to V. Only q
// values in the range [0.0, 1.0] are valid arguments.
double quantile(const struct HistK *h, double q) {
    double t = q * h->totalCount;
    int i = 0;
    double s = 0.0;
    double pv = 0.0;
    for(; i < h->numCentroids; i++) {
        double v = h->cs[i].count / 2.0;
        if (s + v + pv > t) {
            break;
        }
        s += v + pv;
        pv = v;
    }

    struct Centroid ci, cj;
    getBorderingCentroids(h, i, &ci, &cj);

    // Solve for u such that
    // t - s = (ci.count + mu) / 2 * (u-ci.value) / (cj.value - ci.value), where
    // mu = ci.count + (u-ci.value) * (cj.count-ci.count) / (cj.value-ci.value).
    // You can solve for such a u using the quadratic equation as long as
    // ci.count != cj.count, in which case the solution is a little bit simpler.
    // See Algorithm 4 and its description in the Ben-Haim/Tom-Tov paper.
    double d = t - s;
    double a = cj.count - ci.count;
    if (a == 0.0) {
        return ci.value + (cj.value - ci.value) * (d / ci.count);
    }
    double b = 2.0 * ci.count;
    double c = -2.0 * d;
    double z = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);
    return ci.value + (cj.value - ci.value) * z;
}

// Return an estimate of the number of values observed by the sketch that are
// less than or equal to <v>. This is essentially the "Sum" procedure described
// in Ben-Haim and Tom-Tov's paper.
long long countLessThan(const struct HistK *h, double v) {
    if (v >= h->max) {
        return h->totalCount;
    } else if (v < h->min) {
        return 0;
    }
    int i = h->numCentroids - 1;
    for (; i >= 0; i--) {
        if (h->cs[i].value <= v) {
            break;
        }
    }

    struct Centroid ci, cj;
    getBorderingCentroids(h, i + 1, &ci, &cj);

    double s = 0;
    for (int j = 0; j < i; j++) {
        s += h->cs[j].count;
    }
    double x = (v - ci.value) / (cj.value - ci.value);
    double b = ci.count + (cj.count - ci.count) * x;
    double est = s + ci.count / 2.0 + (ci.count + b) * x / 2.0;
    return (long long)round(est);
}

/* HISTK.ADD <KEY> <VALUE1> [<COUNT1>] [<VALUE2> <COUNT2>, ...]
   Add values to the sketch. Returns the total number of values observed by the
   sketch.
*/
int AddCommand(RedisModuleCtx *ctx, RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);

    if (argc < 3) return RedisModule_WrongArity(ctx);
    RedisModuleKey *key = RedisModule_OpenKey(
        ctx, argv[1], REDISMODULE_READ|REDISMODULE_WRITE);
    int keytype = RedisModule_KeyType(key);
    if (keytype != REDISMODULE_KEYTYPE_EMPTY &&
        RedisModule_ModuleTypeGetType(key) != HistKType) {
        return RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
    }

    struct HistK *h;
    if (keytype == REDISMODULE_KEYTYPE_EMPTY) {
        h = createHistK(HISTK_DEFAULT_NUM_CENTROIDS);
        RedisModule_ModuleTypeSetValue(key, HistKType, h);
    } else {
        h = RedisModule_ModuleTypeGetValue(key);
    }

    for (int iarg = 2; iarg < argc;) {
        double value;
        if (RedisModule_StringToDouble(argv[iarg++], &value) !=
            REDISMODULE_OK) {
            return RedisModule_ReplyWithError(ctx,
                                              HISTK_ERRORMSG_VALUENOTDOUBLE);
        }
        long long count = 1;
        if (argc > iarg && RedisModule_StringToLongLong(argv[iarg++], &count) !=
            REDISMODULE_OK) {
            return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_COUNTNOTINT);
        }
        add(h, value, count);
    }

    RedisModule_ReplyWithLongLong(ctx, h->totalCount);
    RedisModule_ReplicateVerbatim(ctx);
    return REDISMODULE_OK;
}

/* HISTK.QUANTILE <KEY> <Q>
   Returns the q-quantile for any 0.0 <= Q <= 1.0
 */
int QuantileCommand(RedisModuleCtx *ctx, RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);
    if (argc != 3) return RedisModule_WrongArity(ctx);
    double q;
    if (RedisModule_StringToDouble(argv[2], &q) == REDISMODULE_ERR) {
        return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_VALUENOTDOUBLE);
    }
    if (q < 0.0 || q > 1.0) {
        return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_BADQUANTILE);
    }

    RedisModuleKey *key = RedisModule_OpenKey(ctx,argv[1],
                                              REDISMODULE_READ);
    int keytype = RedisModule_KeyType(key);
    if (keytype != REDISMODULE_KEYTYPE_EMPTY &&
        RedisModule_ModuleTypeGetType(key) != HistKType) {
        return RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
    }

    if (keytype == REDISMODULE_KEYTYPE_EMPTY) {
        return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_EMPTYSKETCH);
    }
    struct HistK *h = RedisModule_ModuleTypeGetValue(key);
    return RedisModule_ReplyWithDouble(ctx, quantile(h, q));
}

/* HISTK.COUNT <KEY> [<V>]
   Returns an estimate of of the number of values <= V in the sketch. If V is
   omitted, returns the total number of values observed by the sketch so far.
*/
int CountCommand(RedisModuleCtx *ctx, RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);
    if (argc < 2 || argc > 3) return RedisModule_WrongArity(ctx);
    RedisModuleKey *key = RedisModule_OpenKey(ctx,argv[1],
                                              REDISMODULE_READ);
    int keytype = RedisModule_KeyType(key);
    if (keytype != REDISMODULE_KEYTYPE_EMPTY &&
        RedisModule_ModuleTypeGetType(key) != HistKType) {
        return RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
    }

    if (keytype == REDISMODULE_KEYTYPE_EMPTY) {
        return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_EMPTYSKETCH);
    }
    struct HistK *h = RedisModule_ModuleTypeGetValue(key);
    if (argc == 2) {
        return RedisModule_ReplyWithLongLong(ctx, h->totalCount);
    }
    double v;
    if (RedisModule_StringToDouble(argv[2], &v) == REDISMODULE_ERR) {
        return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_VALUENOTDOUBLE);
    }
    return RedisModule_ReplyWithLongLong(ctx, countLessThan(h, v));
}

/* HISTK.MERGESTORE <KEY> HIST1 [HIST2] ... [HISTN]
   Merge HIST1 ... HISTN, store results in KEY. If there's already a histogram
   sketch in KEY before this command is called, the results are merged into that
   sketch.
*/
int MergeStoreCommand(RedisModuleCtx *ctx, RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);
    if (argc < 3) return RedisModule_WrongArity(ctx);

    RedisModuleKey *key = RedisModule_OpenKey(
        ctx, argv[1], REDISMODULE_READ|REDISMODULE_WRITE);
    int keytype = RedisModule_KeyType(key);
    if (keytype != REDISMODULE_KEYTYPE_EMPTY &&
        RedisModule_ModuleTypeGetType(key) != HistKType) {
        return RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
    }

    struct HistK *h;
    if (keytype == REDISMODULE_KEYTYPE_EMPTY) {
        h = createHistK(HISTK_DEFAULT_NUM_CENTROIDS);
        RedisModule_ModuleTypeSetValue(key, HistKType, h);
    } else {
        h = RedisModule_ModuleTypeGetValue(key);
    }

    for (int iarg = 2; iarg < argc; iarg++) {
        RedisModuleKey *akey = RedisModule_OpenKey(ctx, argv[iarg],
                                                   REDISMODULE_READ);
        int akeytype = RedisModule_KeyType(akey);
        if (akeytype == REDISMODULE_KEYTYPE_EMPTY) {
            continue;
        } else if (akeytype != REDISMODULE_KEYTYPE_EMPTY &&
            RedisModule_ModuleTypeGetType(akey) != HistKType) {
            return RedisModule_ReplyWithError(ctx,
                                              REDISMODULE_ERRORMSG_WRONGTYPE);
        }
        struct HistK *ah = RedisModule_ModuleTypeGetValue(akey);
        for (int i = 0; i < ah->numCentroids; i++) {
            add(h, ah->cs[i].value, ah->cs[i].count);
        }
    }

    RedisModule_ReplyWithLongLong(ctx, h->totalCount);
    RedisModule_ReplicateVerbatim(ctx);
    return REDISMODULE_OK;
}

/* HISTK.RESIZE <KEY> <CENTROIDS>
   Resize the sketch to a <CENTROIDS> centroids. In most cases, this should be
   called once before values are added to the sketch. If called on an existing
   sketch with a smaller number of centroids, some centroids will be merged.
*/
int ResizeCommand(RedisModuleCtx *ctx, RedisModuleString **argv, int argc) {
    RedisModule_AutoMemory(ctx);
    if (argc < 3) return RedisModule_WrongArity(ctx);
    long long newSize;
    if (RedisModule_StringToLongLong(argv[2], &newSize) != REDISMODULE_OK) {
        return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_COUNTNOTINT);
    }
    if (newSize > HISTK_MAX_NUM_CENTROIDS) {
        return RedisModule_ReplyWithError(ctx, HISTK_ERRORMSG_CENTROIDLIMIT);
    }

    RedisModuleKey *key = RedisModule_OpenKey(
        ctx, argv[1], REDISMODULE_READ|REDISMODULE_WRITE);
    int keytype = RedisModule_KeyType(key);
    if (keytype != REDISMODULE_KEYTYPE_EMPTY &&
        RedisModule_ModuleTypeGetType(key) != HistKType) {
        return RedisModule_ReplyWithError(ctx,REDISMODULE_ERRORMSG_WRONGTYPE);
    }
    struct HistK *h = createHistK(newSize);
    if (keytype != REDISMODULE_KEYTYPE_EMPTY) {
        struct HistK *oldh = RedisModule_ModuleTypeGetValue(key);
        for (int i = 0; i < oldh->numCentroids; i++) {
            add(h, oldh->cs[i].value, oldh->cs[i].count);
        }
    }
    RedisModule_ModuleTypeSetValue(key, HistKType, h);
    RedisModule_ReplyWithLongLong(ctx, newSize);
    RedisModule_ReplicateVerbatim(ctx);
    return REDISMODULE_OK;
}

void *HistKRdbLoad(RedisModuleIO *rdb, int encver) {
    if (encver > HISTK_ENCODING_VERSION) {
    // TODO: Use RedisModule_Log to log a warning if/when RedisModule_Log exists
        return NULL;
    }
    unsigned int maxCentroids = RedisModule_LoadUnsigned(rdb);
    struct HistK *h = createHistK(maxCentroids);
    h->numCentroids = RedisModule_LoadUnsigned(rdb);
    h->totalCount = RedisModule_LoadUnsigned(rdb);
    for(int i = 0; i < h->numCentroids; i++) {
        h->cs[i].value = RedisModule_LoadDouble(rdb);
        h->cs[i].count = RedisModule_LoadSigned(rdb);
    }
    h->min = RedisModule_LoadDouble(rdb);
    h->max = RedisModule_LoadDouble(rdb);
    return h;
}

void HistKRdbSave(RedisModuleIO *rdb, void *value) {
    struct HistK *h = value;
    RedisModule_SaveUnsigned(rdb, h->maxCentroids);
    RedisModule_SaveUnsigned(rdb, h->numCentroids);
    RedisModule_SaveUnsigned(rdb, h->totalCount);
    for (int i = 0; i < h->numCentroids; i++) {
        RedisModule_SaveDouble(rdb, h->cs[i].value);
        RedisModule_SaveSigned(rdb, h->cs[i].count);
    }
    RedisModule_SaveDouble(rdb, h->min);
    RedisModule_SaveDouble(rdb, h->max);
}

void HistKAofRewrite(RedisModuleIO *aof, RedisModuleString *key, void *value) {
  struct HistK *h = value;
  // Run a RESIZE first to ensure the regenerated histk is the right size.
  RedisModule_EmitAOF(aof, "HISTK.RESIZE", "sl", key, h->numCentroids);
  // http://stackoverflow.com/a/1701272
  char buf[3 + DBL_MANT_DIG - DBL_MIN_EXP];
  for (unsigned int i = 0; i < h->numCentroids; i++) {
    int nbuf = snprintf(buf, sizeof(buf), "%f", h->cs[i].value);
    RedisModule_EmitAOF(aof, "HISTK.ADD", "sbl", key, buf, nbuf,
                        h->cs[i].count);
  }
}

void HistKDigest(RedisModuleDigest *digest, void *value) {
  /* TODO: The DIGEST module interface is yet not implemented. */
    UNUSED(digest);
    UNUSED(value);
}

void HistKFree(void *value) {
    freeHistK(value);
}

/* Registering the module */
int RedisModule_OnLoad(RedisModuleCtx *ctx) {
  if (RedisModule_Init(ctx, "hist", HISTK_MODULE_VERSION, REDISMODULE_APIVER_1)
      == REDISMODULE_ERR) {
      return REDISMODULE_ERR;
    }
    HistKType = RedisModule_CreateDataType(ctx, "aaw-histk",
                                           HISTK_ENCODING_VERSION,
                                           HistKRdbLoad, HistKRdbSave,
                                           HistKAofRewrite, HistKDigest,
                                           HistKFree);
    if (HistKType == NULL) return REDISMODULE_ERR;
    if (RedisModule_CreateCommand(ctx, "histk.add", AddCommand,
                                  "write", 1,1,1) == REDISMODULE_ERR) {
      return REDISMODULE_ERR;
    }
    if (RedisModule_CreateCommand(ctx, "histk.quantile", QuantileCommand,
                                  "readonly", 1,1,1) == REDISMODULE_ERR) {
      return REDISMODULE_ERR;
    }
    if (RedisModule_CreateCommand(ctx, "histk.count", CountCommand,
                                  "readonly", 1,1,1) == REDISMODULE_ERR) {
      return REDISMODULE_ERR;
    }
    if (RedisModule_CreateCommand(ctx, "histk.mergestore", MergeStoreCommand,
                                  "write", 1,1,1) == REDISMODULE_ERR) {
      return REDISMODULE_ERR;
    }
    if (RedisModule_CreateCommand(ctx, "histk.resize", ResizeCommand,
                                  "write", 1,1,1) == REDISMODULE_ERR) {
      return REDISMODULE_ERR;
    }
    return REDISMODULE_OK;
}
