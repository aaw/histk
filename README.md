histk
=====

This is an implementation of the Ben-Haim/Tom-Tov
[streaming histogram sketch](http://www.jmlr.org/papers/volume11/ben-haim10a/ben-haim10a.pdf)
as a Redis module. It's mostly a port of [aaw/histosketch](https://github.com/aaw/histosketch)
to C with the scaffolding needed to run as a module.

The Ben-Haim/Tom-Tov sketch uses a fixed amount of space (about 1KB by default) independent
of the number of values observed and returns good estimates of quantiles and counts
of values observed below a threshold.

The sketch is a list of centroids. Each centroid is a (value, count) pair.
Whenever a new value is added to the sketch, it's added as
a (value, 1) centroid and the two centroids in the sketch with smallest difference in values
are found and merged together by summing their counts and averaging their values. Quantiles and
sums are estimated by viewing two neighboring
centroids (v1, c1) and (v2, c2) as a right trapezoid (along with the additional points
(v1, 0) and (v2, 0)) and interpolating values along the lateral side from (v1, c1) to (v2, c2)
based on the area in the trapezoid.

Commands
--------

*  `HISTK.ADD key value1 [count1] [value2 count2 ...]`:
   Adds values to the sketch. Returns the total number of values observed by the
   sketch so far. When counts are specified, can be used to add multiple observations
   of the same value in one command.

*  `HISTK.QUANTILE key q`:
   Returns an estimate of the q-quantile, the smallest value V observed by the sketch
   such that q times the total number of elements observed were less than or equal to V.
   Only q values in the range [0.0, 1.0] are valid arguments.

* `HISTK.COUNT key [value]`:
   Returns an estimate of of the number of values observed by the sketch that are at most
   the given value. If value is omitted, returns the total number of values observed by
   the sketch so far.

* `HISTK.MERGESTORE key hist1 [hist2] ... [histn]`:
   Merges hist1, hist2, ... histn, storing the results in the given key. If there's
   already a histogram sketch in the key before this command is called, the results are
   merged into that sketch.

* `HISTK.RESIZE key numcentroids`:
   Resize the sketch to numcentroids centroids. In most cases, this should be
   called once before values are added to the sketch. If called on an existing
   sketch with a smaller number of centroids than the sketch currently has, some
   centroids will be merged. The default number of centroids in each sketch is 64
   unless `HISTK.RESIZE` is called.

Trying the module
-----------------

If you have [Docker installed](https://docs.docker.com/engine/installation/), the
easiest way to try out this module is by building the Docker image defined by the
Dockerfile in this repo. Run `make image` from the top-level directory of
this repo, which will create an image named `histk`.

After you've built the image, you can either run a shell in the image:

```docker run -it histk```

and use `redis-cli` to connect inside the container. Alternatively, you can
daemonize the container and expose Redis over a port on your host machine instead:

```docker run -p $LOCAL_PORT:9999 -d histk tail -f /var/log/redis.log```

then use `redis-cli` on your host machine to connect to the Redis instance running
inside the container over `$LOCAL_PORT`.

Using the module
----------------

To use the module in your own Redis instance, you'll need to build it and load it into
a version of Redis that supports modules. Currently, only the unstable branch of Redis
supports modules, so you'll have to build Redis yourself as well.

First, make sure you have C development tools installed (`apt-get install build-essential`).
Clone this repo and run `make` from the top-level directory, which will build `histk.so` in
the src directory.

Next, clone the Redis repo `git clone https://github.com/antirez/redis.git`, make sure
you're on the unstable branch, and run `make` from the top-level directory.

Finally, start Redis and load the module. There are
[a few ways to do this](https://github.com/antirez/redis/blob/unstable/src/modules/INTRO.md#loading-modules),
the easiest is probably running `MODULE LOAD /path/to/histk.so` at a redis-cli prompt.

Testing
-------

Run `make test` from the top level directory to run the tests, which will build
a Docker image with Redis and this module and launch the tests in a container from that
image. You'll need [Docker installed](https://docs.docker.com/engine/installation/) to
run.

Testing on your own data
------------------------

The histk Docker image also comes with a script that lets you analyze how closely histk
approximates quantiles on your own data. You can run this script on any file that
contains one data point per line and it'll load all your samples in a histk and in an
exact in-memory histogram and compare the two results.

histk's tests run against a set of 50,000 pings of facebook.com, so you can use that
file if you don't have one of your own:

```
$ docker run histk ./analyze ping_data
Quantile        Actual                    HISTK.QUANTILE            Relative Error
--------------- ------------------------- ------------------------- -------------------------
0.0001          88.9                      88.85453065363578         0.0005114662133209068
0.001           89.4                      88.9724410677868          0.004782538391646556
0.01            89.7                      89.3453065363578          0.00395421921563218
0.05            90.0                      90.01934248387099         0.00021487030828348433
0.1             90.5                      90.52441067786806         0.00026965851183419844
0.25            91.3                      91.66069646736628         0.003935126845721704
0.5             95.2                      94.77773954944752         0.004435508934374849
0.75            97.6                      97.92153979296629         0.0032836472306922155
0.9             98.8                      101.40473517175185        0.02568652407937506
0.95            100.0                     103.2945520079058         0.031894731560030906
0.99            110.0                     110.23355783785534        0.0021187544195832657
0.999           456.0                     455.4970000000103         0.0011030701754160027
0.9999          772.0                     735.9928000001819         0.046641450776966464
```

You can modify the quantiles tested by passing in your own comma-separated list:

```
$ docker run histk ./analyze ping_data '0.9,0.95,0.99,0.999,0.9999,0.99999'
Quantile        Actual                    HISTK.QUANTILE            Relative Error
--------------- ------------------------- ------------------------- -------------------------
0.9             98.8                      101.40473517175185        0.02568652407937506
0.95            100.0                     103.2945520079058         0.031894731560030906
0.99            110.0                     110.23355783785534        0.0021187544195832657
0.999           456.0                     455.4970000000103         0.0011030701754160027
0.9999          772.0                     735.9928000001819         0.046641450776966464
0.99999         1111.0                    1110.9972100010855        2.5112501480851344e-06
```

If you have your own data file, analyze that by mapping it into the Docker container:

```
$ docker run histk -v /full/path/to/your/file:/tmp/datafile histk ./analyze /tmp/datafile
...
```