#ifndef RAY_H
#define RAY_H

#include "vec3.h"

class ray 
{
    public :
        ray() {}; // if juste writting "ray();" get error, don't find consruction since doesn't have same signature (or thing like that).
        ray(const point3& origin, const vec3& direction, double time = 0.0)
            : orig(origin), dir(direction), tm(time)
        {}

        point3 origin() const {return orig; }
        vec3 direction() const {return dir; }
        double time() const {return tm; }

        point3 at(double t) const {
            return orig + t*dir;
        }

    public :
        point3 orig;
        vec3 dir;
        double tm;
};

#endif
