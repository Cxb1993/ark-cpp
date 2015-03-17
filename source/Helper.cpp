//
// Created by Alex on 17.03.2015.
//

#include <math.h>
#include "../headers/Helper.h"

Helper::Helper(double x0, double y0, double alpha, double betta, double sound, double r0, double ro0, double P0) {
    _x0 = x0;
    _y0 = y0;
    _alpha = alpha;
    _betta = betta;
    _sound = sound;
    _r0 = r0;
    _ro0 = ro0;
    _P0 = P0;
}

double Helper::radius(double x, double y) {
    return sqrt((x - _x0)*(x - _x0) + (y - _y0)*(y - _y0));
}

double Helper::eta(double x, double y) {
    return radius(x, y) / _r0;
}

double Helper::azimuthalVelocity(double x, double y) {
    return _alpha*eta(x, y)*exp(_betta*(1 - eta(x, y)*eta(x, y)));
}

double Helper::pressure(double x, double y) {
    return -_ro0*_alpha*_alpha / (4 * _betta)*exp(2 * _betta*(1 - eta(x, y)*eta(x, y))) + _P0;
}

double Helper::density(double x, double y) {
    return pressure(x, y) / (_sound*_sound) + _ro0;
}

double Helper::getXVelocity(double x, double y) {
    return azimuthalVelocity(x, y) * y/radius(x, y);
}

double Helper::getYVelocity(double x, double y) {
    return azimuthalVelocity(x, y) * x/radius(x, y);
}
