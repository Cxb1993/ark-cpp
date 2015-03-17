//
// Created by Alex on 17.03.2015.
//

#include <cmath>
#include "../headers/Helper.h"

Helper::Helper(double sound, double ro0, double P0) {
    _sound = sound;
    _ro0 = ro0;
    _P0 = P0;
}

double Helper::pressure(double x, double y) {
    return _P0 + 1/16*(cos(2*x)+cos(2*y));
}

double Helper::density(double x, double y) {
    return pressure(x, y) / (_sound*_sound) + _ro0;
}

double Helper::getXVelocity(double x, double y) {
    return sin(x)*cos(y);
}

double Helper::getYVelocity(double x, double y) {
    return -cos(x)*sin(y);
}
