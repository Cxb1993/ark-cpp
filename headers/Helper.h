//
// Created by Alex on 17.03.2015.
//

#ifndef _ARK_CPP_CYLINDRICAL_H_
#define _ARK_CPP_CYLINDRICAL_H_


class Helper {
private:
    double	_x0,
            _y0,
            _alpha,
            _betta,
            _sound,
            _r0,
            _ro0,
            _P0;

    double eta(double x, double y);
    double azimuthalVelocity(double x, double y);
public:
    Helper(double x0, double y0, double alpha, double betta, double sound, double r0, double ro0, double P0);
    double radius(double x, double y);
    double getXVelocity(double x, double y);
    double getYVelocity(double x, double y);
    double pressure(double x, double y);
    double density(double x, double y);
};


#endif //_ARK_CPP_CYLINDRICAL_H_
