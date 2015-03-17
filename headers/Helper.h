//
// Created by Alex on 17.03.2015.
//

#ifndef _ARK_CPP_CYLINDRICAL_H_
#define _ARK_CPP_CYLINDRICAL_H_


class Helper {
private:
    double	_sound,
            _ro0,
            _P0;
public:
    Helper(double sound, double ro0, double P0);
    double getXVelocity(double x, double y);
    double getYVelocity(double x, double y);
    double pressure(double x, double y);
    double density(double x, double y);
};


#endif //_ARK_CPP_CYLINDRICAL_H_
