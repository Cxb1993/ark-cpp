//
// Created by Alex on 15.03.2015.
//

#include "../headers/Arr2d.h"

Arr2d::Arr2d(int n1, int n2) {
    _n1 = n1;
    _n2 = n2;

    _arr = new double[_n1*_n2];
}

double& Arr2d::elem(int i, int j) {
    return _arr[i*_n1 + j];
}


Arr2d::~Arr2d() {
    if (_arr != nullptr) delete _arr;
}

