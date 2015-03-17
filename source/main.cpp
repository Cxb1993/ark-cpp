#include <stdio.h>
#include <cmath>
#include "../headers/constants.h"
#include "../headers/bulk.h"
#include "../headers/forces.h"
#include "../headers/stress.h"
#include "../headers/Helper.h"

void Input();
void InitializeData();
void TimeStepSize();
void Phase1();
void Phase2();
void StressTensor();
void UseForces();
void WriteDataParaView();
void FreeMemory();

double min3d(double, double, double);
double max3d(double, double, double);

int main(int argc, char** argv) {
    Input();
    InitializeData();

    //first output
    WriteDataParaView();

    nStep = 0;
    TIME = 0.0;

    do
    {
        TimeStepSize();
        nStep++;
        TIME += 0.5*dt;
        Phase1();
        StressTensor();
        UseForces();
        Phase2();
        Phase1();
        UseForces();
        TIME += 0.5*dt;

        if (nStep % nPrint == 0)
        {
            WriteDataParaView();
            printf("step: %d dt:%E time:%8.4f\n", nStep, dt, TIME);
        }
    } while (nStep < nStop);
    FreeMemory();
    printf("Done!");
    return 0;
}

void Input() {
    // geometry index
    l = 1;

    // total number of grid nodes along the x1 axis
    n1_g = 128;
    // total number of grid nodes along the X2 axis
    n2_g = 128;
    // total number of grid nodes along the X3 axis
    n3_g = 2;

    // number of grid nodes along the x1 axis to 1 processor
    n1 = n1_g/1;
    // number of grid nodes along the X2 axis to 1 processor
    n2 = n2_g/1;
    // number of grid nodes along the X3 axis to 1 processor
    n3 = n3_g/1;

    // coordinates of west plane
    x1_w = 0.;
    // coordinates of east plane
    x1_e = 2*M_PI;

    // coordinates of south plane
    x2_s = 0.;
    // coordinates of north plane
    x2_n = 2*M_PI;

    // coordinates of bottom plane
    x3_b = 0.;
    // coordinates of top plane
    x3_t = 0.1;

    // total number of steps
    nStop = 1000;
    // print interval
    nPrint = 10;

    // Courant number
    CFL = 0.2;

    // kinematic viscosity
    VIS = 1/100000;  // 0.5
    // initial temperature
    t0 = 1.;

    // initial velocity along the x1 axis
    u10 = 0.;
    // initial velocity along the X2 axis
    u20 = 0.;
    // initial velocity along the X3 axis
    u30 = 0.;
    // sound velocity
    sound = 10.;

    // pressure on the top plane
    pOutlet = 0;
    // velocity on the bottom plane along the X3 axis
    u3Inlet = u30;
    // velocity on the bottom plane along the X2 axis
    u2Inlet = u20;
    // velocity on the bottom plane along the x1 axis
    u1Inlet = u10;
    // temperature on the bottom plane
    tInlet = t0;
    // unperturbed density of the liquid
    ro0_g = 1.;
    // unperturbed density of the borders material
    ro0_s = 100000000.;

    // #####################################################
    // 				block of arrays allocation
    // #####################################################

    // coordinates of grid nodes along the all of the axises
    x1 = new double[n1 + 2];
    x2 = new double[n2 + 2];
    x3 = new double[n3 + 2];

    // variables on the current time step
    roCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    tCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    u1Con = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    u2Con = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    u3Con = new Arr3d(n1 + 1, n2 + 1, n3 + 1);

    // variables on the next time step
    ronCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    tnCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    u1nCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    u2nCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    u3nCon = new Arr3d(n1 + 1, n2 + 1, n3 + 1);

    // forces
    f1 = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    f2 = new Arr3d(n1 + 1, n2 + 1, n3 + 1);
    f3 = new Arr3d(n1 + 1, n2 + 1, n3 + 1);

    // variables perpendicular to the axis x1
    ro1 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    t1 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u11 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u21 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u31 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    p1 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

    // variables perpendicular to the axis X2
    ro2 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    t2 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u12 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u22 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u32 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    p2 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

    // variables perpendicular to the axis X3
    ro3 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    t3 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u13 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u23 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    u33 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    p3 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

    // NMAX = MAX(n1, n2, n3)
    int nmax = (int)max3d(n1, n2, n3);

    // additional buffers for phase 2
    rBuf	= new double[nmax+2];
    qBuf	= new double[nmax+1];
    tfBuf	= new double[nmax+2];
    tbBuf	= new double[nmax+2];
    u2fBuf	= new double[nmax+2];
    u2bBuf	= new double[nmax+2];
    u3fBuf	= new double[nmax+2];
    u3bBuf	= new double[nmax+2];

    // friction stress
    sigm11 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    sigm21 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    sigm31 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

    sigm12 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    sigm22 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    sigm32 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);

    sigm13 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    sigm23 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
    sigm33 = new Arr3d(n1 + 2, n2 + 2, n3 + 2);
}

void InitializeData() {
    // grid step along the X1 axis
    dx1=(x1_e-x1_w)/(n1-1);
    // grid step along the X2 axis
    dx2=(x2_n-x2_s)/(n2-1);
    // grid step along the X3 axis
    dx3=(x3_t-x3_b)/(n3-1);

    x1[0] = x1_w - dx1;
    x2[0] = x2_s - dx2;
    x3[0] = x3_b - dx3;

    // #####################################################
    //			block of arrays initialization
    // #####################################################

    // along the X1 axis
    for (int i = 0; i < n1 + 1 ; ++i) {
        x1[i + 1] = x1[i] + dx1;
    }

    // along the X2 axis
    for (int j = 0; j < n2 + 1; ++j) {
        x2[j + 1] = x2[j] + dx2;
    }

    // along the X3 axis
    for (int k = 0; k < n3 + 1; ++k) {
        x3[k + 1] = x3[k] + dx3;
    }

    Helper helper(sound, ro0_g, 0);

    for (int i = 1; i < n1; ++i) {
        for (int j = 1; j < n2; ++j) {
            for (int k = 1; k < n3; ++k) {
                double x = 0.5*(x1[i] + x1[i + 1]) , y = 0.5*(x2[j] + x2[j+1]);
                u1nCon->elem(i, j, k) = helper.getXVelocity(x, y);
                u2nCon->elem(i, j, k) = helper.getYVelocity(x, y);
                u3nCon->elem(i, j, k) = u30;
                ronCon->elem(i, j, k) = helper.density(x, y);
                tnCon->elem(i, j, k) = t0;
            }
        }
    }

    for (int i = 1; i <= n1; ++i) {
        for (int j = 1; j <= n2; ++j) {
            for (int k = 1; k <= n3; ++k) {
                double x = x1[i], x_c = 0.5*(x1[i] + x1[i + 1]), y = x2[j], y_c = 0.5*(x2[j] + x2[j+1]);
                p1->elem(i, j, k) = helper.pressure(x, y_c);
                p2->elem(i, j, k) = helper.pressure(x_c, y);
                p3->elem(i, j, k)= helper.pressure(x_c, y_c);

                ro1->elem(i, j, k) = helper.density(x, y_c);
                ro2->elem(i, j, k) = helper.density(x_c, y);
                ro3->elem(i, j, k) = helper.density(x_c, y_c);

                u11->elem(i, j, k) = helper.getXVelocity(x, y_c);
                u12->elem(i, j, k) = helper.getXVelocity(x_c, y);
                u13->elem(i, j, k) = helper.getXVelocity(x_c, y_c);

                u21->elem(i, j, k) = helper.getYVelocity(x, y_c);
                u22->elem(i, j, k) = helper.getYVelocity(x_c, y);
                u23->elem(i, j, k) = helper.getYVelocity(x_c, y_c);

                u31->elem(i, j, k) = u32->elem(i, j, k) = u33->elem(i, j, k) = u30;
            }
        }
    }

}

void TimeStepSize() {
    double u1_c, u2_c, u3_c, dtu1, dtu2, dtu3, dtu, dtv1, dtv2, dtv3, dtv;

    dt = pow(10, 8);

    for (int i = 1; i < n1; ++i) {
        for (int j = 1; j < n2; ++j) {
            for (int k = 1; k < n3; ++k) {
                u1_c = u1Con->elem(i, j, k);
                u2_c = u2Con->elem(i, j, k);
                u3_c = u3Con->elem(i, j, k);

                dtu1 = CFL*dx1/(sound + fabs(u1_c));
                dtu2 = CFL*dx2/(sound + fabs(u2_c));
                dtu3 = CFL*dx3/(sound + fabs(u3_c));

                // DTU = MIN(DTU1, DTU2, DTU3)
                dtu = min3d(dtu1, dtu2, dtu3);

                if (VIS > pow(10, -16)) {
                    dtv1 = CFL*dx1*dx1/(2.*VIS);
                    dtv2 = CFL*dx2*dx2/(2.*VIS);
                    dtv3 = CFL*dx3*dx3/(2.*VIS);

                    // DTV = MIN (DTV1, DTV2, DTV3)
                    dtv = min3d(dtv1, dtv2, dtv3);
                } else {
                    dtv = pow(10, 16);
                }

                // DT = MIN(DT, DTU, DTV)
                dt = min3d(dt, dtu, dtv);
            }
        }
    }
}

void Phase1() {
    // initialization
    for (int i = 1; i < n1; i++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int k = 1; k < n3; k++)
            {
                u1Con->elem(i, j, k) = u1nCon->elem(i, j, k);
                u2Con->elem(i, j, k) = u2nCon->elem(i, j, k);
                u3Con->elem(i, j, k) = u3nCon->elem(i, j, k);
                roCon->elem(i, j, k) = ronCon->elem(i, j, k);
                tCon->elem(i, j, k) = tnCon->elem(i, j, k);
            }
        }
    }

    // geometric characteristics of the computational cell
    double ds1, ds2, ds3, dvc;

    // velocity, density, temperature and pressure on the eastern plane
    double	u1_e, u2_e, u3_e, ro_e, t_e, p_e,
    // velocity, density , temperature and pressure on the western plane
            u1_w, u2_w, u3_w, ro_w, t_w, p_w,
    // velocity, density , temperature and pressure on the northern plane
            u1_n, u2_n, u3_n, ro_n, t_n, p_n,
    // velocity, density , temperature and pressure on the southern plane
            u1_s, u2_s, u3_s, ro_s, t_s, p_s,
    // velocity, density , temperature and pressure on the top plane
            u1_t, u2_t, u3_t, ro_t, t_t, p_t,
    // velocity, density , temperature and pressure on the bottom plane
            u1_b, u2_b, u3_b, ro_b, t_b, p_b,
    // velocity, density and temperature in the cell center
            u1_c, u2_c, u3_c, ro_c, t_c,
    // velocity, density and temperature in the cell center on the next time step
            u1_cn, u2_cn, u3_cn, ro_cn, t_cn;

    // plane squares
    ds1 = dx2*dx3;
    ds2 = dx1*dx3;
    ds3 = dx1*dx2;

    // cell volume
    dvc = dx1*dx2*dx3;


    for (int i = 1; i < n1; i++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int k = 1; k < n3; k++)
            {
                // #########################################################
                //	get velocity, density , temperature and pressure values
                // #########################################################

                // east plane
                u1_e = u11->elem(i + 1, j, k);
                u2_e = u21->elem(i + 1, j, k);
                u3_e = u31->elem(i + 1, j, k);
                ro_e = ro1->elem(i + 1, j, k);
                t_e = t1->elem(i + 1, j, k);
                p_e = p1->elem(i + 1, j, k);
                // west plane
                u1_w = u11->elem(i, j, k);
                u2_w = u21->elem(i, j, k);
                u3_w = u31->elem(i, j, k);
                ro_w = ro1->elem(i, j, k);
                t_w = t1->elem(i, j, k);
                p_w = p1->elem(i, j, k);

                // north plane
                u1_n = u12->elem(i, j + 1, k);
                u2_n = u22->elem(i, j + 1, k);
                u3_n = u32->elem(i, j + 1, k);
                ro_n = ro2->elem(i, j + 1, k);
                t_n = t2->elem(i, j + 1, k);
                p_n = p2->elem(i, j + 1, k);
                // south plane
                u1_s = u12->elem(i, j, k);
                u2_s = u22->elem(i, j, k);
                u3_s = u32->elem(i, j, k);
                ro_s = ro2->elem(i, j, k);
                t_s = t2->elem(i, j, k);
                p_s = p2->elem(i, j, k);

                // top plane
                u1_t = u13->elem(i, j, k + 1);
                u2_t = u23->elem(i, j, k + 1);
                u3_t = u33->elem(i, j, k + 1);
                ro_t = ro3->elem(i, j, k + 1);
                t_t = t3->elem(i, j, k + 1);
                p_t = p3->elem(i, j, k + 1);
                // bottom plane
                u1_b = u13->elem(i, j, k);
                u2_b = u23->elem(i, j, k);
                u3_b = u33->elem(i, j, k);
                ro_b = ro3->elem(i, j, k);
                t_b = t3->elem(i, j, k);
                p_b = p3->elem(i, j, k);

                // cell center
                u1_c = u1Con->elem(i, j, k);
                u2_c = u2Con->elem(i, j, k);
                u3_c = u3Con->elem(i, j, k);
                ro_c = roCon->elem(i, j, k);
                t_c = tCon->elem(i, j, k);

                // #####################################################
                //				new values evaluating
                // #####################################################

                // new density
                ro_cn = (ro_c*dvc - 0.5*dt*(
                        (ro_e*u1_e - ro_w*u1_w)*ds1 +
                                (ro_n*u2_n - ro_s*u2_s)*ds2 +
                                (ro_t*u3_t - ro_b*u3_b)*ds3)) / dvc;

                // new conservative velocity along the X1 axis
                u1_cn = (ro_c*u1_c*dvc - 0.5*dt*(
                        ((ro_e*u1_e*u1_e + p_e) - (ro_w*u1_w*u1_w + p_w))*ds1 +
                                (ro_n*u1_n*u2_n - ro_s*u1_s*u2_s)*ds2 +
                                (ro_t*u1_t*u3_t - ro_b*u1_b*u3_b)*ds3)) / (ro_cn*dvc);

                // new conservative velocity along the X2 axis
                u2_cn = (ro_c*u2_c*dvc - 0.5*dt*(
                        (ro_e*u2_e*u1_e - ro_w*u2_w*u1_w)*ds1 +
                                ((ro_n*u2_n*u2_n + p_n) - (ro_s*u2_s*u2_s + p_s))*ds2 +
                                (ro_t*u2_t*u3_t - ro_b*u2_b*u3_b)*ds3)) / (ro_cn*dvc);

                // new conservative velocity along the X3 axis
                u3_cn = (ro_c*u3_c*dvc - 0.5*dt*(
                        (ro_e*u3_e*u1_e - ro_w*u3_w*u1_w)*ds1 +
                                (ro_n*u3_n*u2_n - ro_s*u3_s*u2_s)*ds2 +
                                ((ro_t*u3_t*u3_t + p_t) - (ro_b*u3_b*u3_b + p_b))*ds3)) / (ro_cn*dvc);

                // new temperature
                t_cn = (ro_c*t_c*dvc - 0.5*dt*(
                        (ro_e*t_e*u1_e - ro_w*t_w*u1_w)*ds1 +
                                (ro_n*t_n*u2_n - ro_s*t_s*u2_s)*ds2 +
                                (ro_t*t_t*u3_t - ro_b*t_b*u3_b)*ds3)) / (dvc*ro_cn);

                // finally
                u1nCon->elem(i, j, k) = u1_cn;
                u2nCon->elem(i, j, k) = u2_cn;
                u3nCon->elem(i, j, k) = u3_cn;
                ronCon->elem(i, j, k) = ro_cn;
                tnCon->elem(i, j, k) = t_cn;
            }
        }
    }

    // #####################################################
    // 					boundary conditions
    // #####################################################

    // along the X1 axis
    for (int j = 1; j < n2; ++j)
    {
        for (int k = 1; k < n3; ++k)
        {
            // on the east plane
            u1nCon->elem(n1, j, k) = /*u1nCon->elem(1, j, k);*/0;
            u2nCon->elem(n1, j, k) = /*u2nCon->elem(1, j, k);*/0;
            u3nCon->elem(n1, j, k) = /*u3nCon->elem(1, j, k);*/0;
            ronCon->elem(n1, j, k) = ronCon->elem(n1 - 1, j, k);
            tnCon->elem(n1, j, k) = /*tnCon->elem(1, j, k);*/t0;

            // on the west plane
            u1nCon->elem(0, j, k) = /*u1nCon->elem(n1 - 1, j, k);*/0;
            u2nCon->elem(0, j, k) = /*u2nCon->elem(n1 - 1, j, k);*/0;
            u3nCon->elem(0, j, k) = /*u3nCon->elem(n1 - 1, j, k);*/0;
            ronCon->elem(0, j, k) = ronCon->elem(1, j, k);
            tnCon->elem(0, j, k) = /*tnCon->elem(n1 - 1, j, k);*/t0;
        }
    }

    // along the X2 axis
    for (int i = 1; i < n1; i++)
    {
        for (int k = 1; k < n3; k++)
        {
            // on the north plane
            u1nCon->elem(i, n2, k) = /*u1nCon->elem(i, 1, k);*/0;
            u2nCon->elem(i, n2, k) = /*u2nCon->elem(i, 1, k);*/0;
            u3nCon->elem(i, n2, k) = /*u3nCon->elem(i, 1, k);*/0;
            ronCon->elem(i, n2, k) = ronCon->elem(i, n2 - 1, k);
            tnCon->elem(i, n2, k) = /*tnCon->elem(i, 1, k);*/t0;

            // on the south plane
            u1nCon->elem(i, 0, k) = /*u1nCon->elem(i, n2 - 1, k);*/0;
            u2nCon->elem(i, 0, k) = /*u2nCon->elem(i, n2 - 1, k);*/0;
            u3nCon->elem(i, 0, k) = /*u3nCon->elem(i, n2 - 1, k);*/0;
            ronCon->elem(i, 0, k) = ronCon->elem(i, 1, k);
            tnCon->elem(i, 0, k) = /*tnCon->elem(i, n2 - 1, k);*/t0;
        }
    }

    // along the X3 axis
    for (int i = 1; i < n1; ++i)
    {
        for (int j = 1; j < n2; ++j)
        {
            // on the top plane
            u1nCon->elem(i, j, n3) = u1nCon->elem(i, j, 1);//0;
            u2nCon->elem(i, j, n3) = u2nCon->elem(i, j, 1);//0;
            u3nCon->elem(i, j, n3) = u3nCon->elem(i, j, 1);//0;
            ronCon->elem(i, j, n3) = ronCon->elem(i, j, 1);
            tnCon->elem(i, j, n3) = tnCon->elem(i, j, 1);//t0;

            // on the bottom plane
            u1nCon->elem(i, j, 0) = u1nCon->elem(i, j, n3 - 1);//0;
            u2nCon->elem(i, j, 0) = u2nCon->elem(i, j, n3 - 1);//0;
            u3nCon->elem(i, j, 0) = u3nCon->elem(i, j, n3 - 1);//0;
            ronCon->elem(i, j, 0) = ronCon->elem(i, j, n3 - 1);
            tnCon->elem(i, j, 0) = tnCon->elem(i, j, n3 - 1);//t0;
        }
    }
}

void StressTensor() {
    // initialization of friction stress arrays
    for (int i = 1; i <= n1 + 1; ++i) {
        for (int j = 1; j <= n2 + 1; ++j) {
            for (int k = 1; k <= n3 + 1; ++k) {
                sigm11->elem(i, j, k) = 0.0;
                sigm21->elem(i, j, k) = 0.0;
                sigm31->elem(i, j, k) = 0.0;

                sigm12->elem(i, j, k) = 0.0;
                sigm22->elem(i, j, k) = 0.0;
                sigm32->elem(i, j, k) = 0.0;

                sigm13->elem(i, j, k) = 0.0;
                sigm23->elem(i, j, k) = 0.0;
                sigm33->elem(i, j, k) = 0.0;
            }
        }
    }


    // #####################################################
    // 					boudary conditions
    // #####################################################

    // along the X1 axis
    for (int j = 1; j < n2; ++j)
    {
        for (int k = 1; k < n3; ++k)
        {
            // on the east plane
            u1Con->elem(n1, j, k) = /*u1Con->elem(1, j, k);*/0;
            u2Con->elem(n1, j, k) = /*u2Con->elem(1, j, k);*/0;
            u3Con->elem(n1, j, k) = /*u3Con->elem(1, j, k);*/0;
            roCon->elem(n1, j, k) = roCon->elem(n1 - 1, j, k);
            tCon->elem(n1, j, k) = /*tCon->elem(1, j, k);*/t0;

            // on the west plane
            u1Con->elem(0, j, k) = /*u1Con->elem(n1 - 1, j, k);*/0;
            u2Con->elem(0, j, k) = /*u2Con->elem(n1 - 1, j, k);*/0;
            u3Con->elem(0, j, k) = /*u3Con->elem(n1 - 1, j, k);*/0;
            roCon->elem(0, j, k) = roCon->elem(1, j, k);
            tCon->elem(0, j, k) = /*tCon->elem(n1 - 1, j, k);*/t0;
        }
    }

    // along the X2 axis
    for (int i = 1; i < n1; i++)
    {
        for (int k = 1; k < n3; k++)
        {
            // on the north plane
            u1Con->elem(i, n2, k) = /*u1Con->elem(i, 1, k);*/0;
            u2Con->elem(i, n2, k) = /*u2Con->elem(i, 1, k);*/0;
            u3Con->elem(i, n2, k) = /*u3Con->elem(i, 1, k);*/0;
            roCon->elem(i, n2, k) = roCon->elem(i, n2 - 1, k);
            tCon->elem(i, n2, k) = /*tCon->elem(i, 1, k);*/t0;

            // on the south plane
            u1Con->elem(i, 0, k) = /*u1Con->elem(i, n2 - 1, k);*/0;
            u2Con->elem(i, 0, k) = /*u2Con->elem(i, n2 - 1, k);*/0;
            u3Con->elem(i, 0, k) = /*u3Con->elem(i, n2 - 1, k);*/0;
            roCon->elem(i, 0, k) = roCon->elem(i, 1, k);
            tCon->elem(i, 0, k) = /*tCon->elem(i, n2 - 1, k);*/t0;
        }
    }

    // along the X3 axis
    for (int i = 1; i < n1; ++i)
    {
        for (int j = 1; j < n2; ++j)
        {
            // on the top plane
            u1Con->elem(i, j, n3) = u1Con->elem(i, j, 1);//0;
            u2Con->elem(i, j, n3) = u2Con->elem(i, j, 1);//0;
            u3Con->elem(i, j, n3) = u3Con->elem(i, j, 1);//0;
            roCon->elem(i, j, n3) = roCon->elem(i, j, 1);
            tCon->elem(i, j, n3) = tCon->elem(i, j, 1);//t0;

            // on the bottom plane
            u1Con->elem(i, j, 0) = u1Con->elem(i, j, n3 - 1);//0;
            u2Con->elem(i, j, 0) = u2Con->elem(i, j, n3 - 1);//0;
            u3Con->elem(i, j, 0) = u3Con->elem(i, j, n3 - 1);//0;
            roCon->elem(i, j, 0) = roCon->elem(i, j, n3 - 1);
            tCon->elem(i, j, 0) = tCon->elem(i, j, n3 - 1);//t0;
        }
    }

    // #####################################################
    // 				bypassing along the faces
    // #####################################################

    double u1_c, u1_cw, u2_c, u2_cw, u3_c, u3_cw;

    // bypassing along the face perpendicular to X1
    for (int k = 1; k < n3; ++k) {
        for (int j = 1; j < n2; ++j) {
            for (int i = 1; i <= n1; ++i) {
                // velocity components in cell centers
                u1_c = u1Con->elem(i, j, k);
                u1_cw = u1Con->elem(i - 1, j, k);

                u2_c = u2Con->elem(i, j, k);
                u2_cw = u2Con->elem(i - 1, j, k);

                u3_c = u3Con->elem(i, j, k);
                u3_cw = u3Con->elem(i - 1, j, k);

                // friction stress
                sigm11->elem(i, j, k)=VIS*(u1_c - u1_cw)/dx1;
                sigm21->elem(i, j, k)=VIS*(u2_c - u2_cw)/dx1;
                sigm31->elem(i, j, k)=VIS*(u3_c - u3_cw)/dx1;
            }
        }
    }

    double u1_cs, u2_cs, u3_cs;

    // bypassing along the face perpenditcular to X2
    for (int k = 1; k < n3; ++k) {
        for (int i = 1; i < n1; ++i) {
            for (int j = 1; j <= n2; ++j) {
                // velocity components in cell centers
                u1_c = u1Con->elem(i, j, k);
                u1_cs = u1Con->elem(i, j - 1, k);

                u2_c = u2Con->elem(i, j, k);
                u2_cs = u2Con->elem(i, j - 1, k);

                u3_c = u3Con->elem(i, j, k);
                u3_cs = u3Con->elem(i, j - 1, k);

                // friction stress
                sigm12->elem(i, j, k)=VIS*(u1_c - u1_cs)/dx2;
                sigm22->elem(i, j, k)=VIS*(u2_c - u2_cs)/dx2;
                sigm32->elem(i, j, k)=VIS*(u3_c - u3_cs)/dx2;
            }
        }
    }

    double u1_cb, u2_cb, u3_cb;

    // bypassing along the face perpenditcular to X3
    for (int i = 1; i < n1; ++i) {
        for (int j = 1; j < n2; ++j) {
            for (int k = 1; k <= n3; ++k) {
                // velocity components in the cell centers
                u1_c = u1Con->elem(i, j, k);
                u1_cb = u1Con->elem(i, j, k - 1);

                u2_c = u2Con->elem(i, j, k);
                u2_cb = u2Con->elem(i, j, k - 1);

                u3_c = u3Con->elem(i, j, k);
                u3_cb = u3Con->elem(i, j, k - 1);

                // friction stress
                sigm13->elem(i, j, k)=VIS*(u1_c - u1_cb)/dx3;
                sigm23->elem(i, j, k)=VIS*(u2_c - u2_cb)/dx3;
                sigm33->elem(i, j, k)=VIS*(u3_c - u3_cb)/dx3;
            }
        }
    }

    // #####################################################
    // 				friction forces computation
    // #####################################################

    double ds1, ds2, ds3;

    // area of the face perpendicuar to x1
    ds1 = dx2*dx3;
    ds2 = dx1*dx3;
    ds3 = dx1*dx2;

    for (int i = 1; i < n1; ++i) {
        for (int j = 1; j < n2; ++j) {
            for (int k = 1; k < n3; ++k) {
                // friction forces
                f1->elem(i, j, k) =
                        (sigm11->elem(i + 1, j, k) - sigm11->elem(i, j, k)) * ds1 +
                                (sigm12->elem(i, j + 1, k) - sigm12->elem(i, j, k)) * ds2 +
                                (sigm13->elem(i, j, k + 1) - sigm13->elem(i, j, k)) * ds3;

                f2->elem(i, j, k) =
                        (sigm21->elem(i + 1, j, k) - sigm21->elem(i, j, k)) * ds1 +
                                (sigm22->elem(i, j + 1, k) - sigm22->elem(i, j, k)) * ds2 +
                                (sigm23->elem(i, j, k + 1) - sigm23->elem(i, j, k)) * ds3;

                f3->elem(i, j, k) =
                        (sigm31->elem(i + 1, j, k) - sigm31->elem(i, j, k)) * ds1 +
                                (sigm32->elem(i, j + 1, k) - sigm32->elem(i, j, k)) * ds2 +
                                (sigm33->elem(i, j, k + 1) - sigm33->elem(i, j, k)) * ds3;
            }
        }
    }
}

void UseForces() {
    double dvc, ro_c, ro_cn;

    // cell volume
    dvc = dx1*dx2*dx3;

    for (int i = 1; i < n1; ++i) {
        for (int j = 1; j < n2; ++j) {
            for (int k = 1; k < n3; ++k) {
                ro_c = roCon->elem(i, j, k);
                ro_cn = ronCon->elem(i, j, k);

                u1nCon->elem(i, j, k) = (ro_c*dvc*u1nCon->elem(i, j, k) + 0.5*dt*f1->elem(i, j, k))/(dvc*ro_cn);
                u2nCon->elem(i, j, k) = (ro_c*dvc*u2nCon->elem(i, j, k) + 0.5*dt*f2->elem(i, j, k))/(dvc*ro_cn);
                u3nCon->elem(i, j, k) = (ro_c*dvc*u3nCon->elem(i, j, k) + 0.5*dt*f3->elem(i, j, k))/(dvc*ro_cn);
            }
        }
    }
}

void Phase2() {
    double  u1_f, u1_b, u1_cn, u1_c,
            u2_f, u2_fn, u2_b, u2_bn, u2_cn, u2_c,
            u3_f, u3_fn, u3_b, u3_bn, u3_cn, u3_c,
            ro_cn, ro_c,
            t_f, t_fn, t_b, t_bn, t_cn, t_c,
            p_f, p_b, p_cn, p_c,
            r_f, r_fn, r_b, r_cn, r_c,
            q_f, q_b, q_bn, q_cn, q_c;

    double gr, gt, gu2, gu3, gq;

    double rmax, rmin, qmax, qmin, tmax, tmin, u2_max, u2_min, u3_max, u3_min;

    double qn, pn, rn, ro_n, tn, un, u2_n, u3_n, ucf, ucb;

    // first local invariants for the interior faces puts in the buffer arrays, bypass on the center of the cell
    // then by taking into account the boundary condition calculates extreme elements of the buffers
    // and only then calculates the flow variables

    // flow variables calculation on DS1 faces orthogonal X1 axis

    // bypassing along the X1 axis

    // only interior faces !
    for (int k = 1; k < n3; k++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int i = 1; i < n1; i++)
            {
                u1_f = u11->elem(i + 1, j, k);
                u1_b = u11->elem(i, j, k);
                u1_cn = u1nCon->elem(i, j, k);
                u1_c = u1Con->elem(i, j, k);

                u2_f = u21->elem(i + 1, j, k);
                u2_b = u21->elem(i, j, k);
                u2_cn = u2nCon->elem(i, j, k);
                u2_c = u2Con->elem(i, j, k);

                u3_f = u31->elem(i + 1, j, k);
                u3_b = u31->elem(i, j, k);
                u3_cn = u3nCon->elem(i, j, k);
                u3_c = u3Con->elem(i, j, k);

                ro_cn = ronCon->elem(i, j, k);
                ro_c = roCon->elem(i, j, k);

                t_f = t1->elem(i + 1, j, k);
                t_b = t1->elem(i, j, k);
                t_cn = tnCon->elem(i, j, k);
                t_c = tCon->elem(i, j, k);

                p_f = p1->elem(i + 1, j, k);
                p_b = p1->elem(i, j, k);
                p_cn = sound*sound*(ro_cn - ro0_g);
                p_c = sound*sound*(ro_c - ro0_g);

                // invariant calculation

                r_f = u1_f + p_f / (ro0_g * sound);
                r_b = u1_b + p_b / (ro0_g * sound);
                r_cn = u1_cn + p_cn / (ro0_g * sound);
                r_c = u1_c + p_c / (ro0_g * sound);

                r_fn = 2 * r_cn - r_b;

                q_f = u1_f - p_f / (ro0_g * sound);
                q_b = u1_b - p_b / (ro0_g * sound);
                q_cn = u1_cn - p_cn / (ro0_g * sound);
                q_c = u1_c - p_c / (ro0_g * sound);

                q_bn = 2 * q_cn - q_f;

                t_fn = 2 * t_cn - t_b;
                t_bn = 2 * t_cn - t_f;

                u2_fn = 2 * u2_cn - u2_b;
                u2_bn = 2 * u2_cn - u2_f;

                u3_fn = 2 * u3_cn - u3_b;
                u3_bn = 2 * u3_cn - u3_f;

                // the permissible range of changes
                gr = 2 * (r_cn - r_c) / dt + (u1_cn + sound)*(r_f - r_b) / dx1;
                gq = 2 * (r_cn - r_c) / dt + (u1_cn - sound)*(q_f - q_b) / dx1;

                gt = 2 * (t_cn - t_c) / dt + u1_cn*(t_f - t_b) / dx1;
                gu2 = 2 * (u2_cn - u2_c) / dt + u1_cn*(u2_f - u2_b) / dx1;
                gu3 = 2 * (u3_cn - u3_c) / dt + u1_cn*(u3_f - u3_b) / dx1;

                // RMAX=MAX(RF,RC,RB) +dt*GR
                rmax = max3d(r_f, r_c, r_b) + dt*gr;

                // RMIN=MIN(RF,RC,RB) +dt*GR
                rmin = min3d(r_f, r_c, r_b) + dt*gr;

                // QMAX=MAX(QF,QC,QB) +dt*GQ
                qmax = max3d(q_f, q_c, q_b) + dt*gq;

                // QMIN=MIN(QF,QC,QB) +dt*GQ
                qmin = min3d(q_f, q_c, q_b) + dt*gq;

                // TMAX=MAX(TF,TC,TB) +dt*GT
                tmax = max3d(t_f, t_c, t_b) + dt*gt;

                // TMIN=MIN(TF,TC,TB) +dt*GT
                tmin = min3d(t_f, t_c, t_b) + dt*gt;

                // U2MAX=MAX(U2F,U2C,U2B) +dt*GU2
                u2_max = max3d(u2_f, u2_c, u2_b) + dt*gu2;

                // U2MIN=MIN(U2F,U2C,U2B) +dt*GU2
                u2_min = min3d(u2_f, u2_c, u2_b) + dt*gu2;

                // U3MAX=MAX(U3F,U3C,U3B) +dt*GU3
                u3_max = max3d(u3_f, u3_c, u3_b) + dt*gu3;

                // U3MIN=MIN(U3F,U3C,U3B) +dt*GU3
                u3_min = min3d(u3_f, u3_c, u3_b) + dt*gu3;

                // invariants correction
                if (r_fn > rmax) r_fn = rmax;
                if (r_fn < rmin) r_fn = rmin;

                if (q_bn > qmax) q_bn = qmax;
                if (q_bn < qmin) q_bn = qmin;

                if (t_fn > tmax) t_fn = tmax;
                if (t_fn < tmin) t_fn = tmin;

                if (t_bn > tmax) t_bn = tmax;
                if (t_bn < tmin) t_bn = tmin;

                if (u2_fn > u2_max) u2_fn = u2_max;
                if (u2_fn < u2_min) u2_fn = u2_min;

                if (u2_bn > u2_max) u2_bn = u2_max;
                if (u2_bn < u2_min) u2_bn = u2_min;

                if (u3_fn > u3_max) u3_fn = u3_max;
                if (u3_fn < u3_min) u3_fn = u3_min;

                if (u3_bn > u3_max) u3_bn = u3_max;
                if (u3_bn < u3_min) u3_bn = u3_min;

                // put invariants to buffers
                rBuf[i + 1] = r_fn;
                qBuf[i] = q_bn;
                tfBuf[i + 1] = t_fn;
                tbBuf[i] = t_bn;
                u2fBuf[i + 1] = u2_fn;
                u2bBuf[i] = u2_bn;
                u3fBuf[i + 1] = u3_fn;
                u3bBuf[i] = u3_bn;
            }

            // boundary conditions along the X1 axis
            // assignment of boundary invatiants and add them to the buffer arrays

            // periodicity conditions
//            rBuf[1] = rBuf[n1];
//            tfBuf[1] = tfBuf[n1];
//            u2fBuf[1] = u2fBuf[n1];
//            u3fBuf[1] = u3fBuf[n1];

            // periodicity conditions
//            qBuf[n1] = qBuf[1];
//            tbBuf[n1] = tbBuf[1];
//            u2bBuf[n1] = u2bBuf[1];
//            u3bBuf[n1] = u3bBuf[1];

            // no-slip conditions
            // i == 1
            qn = qBuf[1];
            un = 0;
            pn = -qn*sound*ro0_g;
            ro_n = ro0_g + pn / (sound*sound);

            tn = t0;
            u2_n = u2bBuf[1];
            u3_n = u3bBuf[1];

            p1->elem(1, j, k) = pn;
            u11->elem(1, j, k) = un;
            ro1->elem(1, j, k) = ro_n;
            t1->elem(1, j, k) = tn;
            u21->elem(1, j, k) = u2_n;
            u31->elem(1, j, k) = u3_n;

            // i == n1
            rn = rBuf[n1];

            un = 0;
            pn = rn*sound*ro0_g;
            ro_n = ro0_g + pn / (sound*sound);

            tn = t0;
            u2_n = u2fBuf[n1];
            u3_n = u3fBuf[n1];

            p1->elem(n1, j, k) = pn;
            u11->elem(n1, j, k) = un;
            ro1->elem(n1, j, k) = ro_n;
            t1->elem(n1, j, k) = tn;
            u21->elem(n1, j, k) = u2_n;
            u31->elem(n1, j, k) = u3_n;

            // the flow variables calculations
            for (int i = 2; i < n1; i++)
            {
                rn = rBuf[i];
                qn = qBuf[i];

                pn = (rn - qn)*sound*ro0_g / 2;
                un = (rn + qn) / 2;

                ro_n = ro0_g + pn / (sound*sound);

                ucf = u1nCon->elem(i, j, k);
                ucb = u1nCon->elem(i - 1, j, k);

                if (ucf >= 0 && ucb >= 0)
                {
                    tn = tfBuf[i];
                    u2_n = u2fBuf[i];
                    u3_n = u3fBuf[i];
                }
                else if (ucf <= 0 && ucb <= 0)
                {
                    tn = tbBuf[i];
                    u2_n = u2bBuf[i];
                    u3_n = u3bBuf[i];
                }
                else if (ucb >= 0 && ucf <= 0)
                {
                    if (ucb > -ucf)
                    {
                        tn = tfBuf[i];
                        u2_n = u2fBuf[i];
                        u3_n = u3fBuf[i];
                    }
                    else
                    {
                        tn = tbBuf[i];
                        u2_n = u2bBuf[i];
                        u3_n = u3bBuf[i];
                    }
                }
                else if (ucb <= 0 && ucf >= 0)
                {
                    tn = tnCon->elem(i, j, k) + tnCon->elem(i - 1, j, k) - t1->elem(i, j, k);
                    u2_n = u2nCon->elem(i, j, k) + u2nCon->elem(i - 1, j, k) - u21->elem(i, j, k);
                    u3_n = u3nCon->elem(i, j, k) + u3nCon->elem(i - 1, j, k) - u31->elem(i, j, k);
                }

                p1->elem(i, j, k) = pn;
                u11->elem(i, j, k) = un;
                ro1->elem(i, j, k) = ro_n;
                t1->elem(i, j, k) = tn;
                u21->elem(i, j, k) = u2_n;
                u31->elem(i, j, k) = u3_n;
            }
        }
    }

    // flow variables calculation on DS2 faces orthogonal X2 axis

    // bypassing along the X2 axis

    double u1_fn, u1_bn;

    double gu1;

    double u1_max, u1_min, u1_n;

    for (int k = 1; k < n3; k++)
    {
        for (int i = 1; i < n1; i++)
        {
            for (int j = 1; j < n2; j++)
            {
                u2_f = u22->elem(i, j + 1, k);
                u2_b = u22->elem(i, j, k);
                u2_cn = u2nCon->elem(i, j, k);
                u2_c = u2Con->elem(i, j, k);

                u1_f = u12->elem(i, j + 1, k);
                u1_b = u12->elem(i, j, k);
                u1_cn = u1nCon->elem(i, j, k);
                u1_c = u1Con->elem(i, j, k);

                u3_f = u32->elem(i, j + 1, k);
                u3_b = u32->elem(i, j, k);
                u3_cn = u3nCon->elem(i, j, k);
                u3_c = u3Con->elem(i, j, k);

                ro_cn = ronCon->elem(i, j, k);
                ro_c = roCon->elem(i, j, k);

                t_f = t2->elem(i, j + 1, k);
                t_b = t2->elem(i, j, k);
                t_cn = tnCon->elem(i, j, k);
                t_c = tCon->elem(i, j, k);

                p_f = p2->elem(i, j + 1, k);
                p_b = p2->elem(i, j, k);
                p_cn = sound*sound*(ro_cn - ro0_g);
                p_c = sound*sound*(ro_c - ro0_g);

                // invariant calculation
                r_f = u2_f + p_f / (ro0_g*sound);
                r_b = u2_b + p_b / (ro0_g*sound);
                r_cn = u2_cn + p_cn / (ro0_g*sound);
                r_c = u2_c + p_c / (ro0_g*sound);

                r_fn = 2 * r_cn - r_b;

                q_f = u2_f - p_f / (ro0_g*sound);
                q_b = u2_b - p_b / (ro0_g*sound);
                q_cn = u2_cn - p_cn / (ro0_g*sound);
                q_c = u2_c - p_c / (ro0_g*sound);

                q_bn = 2 * q_cn - q_f;

                t_fn = 2 * t_cn - t_b;
                t_bn = 2 * t_cn - t_f;

                u1_fn = 2 * u1_cn - u1_b;
                u1_bn = 2 * u1_cn - u1_f;

                u3_fn = 2 * u3_cn - u3_b;
                u3_bn = 2 * u3_cn - u3_f;

                // the permissible range of changes
                gr = 2 * (r_cn - r_c) / dt + (u2_cn + sound)*(r_f - r_b) / dx2;
                gq = 2 * (q_cn - q_c) / dt + (u2_cn - sound)*(q_f - q_b) / dx2;
                gt = 2 * (t_cn - t_c) / dt + u2_cn*(t_f - t_b) / dx2;
                gu1 = 2 * (u1_cn - u1_c) / dt + u2_cn*(u1_f - u1_b) / dx2;
                gu3 = 2 * (u3_cn - u3_c) / dt + u2_cn*(u3_f - u3_b) / dx2;


                // RMAX=MAX(RF,RC,RB) +dt*GR
                rmax = max3d(r_f, r_c, r_b) + dt*gr;

                // RMIN=MIN(RF,RC,RB) +dt*GR
                rmin = min3d(r_f, r_c, r_b) + dt*gr;

                // QMAX=MAX(QF,QC,QB) +dt*GQ
                qmax = max3d(q_f, q_c, q_b) + dt*gq;

                // QMIN=MIN(QF,QC,QB) +dt*GQ
                qmin = min3d(q_f, q_c, q_b) + dt*gq;

                // TMAX=MAX(TF,TC,TB) +dt*GT
                tmax = max3d(t_f, t_c, t_b) + dt*gt;

                // TMIN=MIN(TF,TC,TB) +dt*GT
                tmin = min3d(t_f, t_c, t_b) + dt*gt;

                // U1MAX=MAX(U1F,U1C,U1B) +dt*GU1
                u1_max = max3d(u1_f, u1_c, u1_b) + dt*gu1;

                // U1MIN=MIN(U1F,U1C,U1B) +dt*GU1
                u1_min = min3d(u1_f, u1_c, u1_b) + dt*gu1;

                // U3MAX=MAX(U3F,U3C,U3B) +dt*GU3
                u3_max = max3d(u3_f, u3_c, u3_b) + dt*gu3;

                // U3MIN=MIN(U3F,U3C,U3B) +dt*GU3
                u3_min = min3d(u3_f, u3_c, u3_b) + dt*gu3;

                // invariants correction
                if (r_fn > rmax) r_fn = rmax;
                if (r_fn < rmin) r_fn = rmin;

                if (q_bn > qmax) q_bn = qmax;
                if (q_bn < qmin) q_bn = qmin;

                if (t_fn > tmax) t_fn = tmax;
                if (t_fn < tmin) t_fn = tmin;

                if (t_bn > tmax) t_bn = tmax;
                if (t_bn < tmin) t_bn = tmin;

                if (u1_fn > u1_max) u1_fn = u1_max;
                if (u1_fn < u1_min) u1_fn = u1_min;

                if (u1_bn > u1_max) u1_bn = u1_max;
                if (u1_bn < u1_min) u1_bn = u1_min;

                if (u3_fn > u3_max) u3_fn = u3_max;
                if (u3_fn < u3_min) u3_fn = u3_min;

                if (u3_bn > u3_max) u3_bn = u3_max;
                if (u3_bn < u3_min) u3_bn = u3_min;

                // put invariants to buffers
                // ==================================================
                // !!! IMPORTANT !!!
                // ==================================================
                // u2fBuf and u2bBuf are actially the u1fBuf and u1bBuf
                // It's not an error. We do it to save dynamic memory
                // ==================================================
                rBuf[j + 1] = r_fn;
                qBuf[j] = q_bn;
                tfBuf[j + 1] = t_fn;
                tbBuf[j] = t_bn;
                u2fBuf[j + 1] = u1_fn;
                u2bBuf[j] = u1_bn;
                u3fBuf[j + 1] = u3_fn;
                u3bBuf[j] = u3_bn;
            }

            // boundary conditions along the X2 axis
            // assignment of boundary invatiants and add them to the buffer arrays

            // periodicity conditions
//            rBuf[1] = rBuf[n2];
//            tfBuf[1] = tfBuf[n2];
//            u2fBuf[1] = u2fBuf[n2];
//            u3fBuf[1] = u3fBuf[n2];

            // periodicity conditions
//            qBuf[n2] = qBuf[1];
//            tbBuf[n2] = tbBuf[1];
//            u2bBuf[n2] = u2bBuf[1];
//            u3bBuf[n2] = u3bBuf[1];

            // no-slip
            // j == 1
            qn = qBuf[1];

            pn = -qn*sound*ro0_g;
            un = 0.;

            ro_n = ro0_g + pn / (sound*sound);

            tn = t0;
            u1_n = u2bBuf[1];
            u3_n = u3bBuf[1];

            p2->elem(i, 1, k) = pn;
            u22->elem(i, 1, k) = un;
            ro2->elem(i, 1, k) = ro_n;
            t2->elem(i, 1, k) = tn;
            u12->elem(i, 1, k) = u1_n;
            u32->elem(i, 1, k) = u3_n;

            // j == n2
            rn = rBuf[n2];

            pn = rn*sound*ro0_g;
            un = 0.;

            ro_n = ro0_g + pn / (sound*sound);

            tn = t0;
            u1_n = u2fBuf[n2];
            u3_n = u3fBuf[n2];

            p2->elem(i, n2, k) = pn;
            u22->elem(i, n2, k) = un;
            ro2->elem(i, n2, k) = ro_n;
            t2->elem(i, n2, k) = tn;
            u12->elem(i, n2, k) = u1_n;
            u32->elem(i, n2, k) = u3_n;

            // the flow variables calculations
            for (int j = 2; j < n2; j++)
            {
                rn = rBuf[j];
                qn = qBuf[j];

                pn = (rn - qn)*sound*ro0_g / 2;
                un = (rn + qn) / 2;

                ro_n = ro0_g + pn / (sound*sound);

                ucf = u2nCon->elem(i, j, k);
                ucb = u2nCon->elem(i, j - 1, k);

                if (ucf >= 0 && ucb >= 0)
                {
                    tn = tfBuf[j];
                    u1_n = u2fBuf[j];
                    u3_n = u3fBuf[j];
                }
                else if (ucf <= 0 && ucb <= 0)
                {
                    tn = tbBuf[j];
                    u1_n = u2bBuf[j];
                    u3_n = u3bBuf[j];
                }
                else if (ucb >= 0 && ucf <= 0)
                {
                    if (ucb > -ucf)
                    {
                        tn = tfBuf[j];
                        u1_n = u2fBuf[j];
                        u3_n = u3fBuf[j];
                    }
                    else
                    {
                        tn = tbBuf[j];
                        u1_n = u2bBuf[j];
                        u3_n = u3bBuf[j];
                    }
                }
                else
                if (ucb <= 0 && ucf >= 0)
                {
                    tn = tnCon->elem(i, j, k) + tnCon->elem(i, j - 1, k) - t2->elem(i, j, k);
                    u1_n = u1nCon->elem(i, j, k) + u1nCon->elem(i, j - 1, k) - u12->elem(i, j, k);
                    u3_n = u3nCon->elem(i, j, k) + u3nCon->elem(i, j - 1, k) - u32->elem(i, j, k);
                }

                p2->elem(i, j, k) = pn;
                u22->elem(i, j, k) = un;
                ro2->elem(i, j, k) = ro_n;
                t2->elem(i, j, k) = tn;
                u12->elem(i, j, k) = u1_n;
                u32->elem(i, j, k) = u3_n;
            }
        }
    }

    // flow variables calculation on DS3 faces orthogonal X3 axis

    // bypassing along the X3 axis

    for (int i = 1; i < n1; i++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int k = 1; k < n3; k++)
            {
                u3_f = u33->elem(i, j, k + 1);
                u3_b = u33->elem(i, j, k);
                u3_cn = u3nCon->elem(i, j, k);
                u3_c = u3Con->elem(i, j, k);

                u1_f = u13->elem(i, j, k + 1);
                u1_b = u13->elem(i, j, k);
                u1_cn = u1nCon->elem(i, j, k);
                u1_c = u1Con->elem(i, j, k);

                u2_f = u23->elem(i, j, k + 1);
                u2_b = u23->elem(i, j, k);
                u2_cn = u2nCon->elem(i, j, k);
                u2_c = u2Con->elem(i, j, k);

                ro_cn = ronCon->elem(i, j, k);
                ro_c = roCon->elem(i, j, k);

                t_f = t3->elem(i, j, k + 1);
                t_b = t3->elem(i, j, k);
                t_cn = tnCon->elem(i, j, k);
                t_c = tCon->elem(i, j, k);

                p_f = p3->elem(i, j, k + 1);
                p_b = p3->elem(i, j, k);
                p_cn = sound*sound*(ro_cn - ro0_g);
                p_c = sound*sound*(ro_c - ro0_g);

                // invariant calculation
                r_f = u3_f + p_f / (ro0_g*sound);
                r_b = u3_b + p_b / (ro0_g*sound);
                r_cn = u3_cn + p_cn / (ro0_g*sound);
                r_c = u3_c + p_c / (ro0_g*sound);

                r_fn = 2 * r_cn - r_b;

                q_f = u3_f - p_f / (ro0_g*sound);
                q_b = u3_b - p_b / (ro0_g*sound);
                q_cn = u3_cn - p_cn / (ro0_g*sound);
                q_c = u3_c - p_c / (ro0_g*sound);

                q_bn = 2 * q_cn - q_f;

                t_fn = 2 * t_cn - t_b;
                t_bn = 2 * t_cn - t_f;

                u2_fn = 2 * u2_cn - u2_b;
                u2_bn = 2 * u2_cn - u2_f;

                u1_fn = 2 * u1_cn - u1_b;
                u1_bn = 2 * u1_cn - u1_f;

                // the permissible range of changes
                gr = 2 * (r_cn - r_c) / dt + (u3_cn + sound)*(r_f - r_b) / dx3;
                gq = 2 * (r_cn - r_c) / dt + (u3_cn - sound)*(q_f - q_b) / dx3;

                gt = 2 * (t_cn - t_c) / dt + u3_cn*(t_f - t_b) / dx3;
                gu1 = 2 * (u1_cn - u1_c) / dt + u3_cn*(u1_f - u1_b) / dx3;
                gu2 = 2 * (u2_cn - u2_c) / dt + u3_cn*(u2_f - u2_b) / dx3;

                // RMAX=MAX(RF,RC,RB) +dt*GR
                rmax = max3d(r_f, r_c, r_b) + dt*gr;

                // RMIN=MIN(RF,RC,RB) +dt*GR
                rmin = min3d(r_f, r_c, r_b) + dt*gr;

                // QMAX=MAX(QF,QC,QB) +dt*GQ
                qmax = max3d(q_f, q_c, q_b) + dt*gq;

                // QMIN=MIN(QF,QC,QB) +dt*GQ
                qmin = min3d(q_f, q_c, q_b) + dt*gq;

                // TMAX=MAX(TF,TC,TB) +dt*GT
                tmax = max3d(t_f, t_c, t_b) + dt*gt;

                // TMIN=MIN(TF,TC,TB) +dt*GT
                tmin = min3d(t_f, t_c, t_b) + dt*gt;

                // U1MAX=MAX(U1F,U1C,U1B) +dt*GU1
                u1_max = max3d(u1_f, u1_c, u1_b) + dt*gu1;

                // U1MIN=MIN(U1F,U1C,U1B) +dt*GU1
                u1_min = min3d(u1_f, u1_c, u1_b) + dt*gu1;

                // U2MAX=MAX(U2F,U2C,U2B) +dt*GU2
                u2_max = max3d(u2_f, u2_c, u2_b) + dt*gu2;

                // U2MIN=MIN(U2F,U2C,U2B) +dt*GU2
                u2_min = min3d(u2_f, u2_c, u2_b) + dt*gu2;

                // invariants correction
                if (r_fn > rmax) r_fn = rmax;
                if (r_fn < rmin) r_fn = rmin;

                if (q_bn > qmax) q_bn = qmax;
                if (q_bn < qmin) q_bn = qmin;

                if (t_fn > tmax) t_fn = tmax;
                if (t_fn < tmin) t_fn = tmin;

                if (t_bn > tmax) t_bn = tmax;
                if (t_bn < tmin) t_bn = tmin;

                if (u1_fn > u1_max) u1_fn = u1_max;
                if (u1_fn < u1_min) u1_fn = u1_min;

                if (u1_bn > u1_max) u1_bn = u1_max;
                if (u1_bn < u1_min) u1_bn = u1_min;

                if (u2_fn > u2_max) u2_fn = u2_max;
                if (u2_fn < u2_min) u2_fn = u2_min;

                if (u2_bn > u2_max) u2_bn = u2_max;
                if (u2_bn < u2_min) u2_bn = u2_min;

                // put invariants to buffers
                // ====================================================
                // !!! IMPORTANT !!!
                // ====================================================
                // u2fBuf and u2bBuf are actially the u1fBuf and u1bBuf
                // u3fBuf and u3bBuf are actially the u2fBuf and u2bBuf
                // It's not an error. We do it to save dynamic memory
                // ====================================================
                rBuf[k + 1] = r_fn;
                qBuf[k] = q_bn;
                tfBuf[k + 1] = t_fn;
                tbBuf[k] = t_bn;
                u2fBuf[k + 1] = u1_fn;
                u2bBuf[k] = u1_bn;
                u3fBuf[k + 1] = u2_fn;
                u3bBuf[k] = u2_bn;
            }

            // boundary conditions along the X3 axis
            // assignment of boundary invatiants and add them to the buffer arrays

            // periodicity conditions
            rBuf[1] = rBuf[n3];
            tfBuf[1] = tfBuf[n3];
            u2fBuf[1] = u2fBuf[n3];
            u3fBuf[1] = u3fBuf[n3];

            // periodicity conditions
            qBuf[n3] = qBuf[1];
            tbBuf[n3] = tbBuf[1];
            u2bBuf[n3] = u2bBuf[1];
            u3bBuf[n3] = u3bBuf[1];

            // no-slip conditions
            // k == 1;
//			qn = qBuf[1];
//			un = 0.;
//			pn = -qn*sound*ro0_g;
//			ro_n = ro0_g + pn / (sound*sound);
//
//			tn = t0;
//			u1_n = 0.;
//			u2_n = 0.;
//
//			p3->elem(i, j, 1) = pn;
//			u33->elem(i, j, 1) = un;
//			ro3->elem(i, j, 1) = ro_n;
//			t3->elem(i, j, 1) = tn;
//			u13->elem(i, j, 1) = u1_n;
//			u23->elem(i, j, 1) = u2_n;
//
//
//			// k == n3
//			rn = rBuf[n3];
//
//			un = 0.;
//			pn = rn*sound*ro0_g;
//			ro_n = ro0_g + pn / (sound*sound);
//
//			tn = t0;
//			u1_n = 0.;
//			u2_n = 0.;
//
//			p3->elem(i, j, n3) = pn;
//			u33->elem(i, j, n3) = un;
//			ro3->elem(i, j, n3) = ro_n;
//			t1->elem(i, j, n3) = tn;
//			u13->elem(i, j, n3) = u1_n;
//			u23->elem(i, j, n3) = u2_n;

            // the flow variables calculations
            for (int k = 1; k <= n3; k++)
            {
                rn = rBuf[k];
                qn = qBuf[k];

                pn = (rn - qn)*sound*ro0_g / 2;
                un = (rn + qn) / 2;

                ro_n = ro0_g + pn / (sound*sound);

                ucf = u3nCon->elem(i, j, k);
                ucb = u3nCon->elem(i, j, k - 1);

                if (ucf >= 0 && ucb >= 0)
                {
                    tn = tfBuf[k];
                    u1_n = u2fBuf[k];
                    u2_n = u3fBuf[k];
                }
                else if (ucf <= 0 && ucb <= 0)
                {
                    tn = tbBuf[k];
                    u1_n = u2bBuf[k];
                    u2_n = u3bBuf[k];
                }
                else if (ucb >= 0 && ucf <= 0)
                {
                    if (ucb > -ucf)
                    {
                        tn = tfBuf[k];
                        u1_n = u2fBuf[k];
                        u2_n = u3fBuf[k];
                    }
                    else
                    {
                        tn = tbBuf[k];
                        u1_n = u2bBuf[k];
                        u2_n = u3bBuf[k];
                    }
                }
                else
                if (ucb <= 0 && ucf >= 0)
                {
                    tn = tnCon->elem(i, j, k) + tnCon->elem(i, j, k - 1) - t3->elem(i, j, k);
                    u1_n = u1nCon->elem(i, j, k) + u1nCon->elem(i, j, k - 1) - u13->elem(i, j, k);
                    u2_n = u2nCon->elem(i, j, k) + u2nCon->elem(i, j, k - 1) - u23->elem(i, j, k);
                }

                p3->elem(i, j, k) = pn;
                u33->elem(i, j, k) = un;
                ro3->elem(i, j, k) = ro_n;
                t3->elem(i, j, k) = tn;
                u13->elem(i, j, k) = u1_n;
                u23->elem(i, j, k) = u2_n;
            }

            // inlet conditions
            /*qn = qBuf[1];
            rn = u3Inlet + (ronCon->elem(i, j, 1) - ro0_g)*sound / ro0_g;
            un = (rn + qn) / 2;
            pn = (rn - qn)*sound*ro0_g / 2;
            ro_n = ro0_g + pn / sound / sound;
            u2_n = u2Inlet;
            u1_n = u1Inlet;
            tn = tInlet;
            p3->elem(i, j, 1) = pn;
            u33->elem(i, j, 1) = un;
            ro3->elem(i, j, 1) = ro_n;
            t3->elem(i, j, 1) = tn;
            u23->elem(i, j, 1) = u2_n;
            u13->elem(i, j, 1) = u1_n;*/

            // outlet conditions
            /*rn = rBuf[n3];
            pn = pOutlet;
            un = rn - pn / ro0_g / sound;
            tn = tfBuf[n3];
            u2_n = u2fBuf[n3];
            u1_n = u3fBuf[n3];
            p3->elem(i, j, n3) = pn;
            u33->elem(i, j, n3) = un;
            ro3->elem(i, j, n3) = ro_n;
            t3->elem(i, j, n3) = tn;
            u23->elem(i, j, n3) = u2_n;
            u13->elem(i, j, n3) = u1_n;*/
        }
    }
}

void FreeMemory() {
    delete[] x1;
    delete[] x2;
    delete[] x3;

    delete roCon;
    delete u1Con;
    delete u2Con;
    delete u3Con;
    delete tCon;

    delete ronCon;
    delete u1nCon;
    delete u2nCon;
    delete u3nCon;
    delete tnCon;

    delete ro1;
    delete t1;
    delete u11;
    delete u21;
    delete u31;
    delete p1;

    delete ro2;
    delete t2;
    delete u12;
    delete u22;
    delete u32;
    delete p2;

    delete ro3;
    delete t3;
    delete u13;
    delete u23;
    delete u33;
    delete p3;

    delete f1;
    delete f2;
    delete f3;

    delete sigm11;
    delete sigm21;
    delete sigm31;

    delete sigm12;
    delete sigm22;
    delete sigm32;

    delete sigm13;
    delete sigm23;
    delete sigm33;

    delete[] rBuf;
    delete[] qBuf;
    delete[] tfBuf;
    delete[] tbBuf;
    delete[] u2fBuf;
    delete[] u2bBuf;
    delete[] u3fBuf;
    delete[] u3bBuf;
}

void WriteDataParaView() {
    char filename[100];

    sprintf(filename, "C:\\Users\\Alex\\Desktop\\IM\\ark-cpp\\out\\out_%d.vtk", nStep);
    //sprintf_s(filename, 50, "out_%d.vtk", nStep);

    FILE *fd = fopen(filename, "w");
//	FILE *fd;
//	fopen_s(&fd, filename, "w");

    fprintf(fd, "# vtk DataFile Version 3.0\nvtk output\nASCII\n");
    fprintf(fd, "DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d", n1, n2, n3);

    fprintf(fd, "\nX_COORDINATES %d float\n", n1);
    for (int i = 1; i <= n1; i++)
    {
        fprintf(fd, "%f ", x1[i]);
    }

    fprintf(fd, "\nY_COORDINATES %d float\n", n2);
    for (int j = 1; j <= n2; j++)
    {
        fprintf(fd, "%f ", x2[j]);
    }

    fprintf(fd, "\nZ_COORDINATES %d float\n", n3);
    for (int k = 1; k <= n3; k++)
    {
        fprintf(fd, "%f ", x3[k]);
    }

    fprintf(fd, "\nCELL_DATA %d\nscalars U1 float\nLOOKUP_TABLE default\n", (n1-1)*(n2-1)*(n3-1));
    for (int k = 1; k < n3; k++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int i = 1; i < n1; i++)
            {
                fprintf(fd, "%f ", u1nCon->elem(i, j, k));
            }
        }
    }

    fprintf(fd, "\nscalars U2 float\nLOOKUP_TABLE default\n");
    for (int k = 1; k < n3; k++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int i = 1; i < n1; i++)
            {
                fprintf(fd, "%f ", u2nCon->elem(i, j, k));
            }
        }
    }

    fprintf(fd, "\nscalars U3 float\nLOOKUP_TABLE default\n");
    for (int k = 1; k < n3; k++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int i = 1; i < n1; i++)
            {
                fprintf(fd, "%f ", u3nCon->elem(i, j, k));
            }
        }
    }

    fprintf(fd, "\nscalars PC float\nLOOKUP_TABLE default\n");
    for (int k = 1; k < n3; k++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int i = 1; i < n1; i++)
            {
                fprintf(fd, "%f ", ronCon->elem(i, j, k));
            }
        }
    }

    fprintf(fd, "\nscalars TC float\nLOOKUP_TABLE default\n");
    for (int k = 1; k < n3; k++)
    {
        for (int j = 1; j < n2; j++)
        {
            for (int i = 1; i < n1; i++)
            {
                fprintf(fd, "%f ", tnCon->elem(i, j, k));
            }
        }
    }

//	fprintf(fd, "\nscalars R1 float\nLOOKUP_TABLE default\n");
//	for (int k = 1; k < n3; k++)
//	{
//		for (int i = 1; i < n2; i++)
//		{
//			for (int j = 1; j < n1; j++)
//			{
//				double d2u3 = (u3nCon->elem(i, j + 1, k) - u3nCon->elem(i, j, k))/dx2,
//				d3u2 = (u2nCon->elem(i, j, k + 1) - u2nCon->elem(i, j, k))/dx3,
//				r1 = d2u3 - d3u2;
//
//				fprintf(fd, "%f ", r1);
//			}
//		}
//	}
//
//	fprintf(fd, "\nscalars R2 float\nLOOKUP_TABLE default\n");
//	for (int k = 1; k < n3; k++)
//	{
//		for (int i = 1; i < n2; i++)
//		{
//			for (int j = 1; j < n1; j++)
//			{
//				double d3u1 = (u1nCon->elem(i, j, k + 1) - u1nCon->elem(i, j, k))/dx3,
//						d1u3 = (u3nCon->elem(i + 1, j, k) - u3nCon->elem(i, j, k))/dx1,
//						r2 = d3u1 - d1u3;
//
//				fprintf(fd, "%f ", r2);
//			}
//		}
//	}
//
//	fprintf(fd, "\nscalars R3 float\nLOOKUP_TABLE default\n");
//	for (int k = 1; k < n3; k++)
//	{
//		for (int i = 1; i < n2; i++)
//		{
//			for (int j = 1; j < n1; j++)
//			{
//				double d1u2 = (u2nCon->elem(i + 1, j, k) - u2nCon->elem(i, j, k))/dx1,
//						d2u1 = (u1nCon->elem(i, j + 1, k) - u1nCon->elem(i, j, k))/dx2,
//						r3 = d1u2 - d2u1;
//
//				fprintf(fd, "%f ", r3);
//			}
//		}
//	}

    fclose(fd);
}

double min3d(double x1, double x2, double x3) {
    x1 = x1 > x2 ? x2 : x1;
    x1 = x1 > x3 ? x3 : x1;
    return x1;
}

double max3d(double x1, double x2, double x3) {
    x1 = x1 < x2 ? x2 : x1;
    x1 = x1 < x3 ? x3 : x1;
    return x1;
}
