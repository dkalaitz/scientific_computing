#include <stdio.h>
#include <math.h>

#define G 39.47841760435743
#define mA 1.1
#define mB 0.907
#define h 0.01
#define t 100

double xA = 0, yA = 11.29, zA = 0;
double vxA = 0.1, vyA = 0, vzA = 0.295;
double xB = 0, yB = -11.29, zB = 0;
double vxB = -0.061, vyB = 0, vzB = -0.362;

// Derivative functions for body A
double dxA_dt(double vxA) {
    return vxA;
}

double dyA_dt(double vyA) {
    return vyA;
}

double dzA_dt(double vzA) {
    return vzA;
}

// Derivative functions for body B
double dxB_dt(double vxB) {
    return vxB;
}

double dyB_dt(double vyB) {
    return vyB;
}

double dzB_dt(double vzB) {
    return vzB;
}

double fvxA(double xA, double xB, double yA, double yB, double zA, double zB) {
    double r = pow(xA - xB, 2) + pow(yA - yB, 2) + pow(zA - zB, 2);
    return -G * mB * (xA - xB) / pow(r, 1.5);
}

double fvyA(double xA, double xB, double yA, double yB, double zA, double zB) {
    double r = pow(xA - xB, 2) + pow(yA - yB, 2) + pow(zA - zB, 2);
    return -G * mB * (yA - yB) / pow(r, 1.5);
}

double fvzA(double xA, double xB, double yA, double yB, double zA, double zB) {
    double r = pow(xA - xB, 2) + pow(yA - yB, 2) + pow(zA - zB, 2);
    return -G * mB * (zA - zB) / pow(r, 1.5);
}

double fvxB(double xA, double xB, double yA, double yB, double zA, double zB) {
    double r = pow(xA - xB, 2) + pow(yA - yB, 2) + pow(zA - zB, 2);
    return -G * mA * (xB - xA) / pow(r, 1.5);
}

double fvyB(double xA, double xB, double yA, double yB, double zA, double zB) {
    double r = pow(xA - xB, 2) + pow(yA - yB, 2) + pow(zA - zB, 2);
    return -G * mA * (yB - yA) / pow(r, 1.5);
}

double fvzB(double xA, double xB, double yA, double yB, double zA, double zB) {
    double r = pow(xA - xB, 2) + pow(yA - yB, 2) + pow(zA - zB, 2);
    return -G * mA * (zB - zA) / pow(r, 1.5);
}

int main() {
    for (int ti = 1; ti <= t; ti++) {
        // Save current positions and velocities
        double xA_old = xA, yA_old = yA, zA_old = zA;
        double vxA_old = vxA, vyA_old = vyA, vzA_old = vzA;
        double xB_old = xB, yB_old = yB, zB_old = zB;
        double vxB_old = vxB, vyB_old = vyB, vzB_old = vzB;

        // Compute k1
        double k1_xA = h * dxA_dt(vxA_old);
        double k1_yA = h * dyA_dt(vyA_old);
        double k1_zA = h * dzA_dt(vzA_old);
        double k1_vxA = h * fvxA(xA_old, xB_old, yA_old, yB_old, zA_old, zB_old);
        double k1_vyA = h * fvyA(xA_old, xB_old, yA_old, yB_old, zA_old, zB_old);
        double k1_vzA = h * fvzA(xA_old, xB_old, yA_old, yB_old, zA_old, zB_old);

        double k1_xB = h * dxB_dt(vxB_old);
        double k1_yB = h * dyB_dt(vyB_old);
        double k1_zB = h * dzB_dt(vzB_old);
        double k1_vxB = h * fvxB(xA_old, xB_old, yA_old, yB_old, zA_old, zB_old);
        double k1_vyB = h * fvyB(xA_old, xB_old, yA_old, yB_old, zA_old, zB_old);
        double k1_vzB = h * fvzB(xA_old, xB_old, yA_old, yB_old, zA_old, zB_old);

        // Compute k2
        double k2_xA = h * dxA_dt(vxA_old + 0.5 * k1_vxA);
        double k2_yA = h * dyA_dt(vyA_old + 0.5 * k1_vyA);
        double k2_zA = h * dzA_dt(vzA_old + 0.5 * k1_vzA);
        double k2_vxA = h * fvxA(xA_old + 0.5 * k1_xA, xB_old + 0.5 * k1_xB, yA_old + 0.5 * k1_yA, yB_old + 0.5 * k1_yB, zA_old + 0.5 * k1_zA, zB_old + 0.5 * k1_zB);
        double k2_vyA = h * fvyA(xA_old + 0.5 * k1_xA, xB_old + 0.5 * k1_xB, yA_old + 0.5 * k1_yA, yB_old + 0.5 * k1_yB, zA_old + 0.5 * k1_zA, zB_old + 0.5 * k1_zB);
        double k2_vzA = h * fvzA(xA_old + 0.5 * k1_xA, xB_old + 0.5 * k1_xB, yA_old + 0.5 * k1_yA, yB_old + 0.5 * k1_yB, zA_old + 0.5 * k1_zA, zB_old + 0.5 * k1_zB);

        double k2_xB = h * dxB_dt(vxB_old + 0.5 * k1_vxB);
        double k2_yB = h * dyB_dt(vyB_old + 0.5 * k1_vyB);
        double k2_zB = h * dzB_dt(vzB_old + 0.5 * k1_vzB);
        double k2_vxB = h * fvxB(xA_old + 0.5 * k1_xA, xB_old + 0.5 * k1_xB, yA_old + 0.5 * k1_yA, yB_old + 0.5 * k1_yB, zA_old + 0.5 * k1_zA, zB_old + 0.5 * k1_zB);
        double k2_vyB = h * fvyB(xA_old + 0.5 * k1_xA, xB_old + 0.5 * k1_xB, yA_old + 0.5 * k1_yA, yB_old + 0.5 * k1_yB, zA_old + 0.5 * k1_zA, zB_old + 0.5 * k1_zB);
        double k2_vzB = h * fvzB(xA_old + 0.5 * k1_xA, xB_old + 0.5 * k1_xB, yA_old + 0.5 * k1_yA, yB_old + 0.5 * k1_yB, zA_old + 0.5 * k1_zA, zB_old + 0.5 * k1_zB);

        // Compute k3
        double k3_xA = h * dxA_dt(vxA_old + 0.5 * k2_vxA);
        double k3_yA = h * dyA_dt(vyA_old + 0.5 * k2_vyA);
        double k3_zA = h * dzA_dt(vzA_old + 0.5 * k2_vzA);
        double k3_vxA = h * fvxA(xA_old + 0.5 * k2_xA, xB_old + 0.5 * k2_xB, yA_old + 0.5 * k2_yA, yB_old + 0.5 * k2_yB, zA_old + 0.5 * k2_zA, zB_old + 0.5 * k2_zB);
        double k3_vyA = h * fvyA(xA_old + 0.5 * k2_xA, xB_old + 0.5 * k2_xB, yA_old + 0.5 * k2_yA, yB_old + 0.5 * k2_yB, zA_old + 0.5 * k2_zA, zB_old + 0.5 * k2_zB);
        double k3_vzA = h * fvzA(xA_old + 0.5 * k2_xA, xB_old + 0.5 * k2_xB, yA_old + 0.5 * k2_yA, yB_old + 0.5 * k2_yB, zA_old + 0.5 * k2_zA, zB_old + 0.5 * k2_zB);

        double k3_xB = h * dxB_dt(vxB_old + 0.5 * k2_vxB);
        double k3_yB = h * dyB_dt(vyB_old + 0.5 * k2_vyB);
        double k3_zB = h * dzB_dt(vzB_old + 0.5 * k2_vzB);
        double k3_vxB = h * fvxB(xA_old + 0.5 * k2_xA, xB_old + 0.5 * k2_xB, yA_old + 0.5 * k2_yA, yB_old + 0.5 * k2_yB, zA_old + 0.5 * k2_zA, zB_old + 0.5 * k2_zB);
        double k3_vyB = h * fvyB(xA_old + 0.5 * k2_xA, xB_old + 0.5 * k2_xB, yA_old + 0.5 * k2_yA, yB_old + 0.5 * k2_yB, zA_old + 0.5 * k2_zA, zB_old + 0.5 * k2_zB);
        double k3_vzB = h * fvzB(xA_old + 0.5 * k2_xA, xB_old + 0.5 * k2_xB, yA_old + 0.5 * k2_yA, yB_old + 0.5 * k2_yB, zA_old + 0.5 * k2_zA, zB_old + 0.5 * k2_zB);

        // Compute k4
        double k4_xA = h * dxA_dt(vxA_old + k3_vxA);
        double k4_yA = h * dyA_dt(vyA_old + k3_vyA);
        double k4_zA = h * dzA_dt(vzA_old + k3_vzA);
        double k4_vxA = h * fvxA(xA_old + k3_xA, xB_old + k3_xB, yA_old + k3_yA, yB_old + k3_yB, zA_old + k3_zA, zB_old + k3_zB);
        double k4_vyA = h * fvyA(xA_old + k3_xA, xB_old + k3_xB, yA_old + k3_yA, yB_old + k3_yB, zA_old + k3_zA, zB_old + k3_zB);
        double k4_vzA = h * fvzA(xA_old + k3_xA, xB_old + k3_xB, yA_old + k3_yA, yB_old + k3_yB, zA_old + k3_zA, zB_old + k3_zB);

        double k4_xB = h * dxB_dt(vxB_old + k3_vxB);
        double k4_yB = h * dyB_dt(vyB_old + k3_vyB);
        double k4_zB = h * dzB_dt(vzB_old + k3_vzB);
        double k4_vxB = h * fvxB(xA_old + k3_xA, xB_old + k3_xB, yA_old + k3_yA, yB_old + k3_yB, zA_old + k3_zA, zB_old + k3_zB);
        double k4_vyB = h * fvyB(xA_old + k3_xA, xB_old + k3_xB, yA_old + k3_yA, yB_old + k3_yB, zA_old + k3_zA, zB_old + k3_zB);
        double k4_vzB = h * fvzB(xA_old + k3_xA, xB_old + k3_xB, yA_old + k3_yA, yB_old + k3_yB, zA_old + k3_zA, zB_old + k3_zB);

        // Update positions and velocities for body A
        xA += (k1_xA + 2 * k2_xA + 2 * k3_xA + k4_xA) / 6.0;
        yA += (k1_yA + 2 * k2_yA + 2 * k3_yA + k4_yA) / 6.0;
        zA += (k1_zA + 2 * k2_zA + 2 * k3_zA + k4_zA) / 6.0;
        vxA += (k1_vxA + 2 * k2_vxA + 2 * k3_vxA + k4_vxA) / 6.0;
        vyA += (k1_vyA + 2 * k2_vyA + 2 * k3_vyA + k4_vyA) / 6.0;
        vzA += (k1_vzA + 2 * k2_vzA + 2 * k3_vzA + k4_vzA) / 6.0;

        // Update positions and velocities for body B
        xB += (k1_xB + 2 * k2_xB + 2 * k3_xB + k4_xB) / 6.0;
        yB += (k1_yB + 2 * k2_yB + 2 * k3_yB + k4_yB) / 6.0;
        zB += (k1_zB + 2 * k2_zB + 2 * k3_zB + k4_zB) / 6.0;
        vxB += (k1_vxB + 2 * k2_vxB + 2 * k3_vxB + k4_vxB) / 6.0;
        vyB += (k1_vyB + 2 * k2_vyB + 2 * k3_vyB + k4_vyB) / 6.0;
        vzB += (k1_vzB + 2 * k2_vzB + 2 * k3_vzB + k4_vzB) / 6.0;

        // Print the positions for body A and body B at the current time step
        printf("t = %.2f\n", ti * h);
        printf("Body A: (x: %.6f, y: %.6f, z: %.6f)\n", xA, yA, zA);
        printf("Body B: (x: %.6f, y: %.6f, z: %.6f)\n", xB, yB, zB);
    }

    return 0;
}

