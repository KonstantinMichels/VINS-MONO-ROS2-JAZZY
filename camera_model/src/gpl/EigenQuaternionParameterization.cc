#include "camodocal/gpl/EigenQuaternionParameterization.h"

#include <cmath>
#include <Eigen/Geometry>

namespace camodocal {

bool EigenQuaternionParameterization::Plus(const double* x,
                                           const double* delta,
                                           double* x_plus_delta) const {
    const double norm_delta = sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
    if (norm_delta > 0.0) {
        const double sin_delta_by_delta = (sin(norm_delta) / norm_delta);
        double q_delta[4];
        q_delta[0] = sin_delta_by_delta * delta[0];
        q_delta[1] = sin_delta_by_delta * delta[1];
        q_delta[2] = sin_delta_by_delta * delta[2];
        q_delta[3] = cos(norm_delta);
        EigenQuaternionProduct(q_delta, x, x_plus_delta);
    } else {
        for (int i = 0; i < 4; ++i) {
            x_plus_delta[i] = x[i];
        }
    }
    return true;
}

bool EigenQuaternionParameterization::Minus(const double* y, const double* x, double* y_minus_x) const {
    Eigen::Map<const Eigen::Quaterniond> q_y(y);  // Ziel-Quaternion
    Eigen::Map<const Eigen::Quaterniond> q_x(x);  // Ausgangs-Quaternion
    Eigen::Quaterniond q_diff = q_y * q_x.inverse();
    Eigen::AngleAxisd angle_axis(q_diff);
    Eigen::Map<Eigen::Vector3d> delta(y_minus_x);

    delta = angle_axis.axis() * angle_axis.angle();
    return true;
}

bool EigenQuaternionParameterization::PlusJacobian(const double* x, double* jacobian) const {
    // x[0], x[1], x[2], x[3] sind die Komponenten des Quaternions (x, y, z, w)
    const double w = x[3];
    const double x0 = x[0];
    const double x1 = x[1];
    const double x2 = x[2];

    // Füllen der Jacobian-Matrix (4x3)
    jacobian[0] =  w;  jacobian[1] =  x2; jacobian[2] = -x1;
    jacobian[3] = -x2; jacobian[4] =  w;  jacobian[5] =  x0;
    jacobian[6] =  x1; jacobian[7] = -x0; jacobian[8] =  w;
    jacobian[9] = -x0; jacobian[10] = -x1; jacobian[11] = -x2;

    return true;
}


bool EigenQuaternionParameterization::MinusJacobian(const double* x, double* jacobian) const {
    // x[0], x[1], x[2], x[3] sind die Komponenten des Quaternions (x, y, z, w)
    const double w = x[3];
    const double x0 = x[0];
    const double x1 = x[1];
    const double x2 = x[2];

    // Füllen der Jacobian-Matrix (3x4, als Transponierte von PlusJacobian)
    jacobian[0] =  w;  jacobian[1] = -x2; jacobian[2] =  x1; jacobian[3] = -x0;
    jacobian[4] =  x2; jacobian[5] =  w;  jacobian[6] = -x0; jacobian[7] = -x1;
    jacobian[8] = -x1; jacobian[9] =  x0; jacobian[10] =  w; jacobian[11] = -x2;

    return true;
}




void EigenQuaternionParameterization::EigenQuaternionProduct(const double z[4], const double w[4], double zw[4]) const {
    zw[0] = z[3] * w[0] + z[0] * w[3] + z[1] * w[2] - z[2] * w[1];
    zw[1] = z[3] * w[1] - z[0] * w[2] + z[1] * w[3] + z[2] * w[0];
    zw[2] = z[3] * w[2] + z[0] * w[1] - z[1] * w[0] + z[2] * w[3];
    zw[3] = z[3] * w[3] - z[0] * w[0] - z[1] * w[1] - z[2] * w[2];
}

}  // namespace camodocal
