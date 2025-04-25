#pragma once

#include <eigen3/Eigen/Dense>
#include <ceres/ceres.h>
#include "../utility/utility.h"

#include <ceres/ceres.h>

class PoseLocalParameterizationManifold : public ceres::Manifold {
public:
    int AmbientSize() const override { return 7; }
    int TangentSize() const override { return 6; }

    bool Plus(const double* x, const double* delta, double* x_plus_delta) const override;
    bool PlusJacobian(const double* x, double* jacobian) const override;

    bool Minus(const double* y, const double* x, double* y_minus_x) const override {
        Eigen::Map<Eigen::Vector3d> delta_trans(y_minus_x);
        Eigen::Map<const Eigen::Vector3d> trans_x(x);
        Eigen::Map<const Eigen::Vector3d> trans_y(y);
        delta_trans = trans_y - trans_x;

        // Quaternion-Differenz berechnen
        Eigen::Map<const Eigen::Quaterniond> rotation_x(x + 3);
        Eigen::Map<const Eigen::Quaterniond> rotation_y(y + 3);
        Eigen::Quaterniond delta_rotation = rotation_y * rotation_x.inverse();
        Eigen::Map<Eigen::Vector3d> delta_rotation_vec(y_minus_x + 3);
        delta_rotation_vec = delta_rotation.vec();  // Beispiel
        return true;
    }

    bool MinusJacobian(const double* x, double* jacobian) const override {
        Eigen::Map<Eigen::Matrix<double, 6, 6, Eigen::RowMajor>> J(jacobian);
        J.setIdentity();  // Beispiel-Jacobian, spezifisch anpassen
        return true;
    }
};
