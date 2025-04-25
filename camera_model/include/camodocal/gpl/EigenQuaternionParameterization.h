#ifndef EIGENQUATERNIONPARAMETERIZATION_H
#define EIGENQUATERNIONPARAMETERIZATION_H

#include "ceres/manifold.h"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace camodocal {

class EigenQuaternionParameterization : public ceres::Manifold {
public:
    virtual ~EigenQuaternionParameterization() {}

    // Overridden methods from ceres::Manifold
    int AmbientSize() const override { return 4; }  // Quaternion hat 4 Komponenten
    int TangentSize() const override { return 3; }  // 3 Freiheitsgrade für Rotation

    virtual bool Plus(const double* x, const double* delta, double* x_plus_delta) const override;
    virtual bool Minus(const double* y, const double* x, double* y_minus_x) const override;
    virtual bool MinusJacobian(const double* x, double* jacobian) const override;
    virtual bool PlusJacobian(const double* x, double* jacobian) const override;

private:
    void EigenQuaternionProduct(const double z[4], const double w[4], double zw[4]) const;
};


}  // namespace camodocal

#endif