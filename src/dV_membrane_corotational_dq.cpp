#include <dV_membrane_corotational_dq.h>
#include <iostream>
//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dX/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
// Output:
//  f - the per-triangle gradient of the membrane potential energy (the linear model described in the README).
void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX,
                                 Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area,
                                 double mu, double lambda)
{

    // Deformation Gradient
    Eigen::Matrix3d dx; // deformed tangent matrix
    Eigen::Matrix3d U, W, Middledsi_dF;
    Eigen::Vector3d S, x0, x1, x2, X0, X1, X2, n, N, ntiled, x1minusx0, x2minusx0;
    Eigen::Matrix34d LeftF;  // 3*4
    Eigen::Matrix43d RightF; // 4*3
    Eigen::Matrix3d F;
    Eigen::Matrix99d B, dFdq;
    Eigen::Matrix93d N_dFdq;
    Eigen::Matrix39d n_dFdq;
    x0 = q.block(3 * element(0), 0, 3, 1);
    x1 = q.block(3 * element(1), 0, 3, 1);
    x2 = q.block(3 * element(2), 0, 3, 1);
    X0 = V.block(element(0), 0, 1, 3).transpose();
    X1 = V.block(element(1), 0, 1, 3).transpose();
    X2 = V.block(element(2), 0, 1, 3).transpose();
    x1minusx0 = x1 - x0;
    x2minusx0 = x2 - x0;
    ntiled = x1minusx0.cross(x2minusx0);
    N = (X1 - X0).cross(X2 - X0);
    n = ntiled.normalized();
    N.normalize();
    RightF.block(0, 0, 3, 3) = dX;
    RightF.block(3, 0, 1, 3) = N.transpose();
    LeftF.block(0, 0, 3, 1) = x0;
    LeftF.block(0, 1, 3, 1) = x1;
    LeftF.block(0, 2, 3, 1) = x2;
    LeftF.block(0, 3, 3, 1) = n;
    F = LeftF * RightF;

    // TODO: SVD Here
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();
    // Fix for inverted elements (thanks to Danny Kaufman)
    double det = S(0) * S(1);

    if (det <= -1e-10)
    {
        if (S(0) < 0)
            S(0) *= -1;
        if (S(1) < 0)
            S(1) *= -1;
        if (S(2) < 0)
            S(2) *= -1;
    }

    if (U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }

    if (W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }

    B.setZero();
    B.block(0, 0, 3, 1) = dX.block(0, 0, 1, 3).transpose();
    B.block(3, 1, 3, 1) = dX.block(0, 0, 1, 3).transpose();
    B.block(6, 2, 3, 1) = dX.block(0, 0, 1, 3).transpose();

    B.block(0, 3, 3, 1) = dX.block(1, 0, 1, 3).transpose();
    B.block(3, 4, 3, 1) = dX.block(1, 0, 1, 3).transpose();
    B.block(6, 5, 3, 1) = dX.block(1, 0, 1, 3).transpose();

    B.block(0, 6, 3, 1) = dX.block(2, 0, 1, 3).transpose();
    B.block(3, 7, 3, 1) = dX.block(2, 0, 1, 3).transpose();
    B.block(6, 8, 3, 1) = dX.block(2, 0, 1, 3).transpose();

    N_dFdq.setZero();
    N_dFdq.block(0, 0, 3, 1) = N;
    N_dFdq.block(3, 1, 3, 1) = N;
    N_dFdq.block(6, 2, 3, 1) = N;

    Eigen::Matrix3d I;
    Eigen::Matrix39d IOI, IIO;
    Eigen::Matrix3d deltax1_cross, deltax2_cross;
    deltax1_cross.setZero();
    deltax2_cross.setZero();
    I.setIdentity();
    IOI.setZero();
    IIO.setZero();
    IOI.block(0, 0, 3, 3) = -1.0 * I;
    IOI.block(0, 6, 3, 3) = I;

    IIO.block(0, 0, 3, 3) = -1.0 * I;
    IIO.block(0, 3, 3, 3) = I;

    deltax1_cross(0, 1) = -1.0 * x1minusx0.z();
    deltax1_cross(0, 2) = x1minusx0.y();
    deltax1_cross(1, 0) = x1minusx0.z();
    deltax1_cross(1, 2) = -1.0 * x1minusx0.x();
    deltax1_cross(2, 0) = -1.0 * x1minusx0.y();
    deltax1_cross(2, 1) = x1minusx0.x();

    deltax1_cross(0, 1) = -1.0 * x2minusx0.z();
    deltax1_cross(0, 2) = x2minusx0.y();
    deltax1_cross(1, 0) = x2minusx0.z();
    deltax1_cross(1, 2) = -1.0 * x2minusx0.x();
    deltax1_cross(2, 0) = -1.0 * x2minusx0.y();
    deltax1_cross(2, 1) = x2minusx0.x();

    n_dFdq = 1.0 / ntiled.norm() * (I - n * n.transpose()) * (deltax1_cross * IOI - deltax2_cross * IIO);
    dFdq = B + N_dFdq * n_dFdq;
    // TODO: energy model gradient
    Middledsi_dF.setZero();
    Middledsi_dF(0, 0) = 2.0 * mu * (S(0) - 1.0) + lambda * (S(0) + S(1) + S(2) - 3.0);
    Middledsi_dF(1, 1) = 2.0 * mu * (S(1) - 1.0) + lambda * (S(0) + S(1) + S(2) - 3.0);
    Middledsi_dF(2, 2) = 2.0 * mu * (S(2) - 1.0) + lambda * (S(0) + S(1) + S(2) - 3.0);
    dx = U * Middledsi_dF * W.transpose();
    dV.block(0, 0, 3, 1) = dx.block(0, 0, 1, 3).transpose();
    dV.block(3, 0, 3, 1) = dx.block(1, 0, 1, 3).transpose();
    dV.block(6, 0, 3, 1) = dx.block(2, 0, 1, 3).transpose();
    dV = area * dFdq.transpose() * dV;
}