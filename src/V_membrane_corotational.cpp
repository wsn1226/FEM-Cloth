#include <V_membrane_corotational.h>
#include <iostream>
// Allowed to use libigl SVD or Eigen SVD for this part

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
// Output:
//  energy- the per-triangle potential energy (the linear model described in the README).
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX,
                             Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area,
                             double mu, double lambda)
{
    Eigen::Vector3d x0, x1, x2, X0, X1, X2, n, N, S;
    Eigen::Matrix34d LeftF;  // 3*4
    Eigen::Matrix43d RightF; // 4*3
    Eigen::Matrix3d F;
    x0 = q.block(3 * element(0), 0, 3, 1);
    x1 = q.block(3 * element(1), 0, 3, 1);
    x2 = q.block(3 * element(2), 0, 3, 1);
    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));
    LeftF.block(0, 0, 3, 1) = x0;
    LeftF.block(0, 1, 3, 1) = x1;
    LeftF.block(0, 2, 3, 1) = x2;
    n = (x0 - x1).cross(x1 - x2);
    N = (X0 - X1).cross(X1 - X2);
    n.normalize();
    N.normalize();
    RightF.block(0, 0, 3, 3) = dX;
    RightF.block(3, 0, 3, 1) = N.transpose();
    LeftF.block(0, 0, 3, 1) = x0;
    LeftF.block(0, 1, 3, 1) = x1;
    LeftF.block(0, 2, 3, 1) = x2;
    LeftF.block(0, 3, 3, 1) = n;
    F = LeftF * RightF;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    S = svd.singularValues();

    energy = mu * (pow(S[0] - 1.0, 2) + pow(S[1] - 1.0, 2) + pow(S[2] - 1.0, 2)) + 0.5 * lambda * pow(S[0] + S[1] + S[2] - 3.0, 2);
    energy *= area;
}
