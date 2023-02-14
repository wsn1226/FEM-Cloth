#include <dphi_cloth_triangle_dX.h>

// compute 3x3 deformation gradient

// Input:
//   V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//   element - the 1x3 vertex indices for this tetrahedron
//   X - the 3D position in the underformed space at which to compute the gradient
// Output:
//   dphi - the 3x3 gradient of the the basis functions wrt to X. The i'th row stores d phi_i/dX
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X)
{
    Eigen::Matrix32d T;
    Eigen::MatrixXd T_T_inv(2, 3); //(T_t * T).inv * T_t;
    Eigen::Vector3d X0, X1, X2;
    Eigen::RowVector2d ONE;
    ONE << 1.0, 1.0;

    X0 = V.row(element(0));
    X1 = V.row(element(1));
    X2 = V.row(element(2));
    T.setZero();
    T.block(0, 0, 3, 1) = X1 - X0;
    T.block(0, 1, 3, 1) = X2 - X0;
    T_T_inv = (T.transpose() * T).inverse() * T.transpose();
    dphi.block(0, 0, 1, 3) = -1.0 * ONE * T_T_inv;
    dphi.block(1, 0, 2, 3) = T_T_inv;
}