#include <dV_cloth_gravity_dq.h>

void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g)
{
    Eigen::VectorXd g_all;
    g_all.resize(fg.size());
    for (int i = 0; i < fg.size() / 3; i++)
    {
        g_all.segment(i * 3, 3) = g;
    }

    fg = -1.0 * M * g_all;
}
