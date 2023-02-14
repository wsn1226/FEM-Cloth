#include <dV_spring_particle_particle_dq.h>

void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness)
{

    Eigen::Vector3d rel_pos = q1 - q0;
    double current_len = sqrt(rel_pos.transpose() * rel_pos);
    double rel_len = current_len - l0;

    f << -1.0 * rel_pos, rel_pos;
    f *= (stiffness * rel_len / current_len);
}