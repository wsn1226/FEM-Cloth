#include <V_spring_particle_particle.h>

void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0, Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness)
{
    V = 0.0;
    Eigen::Vector3d rel_pos = q1 - q0;
    double rel_len = sqrt(rel_pos.transpose() * rel_pos) - l0;
    V = 0.5 * stiffness * rel_len * rel_len;
}