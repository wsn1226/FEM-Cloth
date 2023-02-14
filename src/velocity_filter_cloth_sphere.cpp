#include <velocity_filter_cloth_sphere.h>

// Input:
//   qdot - the 3nx1 generalized velocities of the cloth mesh
//   index - a list of collision vertex indices from the collision detector
//   normals - a list of collision normals from the collision detector
// Output:
//   qdot- the filtered 3nx1 generalized velocities
void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices,
                                  const std::vector<Eigen::Vector3d> &normals)
{
    Eigen::Vector3d filtered_v, current_v, n;
    for (int i = 0; i < indices.size(); i++)
    {
        n = normals[i];
        unsigned int ind = indices[i];
        current_v = qdot.segment(3 * ind, 3);
        if (n.transpose() * current_v < 0)
        {
            qdot.segment(3 * ind, 3) = current_v - n * (n.transpose() * current_v);
        }
    }
}