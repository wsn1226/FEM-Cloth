#include <collision_detection_cloth_sphere.h>
#include <iostream>

//  q - generalized coordinates for the FEM system
//  center - the position of the sphere center in the world space
//  radius - the radius of the sphere in the world space
// Output:
//  cloth_index - the indices of all vertices currently in contact with the sphere
//  normals - the outward facing contact normals for each contacting vertex.
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius)
{

    cloth_index.clear();
    normals.clear();

    Eigen::Vector3d thisPosition, n;
    double distance;
    for (int i = 0; i < q.rows() / 3; i++)
    {
        thisPosition = q.segment(3 * i, 3);
        distance = (thisPosition - center).norm();
        if (distance <= radius + 0.01)
        {
            n = (thisPosition - center).normalized();
            cloth_index.push_back(i);
            normals.push_back(n);
        }
    }
}