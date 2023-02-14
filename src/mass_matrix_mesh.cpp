
#include <mass_matrix_mesh.h>
// Input:
//   q - generalized coordinates for the FEM system
//   V - the nx3 matrix of undeformed vertex positions
//   F - the mx3 matrix of triangle-vertex indices
//   density - the density of the cloth material
//   areas - the mx1 vector of undeformed triangle areas
// Output:
//   M - sparse mass matrix for the entire mesh
typedef Eigen::Triplet<double> Tri;
void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q,
                      Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                      double density, Eigen::Ref<const Eigen::VectorXd> areas)
{
    std::vector<Tri> M_entries;
    Eigen::RowVectorXi ele(3);
    ele.setZero();
    M.resize(q.size(), q.size());
    M.setZero();

    Eigen::Matrix3d I;
    I.setIdentity();
    I /= 120.0;

    for (int k = 0; k < F.rows(); k++)
    {
        int ind0 = F(k, 0);
        int ind1 = F(k, 1);
        int ind2 = F(k, 2);
        ele << ind0, ind1, ind2;
        Eigen::Matrix99d M0;
        M0 << I * 2.0, I, I,
            I, I * 2.0, I,
            I, I, I * 2.0;
        M0 *= (density * areas(k));
        M0 *= 20.0; // set so that the effect is better
        // std::cout << M0(0, 0) << std::endl;
        //-----------------------------FIRST LINE---------------------------------
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                M_entries.push_back(Tri(3 * ind0 + i, 3 * ind0 + j, M0(i, j)));
                M_entries.push_back(Tri(3 * ind0 + i, 3 * ind1 + j, M0(i, j + 3)));
                M_entries.push_back(Tri(3 * ind0 + i, 3 * ind2 + j, M0(i, j + 6)));

                //-----------------------------SECOND LINE---------------------------------

                M_entries.push_back(Tri(3 * ind1 + i, 3 * ind0 + j, M0(i + 3, j)));
                M_entries.push_back(Tri(3 * ind1 + i, 3 * ind1 + j, M0(i + 3, j + 3)));
                M_entries.push_back(Tri(3 * ind1 + i, 3 * ind2 + j, M0(i + 3, j + 6)));

                //-----------------------------THIRD LINE---------------------------------

                M_entries.push_back(Tri(3 * ind2 + i, 3 * ind0 + j, M0(i + 6, j)));
                M_entries.push_back(Tri(3 * ind2 + i, 3 * ind1 + j, M0(i + 6, j + 3)));
                M_entries.push_back(Tri(3 * ind2 + i, 3 * ind2 + j, M0(i + 6, j + 6)));
            }
        }
    }
    // std::cout << "heyhey" << std::endl;
    M.setFromTriplets(M_entries.begin(), M_entries.end());
}
