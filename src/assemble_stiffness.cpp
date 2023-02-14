#include <assemble_stiffness.h>
#include <iostream>

typedef Eigen::Triplet<double> Tri;
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                        Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                        double mu, double lambda)
{

    std::vector<Tri> K_entries;
    Eigen::RowVectorXi ele(3);
    ele.setZero();
    K.resize(qdot.size(), qdot.size());
    K.setZero();

    for (int k = 0; k < F.rows(); k++)
    {
        Eigen::Matrix<double, 1, 9> tmp_row;
        tmp_row = dX.row(k); // ei is the triangle index.
        Eigen::Matrix3d dphi = Eigen::Map<const Eigen::Matrix3d>(tmp_row.data());
        int ind0 = F(k, 0);
        int ind1 = F(k, 1);
        int ind2 = F(k, 2);
        ele << ind0, ind1, ind2;
        Eigen::Matrix99d K0;
        d2V_membrane_corotational_dq2(K0, q, dphi, V, ele, a0(k), mu, lambda);
        // std::cout << K0.coeff(0, 0) << std::endl;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                //-----------------------------FIRST LINE---------------------------------
                K_entries.push_back(Tri(3 * ind0 + i, 3 * ind0 + j, K0(i, j)));
                K_entries.push_back(Tri(3 * ind0 + i, 3 * ind1 + j, K0(i, j + 3)));
                K_entries.push_back(Tri(3 * ind0 + i, 3 * ind2 + j, K0(i, j + 6)));

                //-----------------------------SECOND LINE---------------------------------

                K_entries.push_back(Tri(3 * ind1 + i, 3 * ind0 + j, K0(i + 3, j)));
                K_entries.push_back(Tri(3 * ind1 + i, 3 * ind1 + j, K0(i + 3, j + 3)));
                K_entries.push_back(Tri(3 * ind1 + i, 3 * ind2 + j, K0(i + 3, j + 6)));

                //-----------------------------THIRD LINE---------------------------------

                K_entries.push_back(Tri(3 * ind2 + i, 3 * ind0 + j, K0(i + 6, j)));
                K_entries.push_back(Tri(3 * ind2 + i, 3 * ind1 + j, K0(i + 6, j + 3)));
                K_entries.push_back(Tri(3 * ind2 + i, 3 * ind2 + j, K0(i + 6, j + 6)));
            }
        }
    }
    // std::cout << "heyhey" << std::endl;
    K.setFromTriplets(K_entries.begin(), K_entries.end());
    K *= -1.0;
};