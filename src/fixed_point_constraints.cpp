#include <fixed_point_constraints.h>
#include <algorithm>

typedef Eigen::Triplet<double> T;
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices)
{

    std::vector<T> P_entries;
    P.resize(q_size - 3 * indices.size(), q_size);
    P.setZero();
    int a = 0;
    for (int i = 0; i < (q_size / 3); i++)
    {
        if (std::find(indices.begin(), indices.end(), i) == indices.end())
        {
            for (int j = 0; j < 3; j++)
            {
                P_entries.push_back(T(3 * a + j, 3 * i + j, 1));
            }
            a += 1;
        }
    }
    P.setFromTriplets(P_entries.begin(), P_entries.end());
}