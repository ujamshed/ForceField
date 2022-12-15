#include <armadillo>

int main()
{
    arma::mat mat_1 = {{2, 0, 6}, {1, -1, 2}};
    arma::mat mat_2 = {{8, -3, 2}, {4, 1, -5}};
    mat_1.print();
    std::cout << arma::dot(mat_1, mat_2) << std::endl;

    // Vectorize
    arma::vec vector_1_col = mat_1.as_col();
    arma::rowvec vector_1_row = mat_1.as_row();

    vector_1_col.print("Vector 1 as a col");
    vector_1_row.print("Vector 1 as a row");
    std::cout << std::endl;
    // Gives me a matrix because 6x1 * 1x6;
    std::cout << vector_1_col*vector_1_row << std::endl;

    // You can recreate a matrix after sending it to a column vector by reshaping
    arma::mat mat_1_recreated = arma::reshape(vector_1_col, 2, 3);
    mat_1_recreated.print("Recreated Matrix 1");
    return 0;
}