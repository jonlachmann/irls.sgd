#include <Rcpp.h>
#include <RcppEigen.h>

//[[Rcpp::depends(RcppEigen)]]
//' QR for least squares
//'
//' @param x The design matrix
//' @param w The weights vector
//' @param z The working response vector
//' @export
//[[Rcpp::export]]
Eigen::VectorXd qrls(const Eigen::Map<Eigen::MatrixXd> &x,
              const Eigen::Map<Eigen::VectorXd> &w,
              const Eigen::Map<Eigen::VectorXd> &z)
    {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> PQR;
    Eigen::VectorXd beta;
    Eigen::VectorXd effects;
    Eigen::MatrixXd Rinv;
    int nvars = x.cols();

    // Decompose the model matrix
    PQR.compute(w.asDiagonal() * x);
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd>::PermutationType Pmat = (PQR.colsPermutation());
    int rank = PQR.rank();

    if (rank == nvars) {
        // Full rank case
        beta = PQR.solve((z.array() * w.array()).matrix());
    } else {
        // Rank deficit case
        Rinv = (PQR.matrixQR().topLeftCorner(rank, rank).triangularView<Eigen::Upper>().
                                              solve(Eigen::MatrixXd::Identity(rank, rank)));
        effects = PQR.householderQ().adjoint() * (z.array() * w.array()).matrix();

        beta = Rinv * effects.head(rank);
        beta = Pmat * beta;

        // Create fitted values from effects
        effects.tail(x.rows() - rank).setZero();
    }
    return beta;
}