#include <Rcpp.h>
#include <RcppEigen.h>

Eigen::MatrixXd getRows(const Eigen::Map<Eigen::MatrixXd> &mat,
                        const Eigen::Map<Eigen::VectorXi> &rows) {
    Eigen::MatrixXd ret(rows.size(), mat.cols());
    for (int i = 0; i < rows.size(); i++) {
        ret.row(i) = mat.row(rows(i)-1);
    }
    return(ret);
}

//[[Rcpp::depends(RcppEigen)]]
//' Logistic gradient
//'
//' @param data A matrix with the data to use
//' @param theta The coefficient vector
//' @param idx The rows to use (base 0 for now)
//' @export
//[[Rcpp::export]]
Eigen::VectorXd logisticG(const Eigen::Map<Eigen::VectorXd> &theta,
                          const Eigen::Map<Eigen::MatrixXd> &data,
                          const Eigen::Map<Eigen::VectorXi> &idx) {
    int n = idx.size();
    Eigen::MatrixXd sub_data = getRows(data, idx);
    Eigen::VectorXd y = sub_data.col(0);
    Eigen::MatrixXd x = sub_data.rightCols(sub_data.cols()-1);
    Eigen::VectorXd eta = 1/(1+exp(-(x*theta).array()));
    return x.transpose()*(eta-y)/n;
}

//[[Rcpp::depends(RcppEigen)]]
//' Gaussian gradient
//'
//' @param data A matrix with the data to use
//' @param theta The coefficient vector
//' @param idx The rows to use (base 0 for now)
//' @export
//[[Rcpp::export]]
Eigen::VectorXd gaussianG(const Eigen::Map<Eigen::VectorXd> &theta,
                          const Eigen::Map<Eigen::MatrixXd> &data,
                          const Eigen::Map<Eigen::VectorXi> &idx) {
    int n = idx.size();
    Eigen::MatrixXd sub_data = getRows(data, idx);
    Eigen::VectorXd y = sub_data.col(0);
    Eigen::MatrixXd x = sub_data.rightCols(sub_data.cols()-1);
    Eigen::VectorXd eta = x*theta;
    return x.transpose()*(eta-y)/n;
}