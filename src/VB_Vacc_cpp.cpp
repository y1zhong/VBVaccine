#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

class BayesLinearReg{
private:
  
  struct LinearRegData{
    int num_reports;
    int num_AEs;
    int num_confounders;
    
    arma::mat A;
    arma::mat X;
    arma::vec V;
    arma::vec nn;
  } dat;
  
  struct HyperParas{
    double a_alpha;
    double b_alpha;
    arma::mat mu_alpha;
    double nu;
    double eta2;
  } hyperparas;
  
  
  struct VBLinearRegParas{
    arma::vec E_beta;
    arma::vec Var_beta;
    arma::vec E_beta_sq;
    
    arma::vec E_inv_sigma2_beta;
    double E_inv_a_beta;
    
    arma::mat E_alpha;
    arma::mat E_inv_sigma2_alpha;
    arma::mat cov_mat_j;
      
    arma::mat Var_alpha;
    arma::mat E_alpha_sq;
    
    arma::mat E_omega;
    arma::mat E_phi;
    
    arma::vec log_det_cov_alpha;
    arma::mat log_cosh_b;
    double ELBO;
    
    arma::mat temp;
  } vb_paras;
  
  
  struct VBProfile{
    arma::vec ELBO;
    arma::mat E_beta_mat;
  } vb_profile;
  
  struct VBControl{
    int max_iter;
    double para_diff_tol;
    int ELBO_stop;
    double ELBO_diff_tol;
    int verbose;
    int save_profile;
    int total_profile;
  } vb_control;
  
  int iter;
  
public:
  void set_method(){
    std::cout << "Mean Field Variational Bayes" << std::endl;
  };
  
  
  void load_data(const arma::mat& in_A, const arma::mat& in_X, const arma::vec& in_V, const arma::vec& in_nn){
    std::cout << "load data" << std::endl;
    dat.A = in_A;
    dat.X = in_X;
    dat.V = in_V; 
    dat.nn = in_nn;

    dat.num_AEs = dat.A.n_cols;
    dat.num_confounders = dat.X.n_cols;
    dat.num_reports = dat.X.n_rows;
  };
  
  void set_hyperparas(const double& in_a_alpha, const double& in_b_alpha,
                      const arma::mat& in_mu_alpha,
                      const double& in_nu,const double& in_eta2){
    hyperparas.a_alpha = in_a_alpha;
    hyperparas.b_alpha = in_b_alpha;
    hyperparas.mu_alpha = in_mu_alpha;
    hyperparas.nu = in_nu;
    hyperparas.eta2 = in_eta2;
  };
  
  void set_paras_initial_values(const arma::vec& in_beta,
                                const arma::mat& in_alpha,
                                const arma::vec& in_inv_sigma2_beta,
                                const arma::mat& in_inv_sigma2_alpha,
                                const double& in_inv_a_beta){
    std::cout << "set initial values" << std::endl;
    
    vb_paras.E_inv_sigma2_beta = in_inv_sigma2_beta;
    vb_paras.Var_beta = 1/vb_paras.E_inv_sigma2_beta;
    vb_paras.E_beta=in_beta;
    vb_paras.E_beta_sq=vb_paras.E_beta%vb_paras.E_beta + vb_paras.Var_beta;
    vb_paras.E_inv_a_beta = in_inv_a_beta;
    
    //vb_paras.E_inv_sigma2_beta.ones(dat.num_AEs);
    
    vb_paras.E_inv_sigma2_alpha = in_inv_sigma2_alpha;
    vb_paras.Var_alpha = 1/in_inv_sigma2_alpha;
    
    vb_paras.E_alpha=in_alpha;
    vb_paras.E_alpha_sq= vb_paras.Var_alpha + vb_paras.E_alpha%vb_paras.E_alpha;
    
    //vb_paras.E_inv_sigma2_alpha.ones(dat.num_confounders, dat.num_AEs);
    
    vb_paras.E_phi.zeros(dat.num_reports, dat.num_AEs);
    update_E_phi();
    //std::cout << find_nonfinite(vb_paras.E_phi)<< std::endl;
    
    vb_paras.E_omega.zeros(dat.num_reports, dat.num_AEs);
    update_E_omega();
    //std::cout << find_nonfinite(vb_paras.E_omega)<< std::endl;
    vb_paras.log_det_cov_alpha.zeros(dat.num_AEs);
  };
  
  
  void set_vb_control(int in_max_iter, 
                      double in_para_diff_tol, 
                      int in_ELBO_stop,
                      double in_ELBO_diff_tol,
                      int in_verbose,
                      int in_save_profile){
    vb_control.max_iter = in_max_iter;
    vb_control.para_diff_tol = in_para_diff_tol;
    vb_control.ELBO_stop = in_ELBO_stop;
    vb_control.ELBO_diff_tol = in_ELBO_diff_tol;
    vb_control.verbose = in_verbose;
    vb_control.save_profile = in_save_profile;
    if(vb_control.save_profile > 0){
      vb_control.total_profile = vb_control.max_iter/vb_control.save_profile;
    } else{
      vb_control.total_profile = 0;
    }
  };
  
  
  
  void update_Var_beta(){
    for(int j=0; j<dat.num_AEs; j++){
      vb_paras.Var_beta(j) = 1/(accu(dat.V % vb_paras.E_omega.col(j)) + vb_paras.E_inv_sigma2_beta(j));
    }
    // vb_paras.Var_beta = 1/(vb_paras.E_omega.t() * dat.V+ vb_paras.E_inv_sigma2_beta);
     //std::cout << vb_paras.Var_beta<< std::endl;
  };
  
  void update_E_beta(){
    
    // for(int j=0; j<dat.num_AEs; j++){
    //   vb_paras.E_beta(j) = accu((dat.A.col(j) - dat.nn/2 -  vb_paras.E_omega.col(j)%(dat.X*vb_paras.E_alpha.col(j)))%dat.V);
    //   vb_paras.E_beta(j) *= vb_paras.Var_beta(j);
    //   vb_paras.E_beta_sq(j) = vb_paras.E_beta(j)*vb_paras.E_beta(j) + vb_paras.Var_beta(j);
    //   //std::cout << vb_paras.E_beta(j) << std::endl;
    // }
    vb_paras.temp = (dat.A.each_col() - dat.nn/2) - vb_paras.E_omega%(dat.X*vb_paras.E_alpha);
    vb_paras.temp = vb_paras.temp.each_col()%dat.V;
    vb_paras.E_beta = conv_to< colvec >::from(sum(vb_paras.temp, 0));
    vb_paras.E_beta %= vb_paras.Var_beta;
    vb_paras.E_beta_sq = vb_paras.E_beta%vb_paras.E_beta + vb_paras.Var_beta;
    // std::cout << vb_paras.E_beta<< std::endl;
  };
  
  
  void update_E_inv_sigma2_beta(){
    // for(int j=0; j<dat.num_AEs; j++){
    //   vb_paras.E_inv_sigma2_beta(j) = 1/(vb_paras.E_beta_sq(j)/2 + vb_paras.E_inv_a_beta);
    // }
    vb_paras.E_inv_sigma2_beta = 1/(vb_paras.E_beta_sq/2 + hyperparas.nu*vb_paras.E_inv_a_beta);
  };
  
  void update_E_inv_sigma2_alpha(){
    // for(int l=0; l<dat.num_confounders; l++){
    //   vb_paras.E_inv_sigma2_alpha(l) = (hyperparas.a_alpha + dat.num_AEs/2) / (hyperparas.b_alpha + accu(vb_paras.E_alpha_sq.row(l))/2);
    // }
    vb_paras.E_inv_sigma2_alpha = (hyperparas.a_alpha + 1/2) / (hyperparas.b_alpha + (vb_paras.E_alpha_sq + hyperparas.mu_alpha%hyperparas.mu_alpha - 2*vb_paras.E_alpha%hyperparas.mu_alpha)/2);
    //std::cout << vb_paras.E_inv_sigma2_alpha<< std::endl;
  };
  
  void update_E_alpha(){
    //arma::vec ind = linspace(0, dat.num_AEs-1, dat.num_AEs);
    for(int j=0; j<dat.num_AEs; j++){
      // vb_paras.cov_mat_j = dat.X.t()*diagmat(vb_paras.E_omega.col(j))*dat.X + diagmat(vb_paras.E_inv_sigma2_alpha.col(j));
      vb_paras.cov_mat_j = dat.X.each_col()%vb_paras.E_omega.col(j);
      vb_paras.cov_mat_j =  vb_paras.cov_mat_j.t()*dat.X;
      vb_paras.cov_mat_j.diag() += vb_paras.E_inv_sigma2_alpha.col(j);
      
      vb_paras.cov_mat_j = inv(vb_paras.cov_mat_j);
      vb_paras.Var_alpha.col(j) = vb_paras.cov_mat_j.diag();
      vb_paras.E_alpha.col(j) = vb_paras.cov_mat_j * (dat.X.t() * (dat.A.col(j) - dat.nn/2 - vb_paras.E_omega.col(j)%dat.V*vb_paras.E_beta(j)) + vb_paras.E_inv_sigma2_alpha.col(j) % hyperparas.mu_alpha.col(j));
    }
  };
  
  
  void update_E_alpha_sq(){
    vb_paras.E_alpha_sq = vb_paras.Var_alpha + vb_paras.E_alpha % vb_paras.E_alpha;
  };
  
  void update_E_inv_a_beta(){
    vb_paras.E_inv_a_beta = ((dat.num_AEs*hyperparas.nu+1)/2)/(hyperparas.nu*accu(vb_paras.E_inv_sigma2_beta)+1/hyperparas.eta2);
  };
  
  void update_E_phi(){
    // for(int i=0; i<dat.num_reports; i++){
    //   for(int j=0; j<dat.num_AEs; j++){
    //     vb_paras.E_phi(i,j) = accu(dat.X.row(i)*vb_paras.E_alpha.col(j)) + dat.V(i)*vb_paras.E_beta(j);
    //     //vb_paras.E_log_cosh_phi(i,j) = log(cosh(vb_paras.E_phi(i,j)/2));
    //   }
    // }
    vb_paras.E_phi = dat.X*vb_paras.E_alpha + dat.V*vb_paras.E_beta.t();
    //std::cout << find_nonfinite(vb_paras.E_alpha)<< std::endl;
    //std::cout << find_nonfinite(vb_paras.E_beta)<< std::endl;
  }
  
  void update_E_omega(){
    // for(int i=0; i<dat.num_reports; i++){
    //   for(int j=0; j<dat.num_AEs; j++){
    //     double value = dat.nn(i)*tanh(vb_paras.E_phi(i,j)/2)/(2*vb_paras.E_phi(i,j));
    //     if(Rcpp::traits::is_nan<REALSXP>(value)){
    //       vb_paras.E_omega(i,j) = 0.25;
    //     } else{
    //       vb_paras.E_omega(i,j) = value;
    //     }
    // 
    //   }
    // }
    vb_paras.temp = tanh(vb_paras.E_phi/2)/(2*vb_paras.E_phi);
    vb_paras.E_omega = vb_paras.temp.each_col()%dat.nn;
    //uvec indices = find_nonfinite(vb_paras.E_omega);
    //vb_paras.E_omega.elem( find_nonfinite(vb_paras.E_omega) ).fill(1000);
  };
  
  void update_log_cosh_b(){
    vb_paras.temp = log(cosh(vb_paras.E_phi/2));
    vb_paras.log_cosh_b = vb_paras.temp.each_col() % dat.nn;
  };
  
  
  void update_log_det_cov_mat(){
    for(int j=0; j<dat.num_AEs; j++){
      // vb_paras.cov_mat_j = dat.X.t()*diagmat(vb_paras.E_omega.col(j))*dat.X + diagmat(vb_paras.E_inv_sigma2_alpha.col(j));
      vb_paras.cov_mat_j = dat.X.each_col()%vb_paras.E_omega.col(j);
      vb_paras.cov_mat_j =  vb_paras.cov_mat_j.t()*dat.X;
      vb_paras.cov_mat_j.diag() += vb_paras.E_inv_sigma2_alpha.col(j);
      
      vb_paras.log_det_cov_alpha(j) = log_det_sympd(vb_paras.cov_mat_j);
      // double sign;
      // bool ok = log_det(vb_paras.log_det_cov_alpha(j), sign, cov_mat_j);
    }
  }
  
  void update_ELBO(){
    vb_paras.ELBO = accu((dat.A.each_col() - dat.nn/2)%vb_paras.E_phi);
    vb_paras.ELBO += accu(log(vb_paras.Var_beta))/2;
    vb_paras.ELBO -= ((1+hyperparas.nu)/2)*accu(log(vb_paras.E_beta_sq/2 + hyperparas.nu*vb_paras.E_inv_a_beta));
    vb_paras.ELBO += hyperparas.nu*vb_paras.E_inv_a_beta * accu(vb_paras.E_inv_sigma2_beta);
    update_log_det_cov_mat();
    vb_paras.ELBO -= accu(vb_paras.log_det_cov_alpha)/2;
    vb_paras.ELBO -= accu(log(hyperparas.b_alpha + (vb_paras.E_alpha_sq + hyperparas.mu_alpha%hyperparas.mu_alpha - 2*vb_paras.E_alpha%hyperparas.mu_alpha)/2)*(hyperparas.a_alpha+0.5));
    vb_paras.ELBO -= (hyperparas.nu*dat.num_AEs+1)/2*log(hyperparas.nu*accu(vb_paras.E_inv_sigma2_beta)+1/hyperparas.eta2);
    //vb_paras.ELBO -= accu(log(cosh(vb_paras.E_phi/2)));
    update_log_cosh_b();
    vb_paras.ELBO -= accu(vb_paras.log_cosh_b);
  };
  
  
  double compute_paras_diff(arma::vec& beta, arma::vec& beta_prev){
    arma::vec temp = beta - beta_prev;
    return accu(temp%temp)/beta.n_elem;
  };
  
  void initialize_vb_profile(){
    if(vb_control.save_profile>0){
      vb_profile.ELBO.zeros(vb_control.total_profile);
      //vb_profile.E_beta_arma::mat.zeros(vb_paras.E_beta.n_elem,vb_control.total_profile);
    }
  }
  
  void save_vb_profile(){
    if(vb_control.save_profile > 0){
      if(iter%vb_control.save_profile==0){
        int profile_iter = iter/vb_control.save_profile;
        if(vb_control.ELBO_stop==0){
          update_ELBO();
        }
        vb_profile.ELBO(profile_iter) = vb_paras.ELBO;
        //vb_profile.E_beta_arma::mat.col(profile_iter) = vb_paras.E_beta;
      }
    }
  }
  
  void monitor_vb(){
    if(vb_control.verbose > 0){
      if(iter%vb_control.verbose==0){
        if(vb_control.ELBO_stop==0){
          update_ELBO();
        }
        if(iter%2000==0){
          std::cout << "iter: " << iter <<  " ELBO: "<< vb_paras.ELBO << std::endl;
        }
      }
    }
  }
  
  void run_mfvb(){
    initialize_vb_profile();
    for(iter=0; iter<vb_control.max_iter; iter++){
      arma::vec E_beta_prev = vb_paras.E_beta;
      update_Var_beta();
      update_E_beta();
      
      update_E_alpha();
      update_E_alpha_sq();
      
      
      update_E_inv_sigma2_alpha();
      //std::cout << vb_paras.E_inv_sigma2_alpha<< std::endl;
      update_E_inv_sigma2_beta();
      //std::cout << vb_paras.E_inv_sigma2_beta<< std::endl;
      
      update_E_inv_a_beta();
      //std::cout << vb_paras.E_inv_a_beta<< std::endl;
      
      update_E_phi();
      update_E_omega();
      
      if(vb_control.ELBO_stop == 0){
        if(compute_paras_diff(vb_paras.E_beta,E_beta_prev) < vb_control.para_diff_tol){
          update_ELBO();
          save_vb_profile();
          monitor_vb();
          break;
        }
      } else {
        double ELBO_prev = vb_paras.ELBO;
        update_ELBO();
        if(abs(vb_paras.ELBO - ELBO_prev) < vb_control.ELBO_diff_tol){
          save_vb_profile();
          monitor_vb();
          break;
        }
      }
      
      save_vb_profile();
      monitor_vb();
    }
  };
  
  List get_vb_data(){
    return List::create(Named("X") = dat.X,
                        Named("V") = dat.V,
                        Named("A") = dat.A,
                        Named("nn") = dat.nn);
  };
  
  List get_vb_hyperparam(){
    return List::create(Named("a_alpha") =hyperparas.a_alpha,
                        Named("b_alpha") = hyperparas.b_alpha,
                        Named("mu_alpha") = hyperparas.mu_alpha,
                        Named("nu") = hyperparas.nu,
                        Named("eta2") = hyperparas.eta2);
  };
  
  List get_vb_post_mean(){
    return List::create(Named("beta") = vb_paras.E_beta,
                        Named("alpha") = vb_paras.E_alpha,
                        Named("inv_sigma2_beta") = vb_paras.E_inv_sigma2_beta,
                        Named("inv_sigma2_alpha") = vb_paras.E_inv_sigma2_alpha,
                        Named("inv_a_beta") = vb_paras.E_inv_a_beta);
  };
  
  
  
  List get_vb_trace(){
    int actual_profile_iter = 1;
    if(iter == 0){
      iter = 1;
    }
    if(vb_control.save_profile>0){
      actual_profile_iter = iter/vb_control.save_profile;
    }
    arma::uvec iters = linspace<uvec>(1,iter,actual_profile_iter);
    return List::create(Named("iters") = iters,
                        Named("ELBO") = vb_profile.ELBO.rows(0,actual_profile_iter-1));
    
    //Named("E_beta_arma::mat") = vb_profile.E_beta_arma::mat) ;
  }
  
  List get_vb_control(){
    return List::create(Named("max_iter")= vb_control.max_iter,
                        Named("para_diff_tol") = vb_control.para_diff_tol,
                        Named("ELBO_stop") = vb_control.ELBO_stop,
                        Named("ELBO_diff_tol") = vb_control.ELBO_diff_tol,
                        Named("verbose") = vb_control.verbose,
                        Named("save_profile") = vb_control.save_profile,
                        Named("total_profile") = vb_control.total_profile);
  };
  
  int get_iter(){
    return iter;
  };
};


// ' @useDynLib VBVaccine, .registration = TRUE
// ' @importFrom Rcpp evalCpp
// ' @export
// [[Rcpp::export]]
List Bayes_Vacc_cpp(arma::mat& A, arma::mat& X, arma::vec& V,  arma::vec& nn, 
                arma::vec initial_beta, 
                arma::mat initial_alpha,
                arma::vec initial_inv_sigma2_beta, 
                arma::mat initial_inv_sigma2_alpha,
                double initial_inv_a_beta,
                arma::mat mu_alpha,
                double a_alpha = 0.1,
                double b_alpha = 0.1,
                double nu = 1,
                double eta2 = 1,
                int max_iter = 1000,
                double paras_diff_tol = 1e-6,
                int ELBO_stop = 1,
                double ELBO_diff_tol = 1e-6,
                int verbose = 0,
                int save_profile = 1){
  
  wall_clock timer;
  timer.tic();
  BayesLinearReg model;
  
  
  model.load_data(A,X,V, nn);
  model.set_hyperparas(a_alpha, b_alpha, mu_alpha, nu, eta2);
  
  
  model.set_vb_control(max_iter,
                       paras_diff_tol,
                       ELBO_stop,
                       ELBO_diff_tol,
                       verbose,
                       save_profile);
  
  //std::cout << "set control done" << std::endl;
  
  model.set_paras_initial_values(initial_beta, initial_alpha,
                                 initial_inv_sigma2_beta,
                                 initial_inv_sigma2_alpha,
                                 initial_inv_a_beta);
  
  //std::cout << "set initial values" << std::endl;
  
  model.run_mfvb(); 
  
  double elapsed = timer.toc();
  
  List output;
  
  output = List::create(Named("post_mean") = model.get_vb_post_mean(),
                        Named("data") = model.get_vb_data(),
                        Named("hyper") = model.get_vb_hyperparam(),
                        Named("iter") = model.get_iter(),
                        Named("trace") = model.get_vb_trace(),
                        Named("vb_control") = model.get_vb_control(),
                        Named("elapsed") = elapsed);
  
  return output;
  
}

