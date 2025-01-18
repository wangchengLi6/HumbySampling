#include <Rcpp.h>
#include <vector>
#include <math.h>
using namespace Rcpp;


int GeneRandomInt(int max_int){
  int toret = R::runif(0,max_int);
  return toret;
}

double runif_c(){
  double toret = R::runif(0,1);
  return toret;
}

double rnorm_c(){
  double toret = R::rnorm(0,1);
  return toret;
}


double runifcos_c(){
  double toret = cos(R::runif(0,1)*3.14159);
  return toret;
}

int RuleIIICheck(int n, double* a1, int infv, double* u,double *v, double* minv, int* p, int* way, int* used) {

  // re-initialize
  for(int i = 0; i<n+1;i++){
    u[i] = 0;
    v[i] = 0;
    minv[i] = 0;
    p[i] = 0;
    way[i] = 0;
    used[i] = 0;
  }


  // 算法模块
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    for (int cs = 0; cs < n + 1; cs++) {
      minv[cs] = infv;
      used[cs] = 0;
    }

    do {
      used[j0] = 1;
      double delta = infv;
      int i0 = p[j0], j1 = 0;
      for (int j = 1; j <= n; ++j) {
        if (!used[j]) {
          double cur = a1[i0*(n+1)+j] - u[i0] - v[j];
          if (cur < minv[j] - 0.000001) {
            minv[j] = cur;
            way[j] = j0;
          }
          if (minv[j] < delta ) {
            delta = minv[j], j1 = j;
          }
        }
      }
      for (int j = 0; j <= n; ++j) {
        if (used[j])
          u[p[j]] += delta, v[j] -= delta;
        else
          minv[j] -= delta;
        j0 = j1;
      }
    } while (p[j0] != 0);


    do {
      int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }


  int toret;
  toret = 1;
  for (int j = 1; j <= n; ++j) {
    if (p[j] != j) {
      toret = 0;
      break;
    }
  }

  return toret;
}

int dis_cal(int M, int n, double* Emat_vec, double* dis_vec) {
  double* disF0T1 = new double[M]();

  for (int rowid = 0; rowid < n; rowid++) {
    double sum_square = 0;
    for (int mth = 0; mth < M; mth++) {
      double mthele = Emat_vec[rowid * M + mth];
      sum_square += pow(mthele,2);
      disF0T1[mth] = pow(1 - mthele, 2) - pow(mthele, 2);
    }
    for (int mth = 0; mth < M; mth++) {
      dis_vec[rowid * M + mth] = sqrt(sum_square + disF0T1[mth]);
    }

  }


  delete[] disF0T1;
  return 0;
}

int iter_Mcate(int* current_point, int* end_point,int dig,int max_dig) {
  if (current_point[dig] < end_point[dig]-1) {
    current_point[dig] += 1;
  }
  else {
    if (dig != max_dig - 1) {
      current_point[dig] = 0;
      iter_Mcate(current_point, end_point, dig + 1, max_dig);
    }
  }
  return 0;
}




// Estmate HUM by check all combinations
double HUMExactRuleOne(int M, int* n_vec, double* emat_vec) {
  // initialization
  int sampid = 0;
  long int negcounter = 0;
  int break_status = 0;
  int max_iter = n_vec[0];
  // allocate the memory and should be deleted on exit.
  int* n_cumsum = new int[M]();
  for (int j = 1; j < M; ++j) {
    n_cumsum[j] = n_vec[j - 1] + n_cumsum[j - 1];
    max_iter *= n_vec[j];
  }
  double* comparison_vec = new double[M]();
  int* start_point = new int[M]();

  // sample part:
  for (int rth = 0; rth < max_iter; rth++) {
    if(rth % 20000 == 0){
      Rcpp::checkUserInterrupt();
    }
    break_status = 0;
    sampid = start_point[0];
    comparison_vec[0] = emat_vec[sampid * M + 0];
    for (int mth = 1; mth < M; mth++) {
      sampid = start_point[mth] + n_cumsum[mth];
      for (int j = 0; j < mth; j++) {
        if (emat_vec[sampid * M + j] >= comparison_vec[j]) {
          negcounter += 1;
          break_status = 1;
          break;
        }
      }
      if (break_status) {
        break;
      }
      else {
        comparison_vec[mth] = emat_vec[sampid * M + mth];
      }
    }
    iter_Mcate(start_point, n_vec, 0, M);
  }


  // deleting and freeing part
  delete[] comparison_vec;
  delete[] n_cumsum;
  delete[] start_point;

  double to_return = 1 - (double)negcounter / (double)max_iter;
  return to_return;
}

// [[Rcpp::export]]
double HUMExactRuleOne_R(int M,  IntegerVector nvec_R, NumericMatrix Ematrix_R){
  int nsum = 0;
  // memory allocation
  int* n_vec = new int[M]();
  for(int i = 0; i < M; i++){
    n_vec[i] = nvec_R[i];
    nsum += nvec_R[i];
  }


  double* Emat_vec = new double[M*nsum];
  // translate R object to c++ object
  for(int i = 0; i < nsum; i++){
    for(int j = 0; j < M; j++){
      Emat_vec[i*M+j] = Ematrix_R(i,j);
    }
  }

  // calculate HUM under rule 1
  double freq_exact = HUMExactRuleOne(M,n_vec,Emat_vec);



  // delete part
  delete[] n_vec;
  delete[] Emat_vec;

  return freq_exact;
}


double HUMExactRuleThree(int M, int* nvec, double* dis_mat_vec) {
  long int pos_counter = 0; // count the number of times the condition is met (positive).
  long int max_iter = nvec[0];
  int sample_id = 0;
  int* n_cumsum_vec = new int[M]();
  // calculate the cumsum of n
  for (int i = 1; i < M; i++) {
    n_cumsum_vec[i] = nvec[i - 1] + n_cumsum_vec[i - 1];
    max_iter *= nvec[i];
  }
  int* start_point = new int[M]();

  // anxilary variables and memory allocation
  double* u = new double[M + 1]();
  double* v = new double[M + 1]();
  double* minv = new double[M + 1]();
  int* p = new int[M + 1]();
  int* way = new int[M + 1]();
  int* used = new int[M + 1]();
  double* sampled_dis_vec_augmented = new double[(M + 1) * (M + 1)]();


  for (int rth = 0; rth < max_iter; rth++) {
    if(rth % 20000 == 0){
      Rcpp::checkUserInterrupt();
    }
    for (int mth = 0; mth < M; mth++) {
      sample_id = start_point[mth] + n_cumsum_vec[mth];
      for (int j = 0; j < M; j++) {
        sampled_dis_vec_augmented[(mth + 1) * (M + 1) + j + 1] = dis_mat_vec[sample_id * M + j];
      }
    }

    pos_counter += RuleIIICheck(M, sampled_dis_vec_augmented, M + 1, u, v, minv, p, way, used);
    iter_Mcate(start_point, nvec, 0, M);
  }



  // delete part
  delete[] sampled_dis_vec_augmented;
  delete[] u;
  delete[] v;
  delete[] minv;
  delete[] way;
  delete[] used;
  delete[] p;

  delete[] start_point;
  delete[] n_cumsum_vec;

  double toret = (double)pos_counter / (double)max_iter;
  return toret;
}

// [[Rcpp::export]]
double HUMExactRuleThree_R(int M,IntegerVector nvec_R, NumericMatrix Ematrix_R) {
  // memory allocation
  int* nvec = new int[M];

  // translate R objects to c++ objects
  int nsum = 0;
  for(int m = 0; m < M; m++){
    nvec[m] = nvec_R[m];
    nsum += nvec_R[m];
  }

  // memory allocation 2 (depend on nsum)
  double* emat =new double [nsum*M];
  double* disv = new double[nsum*M]();
  for(int i = 0; i < nsum ; i++){
    for(int j = 0; j < M ; j++){
      emat[i*M+j] = Ematrix_R(i,j);
    }
  }

  // calculate the distance matrix
  dis_cal(M, nsum, emat, disv);
  // approximate the HUM under rule 3
  double freq_exact = HUMExactRuleThree(M,nvec,disv);


  delete[] nvec;
  delete[] emat;
  delete[] disv;



  return freq_exact;
}





// Approximate HUM by sampling method.
Rcpp::List HUMApproxRuleOne(int R, int M, int n, int* vec_ni, double* mat_E, double eta, int if_std) {
  double var_Jackk = 0; // 计算se，如果不估计se，则默认为0.
  int if_negative = 0;


  // initialization
  int id_sample = 0;
  long int counter_negative = 0;
  // int break_status = 0;

  // allocate the memory and should be deleted on exit.
  int* vec_cumsum_n =new int[M]();
  vec_cumsum_n[0] = 0;
  for (int j = 1; j < M; ++j) {
    vec_cumsum_n[j] = vec_ni[j - 1] + vec_cumsum_n[j - 1];
  }
  double* vec_comparison = new double[M]();


  // module for calculate estimated standard error.
  int* vec_sampled_id = new int[M](); // ****
  double* counter_Mcate_neg = new double[n]();
  for(int mth = 0; mth < M; mth++){
    counter_Mcate_neg[mth] = 0;
  }
  int* vec_expandn = new int[n]();




  for (int rth = 0; rth < R; rth++) {
    if(rth % 20000 == 0){Rcpp::checkUserInterrupt();}
    if_negative = 0;

    // 第一组抽样
    id_sample = GeneRandomInt(vec_ni[0]);
    vec_sampled_id[0] = id_sample;


    // 循环抽样和比较
    if(eta == 0){
      vec_comparison[0] = mat_E[id_sample * M + 0];
      for (int mth = 1; mth < M; mth++) {
        id_sample = GeneRandomInt(vec_ni[mth]) + vec_cumsum_n[mth];

        // 比较
        if(if_negative == 0){
          for (int j = 0; j < mth; j++) {
            if (mat_E[id_sample * M + j] >= vec_comparison[j]) {
              if_negative = 1;
              // break_status = 1;
              break;
            }
          }
          vec_comparison[mth] = mat_E[id_sample * M + mth] ;
        }

        // 记录下抽取的id
        vec_sampled_id[mth] = id_sample;

      }
    }else{
      vec_comparison[0] = mat_E[id_sample * M + 0] + rnorm_c()*eta*0.33;
      for (int mth = 1; mth < M; mth++) {
        id_sample = GeneRandomInt(vec_ni[mth]) + vec_cumsum_n[mth];


        if(if_negative == 0){
          for (int j = 0; j < mth; j++) {
            if ((mat_E[id_sample * M + j] + rnorm_c()*eta*0.33) >= vec_comparison[j]) {
              if_negative = 1;
              // break_status = 1;
              break;
            }
          }
          vec_comparison[mth] = mat_E[id_sample * M + mth] + rnorm_c()*eta*0.33;
        }

        vec_sampled_id[mth] = id_sample;
      }
    }

    counter_negative += if_negative;

    // 如果需要计算std
    if(if_std == 1){
      if(if_negative == 1){
        for(int mth = 0; mth <M ;mth++ ){
          counter_Mcate_neg[vec_sampled_id[mth]] += 1; // 失败计数加一
        }
      }
    }

  }


  // calculate the approximation hum
  double hhum = 1 - (double)counter_negative / (double)R;


  // calculate the standard error
  // if(if_std == 1){
  //   for(int mth = 0; mth < M; mth++){
  //     for(int ith = 0; ith < vec_ni[mth]; ith++){
  //       vec_expandn[vec_cumsum_n[mth]+ith] = vec_ni[mth];
  //     }
  //   }
  //
  //   // 初始化一个用于记录barU的变量
  //   double barU_noinx_mth;
  //   double var_Jackk_mth;
  //
  //   for(int inx = 0; inx < n; inx++){
  //     counter_Mcate_neg[inx] =
  //       ((double)counter_negative - counter_Mcate_neg[inx])/(double)R * (double)vec_expandn[inx] / (double)(vec_expandn[inx] - 1);
  //   }
  //
  //   for(int mth = 0; mth < M; mth++){
  //     barU_noinx_mth = 0;
  //     for(int inx = 0; inx < vec_ni[mth]; inx++){
  //       barU_noinx_mth += counter_Mcate_neg[inx + vec_cumsum_n[mth]];
  //     }
  //     barU_noinx_mth = barU_noinx_mth / vec_ni[mth];
  //
  //     var_Jackk_mth = 0;
  //     for(int inx = 0; inx < vec_ni[mth]; inx++){
  //       var_Jackk_mth += pow(counter_Mcate_neg[inx + vec_cumsum_n[mth]] - barU_noinx_mth,2);
  //     }
  //     var_Jackk_mth = var_Jackk_mth * (vec_ni[mth] - 1) / vec_ni[mth];
  //
  //     var_Jackk += var_Jackk_mth;
  //   }
  // }



  if(if_std == 1){
    // 初始化一个用于记录barU的变量
    double barU_noinx = 0;

    for(int mth = 0; mth < M; mth++){
      for(int ith = 0; ith < vec_ni[mth]; ith++){
        vec_expandn[vec_cumsum_n[mth]+ith] = vec_ni[mth];
      }
    }

    for(int inx = 0; inx < n; inx++){
      counter_Mcate_neg[inx] =
        ((double)counter_negative - counter_Mcate_neg[inx])/(double)R * (double)vec_expandn[inx] / (double)(vec_expandn[inx] - 1);
      barU_noinx += counter_Mcate_neg[inx];
    }
    barU_noinx = barU_noinx/n;
    for(int inx = 0; inx < n; inx++){
      var_Jackk += pow(counter_Mcate_neg[inx] - barU_noinx ,2);
    }
    var_Jackk = var_Jackk * (n - 1) / n;
  }



  // deleting and freeing part
  delete[] vec_comparison;
  delete[] vec_cumsum_n;

  delete[] vec_sampled_id;
  delete[] counter_Mcate_neg;
  delete[] vec_expandn;

  Rcpp::List toret;
  toret["hhum"] = hhum;
  if(if_std == 1){
    toret["hstd"] = sqrt(var_Jackk);
  }
  return toret;
}





Rcpp::List HUMApproxRuleThree(long int R, int M,int n, int* vec_ni, double* mat_dis_to_v, double eta, int if_std) {
  long int counter_postive = 0; // count the number of times the condition is met (positive).
  double var_Jackk = 0; // 计算se，如果不估计se，则默认为0.
  int if_positive = 0;

  // calculate the cumsum of n
  int* vec_cumsum_n = new int[M]();
  vec_cumsum_n[0] = 0;
  for (int i = 1; i < M; i++) {
    vec_cumsum_n[i] = vec_ni[i - 1] + vec_cumsum_n[i - 1];
  }


  // anxilary variables and memory allocation
  double* u = new double[M + 1]();
  double* v = new double[M + 1]();
  double* minv = new double[M + 1]();
  int* p = new int[M + 1]();
  int* way = new int[M + 1]();
  int* used = new int[M + 1]();
  double* sampled_dis_vec_augmented = new double[(M+1)*(M+1)]();


  // module for calculate estimated standard error.
  int* vec_sampled_id = new int[M](); // ****
  double* counter_Mcate_pos = new double[n]();
  for(int mth = 0; mth < M; mth++){
    counter_Mcate_pos[mth] = 0;
  }
  int* vec_expandn = new int[n]();


  for (int rth = 0; rth < R; rth++) {
    // 允许中途打断
    if(rth % 20000 == 0){Rcpp::checkUserInterrupt();}

    // 抽样并赋值
    for (int mth = 0; mth < M; mth++) {
      int sample_id = 0;
      sample_id = GeneRandomInt(vec_ni[mth]) + vec_cumsum_n[mth];
      // 把抽到的id记录下来
      vec_sampled_id[mth] = sample_id;

      // 如果不需要扰动
      if(eta == 0){
        for (int j = 0; j < M; j++) {
          sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = mat_dis_to_v[sample_id * M + j];
        }
      }else{
        // 否则增加一个随机扰动项
        for (int j = 0; j < M; j++) {
          sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = mat_dis_to_v[sample_id * M + j] + runifcos_c()*eta;
        }
      }

    }

    // 计算是否为真
    if_positive = RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
    counter_postive += if_positive;


    // 如果需要计算std
    if(if_std == 1){
      if(if_positive == 1){
        for(int mth = 0; mth <M ;mth++ ){
          counter_Mcate_pos[vec_sampled_id[mth]] += 1; // 成功计数加一
        }
      }
    }



  }


  // calculate the approximation hum
  double hhum = (double)counter_postive / (double)R;


  // calculate the standard error
  if(if_std == 1){
    // 初始化一个用于记录barU的变量
    double barU_noinx = 0;

    for(int mth = 0; mth < M; mth++){
      for(int ith = 0; ith < vec_ni[mth]; ith++){
        vec_expandn[vec_cumsum_n[mth]+ith] = vec_ni[mth];
      }
    }

    for(int inx = 0; inx < n; inx++){
      counter_Mcate_pos[inx] =
        ((double)counter_postive - counter_Mcate_pos[inx])/(double)R * (double)vec_expandn[inx] / (double)(vec_expandn[inx] - 1);
      barU_noinx += counter_Mcate_pos[inx];
    }
    barU_noinx = barU_noinx/n;
    for(int inx = 0; inx < n; inx++){
      var_Jackk += pow(counter_Mcate_pos[inx] - barU_noinx ,2);
    }
    var_Jackk = var_Jackk * (n - 1) / n;
  }




  // delete part
  delete[] sampled_dis_vec_augmented;
  delete[] u;
  delete[] v;
  delete[] minv;
  delete[] way;
  delete[] used;
  delete[] p;

  delete[] vec_cumsum_n;

  delete[] vec_sampled_id;
  delete[] counter_Mcate_pos;
  delete[] vec_expandn;
  // gsl_rng_free(sampler);


  Rcpp::List toret;
  toret["hhum"] = hhum;
  if(if_std == 1){
    toret["hstd"] = sqrt(var_Jackk);
  }
  return toret;
}





// [[Rcpp::export]]
Rcpp::List HUMApproxRuleAll_R(int R, int M,IntegerVector vec_ni_R, NumericMatrix Ematrix_R, double eta=0, int if_std = 0,int rule_type = 4) {


  // translate R objects to c++ objects
  int n = 0;
  int* vec_ni = new int[M]; // memory allocation
  for(int mth = 0; mth < M; mth++){
    vec_ni[mth] = vec_ni_R[mth];
    n += vec_ni_R[mth];
  }

  // memory allocation (depend on n)
  double* mat_E =new double [n*M];
  for(int inx = 0; inx < n ; inx++){
    for(int jnx = 0; jnx < M ; jnx++){
      mat_E[inx*M+jnx] = Ematrix_R(inx,jnx);
    }
  }
  // calculate the distance matrix
  double* mat_dis_to_v = new double[n*M]();
  dis_cal(M, n, mat_E, mat_dis_to_v);


  Rcpp::List toret;
  // approximate the HUM under rule 3
  if(rule_type == 3){
    toret["Rule3"] = HUMApproxRuleThree(R,M,n,vec_ni,mat_dis_to_v,eta,if_std);
  }
  if(rule_type == 1){
    toret["Rule1"] =   HUMApproxRuleOne(R,M,n,vec_ni,mat_E,       eta,if_std);
  }
  if(rule_type == 4){
    toret["Rule1"] =   HUMApproxRuleOne(R,M,n,vec_ni,mat_E,       eta,if_std);
    toret["Rule3"] = HUMApproxRuleThree(R,M,n,vec_ni,mat_dis_to_v,eta,if_std);
  }


  delete[] vec_ni;
  delete[] mat_E;
  delete[] mat_dis_to_v;



  return toret;
}


