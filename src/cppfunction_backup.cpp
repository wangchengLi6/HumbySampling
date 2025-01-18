// #include <Rcpp.h>
// #include <vector>
// #include <math.h>
// using namespace Rcpp;
//
//
// int GeneRandomInt(int max_int){
//   int toret = R::runif(0,max_int);
//   return toret;
// }
//
// double runif_c(){
//   double toret = R::runif(0,1);
//   return toret;
// }
//
// double rnorm_c(){
//   double toret = R::rnorm(0,1);
//   return toret;
// }
//
//
// double runifcos_c(){
//   double toret = cos(R::runif(0,1)*3.14159);
//   return toret;
// }
//
// int RuleIIICheck(int n, float* a1, int infv, float* u,float *v, float* minv, int* p, int* way, int* used) {
//
//   // re-initialize
//   for(int i = 0; i<n+1;i++){
//     u[i] = 0;
//     v[i] = 0;
//     minv[i] = 0;
//     p[i] = 0;
//     way[i] = 0;
//     used[i] = 0;
//   }
//
//
//   // 算法模块
//   for (int i = 1; i <= n; ++i) {
//     p[0] = i;
//     int j0 = 0;
//     for (int cs = 0; cs < n + 1; cs++) {
//       minv[cs] = infv;
//       used[cs] = 0;
//     }
//
//     do {
//       used[j0] = 1;
//       float delta = infv;
//       int i0 = p[j0], j1 = 0;
//       for (int j = 1; j <= n; ++j) {
//         if (!used[j]) {
//           float cur = a1[i0*(n+1)+j] - u[i0] - v[j];
//           if (cur < minv[j] - 0.000001) {
//             minv[j] = cur;
//             way[j] = j0;
//           }
//           if (minv[j] < delta ) {
//             delta = minv[j], j1 = j;
//           }
//         }
//       }
//       for (int j = 0; j <= n; ++j) {
//         if (used[j])
//           u[p[j]] += delta, v[j] -= delta;
//         else
//           minv[j] -= delta;
//         j0 = j1;
//       }
//     } while (p[j0] != 0);
//
//
//     do {
//       int j1 = way[j0];
//       p[j0] = p[j1];
//       j0 = j1;
//     } while (j0);
//   }
//
//
//   int toret;
//   toret = 1;
//   for (int j = 1; j <= n; ++j) {
//     if (p[j] != j) {
//       toret = 0;
//       break;
//     }
//   }
//
//   return toret;
// }
//
// int dis_cal(int M, int n, float* Emat_vec, float* dis_vec) {
//   float* disF0T1 = new float[M]();
//
//   for (int rowid = 0; rowid < n; rowid++) {
//     float sum_square = 0;
//     for (int mth = 0; mth < M; mth++) {
//       float mthele = Emat_vec[rowid * M + mth];
//       sum_square += pow(mthele,2);
//       disF0T1[mth] = pow(1 - mthele, 2) - pow(mthele, 2);
//     }
//     for (int mth = 0; mth < M; mth++) {
//       dis_vec[rowid * M + mth] = sqrt(sum_square + disF0T1[mth]);
//     }
//
//   }
//
//
//   delete[] disF0T1;
//   return 0;
// }
//
// int iter_Mcate(int* current_point, int* end_point,int dig,int max_dig) {
//   if (current_point[dig] < end_point[dig]-1) {
//     current_point[dig] += 1;
//   }
//   else {
//     if (dig != max_dig - 1) {
//       current_point[dig] = 0;
//       iter_Mcate(current_point, end_point, dig + 1, max_dig);
//     }
//   }
//   return 0;
// }
//
//
// Rcpp::List HUMApproxRuleOne(int R, int M, int nsum, int* n_vec, float* emat_vec, double eta, int estse) {
//
//   float var_Jackk = 0; // 计算se，如果不估计se，则默认为0.
//   int nega = 0;
//
//
//   // initialization
//   int sampid = 0;
//   long int negcounter = 0;
//   // int break_status = 0;
//
//   // allocate the memory and should be deleted on exit.
//   int* n_cumsum =new int[M]();
//   n_cumsum[0] = 0;
//   for (int j = 1; j < M; ++j) {
//     n_cumsum[j] = n_vec[j - 1] + n_cumsum[j - 1];
//   }
//   float* comparison_vec = new float[M]();
//
//
//   // module for calculate estimated standard error.
//   int* selectid = new int[M](); // ****
//   float* counter_Mcate_neg = new float[nsum]();
//   for(int mth = 0; mth < M; mth++){
//     counter_Mcate_neg[mth] = 0;
//   }
//   int* nvec_expand = new int[nsum]();
//
//
//
//   // sample part:
//   // for (int rth = 0; rth < R; rth++) {
//   //   // 允许中途打断
//   //   if(rth % 20000 == 0){Rcpp::checkUserInterrupt();}
//   //
//   //   // 抽样并赋值
//   //   for (int mth = 0; mth < M; mth++) {
//   //     int sample_id = 0;
//   //     sample_id = GeneRandomInt(nvec[mth]) + n_cumsum_vec[mth];
//   //     // 把抽到的id记录下来
//   //     selectid[mth] = sample_id;
//   //
//   //     // 如果不需要扰动
//   //     if(eta == 0){
//   //       for (int j = 0; j < M; j++) {
//   //         sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j];
//   //       }
//   //     }else{
//   //       // 否则增加一个随机扰动项
//   //       for (int j = 0; j < M; j++) {
//   //         sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j] + runifcos_c()*eta;
//   //       }
//   //     }
//   //
//   //   }
//   //
//   //   // 计算是否为真
//   //   post = RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
//   //   pos_counter += post;
//   //
//   //
//   //   // 如果需要计算std
//   //   if(estse == 1){
//   //     if(post == 1){
//   //       for(int mth = 0; mth <M ;mth++ ){
//   //         counter_Mcate_pos[selectid[mth]] += 1; // 成功计数加一
//   //       }
//   //     }
//   //   }
//   //
//   //
//   //
//   // }
//
//   for (int rth = 0; rth < R; rth++) {
//     if(rth % 20000 == 0){Rcpp::checkUserInterrupt();}
//     nega = 0;
//
//     // 第一组抽样
//     sampid = GeneRandomInt(n_vec[0]);
//     selectid[0] = sampid;
//
//
//     // 循环抽样和比较
//     if(eta == 0){
//       comparison_vec[0] = emat_vec[sampid * M + 0];
//       for (int mth = 1; mth < M; mth++) {
//         sampid = GeneRandomInt(n_vec[mth]) + n_cumsum[mth];
//
//         // 比较
//         if(nega == 0){
//           for (int j = 0; j < mth; j++) {
//             if (emat_vec[sampid * M + j] >= comparison_vec[j]) {
//               nega = 1;
//               // break_status = 1;
//               break;
//             }
//           }
//           comparison_vec[mth] = emat_vec[sampid * M + mth] ;
//         }
//
//         // 记录下抽取的id
//         selectid[mth] = sampid;
//
//       }
//     }else{
//       comparison_vec[0] = emat_vec[sampid * M + 0] + rnorm_c()*eta*0.33;
//       for (int mth = 1; mth < M; mth++) {
//         sampid = GeneRandomInt(n_vec[mth]) + n_cumsum[mth];
//
//
//         if(nega == 0){
//           for (int j = 0; j < mth; j++) {
//             if ((emat_vec[sampid * M + j] + rnorm_c()*eta*0.33) >= comparison_vec[j]) {
//               nega = 1;
//               // break_status = 1;
//               break;
//             }
//           }
//           comparison_vec[mth] = emat_vec[sampid * M + mth] + rnorm_c()*eta*0.33;
//         }
//
//         selectid[mth] = sampid;
//       }
//     }
//
//     negcounter += nega;
//
//
//
//     // 如果需要计算std
//     if(estse == 1){
//       if(nega == 1){
//         for(int mth = 0; mth <M ;mth++ ){
//           counter_Mcate_neg[selectid[mth]] += 1; // 失败计数加一
//         }
//       }
//     }
//   }
//
//
//   // calculate the approximation hum
//   double hhum = 1 - (float)negcounter / (float)R;
//
//
//   // calculate the standard error
//   if(estse == 1){
//     // 初始化一个用于记录barU的变量
//     float barU_noinx = 0;
//
//     for(int mth = 0; mth < M; mth++){
//       for(int ith = 0; ith < n_vec[mth]; ith++){
//         nvec_expand[n_cumsum[mth]+ith] = n_vec[mth];
//       }
//     }
//
//     for(int inx = 0; inx < nsum; inx++){
//       counter_Mcate_neg[inx] =
//         ((float)negcounter - counter_Mcate_neg[inx])/(float)R * (float)nvec_expand[inx] / (float)(nvec_expand[inx] - 1);
//       barU_noinx += counter_Mcate_neg[inx];
//     }
//     barU_noinx = barU_noinx/nsum;
//     for(int inx = 0; inx < nsum; inx++){
//       var_Jackk += pow(counter_Mcate_neg[inx] - barU_noinx ,2);
//     }
//     var_Jackk = var_Jackk * (nsum - 1) / nsum;
//   }
//
//
//   // if(eta == 0){
//   //   for (int rth = 0; rth < R; rth++) {
//   //     if(rth % 20000 == 0){Rcpp::checkUserInterrupt();}
//   //
//   //
//   //     break_status = 0;
//   //     sampid = GeneRandomInt(n_vec[0]);
//   //     comparison_vec[0] = emat_vec[sampid * M + 0];
//   //     for (int mth = 1; mth < M; mth++) {
//   //       // sampid = gsl_rng_uniform_int(sampler, n_vec[mth]) + n_cumsum[mth];
//   //       sampid = GeneRandomInt(n_vec[mth]) + n_cumsum[mth];
//   //       for (int j = 0; j < mth; j++) {
//   //         if (emat_vec[sampid * M + j] >= comparison_vec[j]) {
//   //           negcounter += 1;
//   //           break_status = 1;
//   //           break;
//   //         }
//   //       }
//   //       if (break_status) {
//   //         break;
//   //       }
//   //       else {
//   //         comparison_vec[mth] = emat_vec[sampid * M + mth] ;
//   //       }
//   //     }
//   //   }
//   // }else{
//   //   for (int rth = 0; rth < R; rth++) {
//   //     if(rth % 20000 == 0){Rcpp::checkUserInterrupt();}
//   //
//   //
//   //     break_status = 0;
//   //     // sampid = gsl_rng_uniform_int(sampler, n_vec[0]);
//   //     sampid = GeneRandomInt(n_vec[0]);
//   //     comparison_vec[0] = emat_vec[sampid * M + 0] + rnorm_c()*eta*0.33; // permutation by rnorm, rescale by eta and 0.33
//   //     for (int mth = 1; mth < M; mth++) {
//   //       // sampid = gsl_rng_uniform_int(sampler, n_vec[mth]) + n_cumsum[mth];
//   //       sampid = GeneRandomInt(n_vec[mth]) + n_cumsum[mth];
//   //       for (int j = 0; j < mth; j++) {
//   //         if ((emat_vec[sampid * M + j] + rnorm_c()*eta*0.33) >= comparison_vec[j]) {
//   //           negcounter += 1;
//   //           break_status = 1;
//   //           break;
//   //         }
//   //       }
//   //       if (break_status) {
//   //         break;
//   //       }
//   //       else {
//   //         comparison_vec[mth] = emat_vec[sampid * M + mth] + rnorm_c()*eta*0.33;
//   //       }
//   //     }
//   //   }
//   // }
//
//
//
//
//
//
//
//
//
//
//
//   // deleting and freeing part
//   delete[] comparison_vec;
//   delete[] n_cumsum;
//   // gsl_rng_free(sampler);
//
//   // double to_return = 1 - (float)negcounter / (float)R;
//   // return to_return;
//
//   delete[] selectid;
//   delete[] counter_Mcate_neg;
//   delete[] nvec_expand;
//
//   Rcpp::List toret;
//   toret["hhum"] = hhum;
//   toret["hstd"] = sqrt(var_Jackk);
//   return toret;
// }
//
//
//
// // [[Rcpp::export]]
// Rcpp::List HUMApproxRuleOne_R(int R,int M, IntegerVector nvec_R, NumericMatrix Ematrix_R, double eta = 0, int estse = 0){
//   int nsum = 0;
//   // memory allocation
//   int* n_vec = new int[M]();
//   for(int i = 0; i < M; i++){
//     n_vec[i] = nvec_R[i];
//     nsum += nvec_R[i];
//   }
//
//
//   float* Emat_vec = new float[M*nsum];
//   // translate R object to c++ object
//   for(int i = 0; i < nsum; i++){
//     for(int j = 0; j < M; j++){
//       Emat_vec[i*M+j] = Ematrix_R(i,j);
//     }
//   }
//
//   // calculate HUM under rule 1
//   Rcpp::List freq_approx = HUMApproxRuleOne(R,M,nsum,n_vec,Emat_vec,eta,estse);
//
//
//
//   // delete part
//   delete[] n_vec;
//   delete[] Emat_vec;
//
//   return freq_approx;
// }
//
//
//
// Rcpp::List HUMApproxRuleThree(long int R, int M,int nsum, int* nvec, float* dis_mat_vec, double eta, int estse) {
//   long int pos_counter = 0; // count the number of times the condition is met (positive).
//   float var_Jackk = 0; // 计算se，如果不估计se，则默认为0.
//   int post = 0;
//
//   // calculate the cumsum of n
//   int* n_cumsum_vec = new int[M]();
//   n_cumsum_vec[0] = 0;
//   for (int i = 1; i < M; i++) {
//     n_cumsum_vec[i] = nvec[i - 1] + n_cumsum_vec[i - 1];
//   }
//
//
//   // anxilary variables and memory allocation
//   float* u = new float[M + 1]();
//   float* v = new float[M + 1]();
//   float* minv = new float[M + 1]();
//   int* p = new int[M + 1]();
//   int* way = new int[M + 1]();
//   int* used = new int[M + 1]();
//   float* sampled_dis_vec_augmented = new float[(M+1)*(M+1)]();
//
//
//   // module for calculate estimated standard error.
//   int* selectid = new int[M](); // ****
//   float* counter_Mcate_pos = new float[nsum]();
//   for(int mth = 0; mth < M; mth++){
//     counter_Mcate_pos[mth] = 0;
//   }
//   int* nvec_expand = new int[nsum]();
//
//
//   for (int rth = 0; rth < R; rth++) {
//     // 允许中途打断
//     if(rth % 20000 == 0){Rcpp::checkUserInterrupt();}
//
//     // 抽样并赋值
//     for (int mth = 0; mth < M; mth++) {
//       int sample_id = 0;
//       sample_id = GeneRandomInt(nvec[mth]) + n_cumsum_vec[mth];
//       // 把抽到的id记录下来
//       selectid[mth] = sample_id;
//
//       // 如果不需要扰动
//       if(eta == 0){
//         for (int j = 0; j < M; j++) {
//           sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j];
//         }
//       }else{
//         // 否则增加一个随机扰动项
//         for (int j = 0; j < M; j++) {
//           sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j] + runifcos_c()*eta;
//         }
//       }
//
//     }
//
//     // 计算是否为真
//     post = RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
//     pos_counter += post;
//
//
//     // 如果需要计算std
//     if(estse == 1){
//       if(post == 1){
//         for(int mth = 0; mth <M ;mth++ ){
//           counter_Mcate_pos[selectid[mth]] += 1; // 成功计数加一
//         }
//       }
//     }
//
//
//
//   }
//
//
//   // calculate the approximation hum
//   double hhum = (float)pos_counter / (float)R;
//
//
//   // calculate the standard error
//   if(estse == 1){
//     // 初始化一个用于记录barU的变量
//     float barU_noinx = 0;
//
//     for(int mth = 0; mth < M; mth++){
//       for(int ith = 0; ith < nvec[mth]; ith++){
//         nvec_expand[n_cumsum_vec[mth]+ith] = nvec[mth];
//       }
//     }
//
//     for(int inx = 0; inx < nsum; inx++){
//       counter_Mcate_pos[inx] =
//         ((float)pos_counter - counter_Mcate_pos[inx])/(float)R * (float)nvec_expand[inx] / (float)(nvec_expand[inx] - 1);
//       barU_noinx += counter_Mcate_pos[inx];
//     }
//     barU_noinx = barU_noinx/nsum;
//     for(int inx = 0; inx < nsum; inx++){
//       var_Jackk += pow(counter_Mcate_pos[inx] - barU_noinx ,2);
//     }
//     var_Jackk = var_Jackk * (nsum - 1) / nsum;
//   }
//
//
//
//   // if(estse == 0){
//   //   if(eta == 0){
//   //     for (int rth = 0; rth < R; rth++) {
//   //       if(rth % 20000 == 0){
//   //         Rcpp::checkUserInterrupt();
//   //       }
//   //       for (int mth = 0; mth < M; mth++) {
//   //         int sample_id = 0;
//   //         sample_id = GeneRandomInt(nvec[mth]) + n_cumsum_vec[mth];
//   //         // sample_id = gsl_rng_uniform_int(sampler, nvec[mth]) + n_cumsum_vec[mth];
//   //         for (int j = 0; j < M; j++) {
//   //           sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j];
//   //         }
//   //       }
//   //       pos_counter += RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
//   //     }
//   //   }else{
//   //     for (int rth = 0; rth < R; rth++) {
//   //       if(rth % 20000 == 0){
//   //         Rcpp::checkUserInterrupt();
//   //       }
//   //       for (int mth = 0; mth < M; mth++) {
//   //         int sample_id = 0;
//   //         sample_id = GeneRandomInt(nvec[mth]) + n_cumsum_vec[mth];
//   //         // sample_id = gsl_rng_uniform_int(sampler, nvec[mth]) + n_cumsum_vec[mth];
//   //         for (int j = 0; j < M; j++) {
//   //           // 增加扰动
//   //           sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j] + runifcos_c()*eta;
//   //         }
//   //       }
//   //       pos_counter += RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
//   //     }
//   //   }
//   // }else{
//   //   // int* selectid = new int[M](); // ****
//   //   // float* counter_Mcate_pos = new float[nsum]();
//   //   // for(int mth = 0; mth < M; mth++){
//   //   //   counter_Mcate_pos[mth] = 0;
//   //   // }
//   //
//   //   // Rcpp::Rcout << "over1" << std::endl;
//   //   if(eta == 0){
//   //     for (int rth = 0; rth < R; rth++) {
//   //       if(rth % 20000 == 0){
//   //         Rcpp::checkUserInterrupt();
//   //       }
//   //       for (int mth = 0; mth < M; mth++) {
//   //         int sample_id = 0;
//   //         sample_id = GeneRandomInt(nvec[mth]) + n_cumsum_vec[mth];
//   //         // sample_id = gsl_rng_uniform_int(sampler, nvec[mth]) + n_cumsum_vec[mth];
//   //
//   //         for (int j = 0; j < M; j++) {
//   //           sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j];
//   //         }
//   //
//   //         // 把抽到的id记录下来
//   //         selectid[mth] = sample_id; // ****
//   //       }
//   //       int indx = RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
//   //       if(indx == 1){
//   //         pos_counter+=1;
//   //         for(int mth = 0; mth <M ;mth++ ){
//   //           int index_mth = selectid[mth];
//   //           counter_Mcate_pos[index_mth] += 1; // 成功计数加一
//   //         }
//   //       }
//   //     }
//   //   }else{
//   //     for (int rth = 0; rth < R; rth++) {
//   //       if(rth % 20000 == 0){
//   //         Rcpp::checkUserInterrupt();
//   //       }
//   //       for (int mth = 0; mth < M; mth++) {
//   //         int sample_id = 0;
//   //         sample_id = GeneRandomInt(nvec[mth]) + n_cumsum_vec[mth];
//   //         // sample_id = gsl_rng_uniform_int(sampler, nvec[mth]) + n_cumsum_vec[mth];
//   //
//   //         for (int j = 0; j < M; j++) {
//   //           // 增加扰动
//   //           sampled_dis_vec_augmented[(mth+1) * (M+1) + j+1] = dis_mat_vec[sample_id * M + j] + runifcos_c()*eta;
//   //         }
//   //
//   //         // 把抽到的id记录下来
//   //         selectid[mth] = sample_id; // ****
//   //       }
//   //       // pos_counter += RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
//   //       // ****
//   //       int indx = RuleIIICheck(M,sampled_dis_vec_augmented,M+1,u,v,minv,p,way,used);
//   //       if(indx==1){
//   //         pos_counter+=1;
//   //         for(int mth = 0; mth <M ;mth++ ){
//   //           int index_mth = selectid[mth];
//   //           counter_Mcate_pos[index_mth] += 1; // 成功计数加一
//   //         }
//   //       }
//   //     }
//   //   }
//   //
//   //
//   //   // Rcpp::Rcout << "over2" << std::endl;
//   //   float barU_noinx = 0;
//   //   int icate;
//   //
//   //   for(int inx = 0; inx < nsum; inx++){
//   //     icate = 0;
//   //     while(inx >= n_cumsum_vec[icate]){
//   //       icate+=1;
//   //       if(icate == M){
//   //         break;
//   //       }
//   //     }
//   //     icate -= 1;
//   //     counter_Mcate_pos[inx] =
//   //       ((float)pos_counter - counter_Mcate_pos[inx])/(float)R * (float)nvec[icate] / (float)(nvec[icate] - 1);
//   //     barU_noinx += counter_Mcate_pos[inx];
//   //     // vecret[inx] =  counter_Mcate_pos[inx];
//   //     // Rcpp::Rcout << counter_Mcate_pos[inx] << std::endl;
//   //   }
//   //   barU_noinx = barU_noinx/nsum;
//   //   for(int inx = 0; inx < nsum; inx++){
//   //     var_Jackk += pow(counter_Mcate_pos[inx] - barU_noinx ,2);
//   //   }
//   //   var_Jackk = var_Jackk * (nsum - 1) / nsum;
//   //
//   //
//   //   delete[] selectid;
//   //   delete[] counter_Mcate_pos;
//   // }
//
// //
// //
// //   // calculate the approximation hum
// //   double hhum = (float)pos_counter / (float)R;
// //
// //
// //   // calculate the standard error
//
//
//
//
//   // delete part
//   delete[] sampled_dis_vec_augmented;
//   delete[] u;
//   delete[] v;
//   delete[] minv;
//   delete[] way;
//   delete[] used;
//   delete[] p;
//
//   delete[] n_cumsum_vec;
//
//   delete[] selectid;
//   delete[] counter_Mcate_pos;
//   delete[] nvec_expand;
//   // gsl_rng_free(sampler);
//
//
//   Rcpp::List toret;
//   toret["hhum"] = hhum;
//   toret["hstd"] = sqrt(var_Jackk);
//   return toret;
// }
//
// // [[Rcpp::export]]
// Rcpp::List HUMApproxRuleThree_R(int R, int M,IntegerVector nvec_R, NumericMatrix Ematrix_R, double eta=0, int estse = 0) {
//   // memory allocation
//   int* nvec = new int[M];
//
//   // translate R objects to c++ objects
//   int nsum = 0;
//   for(int m = 0; m < M; m++){
//     nvec[m] = nvec_R[m];
//     nsum += nvec_R[m];
//   }
//
//   // memory allocation 2 (depend on nsum)
//   float* emat =new float [nsum*M];
//   float* disv = new float[nsum*M]();
//   for(int i = 0; i < nsum ; i++){
//     for(int j = 0; j < M ; j++){
//       emat[i*M+j] = Ematrix_R(i,j);
//     }
//   }
//
//   // calculate the distance matrix
//   dis_cal(M, nsum, emat, disv);
//   // approximate the HUM under rule 3
//   Rcpp::List toret = HUMApproxRuleThree(R,M,nsum,nvec,disv,eta,estse);
//
//
//   delete[] nvec;
//   delete[] emat;
//   delete[] disv;
//
//
//
//   return toret;
// }
//
//
//
// double HUMExactRuleOne(int M, int* n_vec, float* emat_vec) {
//   // initialization
//   int sampid = 0;
//   long int negcounter = 0;
//   int break_status = 0;
//   int max_iter = n_vec[0];
//   // allocate the memory and should be deleted on exit.
//   int* n_cumsum = new int[M]();
//   for (int j = 1; j < M; ++j) {
//     n_cumsum[j] = n_vec[j - 1] + n_cumsum[j - 1];
//     max_iter *= n_vec[j];
//   }
//   float* comparison_vec = new float[M]();
//   int* start_point = new int[M]();
//
//   // sample part:
//   for (int rth = 0; rth < max_iter; rth++) {
//     if(rth % 20000 == 0){
//       Rcpp::checkUserInterrupt();
//     }
//     break_status = 0;
//     sampid = start_point[0];
//     comparison_vec[0] = emat_vec[sampid * M + 0];
//     for (int mth = 1; mth < M; mth++) {
//       sampid = start_point[mth] + n_cumsum[mth];
//       for (int j = 0; j < mth; j++) {
//         if (emat_vec[sampid * M + j] >= comparison_vec[j]) {
//           negcounter += 1;
//           break_status = 1;
//           break;
//         }
//       }
//       if (break_status) {
//         break;
//       }
//       else {
//         comparison_vec[mth] = emat_vec[sampid * M + mth];
//       }
//     }
//     iter_Mcate(start_point, n_vec, 0, M);
//   }
//
//
//   // deleting and freeing part
//   delete[] comparison_vec;
//   delete[] n_cumsum;
//   delete[] start_point;
//
//   double to_return = 1 - (double)negcounter / (double)max_iter;
//   return to_return;
// }
//
// // [[Rcpp::export]]
// double HUMExactRuleOne_R(int M,  IntegerVector nvec_R, NumericMatrix Ematrix_R){
//   int nsum = 0;
//   // memory allocation
//   int* n_vec = new int[M]();
//   for(int i = 0; i < M; i++){
//     n_vec[i] = nvec_R[i];
//     nsum += nvec_R[i];
//   }
//
//
//   float* Emat_vec = new float[M*nsum];
//   // translate R object to c++ object
//   for(int i = 0; i < nsum; i++){
//     for(int j = 0; j < M; j++){
//       Emat_vec[i*M+j] = Ematrix_R(i,j);
//     }
//   }
//
//   // calculate HUM under rule 1
//   double freq_exact = HUMExactRuleOne(M,n_vec,Emat_vec);
//
//
//
//   // delete part
//   delete[] n_vec;
//   delete[] Emat_vec;
//
//   return freq_exact;
// }
//
//
// double HUMExactRuleThree(int M, int* nvec, float* dis_mat_vec) {
//   long int pos_counter = 0; // count the number of times the condition is met (positive).
//   long int max_iter = nvec[0];
//   int sample_id = 0;
//   int* n_cumsum_vec = new int[M]();
//   // calculate the cumsum of n
//   for (int i = 1; i < M; i++) {
//     n_cumsum_vec[i] = nvec[i - 1] + n_cumsum_vec[i - 1];
//     max_iter *= nvec[i];
//   }
//   int* start_point = new int[M]();
//
//   // anxilary variables and memory allocation
//   float* u = new float[M + 1]();
//   float* v = new float[M + 1]();
//   float* minv = new float[M + 1]();
//   int* p = new int[M + 1]();
//   int* way = new int[M + 1]();
//   int* used = new int[M + 1]();
//   float* sampled_dis_vec_augmented = new float[(M + 1) * (M + 1)]();
//
//
//   for (int rth = 0; rth < max_iter; rth++) {
//     if(rth % 20000 == 0){
//       Rcpp::checkUserInterrupt();
//     }
//     for (int mth = 0; mth < M; mth++) {
//       sample_id = start_point[mth] + n_cumsum_vec[mth];
//       for (int j = 0; j < M; j++) {
//         sampled_dis_vec_augmented[(mth + 1) * (M + 1) + j + 1] = dis_mat_vec[sample_id * M + j];
//       }
//     }
//
//     pos_counter += RuleIIICheck(M, sampled_dis_vec_augmented, M + 1, u, v, minv, p, way, used);
//     iter_Mcate(start_point, nvec, 0, M);
//   }
//
//
//
//   // delete part
//   delete[] sampled_dis_vec_augmented;
//   delete[] u;
//   delete[] v;
//   delete[] minv;
//   delete[] way;
//   delete[] used;
//   delete[] p;
//
//   delete[] start_point;
//   delete[] n_cumsum_vec;
//
//   double toret = (double)pos_counter / (double)max_iter;
//   return toret;
// }
//
// // [[Rcpp::export]]
// double HUMExactRuleThree_R(int M,IntegerVector nvec_R, NumericMatrix Ematrix_R) {
//   // memory allocation
//   int* nvec = new int[M];
//
//   // translate R objects to c++ objects
//   int nsum = 0;
//   for(int m = 0; m < M; m++){
//     nvec[m] = nvec_R[m];
//     nsum += nvec_R[m];
//   }
//
//   // memory allocation 2 (depend on nsum)
//   float* emat =new float [nsum*M];
//   float* disv = new float[nsum*M]();
//   for(int i = 0; i < nsum ; i++){
//     for(int j = 0; j < M ; j++){
//       emat[i*M+j] = Ematrix_R(i,j);
//     }
//   }
//
//   // calculate the distance matrix
//   dis_cal(M, nsum, emat, disv);
//   // approximate the HUM under rule 3
//   double freq_exact = HUMExactRuleThree(M,nvec,disv);
//
//
//   delete[] nvec;
//   delete[] emat;
//   delete[] disv;
//
//
//
//   return freq_exact;
// }
