#' A function to calculate HUM
#' @export
CalculateHUM_inputMatrix = function(in_ecmat,...){
  in_emat = in_ecmat[,-ncol(in_ecmat)]
  in_class = in_ecmat[,ncol(in_ecmat)]
  t1 = sort(unique(in_class))
  num_cates = length(t1)

  if(sum(abs(t1 - (1:num_cates))) != 0){
    stop("Please check the index of the classes, they must be like 1,...,M .")
  }

  emat_list = list()
  for(i in 1:num_cates){
    emat_list[[i]] = in_emat[in_class == i,,drop = F]
  }
  # return(emat_list)
  return(CalculateHUM(emat_list,...))
}

#' The main function to calculate HUM with specified method
#'
#' We offer four method here.
#' ACC: All combinations checking
#' AP: Adaptive precision with parameter kappa (default as 5)
#' FP: Fixed precision with parameter precision_type (i.e. the required precision as 10^(-precision_type))
#'
#' @examples
#' dt1 = gene_data(50,4,0.7,1)
#' CalculateHUM(dt1,para_type = c(1,3),para_exact = T)
#'
#' @export
CalculateHUM = function(
    in_emat_list,
    para_exact = F,para_type = c(1,3),
    para_precision_type = 3,
    para_kappa2 = 25,
    para_prob_bound = 0.05,
    para_R,para_random_seed = 1,eta = "auto",estse = 0){
  # in_emat_list = elist
  # para_random_seed = 1

  # in_emat_list = dt_ri
  num_cates = length(in_emat_list)
  n_vec = c()
  M_vec = c()
  for(item in in_emat_list){
    n_vec = c(n_vec,nrow(item))
    M_vec = c(M_vec,ncol(item))
  }
  if(sum(M_vec == num_cates) != num_cates){
    stop("Check the input E matrix list, they correspond to distinct categories.")
  }

  E_matrix = do.call(rbind,in_emat_list)


  toret = list()
  if(para_exact){
    total_iter = prod(n_vec)
    # if(total_iter >= 1e10){stop("over 1e10, please ignore be true")}
    if(total_iter >= 1e10){
      cat(total_iter)
      message("over 1e10, please ignore be true")
      if(1 %in% para_type){toret[["rule1_exact"]] = 0}
      if(3 %in% para_type){toret[["rule3_exact"]] = 0}
    }else{
      if(1 %in% para_type){toret[["rule1_exact"]] = HUMExactRuleOne_R(M = num_cates,nvec_R = n_vec,Ematrix_R = E_matrix)}
      if(3 %in% para_type){toret[["rule3_exact"]] = HUMExactRuleThree_R(M = num_cates,nvec_R = n_vec,Ematrix_R = E_matrix)}
    }

  }
  if(!para_exact){
    if(missing(para_R)){
      if(para_precision_type > 0){
        final_R = ceiling(log(2/para_prob_bound,exp(1))/(2*10**((-para_precision_type)*2)))
      }else{
        n_sum = sum(n_vec)
        final_R = ceiling(n_sum * log(n_sum) * para_kappa2)
      }
    }else{
      final_R = para_R
    }
    # to_msg = paste0("The final sampling times is ", final_R)
    # message(to_msg)
    # if(1 %in% para_type){
    #   set.seed(para_random_seed)
    #   toret[["rule1_approx"]] = HUMApproxRuleOne_R(R = final_R,M = num_cates,nvec_R = n_vec,Ematrix_R = E_matrix,eta = eta,estse = estse)
    # }
    # if(3 %in% para_type){
    #   set.seed(para_random_seed)
    #   toret[["rule3_approx"]] = HUMApproxRuleThree_R(R = final_R,M = num_cates,nvec_R = n_vec,Ematrix_R = E_matrix,eta = eta,estse = estse)
    # }
    if(eta == "auto"){
      a1 = sort(unlist(E_matrix))
      b1 = abs(a1[-1] - a1[-length(a1)])
      if(sum(b1 <= 1e-5) > 0){
        eta = max(min(abs(b1[b1 > 1e-5]))/(1e4),1e-8)
      }else{
        eta = 0
      }
      message(paste0("eta:",eta))
    }

    toret = HUMApproxRuleAll_R(
      R = final_R,
      M = num_cates,
      vec_ni_R = n_vec,
      Ematrix_R = E_matrix,
      eta = eta,
      if_std = estse,
      rule_type = sum(para_type)
    )
  }

  return(toret)
}

#' @export
gene_data = function(para_ni,para_M,para_tau,para_seed = 1){
  if(length(para_ni) == 1){
    ni_vec = rep(para_ni, para_M)
  }else if(length(para_ni) != para_M){
    stop("The length of para_ni should be the same with the value of para_M.")
  }else{
    ni_vec = para_ni
  }

  set.seed(para_seed)
  dt_list = list()
  for(mi in 1:para_M){
    ni = ni_vec[mi]
    tep = matrix(runif(ni*para_M),nrow = ni)
    tep[,mi] = tep[,mi] + para_tau
    tep = t(apply(tep,1,function(x){x/sum(x)}))
    dt_list[[mi]] = tep
  }
  return(dt_list)
}

#' @export
gene_data_withtie = function(para_ni,para_M,para_tau,para_seed = 1){
  if(length(para_ni) == 1){
    ni_vec = rep(para_ni, para_M)
  }else if(length(para_ni) != para_M){
    stop("The length of para_ni should be the same with the value of para_M.")
  }else{
    ni_vec = para_ni
  }

  set.seed(para_seed)
  dt_list = list()
  for(mi in 1:para_M){
    ni = ni_vec[mi]
    tep = matrix(runif(ni*para_M),nrow = ni)
    tep[,mi] = tep[,mi] + para_tau
    tep = t(apply(tep,1,function(x){
      xtep = round(x/sum(x),1)
      xtep[which.max(xtep)] = xtep[which.max(xtep)] - sum(xtep) + 1
      xtep
    }))
    dt_list[[mi]] = tep
  }
  return(dt_list)
}
# TEST PART

