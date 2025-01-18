#' @export
EstimateHUM = function(X,Class,id_train,train_ratio = 0.75,random_seed = 1,method = "logistic",para_exact = F,para_rule = 1,para_precision_type = 3,para_print = T,balance = F,...){
  # X = mydata[,1:2]
  # Class = mydata[,3]
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)
  if(n != length(Class)){
    stop("The number of rows of matrix X should be the same as the length of vector of Class!")
  }
  if(missing(id_train)){
    if(para_print){
      message("Missing the index of training data, automatically generate it with random seed set as random_seed!")
    }


    set.seed(random_seed)
    if(balance){
      id_train = list()
      for(cnx in 1:length(table(Class))){
        id.cnx = which(Class == cnx)
        n.cnx = length(id.cnx)
        id.train.cnx = sample(1:n.cnx,floor(n.cnx*train_ratio))

        id_train[[cnx]] = id.cnx[id.train.cnx]
      }
      id_train = unlist(id_train)
    }else{
      id_train = sort(sample(1:n,floor(n*train_ratio)))
    }






    if(para_print){
      tomsg = paste0("The size of training set is ",length(id_train),".")
      message(tomsg)
    }
  }

  fulldata = as.data.frame(cbind(X,Class))
  colnames(fulldata) = c(paste0("X",1:p),"Class")
  fulldata$Class = factor(fulldata$Class)

  if(method == "logistic"){
    if(!para_print){
      caop = capture.output({md_train = nnet::multinom(Class~.,data = fulldata[id_train,])})
    }else{
      md_train = nnet::multinom(Class~.,data = fulldata[id_train,])
    }

    hatEC = cbind(
      predict(md_train,newdata = fulldata,type = "prob"),
      fulldata$Class
    )
  }

  hum_train = CalculateHUM_inputMatrix(hatEC[id_train,],para_exact = para_exact,para_precision_type = para_precision_type,para_type = para_rule,para_random_seed = random_seed,...)

  if(length(id_train) == n){
    message("No test dataset")
    hum_test = c(-1,-1)
  }else{
    hum_test = CalculateHUM_inputMatrix(hatEC[-id_train,],para_exact = para_exact,para_precision_type = para_precision_type,para_type = para_rule,para_random_seed = random_seed,...)
  }
  return(list(
    "Train" = hum_train,
    "Test" = hum_test))
}


#' @export
SeHUM = function(X,Class,train_ratio = 0.75,random_seed = 1,method = "logistic",para_B = 64,para_rule = 1,para_precision_type = 3,para_print = T,...){
  X = as.matrix(X)
  n = nrow(X)
  p = ncol(X)

  fulldata = as.data.frame(cbind(X,Class))
  colnames(fulldata) = c(paste0("X",1:p),"Class")
  set.seed(random_seed)
  random_seed_vec = sample(1:1e6,para_B)

  bootstrap_list = list("Train" = list(),"Test" = list())
  for(bt in 1:para_B){
    if(para_print){message(bt)}

    set.seed(bt)
    id_boot = sample(1:n,n,replace = T)
    hum_boot = EstimateHUM(X[id_boot,],Class[id_boot],
                           train_ratio = train_ratio,
                           random_seed = random_seed_vec[bt],
                           method = method,
                           para_rule = para_rule,para_precision_type = para_precision_type,
                           para_print = F,...)
    bootstrap_list[["Train"]][[bt]] = unlist(hum_boot[["Train"]])
    bootstrap_list[["Test"]][[bt]] = unlist(hum_boot[["Test"]])
  }
  bootstrap_matrix_Train = do.call(rbind,bootstrap_list[["Train"]])
  bootstrap_matrix_Test =  do.call(rbind,bootstrap_list[["Test"]])
  se_Train = sqrt(apply(bootstrap_matrix_Train,2,var))
  se_Test =  sqrt(apply(bootstrap_matrix_Test,2,var))
  return(list(
    "se_Train" = se_Train,
    "se_Test" = se_Test
  ))
}
