# NMF Regression functions
# Mar 29, 2023
# Andy Zeng 

library(tidyverse)
library(ggpubr)
library(glmnet)
library(rsample)


## Train LASSO ##########################################################################
# Train on patient survival using various parameters
#
train_LASSO <- function(x_train, y_train, alpha = 1){
  
  train_y <- y_train$NMFscore
  
  # Perform Lasso regression with LOOCV 
  model <- cv.glmnet(x = x_train, y = train_y, nfold = dim(x_train)[1], family = 'gaussian', alpha = alpha, maxit=1000000, standardize=FALSE)
  #plot(model)
  
  return(model)
}

### Evaluate correlation ##########################################################################
# Evaluate correlation on validation set
#
evaluate_model <- function(model, x_val, anno_val, lambda, 
                           feature_name, iteration, foldname){

  # Create score classification with survival and get covariates
  pred_y <- predict(model, x_val, s = lambda) %>% data.frame()
  colnames(pred_y) <- 'PredScore'
  pred_y <- pred_y %>% rownames_to_column('Patient') %>% 
    # add anno to get covariates
    left_join(anno_val, by = 'Patient')
  
  # Calculate correlation in validation set
  pearson <- cor(pred_y$PredScore, pred_y$NMFscore, method = 'pearson')
  spearman <- cor(pred_y$PredScore, pred_y$NMFscore, method = 'spearman')
  
  # Summary Metrics
  summary_metrics <- data.frame(
    'model_id' = paste0(feature_name, '_iter', iteration, '_', foldname),
    'lambda' = lambda,
    'model_size' = sum(coef(model, s = lambda)!=0),
    'pearson' = pearson,
    'spearman' = spearman,
    'features' = feature_name,
    'iteration' = iteration,
    'foldname' = foldname
  )
  return(summary_metrics)
}


## GridSearch Lasso ##########################################################################
# Perform lasso model gridsearch given a specified feature set 
#
gridsearch_lasso <- function(expr_train, expr_val, anno_train, anno_val, features, feature_name,
                             iteration, foldname, summary_metrics){
  
  # Filter expr matrix for feature set
  x_train <- expr_train[, colnames(expr_train) %in% features]
  x_val <- expr_val[, colnames(expr_val) %in% features]

  # Train LASSO 
  model <- train_LASSO(x_train, anno_train)

  # Get summary metrics for lambda.min and lambda.1se
  for(lambda in c('lambda.min', 'lambda.1se')){
    summary_metrics <- summary_metrics %>% rbind(
      evaluate_model(model = model, x_val = x_val, anno_val = anno_val, lambda = lambda, 
                     feature_name = feature_name, iteration = iteration, foldname = foldname))
  }
  
  return(summary_metrics)
}


## Nested CV regression ##########################################################################
# Run survival analysis with nested cross validation 
### 5-fold outer cross validation for each iteration
### Grid search for each feature set and hyperparameter
#
nestedCV_regression <- function(train_anno, train_expr, iteration, feature_sets, summary_metrics){
  # set up random seed and shuffle data 
  set.seed(iteration)
  train_anno <- train_anno[sample(nrow(train_anno)),]
  train_expr <- train_expr[sample(nrow(train_expr)),]
  
  ## 5-fold outer cross validation
  folds <- rsample::vfold_cv(train_anno, 5)
  for(outer_cv in 1:5){
    # fold ID
    foldname <- folds$id[[outer_cv]]
    # get anno splits
    anno_train <- analysis(folds$splits[[outer_cv]])
    anno_val <- assessment(folds$splits[[outer_cv]])
    # get expr splits
    expr_train <- train_expr[anno_train$Patient,]
    expr_val <- train_expr[anno_val$Patient,]
    
    # Iterate through feature set and run gridsearch to train survival functions
    for(feature_name in names(feature_sets)){
      # get feature list
      features <- feature_sets[[feature_name]]
      # run gridsearch and get results
      summary_metrics <- gridsearch_lasso(expr_train = expr_train, expr_val = expr_val, anno_train = anno_train, anno_val = anno_val, 
                                          features = features, feature_name = feature_name, iteration = iteration, foldname = foldname,
                                          summary_metrics = summary_metrics)
    }
  }
  return(summary_metrics)
}




