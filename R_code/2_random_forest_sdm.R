####random forest species distribution models
library(tidyr); library(dplyr); library(spatialsample); library(sf); library(ranger)

#------------------------------------------------------------------------#
#read in data & merge by row id  (or merge from output of above code)    #
#------------------------------------------------------------------------#

lulc <- read.csv("data/m_noctivagus_ter_mammals_lulc_cleaned_Oct2022.csv")
amazon_basin_pnts <-  read.csv("data/m_noctivagus_ter_mammals_amazon_thinned_Oct22.csv")

data0 <- left_join(amazon_basin_pnts, lulc, by = "row_code")

#------------------------------------#
#get fold id by k-means clustering   #
#------------------------------------#

#convert analysis_data to a spatial object
data0_sf <- st_as_sf(x = data0, 
                             coords = c("lon", "lat"),
                             crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

#identify groups of 5 clusters using the spatialsample package
set.seed(909) #set seed to get same split each time
clusters <- spatial_clustering_cv(data0_sf, v = 3) #k-means clustering to identify 5 cross-validation folds

#for loop to create a dataframe that assigns a fold number to each data point
splits_df <- c()
for(i in 1:3){
  new_df <- assessment(clusters$splits[[i]]) #extract points in fold number i
  new_df$fold <- i
  new_df <- new_df[,c("row_code", "fold")]
  splits_df <- rbind(splits_df, new_df) #bind all points x fold id together
}

splits_df <- st_drop_geometry(splits_df) #drop shapefiles

#final data - merge cluster id to final dataset for analysis
analysis_data <- merge(data0, splits_df, by = "row_code")

#sanity check: check how many data points are in each fold
table(analysis_data$fold)

#write df to save for later
write.csv(analysis_data, "data/m_noctivagus_ter_mammals_finalData_Oct22.csv")


#------------------------------------------------#
#run spatial cv to evaluate model performance    #
#------------------------------------------------#



#------------------------------------#
#tune, train model,                  #
#------------------------------------#

#create empty dataframe to for loop to store results, one row for each fold
rf_performance <- data.frame(model = rep("RF", 3),
                             fold_id = 1:3,
                             auc = rep(NA, 3),
                             sensitivity = rep(NA, 3),
                             specificity = rep(NA, 3),
                             oob_error = rep(NA, 3),
                             presence = rep(NA, 3), #number of presence points in the fold
                             background = rep(NA, 3)) #number of bkg points in the fold

#create empty dataframe to store parameters used to train each model
hypergrid_final <- data.frame(mtry = rep(NA, 3), #the number of variables to randomly sample as candidates at each split
                              node_size  = rep(NA, 3), #minimum number of samples within the terminal nodes
                              sampe_size = rep(NA, 3), #the number of samples to train on
                              num.trees  = rep(NA, 3)) #number of trees
 
  
  
for(i in 1:3){ # run one iteration per fold
  
  train <- analysis_data[analysis_data$fold != i, ]
  test <- analysis_data[analysis_data$fold == i, ]
  
  #remove any rows with NAs bc RF can't handle missing data
  train_complete <- train[complete.cases(train),]
  test_complete <- test[complete.cases(test),]
  
  
  #------------------------------------#
  #define the grid to search over      #
  #------------------------------------#
  # hyperparameter grid search - keeping it small to save time
  # the function below creates a grid with all combinations of parameters
  
  hyper_grid <- expand.grid(
    mtry       = seq(1, 4, by = 1), #the number of variables to randomly sample as candidates at each split
    node_size  = seq(3, 12, by = 3), #minimum number of samples within the terminal nodes
    sampe_size = c(.6, .70, .80), #the number of samples to train on
    num.trees  = c(500, 1000), #number of trees
    OOB_RMSE   = 0
  )
  
  #tune model
  for(j in 1:nrow(hyper_grid)) {
    
    # train model
    model <- ranger(
      formula = presence ~ farming + urban + flooded_forest + forest_formation + river_lake_ocean, 
      data = train_complete, 
      num.trees = hyper_grid$num.trees[j],
      mtry = hyper_grid$mtry[j],
      min.node.size = hyper_grid$node_size[j],
      sample.fraction = hyper_grid$sampe_size[j],
      probability = TRUE,
      seed = 123
    )
    
    # add OOB error to grid
    hyper_grid$OOB_RMSE[i] <- sqrt(model$prediction.error)
  }
  
  #arrange the hypergrid so the lowest out-of-bag error (best performing set of parameters) is in the first row
  hyper_grid2 <- hyper_grid %>% 
    dplyr::arrange(OOB_RMSE)
  
  #train model
  train_model <- ranger(
    formula = presence ~ farming + urban + flooded_forest + forest_formation + river_lake_ocean, 
    data = train_complete, 
    #use the first row of the grid as model parameters
    num.trees = hyper_grid2$num.trees[1],
    mtry = hyper_grid2$mtry[1],
    min.node.size = hyper_grid2$node_size[1],
    sample.fraction = hyper_grid2$sampe_size[1],
    probability = TRUE,
    seed = 123)
  
  #save model performance results
  pred0 <- predict(train_model, data=test_complete); pred <- pred0$predictions[,1]
  auc <- pROC::roc(response=test_complete[,"presence"], predictor=pred, levels=c(0,1), auc = TRUE)
  rf_performance[i, "auc"] <- auc$auc
  best.threshold <- pROC::coords(auc, "best", ret = "threshold")
  metrica.format <- data.frame(cbind(ifelse(test_complete[,"presence"]==1,1,0)),ifelse(pred >= best.threshold[1,1],1,0)); colnames(metrica.format) <- c("labels","predictions"); rownames(metrica.format) <- 1:dim(metrica.format)[1]
  sensitivity <- metrica::recall(data = metrica.format, obs = labels, pred = predictions)$recall 
  rf_performance[i, "sensitivity"] <- sensitivity 
  specificity <- metrica::specificity(data = metrica.format, obs = labels, pred = predictions)$spec
  rf_performance[i, "specificity"] <- specificity
  rf_performance[i, "oob_error"] <- train_model$prediction.error
  rf_performance[i, "presence"] <- nrow(subset(test, presence == 1))
  rf_performance[i, "background"] <- nrow(subset(test, presence == 0))
  
  hypergrid_final[i, "mtry"] <- hyper_grid2$mtry[1]
  hypergrid_final[i, "node_size"] <- hyper_grid2$node_size[1]
  hypergrid_final[i, "sampe_size"] <- hyper_grid2$sampe_size[1]
  hypergrid_final[i, "num.trees"] <- hyper_grid2$num.trees[1]
  
}

#------------------------------------------------------------#
#calculate average out of sample performance                 #
#------------------------------------------------------------#

model_performance <- data.frame(metric = names(rf_performance)[2:ncol(rf_performance)],
                                mean_metric = colMeans(rf_performance[2:ncol(rf_performance)]))

#------------------------------------------------------------#
#train final model                                           #
#------------------------------------------------------------#

final_model <- ranger(
  formula = presence ~ farming + urban + flooded_forest + forest_formation + river_lake_ocean, 
  data = analysis_data[complete.cases(analysis_data),], 
  num.trees = 500,
  mtry = 2,
  min.node.size = 3,
  sample.fraction = 0.6,
  probability = TRUE,
  importance = 'permutation', #specify this to get variable importance in the next step
  seed = 123)

#------------------------------------------------------------#
#variable importance                                         #
#------------------------------------------------------------#

#extract model results to get permutation importance
permutation_importance <- data.frame(variable = rownames(as.data.frame(final_model$variable.importance)),
                                     importance = as.vector(final_model$variable.importance))
  
#plot importance
ggplot(permutation_importance, aes(x = variable, y = importance)) +
  geom_bar(stat="identity") +
  ggtitle("permutation importance") +
  theme_classic()




#------------------------------------------------------------#
#pdps                                                        #
#------------------------------------------------------------#










