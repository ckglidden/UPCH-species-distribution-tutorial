####random  forest  species  distribution  models

#load packages
library(tidyr);library(dplyr);library(spatialsample);library(sf);library(ranger);library(ggplot2)

analysis_data <- read.csv("data/a_chamek_ter_mammals_finalData_Oct22.csv")

#-------------------------------------------------------------------------#
#CHOICE OF STATISTICAL MODEL -- RANDOM FOREST                             #
#run model using spatial  cv  to  evaluate  model  performance            #
#-------------------------------------------------------------------------#

#first  reduce  data  down  to  covariates  of  interest  (or  you  could  specify  it  in  the  formula  below)

analysis_data_v2  <-  analysis_data[  c("presence", "fold", "bio13_precip_wettest_month", "cmi_min",
                                        "mean_forest",  "mean_farming", "diff_forest_formation")]

#for many ML algorithms  you should make sure your response variables are not grouped in the data-frame (e.g.  all 1s and 0s next to each other)
analysis_data_v2 <- analysis_data_v2[sample(1:nrow(analysis_data_v2)), ]

#------------------------------------#
#tune, train, model                  #
#------------------------------------#

#create  empty  dataframe  to  for  loop  to  store  results    one  row  for  each  fold
rf_performance  <-  data.frame(model  =  rep("RF", 3),
                               fold_id  =  1:3,
                               auc  =  rep(NA, 3),
                               presence  =  rep(NA, 3),    #number  of  presence  points  in  the  fold
                               background  =  rep(NA, 3))  #number  of  bkg  points  in  the  fold

#create  empty  dataframe  to  store  parameters  used  to  train  each  model
hypergrid_final  <-  data.frame(mtry  =  rep(NA,3),    #the  number  of  variables  to  randomly  sample  as  candidates  at  each  split
                                                            node_size    =  rep(NA, 3),    #minimum  number  of  samples  within  the  terminal  nodes
                                                            sampe_size  =  rep(NA, 3))    #the  number  of  samples  to  train  on

    
    
for(i  in  1:3){  #  run  one  iteration  per  fold
    
    train  <-  analysis_data_v2[analysis_data_v2$fold !=  i, ];  train  <-  train[-2]
    test  <-  analysis_data_v2[analysis_data_v2$fold ==  i, ];  test  <-  test[-2]
    
    #remove  any  rows  with  NAs  bc  RF  can't  handle  missing  data
    train_complete  <-  train[complete.cases(train), ]
    test_complete  <-  test[complete.cases(test), ]
    
    
    #-----------------------------------------------#
    #define  the  grid  to  search  over            #
    #-----------------------------------------------#
    #  the  function  below  creates  a  grid  with  all  combinations  of  parameters
    
    hyper_grid  <-  expand.grid(
        mtry =  seq(1, 3, by  =  1),    #the  number  of  variables  to  randomly  sample  as  candidates  at  each  split
        node_size =  seq(1,4, by  =  1),    #shallow trees
        sampe_size  =  c(.6, .70, .80),    #the  number  of  samples  to  train  on
        OOB_RMSE =  0
    )
    
    #tune  model
    for(j  in  1:nrow(hyper_grid)){
        
        #  train  model
        model  <-  ranger(
            formula  =  presence  ~  .,    
            data  =  train_complete,    
            num.trees  =  2000,  
            mtry  =  hyper_grid$mtry[j],
            min.node.size  =  hyper_grid$node_size[j],  
            sample.fraction  =  hyper_grid$sampe_size[j], 
            probability  =  TRUE, 
            replace = TRUE,
            splitrule = "hellinger",
            seed  =  123
            )
        
        #  add  OOB  error  to  grid
        hyper_grid$OOB_RMSE[j]  <-  sqrt(model$prediction.error)
    }
    
    #arrange  the  hypergrid  so  the  lowest  out-of-bag  error  (best  performing  set  of  parameters)  is  in  the  first  row
    hyper_grid2  <-  hyper_grid  %>%  
        dplyr::arrange(OOB_RMSE)
    
    #train  model
    train_model  <-  ranger(
        formula  =  presence  ~  .,    
        data  =  train_complete,   
        #use  the  first  row  of  the  grid  as  model  parameters
        num.trees  =  2000,  
        mtry  =  hyper_grid2$mtry[1], 
        min.node.size  =  hyper_grid2$node_size[1], 
        sample.fraction  =  hyper_grid2$sampe_size[1],
        probability  =  TRUE, 
        replace = TRUE,
        splitrule = "hellinger",
        seed  =  123)
    
    #save  model  performance  results
    pred0  <-  predict(train_model, data=test_complete);  pred  <-  pred0$predictions[,1]
    auc  <-  pROC::roc(response=test_complete[,"presence"], predictor=pred, levels=c(0, 1), auc  =  TRUE)
    rf_performance[i, "auc"]  <-  auc$auc
    rf_performance[i, "presence"]  <-  nrow(subset(test, presence  ==  1))
    rf_performance[i, "background"]  <-  nrow(subset(test, presence  ==  0))
    
    #save hypergrid results to use for final model
    hypergrid_final[i, "mtry"]  <-  hyper_grid2$mtry[1]
    hypergrid_final[i, "node_size"]  <-  hyper_grid2$node_size[1]
    hypergrid_final[i, "sampe_size"]  <-  hyper_grid2$sampe_size[1]

}

#mean auc = 0.71
mean(rf_performance$auc)

#------------------------------------------------------------#
#train  final  model                                                                                      #
#------------------------------------------------------------#

final_model  <-  ranger(
    formula  =  presence  ~  .,    
    data  =  analysis_data_v2[complete.cases(analysis_data_v2), -2],    #complete case dataset without fold column
    #parameters used here are the averages from hypergrid_final
    num.trees  =  2000 ,
    mtry  =  2,  
    min.node.size  =  1,  
    sample.fraction  =  0.6,
    probability  =  TRUE,  
    replace = TRUE,
    splitrule = "hellinger",
    importance  =  'permutation',    #specify  this  to  get  variable  importance  in  the  next  step
    seed  =  123)

#check in-sample auc - if the in-sample AUC is very high but the out-of-sample AUC is low, then your model is likely overfit
pred0  <-  predict(final_model, data=analysis_data_v2[complete.cases(analysis_data_v2), -2]);  pred  <-  pred0$predictions[,1]
auc  <-  pROC::roc(response=analysis_data_v2[complete.cases(analysis_data_v2), -2][,"presence"], predictor=pred, levels=c(0, 1), auc  =  TRUE)
auc$auc

#------------------------------------------------------------#
#MODEL INTERPRETATION AND EVALUATION                         #
#variable  importance                                        #
#------------------------------------------------------------#

#extract  model  results  to  get  permutation  importance
permutation_importance  <-  data.frame(variable  =  rownames(as.data.frame(final_model$variable.importance)),  
                                                                          importance  =  as.vector(final_model$variable.importance))
    
#plot  importance
ggplot(permutation_importance, aes(x  =  variable, y  =  importance))  +
    geom_bar(stat="identity")  +
    ggtitle("permutation  importance")  +
    coord_flip()  +
    theme_classic(base_size = 14)



#------------------------------------------------------------#
#pdps                                                                                                                #
#------------------------------------------------------------#

#try plotting a PDP for just one variable
pdp::partial(final_model, pred.var  =  "mean_forest", prob  =  TRUE, plot = TRUE, train  =  analysis_data_v2[complete.cases(analysis_data_v2), -2])
#train = data without NAs & without "fold" column

###now run a for loop to get plotting data for all variables in the model (or in the 'var_names' list
#list of covariate names to generate pdps for and loop through
var_names  <-  names(analysis_data_v2[complete.cases(analysis_data_v2),  -c(1, 2)]) #analysis dataset exlucing 'presence' and 'fold' column

#dataframe  to  make  partial  dependence  plots
pd_df  =  data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value', 'yhat'))),  
                     row.names  =  NULL, stringsAsFactors=F)

#loop  through  each  variable
for  (j  in  1:length(var_names))  { 
    
    output  <-  as.data.frame(pdp::partial(final_model, pred.var  =  var_names[j], prob  =  TRUE, train  = analysis_data_v2[complete.cases(analysis_data_v2), -2]))
    
    loop_df  <-  data.frame(variable = rep(var_names[j], length(output[[1]])),
                            value  =  output[[1]],
                            yhat  =  output[[2]])
    
    pd_df  <-  rbind(pd_df, loop_df)  
}


#plot  pdps  for  each  variable
ggplot(pd_df, aes(x  =  value, y=  yhat))  +
    geom_smooth()  +
    ylab("probability")  +
    facet_wrap(~variable, scales  =  "free")  +
    theme_bw(base_size = 14)



#------------------------------------------------------------#
#prediction  maps                                            #                                              
#------------------------------------------------------------#

#path for rasters of each covariate in Madre de Dios (we will just plot MDD for now to reduce computational time and file size)
env_data <- list.files(path="env_data", pattern="tif", all.files=FALSE, full.names=TRUE,recursive=TRUE)
e <- raster::stack(env_data)

prediction_df <- as.data.frame(rasterToPoints(e)) ##this gets raster value for every grid cell of MDD

#reduce dataset to complete cases
prediction_df_complete <- prediction_df[complete.cases(prediction_df), ]

#predict probability of species occurrence each 1km grid cell of the area of interest
predictions <- predict(final_model,
                       data=prediction_df_complete)

prediction_df_full <- cbind(prediction_df_complete, as.data.frame(predictions$predictions)[,2])
names(prediction_df_full)[ncol(prediction_df_full)] <- "probability"

#reduce dataset to only the long(x), lat (y), and variable of interest (probability)
rf_tiff_df <- prediction_df_full[, c("x", "y", "probability")]

rf_sdm_raster <- rasterFromXYZ(rf_tiff_df)
outfile <- writeRaster(rf_sdm_raster, filename='final_figures/rf_sdm_example_predictions.tif', format="GTiff",options=c("INTERLEAVE=BAND","COMPRESS=LZW"), overwrite=TRUE)

