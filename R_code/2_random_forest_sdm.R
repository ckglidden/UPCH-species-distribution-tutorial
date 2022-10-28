####random  forest  species  distribution  models
library(tidyr);library(dplyr);library(spatialsample);library(sf);library(ranger)

#------------------------------------------------------------------------#
#read  in  data  &  merge  by  row  id    (or  merge  from  output  of  above  code)        #
#------------------------------------------------------------------------#

lulc <- read.csv("data/a_chamek_ter_mammals_lulc_cleaned_Oct2022.csv")
climate <- read.csv("data/a_chamek_ter_mammals_climate_Oct2022.csv"); climate <- climate[  c("row_code", "bio10_temp_warmest_qt", "bio13_precip_wettest_month", "cmi_min")]
amazon_basin_pnts <-  read.csv("data/a_chamek_ter_mammals_amazon_thinned_Oct22.csv")

covariates <- left_join(climate, lulc, by = "row_code")
data0 <- left_join(amazon_basin_pnts, covariates, by = "row_code")

#for many ML algorithms  you should make sure your reponse variables are not grouped in the dataframe (e.g.  all 1s and 0s next to each other)
data0 <- data0[sample(1:nrow(data0)), ]

#--------------------------------------------#
#get  fold  id  by  k-means  clustering      #
#--------------------------------------------#

#convert  analysis_data  to  a  spatial  object
data0_sf  <-  st_as_sf(x  =  data0,
                       coords  =  c("lon", "lat"),
                       crs  =  "+proj=longlat +datum=WGS84 +ellps=WGS84")

#identify  groups  of  5  clusters  using  the  spatialsample  package
set.seed(99)  #set  seed  to  get  same  split  each  time
clusters  <-  spatial_block_cv(data0_sf, 
                               method = "random", n = 30, #method for how blocks are oriented in space & number of blocks
                               relevant_only = TRUE,  v = 3)  #k-means  clustering  to  identify  cross-validation  folds  (3  is  too  few  to  be  robust  but  using  here  to  save  time)

#for  loop  to  create  a  dataframe  that  assigns  a  fold  number  to  each  data  point
splits_df  <-  c()
for(i  in  1:3){
    new_df  <-  assessment(clusters$splits[[i]])  #extract  points  in  fold  number  i
    new_df$fold  <-  i
    new_df  <-  new_df[  c("row_code", "fold")]
    splits_df  <-  rbind(splits_df, new_df)  #bind  all  points  x  fold  id  together
}

splits_df  <-  st_drop_geometry(splits_df)  #drop  shapefiles

#final  data  -  merge  cluster  id  to  final  dataset  for  analysis
analysis_data  <-  merge(data0, splits_df, by  =  "row_code")

#sanity  check:  check  how  many  data  points  are  in  each  fold
table(analysis_data$fold)

#write  df  to  save  for  later
write.csv(analysis_data, "data/a_chamek_ter_mammals_finalData_Oct22.csv")


#------------------------------------------------#
#run  spatial  cv  to  evaluate  model  performance        #
#------------------------------------------------#

#first  reduce  data  down  to  covariates  of  interest  (or  you  could  specify  it  in  the  formula  below)
analysis_data_v2  <-  analysis_data[  c("presence", "fold", "bio10_temp_warmest_qt", "bio13_precip_wettest_month", "cmi_min",
                                        "mean_forest",  "mean_farming", "mean_urban",
                                        "diff_forest_formation")]

#------------------------------------#
#tune, train, model                  # fix commas and spaces
#------------------------------------#

#create  empty  dataframe  to  for  loop  to  store  results    one  row  for  each  fold
rf_performance  <-  data.frame(model  =  rep("RF", 3),
                               fold_id  =  1:3,
                               auc  =  rep(NA, 3),
                               sensitivity  =  rep(NA, 3),
                               specificity  =  rep(NA, 3),
                               oob_error  =  rep(NA, 3),
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
        mtry =  seq(1, 6, by  =  2),    #the  number  of  variables  to  randomly  sample  as  candidates  at  each  split
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
        hyper_grid$OOB_RMSE[i]  <-  sqrt(model$prediction.error)
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
    best.threshold  <-  pROC::coords(auc, "best", ret  =  "threshold")
    metrica.format  <-  data.frame(cbind(ifelse(test_complete[,"presence"]==1, 1, 0)),  ifelse(pred  >=  best.threshold[1, 1], 1, 0));  colnames(metrica.format)  <-  c("labels", "predictions");  rownames(metrica.format)  <-  1:dim(metrica.format)[1]
    sensitivity  <-  metrica::recall(data  =  metrica.format, obs  =  labels, pred  =  predictions)$recall  
    rf_performance[i, "sensitivity"]  <-  sensitivity  
    specificity  <-  metrica::specificity(data  =  metrica.format, obs  =  labels, pred  =  predictions)$spec
    rf_performance[i, "specificity"]  <-  specificity
    rf_performance[i, "oob_error"]  <-  train_model$prediction.error
    rf_performance[i, "presence"]  <-  nrow(subset(test, presence  ==  1))
    rf_performance[i, "background"]  <-  nrow(subset(test, presence  ==  0))
    
    hypergrid_final[i, "mtry"]  <-  hyper_grid2$mtry[1]
    hypergrid_final[i, "node_size"]  <-  hyper_grid2$node_size[1]
    hypergrid_final[i, "sampe_size"]  <-  hyper_grid2$sampe_size[1]

}

#------------------------------------------------------------#
#calculate  average  out  of  sample  performance           #
#------------------------------------------------------------#

model_performance  <-  data.frame(metric  =  names(rf_performance)[2:ncol(rf_performance)],  
                                                                mean_metric  =  colMeans(rf_performance[2:ncol(rf_performance)]))

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

#------------------------------------------------------------#
#variable  importance                                                                                  #
#------------------------------------------------------------#

#extract  model  results  to  get  permutation  importance
permutation_importance  <-  data.frame(variable  =  rownames(as.data.frame(final_model$variable.importance)),  
                                                                          importance  =  as.vector(final_model$variable.importance))
    
#plot  importance
ggplot(permutation_importance, aes(x  =  variable, y  =  importance))  +
    geom_bar(stat="identity")  +
    ggtitle("permutation  importance")  +
    coord_flip()  +
    theme_classic()




#------------------------------------------------------------#
#pdps                                                                                                                #
#------------------------------------------------------------#

var_names  <-  names(analysis_data_v2[complete.cases(analysis_data_v2),  -c(1, 2)])

#dataframe  to  make  partial  dependence  plots
pd_df  =  data.frame(matrix(vector(), 0, 3, dimnames=list(c(), c('variable', 'value', 'yhat'))),  
                                      row.names  =  NULL, stringsAsFactors=F)

for  (j  in  1:length(var_names))  {  #loop  through  each  variable
    
    output  <-  as.data.frame(pdp::partial(final_model, pred.var  =  var_names[j], prob  =  TRUE, train  =  analysis_data_v2[complete.cases(analysis_data_v2), -2]))
    
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
    theme_bw()



#------------------------------------------------------------#
#prediction  maps                                                                                          #
#------------------------------------------------------------#



