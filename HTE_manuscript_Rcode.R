
#####################
#####################
## Bellavia A et al. Manuscript under review.
## R code to conduct cluster-based (phenotype-based) HTE assessment
#####################
#####################

Packages <- c("table1","survival","eventglm",
              "ggradar","radiant.data","VarSelLCM")
lapply(Packages, library, character.only = TRUE)

data <- "open your data"

set.seed(123)

#####################
# Step 1. Define your covariates 
#####################

## Here selecting the first 7 columns of data

covariates<-data[,c(1:7)]

#####################
# Step 2. Conduct model-based clustering, as presented in https://cran.r-project.org/web/packages/VarSelLCM/vignettes/VarSelLCM.pdf 
#####################


## Here comparing results with 3, 4, and 5 cluster, and conducting variable selection

res_3 <- VarSelCluster(covariates, 3, nbcores = 2, initModel=40, crit.varsel = "BIC")
res_4 <- VarSelCluster(covariates, 4, nbcores = 2, initModel=40, crit.varsel = "BIC")
res_5 <- VarSelCluster(covariates, 5, nbcores = 2, initModel=40, crit.varsel = "BIC")

## Using BIC to compare models

BIC(res_3)
BIC(res_4)
BIC(res_5)

## Selecting the model with 3 groups, print summary

print(res_3)

#####################
# Step 3. Clusering diagnostic
#####################


## Probability of Misclassification (Figure 1a of poster)

plot(res_3, type="probs-class")

## See other diagnostic tools in varselLCM vignette

#####################
# Step 4. Clustering Interpretation
#####################

## Variables importance (Figure 1b in Poster)

plot(res_3)


## merge clusters with original data

datanew<-cbind(data,
               fitted(res_3))


## Covariates by cluster (Table 1 in poster)

table1( ~var1+var2|`fitted(res_3)`,data=datanew )

## Radar plot (package ggradar)

varlist <- c(" list your variables")


r1 <- apply( (datanew%>%filter(`fitted(res_3)`==1))[,varlist], 2, function(x) { round((table(x)[2]/(table(x)[1]+table(x)[2]))*100,0)  })
r2 <- apply( (datanew%>%filter(`fitted(res_3)`==2))[,varlist], 2, function(x) { round((table(x)[2]/(table(x)[1]+table(x)[2]))*100,0)  })
r3 <- apply( (datanew%>%filter(`fitted(res_3)`==3))[,varlist], 2, function(x) { round((table(x)[2]/(table(x)[1]+table(x)[2]))*100,0)  })

t_cat_3 <- rbind.data.frame(r1,r2,r3)
colnames(t_cat_3) <- varlist

group <- c("your labels")
t_cat_3 <- cbind.data.frame(group,t_cat_3)




lcols <- c("darkorange1", "red", "blue")

ggradar(t_cat_3, grid.min=0, grid.mid = 50, grid.max = 100, legend.text.size = 13, 
        legend.position = "bottom",
        values.radar = c("0%", "50%", "100%"),
        axis.label.size = 5, 
        background.circle.colour = "white",
        axis.line.colour = "gray60",
        gridline.min.colour = "gray60",
        gridline.mid.colour = "gray60",
        gridline.max.colour = "gray60",
        group.colours = lcols,
        plot.extent.x.sf=1.6)+ 
  guides(color = guide_legend(nrow = 2))




#####################
# Step 5. TRT effect by clustering (Figure 3 in poster)
#####################

## Interaction test

model_ov <- coxph(Surv(time,status)~trt,
                      data=datanew,ties="breslow",x=TRUE,y=TRUE)

model_ov_int <- coxph(Surv(time,status)~trt*as.factor(`fitted(res_3)`),
                      data=datanew,ties="breslow",x=TRUE,y=TRUE)


anova(model_ov,model_ov_int)

## Stratified effects (both HR and ARD, here at 3 years)

model_cl1 <- coxph(Surv(time,status)~trt,
                   data=datanew[datanew$`fitted(res_3)`==1,],ties="breslow",x=TRUE,y=TRUE)

model_cl1_3risk <- cumincglm(Surv(time,status)~trt,
                             data=datanew[datanew$`fitted(res_3)`==1,],
                             time=1096, model.censoring = "coxph")

model_cl2 <- coxph(Surv(time,status)~trt,
                   data=datanew[datanew$`fitted(res_3)`==2,],ties="breslow",x=TRUE,y=TRUE)

model_cl2_3risk <- cumincglm(Surv(time,status)~trt,
                             data=datanew[datanew$`fitted(res_3)`==2,],
                             time=1096, model.censoring = "coxph")

model_cl3 <- coxph(Surv(time,status)~trt,
                   data=datanew[datanew$`fitted(res_3)`==3,],ties="breslow",x=TRUE,y=TRUE)

model_cl3_3risk <- cumincglm(Surv(time,status)~trt,
                             data=datanew[datanew$`fitted(res_3)`==3,],
                             time=1096, model.censoring = "coxph")

## For model-based and risk-based stratification, refer to https://grf-labs.github.io/grf/articles/survival.html 

#####################
# Step 6: Validation
#####################

## Training validation split (2/3 vs 1/3)

set.seed(123)
dt = sort(sample(nrow(data), nrow(data)*.6666))
data_tr<-data[dt,]
data_val<-data[-dt,]
covariates_tr<-data_tr[,c(1:7)]
covariates_val<-data_val[,c(1:7)]


res_3_trval <- VarSelCluster(covariates_tr, 3, nbcores = 2, initModel=40, crit.varsel = "BIC")

## predict individual probablities in validation set
g<-as.data.frame(
  cbind(c(1,2,3),predict(res_3_trval,newdata=data_val))
)
colnames(g)<-c("Individual","p1","p2","p3")


data_new_val<-cbind(data_val,g)

## Assing individuals to the group with the highest probability

data_new_val$group <- which.pmax(data_new_val$p1, data_new_val$p2,data_new_val$p3)

table(data_new_val$group)

## Table 1 in validation

table1( ~var1+var2|group,data=data_new_val )

## TRT effects in validation sample


model_cl1 <- coxph(Surv(dtime, status)~trt,
                   data=data_new_val[data_new_val$group==1,],ties="breslow",x=TRUE,y=TRUE)

model_cl2 <- coxph(Surv(dtime, status)~trt,
                   data=data_new_val[data_new_val$group==2,],ties="breslow",x=TRUE,y=TRUE)

model_cl3 <- coxph(Surv(dtime, status)~trt,
                   data=data_new_val[data_new_val$group==3,],ties="breslow",x=TRUE,y=TRUE)


## Repeat more times for several splits. Do not integrate within a single bootstrap function as labels for the same group
## will randomly change between iterations
