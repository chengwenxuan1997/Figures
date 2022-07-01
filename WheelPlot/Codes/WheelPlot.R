workdir <- "I:/genomicdata/External/Xlu/WheelPlot"; setwd(workdir)
fig.path <- file.path(workdir,"Figures")
res.path <- file.path(workdir,"Results")
data.path <- file.path(workdir,"InputData")

invisible(lapply(ls()[grep("path", ls())], function(x){
  if (!dir.exists(get(x))) dir.create(get(x))
}))


library(e1071)
library(caret)
library(data.table)
library(Matrix)
library(matrixStats)
library(glmnet)


# function needed ---------------------------------------------------------

## read cef file
cef2df <- function(filename){
  cat("reading cef file \n")
  full.cef <- as.data.frame(fread(file = filename, skip = 1))
  
  cat("prepare expression matrix\n")
  expr.cef <- full.cef[-c(2, 3, 4), ]
  rownames(expr.cef) <- expr.cef[, 1]
  colnames(expr.cef) <- expr.cef[1, ]
  expr.cef <- expr.cef[-1, -c(1, 2)]
  expr.cef <- as.data.frame(lapply(expr.cef, as.numeric),
                            row.names = rownames(expr.cef), 
                            check.names = F)
  expr.cef <- as.matrix(expr.cef)
  
  cat("prepare sample information\n")
  info.cef <- as.data.frame(t(full.cef[1:3, ]))
  colnames(info.cef) <- info.cef[2, ]
  info.cef <- info.cef[-c(1, 2), ]
  
  ## output a list which consist of,
  ## expression matrix, (a interger matrix)
  ## sample info, (a data frame, CeLL_ID, Cell_type, Timepoint)
  return(list(expr.mat = expr.cef, sample.info = info.cef))
}

## fit the CV
fit_CV <- function(mu, cv, fit_method="SVR", svr_gamma=0.005, x0=c(0.5, 0.5), verbose=FALSE){

  ## in the example jupyter-notebook, only SVR method was applied
  ## so, only the SVR method is included
  log2_m = log2(mu)
  log2_cv = log2(cv)
  
  if (fit_method == "SVR"){
    model = svm(x = log2_m, 
                y = log2_cv, 
                gamma = svr_gamma,
                scale = F)
    fitted_fun = predict(object = model, newdata = log2_m)
    score = log2_cv - fitted_fun
    mu_linspace = seq(min(log2_m), max(log2_m), by = (max(log2_m) - min(log2_m))/49)
    cv_fit = predict(object = model, newdata = mu_linspace)
    return(list("score" = score, 
                "mu_linspace" = mu_linspace, 
                "cv_fit" = cv_fit, 
                "params" = "None"))
  }
}

py_seq <- function(from, to, n){
  seq(from = from, to = to, by = (to-from)/(n-1))
}


# read data, merge matrix --------------------------------------


he.cef <- cef2df(file.path(data.path, "Human_Embryo_fulldataset.cef"))

## remove genes related to cellcycle, blood and SNAR
panther_cellcycle <- read.table(file.path(data.path, "PANTHER_cell_cycle_genes.txt"))$V1
blood_genes = c('HBG1','HBA1','HBA2','HBE1','HBZ','BLVRB','S100A6')
SNAR_genes = c('SNAR-E', 'SNAR-A13_loc1', 'SNAR-C1_loc1', 'SNAR-A1_loc2',
               'SNAR-A8_loc1', 'SNAR-C1_loc2', 'SNAR-A2_loc2', 'SNAR-C4', 
               'SNAR-A12_loc1', 'SNAR-C3', 'SNAR-C1_loc3', 'SNAR-G2', 'SNAR-G1',
               'SNAR-A11_loc9', 'SNAR-A6_loc3', 'SNAR-A14_loc7', 'SNAR-A6_loc5', 
               'SNAR-A10_loc6', 'SNAR-A5_loc9', 'SNAR-A14_loc3', 'SNAR-A9_loc9', 
               'SNAR-A11_loc7', 'SNAR-B1_loc1', 'SNAR-B1_loc2', 'SNAR-D', 'SNAR-F')

gene.to.remove <- unique(c(panther_cellcycle, blood_genes, SNAR_genes))
index.to.remove <- match(gene.to.remove, rownames(he.cef$expr.mat))
index.to.remove <- index.to.remove[!is.na(index.to.remove)]
df <- he.cef$expr.mat[-index.to.remove,]

## load iPSC and ES dataset
ips.cef <- cef2df(file.path(data.path, "iPSC_fulldataset.cef"))
es.cef <- cef2df(file.path(data.path, "ES_fulldataset.cef"))

gene.to.keep <- intersect(rownames(df), rownames(es.cef$expr.mat))

## Add day0 ES to the Reference Dataset
df.merge <- cbind(df[gene.to.keep, ],
                  es.cef$expr.mat[gene.to.keep, es.cef$sample.info$Cell_type %in% c("eSCa", "eSCb", "eSCc")])



# Filter matrix -----------------------------------------------------------

## Filter for Development genes
thrs = 4500
df_f <- df.merge
df_f <- df_f[rowSums(df_f >= 1) >= 5, ]
df_f <- df_f[rowSums(df_f >= 2) >= 2, ]
df_f <- df_f[rowSums(df_f >= 3) >= 1, ]

mu <- rowMeans(df_f)
sigma <- sqrt(rowVars(df_f))
cv <- sigma/mu
para <- fit_CV(mu, cv)

df_dev <- df_f[order(para$score, decreasing = T)[1:thrs], ]

## Filter variable genes in the whole ES dataset
thrs = 800
df_f = es.cef$expr.mat
df_f <- df_f[rowSums(df_f >= 1) >= 7, ]
df_f <- df_f[rowSums(df_f >= 2) >= 3, ]
df_f <- df_f[rowSums(df_f >= 3) >= 1, ]

mu <- rowMeans(df_f)
sigma <- sqrt(rowVars(df_f))
cv <- sigma/mu
para <- fit_CV(mu, cv)

variable_es <- rownames(df_f[order(para$score, decreasing = T)[1:thrs], ])

## Avoid to select genes that are cell-culture specific 

thrs = 4500
df_f = df
df_f <- df_f[rowSums(df_f >= 1) >= 5, ]
df_f <- df_f[rowSums(df_f >= 2) >= 2, ]
df_f <- df_f[rowSums(df_f >= 3) >= 1, ]

mu <- rowMeans(df_f)
sigma <- sqrt(rowVars(df_f))
cv <- sigma/mu
para <- fit_CV(mu, cv)

variable_emb <- rownames(df_f[order(para$score, decreasing = T)[1:thrs], ])

## make final preparation
es_batch_gene <- setdiff(setdiff(rownames(df_dev), variable_emb), variable_es)
df_dev <- df_dev[-match(es_batch_gene, rownames(df_dev)),]



# Calculate the probability -----------------------------------------------

# prototype dictionary

proto = list(
  'hEndo'    = 'Endo',
  'hPeric'   = 'Peric',
  'hMgl'     = 'Mgl',
  'hDA1'     = 'DA',
  'hDA2'     = 'DA',
  'hDA0'     = 'DA',
  'hSert'    = 'Sert',
  'hOMTN'    = 'OMTN',
  'hRgl1'    = 'Rgl',
  'hRgl3'    = 'Rgl',
  'hRgl2c'   = 'Rgl',
  'hRgl2b'   = 'Rgl',
  'hRgl2a'   = 'Rgl',
  'hOPC'     = 'OPC',
  'hProgFPM' = 'ProgFP',
  'hProgFPL' = 'ProgFP',
  'hProgM'   = 'ProgFP',
  'hProgBP'  = 'ProgBP',
  'hNbML5'   =  'Gaba',
  'hGaba'    =  'Gaba',
  'hNbGaba'  =  'Gaba',
  'hNbML1'   =  'NbML1',
  'hNProg'   =  'NProg',
  'hNbM'     =  'NbM',
  'hRN'      =  'RN',
  'eSCa'     =  'eES',
  'eSCb'     =  'eES',
  'eSCc'     =  'eES',
  'Unk'      =  'none'
)
proto <- unlist(proto)

cols_annot_all <- rbind(he.cef$sample.info,
                        es.cef$sample.info,
                        ips.cef$sample.info)
ct_dev <- cols_annot_all$Cell_type[match(colnames(df_dev), cols_annot_all$Cell_ID) ]
protogroup <- proto[ct_dev]
# write.table(data.frame("sample" = colnames(df_dev), "name" = names(protogroup), "value" = protogroup),
#             file.path(res.path, "protogroup.txt"),
#             row.names = F, col.names = T, quote = F)


## prepare the reference

df_dev_log = log2(df_dev + 1)
df_ips_log = log2(ips.cef$expr.mat + 1)
df_es_log = log2(es.cef$expr.mat + 1 )

ct_dev <- cols_annot_all$Cell_type[match(colnames(df_dev), cols_annot_all$Cell_ID) ]
protogroup <- proto[ct_dev]

bool1 <- protogroup != "none"
# classes_names = unique(protogroup[bool1])
classes_names <- c('DA', 'Endo', 'Gaba', 'Mgl', 'NProg', 'NbM', 'NbML1', 'OMTN',
                   'OPC', 'Peric', 'ProgBP', 'ProgFP', 'RN', 'Rgl', 'Sert', 'eES')
classes_index <- sapply(protogroup[bool1], function(i) which(i == classes_names))

train_index <- classes_index
df_train_set <- df_dev_log[, bool1]

## training model

normalizer <- 0.9 * rowMaxs(df_train_set)
class_weight <- 1/(table(train_index))
CrossValidation = F
if (CrossValidation == T){
  fit <- cv.glmnet(x = t(df_train_set/normalizer),# or res.path/LR.InputData.csv
                   y = train_index,# or res.path/LR.label.csv
                   standardize = F,
                   intercept = F,
                   weights = c(class_weight[train_index]),
                   family = "multinomial",
                   type.multinomial = "grouped",
                   lambda = 10^(py_seq(-3.25, 0.8, 30)), # given by source code
                   nfolds = 10,
                   alpha = 0)
  View(coef(fit))
  plot(fit, label = T)
  print(fit)
  plot(fit, xvar = "lambda", label =TRUE)
  plot(fit, xvar = "dev", label =TRUE)
}

# Wheel/Polygonal plot ----------------------------------------------------

## calculate

wanted_order = c('NProg', 'ProgFP','eES', 'ProgBP', 'Rgl', 
                 'Gaba', 'Sert','DA','OMTN','RN', 'NbM','NbML1')
reorder_ix = match(wanted_order, classes_names)
bool00 = classes_names[classes_index] %in% wanted_order

final.fit = glmnet(x = t(df_train_set/normalizer), # or res.path/LR.InputData.csv
                   y = train_index, # or res.path/LR.label.csv
                   standardize = F,
                   intercept = F,
                   weights = c(class_weight[train_index]),
                   family = "multinomial",
                   type.multinomial = "grouped",
                   lambda = 0.2, # given by source code
                   alpha = 0)


PlotReproduce <- T # if want to reproduce the plot in jupyter, set True
if (PlotReproduce == T){
  predicted.data = as.matrix(read.csv(file.path(res.path, "PredictedProbability.csv"), colClasses = "numeric")[, -1])
}else{
  predicted.data = predict(object = final.fit,
                           newx = t(df_train_set/rowMaxs(df_train_set)),
                           type = "response",
                           exact = T)[, , 1][, reorder_ix]
}
summary(predicted.data)

sides = length(reorder_ix)
start_angle=90
basis = t(sapply((1:sides)-1, function(i){
  c(cos(2*i*pi/sides + start_angle*pi/180),
    sin(2*i*pi/sides + start_angle*pi/180))
}))
data = predicted.data %*% basis

## prepare the plot data
### data for point
plot.data <- data.frame(
  "x" = data[bool00, 1],
  "y" = data[bool00, 2],
  "group" = names(classes_index)[bool00]
)
### data for the polygonal
basis <- data.frame(
  "x" = basis[, 1],
  "y" = basis[, 2]
)
### data for the class label
label <- data.frame(
  "x" = 1.1 * basis[, 1],
  "y" = 1.1 * basis[, 2],
  "label" = wanted_order
)

### color for the plot
color.palette <- list(
  'hEndo'=   c(190,  10,  10),'hPeric'= c(225, 160,  30),'hMgl'=    c(217, 245,   7),
  'hDA1'=    c(170, 180, 170),'hDA2'=   c(130, 140, 140),'hNbM'=    c(180, 140, 130),
  'hNbML1'=  c(100, 100, 240),'hProgM'= c( 80, 235, 255),'hProgFPM'=c(190, 235, 255),
  'hProgFPL'=c(210, 255, 215),'hProgBP'=c(230, 140, 120),'hNProg'=  c(255, 195,  28),
  'hNbML5'=  c(139, 101, 100),'hRgl1'=  c(252, 183,  26),'hRgl3'=   c(214, 194,  39),
  'hRgl2c'=  c(255, 120, 155),'hRgl2b'= c(250, 145,  45),'hRgl2a'=  c(250, 125,  25),
  'hDA0'=    c(190, 200, 190),'hOPC'=   c(255,  35, 155),'hRN'=     c(199, 121,  41),
  'hNbGaba'= c( 40,  55, 130),'hGaba'=  c(  7,  121, 61),'hOMTN'=   c( 95, 186,  70),
  'hSert'=   c( 50, 180, 180),'eSCa'=   c(245, 205, 170),'eSCb'=    c(205, 245, 170),
  'eSCc'=    c(205, 205, 220)
)
color.palette <- lapply(color.palette, function(x){
  rgb(red = x[1]/255, green = x[2]/255, blue = x[3]/255)
})
color.palette <- unlist(color.palette)

## plot 
ggplot(plot.data, aes(x = x, y = y)) +
  geom_point(aes(col = group)) +
  geom_polygon(mapping = aes(x = x, y = y), data = basis,
               fill = "transparent", color = "black", size = 1) +
  geom_text(mapping = aes(x = x, y = y, label = label), data = label) +
  scale_color_manual(values = color.palette) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.background = )
if (PlotReproduce == T)
  ggsave(file.path(fig.path, "Example_using_Raw_Output.pdf"), width = 10, height = 8)
if (PlotReproduce == F)
  ggsave(file.path(fig.path, "Example_using_R_Output.pdf"), width = 10, height = 8)
