workdir <- "I:/genomicdata/External/Xlu/WheelPlot"; setwd(workdir)
fig.path <- file.path(workdir,"Figures")
res.path <- file.path(workdir,"Results")
data.path <- file.path(workdir,"InputData")

invisible(lapply(ls()[grep("path", ls())], function(x){
  if (!dir.exists(get(x))) dir.create(get(x))
}))


## read table 

raw.data <- read.table(file.path(res.path, "raw.data.txt"), header = T, stringsAsFactors = T)

## calculate the basis of axis switch
sides = ncol(raw.data)-1
start_angle=90
basis = t(sapply((1:sides)-1, function(i){
  c(cos(2*i*pi/sides + start_angle*pi/180),
    sin(2*i*pi/sides + start_angle*pi/180))
}))
## cat the 12 dimensional coordinate to 2 dimensional coordinate
data = as.matrix(raw.data[, -match("celltype", colnames(raw.data))]) %*% basis

## prepare the plot data
### data for point
plot.data <- data.frame(
  "x" = data[, 1],
  "y" = data[, 2],
  "group" = raw.data$celltype
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
  "label" = colnames(raw.data)[-13]
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
        axis.text = element_blank())

ggsave(file.path(fig.path, "Example_using_Raw_Output.pdf"), width = 10, height = 8)

