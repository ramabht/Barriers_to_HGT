GrowthProfiler960
================
RP Bhatia
16 September 2022

# data analysis performed in paper "Evolutionary barriers to horizontal gene transfer in macrophage associated Salmonella" doi: <https://doi.org/10.1101/2022.04.01.486712>

# library load

``` r
library(knitr)
```

``` r
library(data.table)
library(fs)
```

``` r
library(dplyr)
```

# Calculating growth rates from log linear part of the growth curve using linear regression

``` r
file_path = fs::dir_ls("C:\\Users\\ramab\\Documents\\R\\gene_dosage_new\\GP960\\input_files")
file_path
```

    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b0119-b0423_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b0440-b0695_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b0780-b0882_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b0891-b1094_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b1763-b1913_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b2341-b3006_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b3071-b3590_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b3602-b4043_inspi2.csv
    ## C:/Users/ramab/Documents/R/gene_dosage_new/GP960/input_files/b4172-b4373_inspi2.csv

``` r
file_contents <- list()
samples <- list()
for(i in seq_along(file_path)){
  file_contents[[i]] <- read.csv(file = file_path[[i]], header = T, check.names = F)
  names(file_contents[[i]]) <- make.unique(names(file_contents[[i]]), sep="_")
  samples[[i]] <- file_contents[[i]][,-1]
  samples[[i]][samples[[i]]<= 0] <- 0.001
}

file_contents <- setNames(file_contents, file_path)
samples <- setNames(samples, file_path)


names(file_contents) <- basename(file_path) 
# extract file name from file path
names(samples) <- basename(file_path)

names(file_contents) <- tools::file_path_sans_ext(names(file_contents))
names(samples) <- tools::file_path_sans_ext(names(samples))

file_contents2 <- file_contents[-9]
samples2 <- samples[-9]
```

# by using the growth rate best line equation ax+b=y, a is the slope = growth rate, b is the intercept, and R is the regression coefficient

``` r
line=list()
GR_a = rep(0,120)
GR_b = rep(0,120)
GR_R = rep(0,120)

for (i in 1:length(file_contents2)) {
  name = names(file_contents2[i])
  file <- paste0(name, ".pdf")
  pdf(file , paper = "a4r", width = 14, height = 18)
  par(mfrow=c(12,10), mar=c(0,0,0,0), oma=c(5,5,5,5))
  time <- file_contents2[[i]]$Time_mins
  for (j in 1:120) {
    OD=file_contents2[[i]][j+1] #take OD values for each well into OD
    lnOD = log(OD[(OD>=0.05) & (OD<=0.18)]) #select the points between the OD range, take the ln and write them into vector lnOD
    if (length(lnOD)==0){
      GR_a[j] = 0
      GR_b[j] = 0
      GR_R[j] = 0
    }
    else {
      x = (OD>=0.05)&(OD<=0.18)
      mins = c(((which(x==1))*20)-20)
      line = lm(lnOD~mins) #applies regression analysis
      GR_a[i] = summary(line)$coefficients[2,1] # write slope in a
      GR_b[j] = summary(line)$coefficients[1,1] # write intercept in b
      GR_R[j] = summary(line)[8]# write regression coefficient in R
      if (dim(summary(line)$coefficients)[[1]]==1){
        GR_a[j]<- NA
      } 
      else {
        GR_a[j] = summary(line)$coefficients[2,1]
      }
    }
    plot(time, log(file_contents2[[i]][,j+1]), type = "l", 
         xlab="", ylab="", xaxt="n", yaxt="n")
    abline(GR_b[j],GR_a[j], col='red')
    SIR = matrix(nrow = 121, ncol= 3)
    SIR[,1] = c("Well ID", colnames(file_contents2[[i]][-1]))
    SIR[,2] = c("GrowthRate", GR_a)
    SIR[,3] = c("R", unlist(GR_R))
    filename <- paste("gr_", name, ".csv")
    write.table(SIR,filename,row.names= FALSE, col.names= FALSE, 
                sep=",", quote = FALSE)
  }
  dev.off()
}
```


# For Plate 9 only

``` r
GR_a = rep(0,100)
GR_b = rep(0,100)
GR_R = rep(0,100)

for (i in 1:length(file_contents)) {
  name = names(samples[9])
  file <- paste0(name, ".pdf")
  pdf(file , paper = "a4r", width = 14, height = 18)
  par(mfrow=c(10,10), mar=c(0,0,0,0), oma=c(5,5,5,5))
  time <- file_contents[[9]]$Time_mins
  for (j in 1:100) {
    OD=file_contents[[9]][j+1] #take OD values of growth for one well into OD
    lnOD = log(OD[(OD>=0.05) & (OD<=0.18)]) #select the points between the OD range, take the ln and write them into vector lnOD
    if (length(lnOD)==0){
      GR_a[j] = 0
      GR_b[j] = 0
      GR_R[j] = 0
    }
    else {
      x = (OD>=0.05)&(OD<=0.18)
      mins = c(((which(x==1))*20)-20)
      line = lm(lnOD~mins) #applies regression analysis
      GR_a[i] = summary(line)$coefficients[2,1] # write slope in a
      GR_b[j] = summary(line)$coefficients[1,1] # write intercept in b
      GR_R[j] = summary(line)[8]# write regression coefficient in R
      if (dim(summary(line)$coefficients)[[1]]==1){
        GR_a[j]<- NA
      } 
      else {
        GR_a[j] = summary(line)$coefficients[2,1]
      }
    }
        plot(time, log(file_contents[[9]][,j+1]), type = "l", 
             xlab="", ylab="",
             xaxt="n", yaxt="n")
        abline(GR_b[j],GR_a[j], col='red')
        SIR = matrix(nrow = 101, ncol= 3)
        SIR[,1] = c("Well ID", colnames(file_contents[[9]][-1]))
        SIR[,2] = c("GrowthRate", GR_a)
        SIR[,3] = c("R", unlist(GR_R))
        filename <- paste("gr_", name, ".csv")
        write.table(SIR,filename,row.names= FALSE, col.names= FALSE, 
                    sep=",", quote = FALSE)
  }
  dev.off()
}
```

