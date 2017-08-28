
# One-Way S3 Implementation


In this tutorial we try to develop an S3 class to estimate the relevant parameters in a one-way AOV with possibly unequal sample sizes (e.g., see Dowdy and Wearden, Statistics for Research, Wiley). This code is part of my project for Statistical Computing with Prof. Harner.

#### 1. Develop a default method, **oneway.default**, for the generic function **oneway**.

``` r
oneway <- function(z, ...) UseMethod("oneway")  

  oneway.default <- function(z, ...) {
      if (!is.data.frame(z)) {  
       df <- data.frame(matrix(unlist(z), nrow=length(unlist(z)), byrow=T));df
      ##first one
       
      quant <- lapply(z[] ,quantile)
      
       remove(f)
       #remove(l)
       f<-numeric()
       l <- numeric()
       for (i in seq(1,length(z))){
        
         f <-  c(f,replicate(length(z[[i]]),paste("fact",i) ))
         l[i] <- length(z[[i]]); #print(l)
         
      }
      f
      unlist(f)
       df2 <- data.frame(matrix(unlist(f), nrow=length(unlist(f)), byrow=T));df2
      ### Second one
      
      print(df)
      z <- cbind(df,df2)
      }
      
    ###################  
    
    
      l <- length(levels(z[,1])) # If z[,1] is not the column of labels then its length will be zero so it   means that this column is assigned for the data and the next column which is z[,2] is assigned for the labels.
      
      if (l == 0)   {lev <- levels(z[, 2]) ; n=2;m=1 } else {lev <- levels(z[, 1]); n=1;m=2}
       lev ;n;m; # This line will find which column of the input z is assigend for the labels. as a consequence the other column contains the data.
      varOfGroups=0 
       meanOfGroups=0
      remove(quant)
       quant <- numeric()
      
      for (i in 1:length(lev)) {
        meanOfGroups[i] <- mean(z[which(z[,n]== lev[i]),m]) # this line will find the mean for the each group ; # More detials of the line is in the following line
        # which(z[,2]== lev[i] this part of code will select the groups that have the same lable as lev[i]  then by assinging the number of those rows to this code z[which(z[,2]== lev[i]),1] the code will select those rows but from the other column which has the data. Finally by apllying mean the mean of different groups will be found.
      l[i] <- length(which(z[,n]== lev[i])); #print("inja"); print(l)
      quant <-  c(quant, quantile(c(z[which(z[,2]== lev[i]),1])));quant
      }
       print(quant)
      meanOfGroups ; 
      #meanOfMeans <- mean(meanOfGroups) ; meanOfMeans 
      meanOfMeans <- mean(z[,m])
      
      
      S_BG=0;
      
      for(i in 1:length(lev))
      { S_BG[i] <- length(which(z[,n] == lev[i])) *(meanOfGroups[i] - meanOfMeans)^2}
      LSM <- sum(meanOfGroups)/length(lev) # least squares mean
      SS_BG = sum(S_BG) ; SS_BG  #Sum of Squered Between Group
      df_BG <- (length(lev)-1) ;df_BG #Degree of Freedon Between Group
      # meanOfSS_BG <- (SST_BG / df_BG); meanOfSS_BG; # Mean of SS Between Groups
      S_WG=0; # Squar within Groups
      for (i in 1:length(lev)) {
        S_WG[i] <- sum((z[which(z[,n]== lev[i]),m]-meanOfGroups[i])^2);
      }
      SS_WG <- sum(S_WG) ; SS_WG # Sum of Squeared within Groups
      df_WG <- length(z[,n]) - length(lev); df_WG   # Degree of Freedom in between N-K
      ## meanOfSS_WG <- SS_WG/df_WG ; meanOfSS_WG
      ## F_test <- meanOfSS_BG / meanOfSS_WG ;F_test
      ## oneway <- list(SS_BG,SS_WG,meanOfSS_BG, df_BG,df_WG, meanOfSS_WG ,F_test); return(oneway);
      l <- matrix(l,nrow=2, byrow=T); print(l)
      oneway <- list(df_BG,df_WG,SS_BG,SS_WG,LSM,l); oneway; print("oneway");print(oneway)
    # print(l)
      # return(oneway) 
     #print(oneway)
      
      oneway <- oneway.table(oneway)
      #print(oneway)
      print("meanOfGRoups_itself"); print(meanOfGroups)
      meanOfGroups <- matrix(meanOfGroups, nrow=2 , byrow=T); 
      print("table_meanofgroups") ; print(meanOfGroups)
      oneway <- cbind(oneway,meanOfGroups)
      print("cbind_oneway_meanofGroups", oneway)
      quant <- matrix(quant, nrow=2 , byrow=F)
      oneway <- cbind(oneway,quant)
      print("final");print(oneway)
      
      return(oneway)
  }

  oneway.table <- function(x) {
  #Since the first Question didn't talk about the F test and some other values in ANOVA table, I assume that those values should be calculated in this part.   #If this assumpton is not true then the comment lines that have ## should be uncommented in the first question.
  
  x <- unlist(x); x; print("x");print(x)
  NOb <- numeric() 
  l_2 <- 6+x[1]; print("l_2"); print(l_2) ; print(6:l_2)
  NOb <- x[6:l_2] # Number of Observations which is needed for calculation Fisher
  print("NOb"); print(NOb)
 
  meanOfSS_WG=0;
  meanOfSS_WG <- x[4]/x[2];  meanOfSS_WG # In this line of the code the meanOfSS_WG is calculated by this division of SS_WG/df_WG
  
  meanOfSS_BG=0;
  meanOfSS_BG <- x[3]/x[1];  meanOfSS_BG # In this line of the code the meanOfSS_BG is calculated by this division of SS_BG/df_BG
  
  meanT <- c(meanOfSS_BG, meanOfSS_WG); meanT
  dim(meanT) <- c(2,1); meanT
  
  #remove(F_test)
  F_test <- meanOfSS_BG / meanOfSS_WG ;F_test # Now F test is calculated here
  F_test <- matrix(F_test,2,1)
  dim(F_test) <- c(2,1); F_test
  
  
  #remove(p_Value)
  p_Value <- 1 - pf(F_test, x[1] ,x[2] ); p_Value; # Since I saw that builtin function of AOV in the R is calculating this parameter so I have decided to calculate this paramter too.
  #p_Value <- c(p_Value,p_Value) ; p_Value
  dim(p_Value) <- c(2,1) ; p_Value
  
  LSM <- x[5]
  x <- x[1:4]
  dim(x) <- c(2,2); x
  
  NOb <- matrix(NOb, nrow=2 , byrow=TRUE ) 
  print(NOb)
  
  print("cbind"); print(cbind(p_Value,NOb)) ;
 

   output.table<- cbind(x,meanT,F_test,p_Value,LSM)
  dimnames(output.table) <- list(c("diet","Residuals"),c( "Df","Sum Sq","Mean sq","F value","Pr(>F)","LSM")) ;output.table ; print(output.table)   
  
  output.table<- cbind(x,meanT,F_test,p_Value,LSM,NOb)
  print("table"); print(output.table)
  print(output.table)
  #print(output.table)
  ##output <- as.table(output) # if I want to have the output strictly in table format I can use this command
  
  return(output.table)
  }

  
  ## can check with 
z <- list(A=c(62, 60, 63, 59), B=c( 63 ,67 ,71, 64, 65, 66 ), C=c(68, 66, 71, 67, 68, 68),D=c(56 ,62,60, 61, 63, 64, 63, 59))
  oneway.default(z)
```

    ## Warning in remove(f): object 'f' not found

    ##    matrix.unlist.z...nrow...length.unlist.z....byrow...T.
    ## 1                                                      62
    ## 2                                                      60
    ## 3                                                      63
    ## 4                                                      59
    ## 5                                                      63
    ## 6                                                      67
    ## 7                                                      71
    ## 8                                                      64
    ## 9                                                      65
    ## 10                                                     66
    ## 11                                                     68
    ## 12                                                     66
    ## 13                                                     71
    ## 14                                                     67
    ## 15                                                     68
    ## 16                                                     68
    ## 17                                                     56
    ## 18                                                     62
    ## 19                                                     60
    ## 20                                                     61
    ## 21                                                     63
    ## 22                                                     64
    ## 23                                                     63
    ## 24                                                     59
    ##    0%   25%   50%   75%  100%    0%   25%   50%   75%  100%    0%   25% 
    ## 59.00 59.75 61.00 62.25 63.00 63.00 64.25 65.50 66.75 71.00 66.00 67.25 
    ##   50%   75%  100%    0%   25%   50%   75%  100% 
    ## 68.00 68.00 71.00 56.00 59.75 61.50 63.00 64.00 
    ##      [,1] [,2]
    ## [1,]    4    6
    ## [2,]    6    8
    ## [1] "oneway"
    ## [[1]]
    ## [1] 3
    ## 
    ## [[2]]
    ## [1] 20
    ## 
    ## [[3]]
    ## [1] 228
    ## 
    ## [[4]]
    ## [1] 112
    ## 
    ## [[5]]
    ## [1] 64
    ## 
    ## [[6]]
    ##      [,1] [,2]
    ## [1,]    4    6
    ## [2,]    6    8
    ## 
    ## [1] "x"
    ## [1]   3  20 228 112  64   4   6   6   8
    ## [1] "l_2"
    ## [1] 9
    ## [1] 6 7 8 9
    ## [1] "NOb"
    ## [1] 4 6 6 8
    ##      [,1] [,2]
    ## [1,]    4    6
    ## [2,]    6    8
    ## [1] "cbind"
    ##              [,1] [,2] [,3]
    ## [1,] 4.658471e-05    4    6
    ## [2,] 4.658471e-05    6    8
    ##           Df Sum Sq Mean sq  F value       Pr(>F) LSM
    ## diet       3    228    76.0 13.57143 4.658471e-05  64
    ## Residuals 20    112     5.6 13.57143 4.658471e-05  64
    ## [1] "table"
    ##                                        LSM    
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8
    ##                                        LSM    
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8
    ## [1] "meanOfGRoups_itself"
    ## [1] 61 66 68 61
    ## [1] "table_meanofgroups"
    ##      [,1] [,2]
    ## [1,]   61   66
    ## [2,]   68   61
    ## [1] "cbind_oneway_meanofGroups"
    ## [1] "final"
    ##                                        LSM                               
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6 61 66 59.00 61.00 63 64.25
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8 68 61 59.75 62.25 63 65.50
    ##                                
    ## [1,] 66.75 66.00 68 71 59.75 63
    ## [2,] 71.00 67.25 68 56 61.50 64

    ##                                        LSM                               
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6 61 66 59.00 61.00 63 64.25
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8 68 61 59.75 62.25 63 65.50
    ##                                
    ## [1,] 66.75 66.00 68 71 59.75 63
    ## [2,] 71.00 67.25 68 56 61.50 64

The *z* argument for **oneway.default** should be a list of (possibly) named components, one for each sample. The computations for the one-way ANOVA should be done in **oneway.default**.

#### 2. This method uses the more standard input of a factor representing groups (or samples) and a numeric response.

``` r
oneway.factor <- function(z, y, ...) {
## Your code here
  
      df2 <- data.frame(matrix(unlist(z), nrow=length(unlist(z)), byrow=T));df2 ##factors
      ### Second one
       df <- data.frame(matrix(unlist(y), nrow=length(unlist(y)), byrow=T));df ## values
       
      z <- cbind(df,df2) # combine values and factors
    
    ###################
  o <- oneway.default(z)

return(o)
}
  
  ##can check with 
  #factor = c(rep("A",4),rep("B",6),rep("C",6),rep("D",8))
  #values = c( 62, 60, 63, 59,63 ,67 ,71, 64, 65, 66 ,68, 66, 71, 67, 68, 68 ,56 ,62,60, 61, 63, 64, 63, 59)
  #oneway.factor(factor,,values)
```

The *z* argument should be the factor with levels representing samples. The factor should be used to deconstruct *y*, the response, into a list as in the default.

#### 3. The model formula is the standard for R models, but do not use **model.matrix** to implement **oneway**.

``` r
oneway.formula <- function(formula, data=list(), ...) {
## Your code here
  out <- model.frame(formula,data) 
  fac <- out[,1]
  values <- out[,2]
  o2 <- oneway.factor(fac,values)
  return(o2)
}
  
  ## to check we can run
  z <- list( factor= c(rep("A",4),rep("B",2)) , factor2= c(rep("F",4)) , values = c( 63 ,67 ,71, 64, 65, 66), values2 =c(23,45,67) ) 
  
 #oneway.formula(factor ~ values , z)
```

You might want to extract the factor term and the response from the **model.frame** and then call **oneway.factor**, which in turn calls **oneway.default**. \#\#\#???? shall I use oneway.factor insides oneway.default?

#### 4. The default **print** method should be short and provide essential information.

``` r
print.oneway <- function(x, ...) {
## Your code here
  o <- x
  o <- o[1:4]
  dim(o) <- c(2,2)
  dimnames(o) <- list(c("diet","Residuals"),c( "Df","Sum Sq")) 
  print(o)
  
}
  
  # to check run
  # z <- list(A=c(62, 60, 63, 59), B=c( 63 ,67 ,71, 64, 65, 66 ), C=c(68, 66, 71, 67, 68, 68),D=c(56 ,62,60, 61, 63, 64, 63, 59))
#  print.oneway(oneway.default(z))
```

#### 5. The summary method should create a summary object---not print directly.

``` r
summary.oneway <- function(object, ...) {
## Your code here
   # o <- oneway.default(object)
  o <- object
  o <- o[1:10]
  dim(o) <- c(2,5)
  dimnames(o) <- list(c("diet","Residuals"),c( "Df","Sum Sq","Mean sq", "F value", "Pr(>F)")) 
  #o
 # return o ##???   ### inherits? how and why?
}
 
 # to check 
  z <- list(A=c(62, 60, 63, 59), B=c( 63 ,67 ,71, 64, 65, 66 ), C=c(68, 66, 71, 67, 68, 68),D=c(56 ,62,60, 61, 63, 64, 63, 59))
 #  summary.oneway(oneway.default(z))   
```

The argument is a **oneway** object. The summary object should include a component representing an AOV table, e.g., see Dowdy and Wearden. You might want to decide whether objects of class *summary.oneway* inherit from the class *oneway*.

#### 6. The print method for the summary object should provide more detailed information about the summary object.

``` r
print.summary.oneway <- function(x, ...) {
## Your code here

  o <- x
  o <- o[1:12]
  dim(o) <- c(2,6)
  dimnames(o) <- list(c("diet","Residuals"),c( "Df","Sum Sq","Mean sq", "F value", "Pr(>F)","LSM")) 
  print(o)
 }
 
 # to check 
  #z <- list(A=c(62, 60, 63, 59), B=c( 63 ,67 ,71, 64, 65, 66 ), C=c(68, 66, 71, 67, 68, 68),D=c(56 ,62,60, 61, 63, 64, 63, 59))
 #print.summary.oneway(oneway.default(z))  
   ## I add LSM to the AOV table. It could be printed seperately too but i prefered this format.
```

The AOV table should be formatted nicely. The least squares means should also be formated and printed.

#### 7. Implement Fisher's LSD multiple comparison procedure for your oneway.

``` r
lsmeans.oneway <- function(object, ...) {
## Your code here
   o <- unlist(object)
   o[2];  print(o[2]) #DFW
   o[6];  print(o[6])   #MSW 
   
   o[1] ;
   ob <- numeric()  # ??????????????????????????? rooye 13 deghat konam
   ob <- o[13:(13+o[1])]; # vector inculding the Number of Observations in each group
   print(o[13:(13+o[1])])# DFB 
   print("o"); print(o);
   meanofGroups <- o[ (12+o[1]+2) : (12 + o[1]+2 + o[1])]
   print("meanofgroups222"); print(meanofGroups)
   
   LSD <- numeric();
   sig_dif <- numeric()
   
    
   LSD <- matrix(LSD,nrow= (o[1]+1) ,ncol= (o[1]+1) )
   sig_dif <- matrix(sig_dif,nrow= (o[1]+1) ,ncol= (o[1]+1) )
   
   
   a <- numeric
   for ( i in seq(1,(o[1])+1)) 
     for ( j in seq(i+1,(o[1])+1)) {
      if (i!=o[1]+1) {
      LSD[i,j] <- ( qt(0.975,o[6]) * sqrt( o[6] * ( ( 1/ ob[i]) + (1/ ob[j]) )) ) ; print("result"); print(LSD)
      a <- (abs(meanofGroups[i]-meanofGroups[j]))
      print("a"); print(a)
      if (a < LSD[i,j]) 
      {sig_dif[i,j] <- T ; 
       cat('LSD between group =\"',i,'\" and group =\"',j , '\"' , 'is' , '\"', LSD[i,j], '\"') ;
       cat('No Significance is found on Mean difference between Group' ,i, '& Group', j , 'by knowing that  the difference between their mean is', a, '\n');
      }
      else 
      {sig_dif[i,j] <- F ; 
      cat('LSD between group =\"',i,'\" and group =\"',j , '\"' , 'is' , '\"', LSD[i,j], '\"') ;
      cat('Significant difference is found on Mean difference between Group' ,i, '& Group', j , 'by knowing that  the difference between their mean is', a, '\n');
      }
     }
   
     }
   print(sig_dif)
   print(LSD)
}
   #to check it 
    z <- list(A=c(62, 60, 63, 59), B=c( 63 ,67 ,71, 64, 65, 66 ), C=c(68, 66, 71, 67, 68, 68),D=c(56 ,62,60, 61, 63, 64, 63, 59))
  lsmeans.oneway(oneway.default(z))
```

    ## Warning in remove(f): object 'f' not found

    ##    matrix.unlist.z...nrow...length.unlist.z....byrow...T.
    ## 1                                                      62
    ## 2                                                      60
    ## 3                                                      63
    ## 4                                                      59
    ## 5                                                      63
    ## 6                                                      67
    ## 7                                                      71
    ## 8                                                      64
    ## 9                                                      65
    ## 10                                                     66
    ## 11                                                     68
    ## 12                                                     66
    ## 13                                                     71
    ## 14                                                     67
    ## 15                                                     68
    ## 16                                                     68
    ## 17                                                     56
    ## 18                                                     62
    ## 19                                                     60
    ## 20                                                     61
    ## 21                                                     63
    ## 22                                                     64
    ## 23                                                     63
    ## 24                                                     59
    ##    0%   25%   50%   75%  100%    0%   25%   50%   75%  100%    0%   25% 
    ## 59.00 59.75 61.00 62.25 63.00 63.00 64.25 65.50 66.75 71.00 66.00 67.25 
    ##   50%   75%  100%    0%   25%   50%   75%  100% 
    ## 68.00 68.00 71.00 56.00 59.75 61.50 63.00 64.00 
    ##      [,1] [,2]
    ## [1,]    4    6
    ## [2,]    6    8
    ## [1] "oneway"
    ## [[1]]
    ## [1] 3
    ## 
    ## [[2]]
    ## [1] 20
    ## 
    ## [[3]]
    ## [1] 228
    ## 
    ## [[4]]
    ## [1] 112
    ## 
    ## [[5]]
    ## [1] 64
    ## 
    ## [[6]]
    ##      [,1] [,2]
    ## [1,]    4    6
    ## [2,]    6    8
    ## 
    ## [1] "x"
    ## [1]   3  20 228 112  64   4   6   6   8
    ## [1] "l_2"
    ## [1] 9
    ## [1] 6 7 8 9
    ## [1] "NOb"
    ## [1] 4 6 6 8
    ##      [,1] [,2]
    ## [1,]    4    6
    ## [2,]    6    8
    ## [1] "cbind"
    ##              [,1] [,2] [,3]
    ## [1,] 4.658471e-05    4    6
    ## [2,] 4.658471e-05    6    8
    ##           Df Sum Sq Mean sq  F value       Pr(>F) LSM
    ## diet       3    228    76.0 13.57143 4.658471e-05  64
    ## Residuals 20    112     5.6 13.57143 4.658471e-05  64
    ## [1] "table"
    ##                                        LSM    
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8
    ##                                        LSM    
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8
    ## [1] "meanOfGRoups_itself"
    ## [1] 61 66 68 61
    ## [1] "table_meanofgroups"
    ##      [,1] [,2]
    ## [1,]   61   66
    ## [2,]   68   61
    ## [1] "cbind_oneway_meanofGroups"
    ## [1] "final"
    ##                                        LSM                               
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6 61 66 59.00 61.00 63 64.25
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8 68 61 59.75 62.25 63 65.50
    ##                                
    ## [1,] 66.75 66.00 68 71 59.75 63
    ## [2,] 71.00 67.25 68 56 61.50 64
    ## [1] 20
    ## [1] 5.6
    ## [1] 4 6 6 8
    ## [1] "o"
    ##                                        LSM                               
    ## [1,]  3 228 76.0 13.57143 4.658471e-05  64 4 6 61 66 59.00 61.00 63 64.25
    ## [2,] 20 112  5.6 13.57143 4.658471e-05  64 6 8 68 61 59.75 62.25 63 65.50
    ##                                
    ## [1,] 66.75 66.00 68 71 59.75 63
    ## [2,] 71.00 67.25 68 56 61.50 64
    ## [1] "meanofgroups222"
    ## [1] 61 68 66 61
    ## [1] "result"
    ##      [,1]   [,2] [,3] [,4]
    ## [1,]   NA 3.8034   NA   NA
    ## [2,]   NA     NA   NA   NA
    ## [3,]   NA     NA   NA   NA
    ## [4,]   NA     NA   NA   NA
    ## [1] "a"
    ## [1] 7
    ## LSD between group =" 1 " and group =" 2 " is " 3.8034 "Significant difference is found on Mean difference between Group 1 & Group 2 by knowing that  the difference between their mean is 7 
    ## [1] "result"
    ##      [,1]   [,2]   [,3] [,4]
    ## [1,]   NA 3.8034 3.8034   NA
    ## [2,]   NA     NA     NA   NA
    ## [3,]   NA     NA     NA   NA
    ## [4,]   NA     NA     NA   NA
    ## [1] "a"
    ## [1] 5
    ## LSD between group =" 1 " and group =" 3 " is " 3.8034 "Significant difference is found on Mean difference between Group 1 & Group 3 by knowing that  the difference between their mean is 5 
    ## [1] "result"
    ##      [,1]   [,2]   [,3]     [,4]
    ## [1,]   NA 3.8034 3.8034 3.608222
    ## [2,]   NA     NA     NA       NA
    ## [3,]   NA     NA     NA       NA
    ## [4,]   NA     NA     NA       NA
    ## [1] "a"
    ## [1] 0
    ## LSD between group =" 1 " and group =" 4 " is " 3.608222 "No Significance is found on Mean difference between Group 1 & Group 4 by knowing that  the difference between their mean is 0 
    ## [1] "result"
    ##      [,1]   [,2]     [,3]     [,4]
    ## [1,]   NA 3.8034 3.803400 3.608222
    ## [2,]   NA     NA 3.401865       NA
    ## [3,]   NA     NA       NA       NA
    ## [4,]   NA     NA       NA       NA
    ## [1] "a"
    ## [1] 2
    ## LSD between group =" 2 " and group =" 3 " is " 3.401865 "No Significance is found on Mean difference between Group 2 & Group 3 by knowing that  the difference between their mean is 2 
    ## [1] "result"
    ##      [,1]   [,2]     [,3]     [,4]
    ## [1,]   NA 3.8034 3.803400 3.608222
    ## [2,]   NA     NA 3.401865 3.182153
    ## [3,]   NA     NA       NA       NA
    ## [4,]   NA     NA       NA       NA
    ## [1] "a"
    ## [1] 7
    ## LSD between group =" 2 " and group =" 4 " is " 3.182153 "Significant difference is found on Mean difference between Group 2 & Group 4 by knowing that  the difference between their mean is 7 
    ## [1] "result"
    ##      [,1]   [,2]     [,3]     [,4]
    ## [1,]   NA 3.8034 3.803400 3.608222
    ## [2,]   NA     NA 3.401865 3.182153
    ## [3,]   NA     NA       NA 3.182153
    ## [4,]   NA     NA       NA       NA
    ## [1] "a"
    ## [1] 5
    ## LSD between group =" 3 " and group =" 4 " is " 3.182153 "Significant difference is found on Mean difference between Group 3 & Group 4 by knowing that  the difference between their mean is 5 
    ##      [,1] [,2] [,3] [,4]
    ## [1,]   NA    0    0    1
    ## [2,]   NA   NA    1    0
    ## [3,]   NA   NA   NA    0
    ## [4,]   NA   NA   NA   NA
    ##      [,1]   [,2]     [,3]     [,4]
    ## [1,]   NA 3.8034 3.803400 3.608222
    ## [2,]   NA     NA 3.401865 3.182153
    ## [3,]   NA     NA       NA 3.182153
    ## [4,]   NA     NA       NA       NA

The argument is a *oneway* object, which should include the least-squares means as a component. Fisher's LSD should be computed and formatted nicely.

#### 8. A plot generic function should be implemented for *oneway* objects.

``` r
plot.oneway <- function(x, ...) {
  ## Your code here
#   o <- unlist(x); print("ooooo");print(o)
#   l1<- 2*(o[1]+1)+12 +1
#   print(l1)
#   l2<- l1+4
#   print(l2)
#   box <- numeric()
#   
#   print("bale")
#   print(seq(1:(o[1]+1)))
#   print(o[l1:l2])
#   
#   for ( i in seq(1:(o[1]+1))){
#     box[i] <-  o[l1:l2] ; print(box)
#     #box <- c(box,o[l1:l2])
#     print(box)
#     #b<- c( o[],   2 ,   3 ,   8 ,  13 )
#     l1 <- l2+1
#     l2 <- l1+4
#    ## names(box[i]) <- c("0%" , "25%" , "50%" , "75%" , "100%"); print(box)
#   }
#   
#   
#   #names(box) <- rep(c("0%" , "25%" , "50%" , "75%" , "100%" ),(o[1]+1))
#   print(box)
#   
#   for ( i in seq(1:(o[1]+1))){
#   boxplot(box[i:i+5])}
#   
 
}
  plot.oneway(oneway.default(z))
```

    ## NULL

The plot should compare the distributions of the groups in a side-by-side manner.

#### 9. Your S3 class implementation should be illustrated with the *coagulation* data set. The data consists of blood coagulation times for 24 animals randomly assigned to four different diets.

``` r
library(faraway)
data(coagulation)
coagulation[1:4,]
```

    ##   coag diet
    ## 1   62    A
    ## 2   60    A
    ## 3   63    A
    ## 4   59    A

``` r
## Your implementation code here
```

You should provide brief explanations of the output along with the output, which implies that you may want multiple chucks of R code interspersed with markdown.
