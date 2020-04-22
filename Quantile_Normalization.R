quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

rownames(df) <- make.names(df[,1], unique = TRUE)
df<-df[,-1]

#Example
df <- data.frame(one=c(5,2,3,4),
                 two=c(4,1,4,2),
                 three=c(3,4,6,8)
)

#Function Calling

quantile_normalisation(df)
#       one      two    three
#A 5.666667 4.666667 2.000000
#B 2.000000 2.000000 3.000000
#C 3.000000 4.666667 4.666667
#D 4.666667 3.000000 5.666667

