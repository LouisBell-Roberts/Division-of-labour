grand.mean <- function(M, N) {weighted.mean(M, N)}
grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                        weighted.mean(M, N)^2)}

#Create example dataframe:
spec_names<-c("a", "b", "c", "d")
head_width_mean<-c(1, 2, 1, 1.5)
sample_size<-c(5,2,3,7)
hw_2<-c(1.5)
ss_2<-c(4)
hw_2<-c(1.5, 1.5, 1.5, 0)
ss_2<-c(4,4,4,0)
SD<-c(1,3,2,4)
SD_2<-c(4,2,3,0)
df.3<-data.frame(spec_names, head_width_mean, sample_size, hw_2, ss_2, SD, SD_2)

#########Calculates the grand standard deviation for multiple groups. Missing input data must be assigned as zeros to work
#The function
grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                        weighted.mean(M, N)^2)}
##Create empty vector to assign the results of the for loop to
fj<-c()
#Standard form of the for loop - doesn't assign to the vector
for (i in 1:4) {
  print(grand.sd(c(df.3[i,6], df.3[i,7]), c(df.3[i,2], df.3[i,4]), c(df.3[i,3], df.3[i,5])))
}
#The working for loop - assigns values to the empty vector
for (i in 1:4) {
  fj<-append(fj,print(grand.sd(c(df.3[i,6], df.3[i,7]), c(df.3[i,2], df.3[i,4]), c(df.3[i,3], df.3[i,5]))))
}
##########################

#########Calculates weighted means. Missing input data must be assigned as zeros to work
grand.mean <- function(M, N) {weighted.mean(M, N)}
for (i in 1:4) {
  print(grand.mean(c(df.3[i,2], df.3[i,4]), c(df.3[i,3], df.3[i,5])))
}
##########################


##Stuff I was practicing with
grand.mean(c(df.3[1,2], df.3[1,4]), c(df.3[1,3], df.3[1,5]))

grand.sd(c(), c(), c())

g<-df.3$head_width_mean
df.3[1,]

spec_names<-c("a", "b", "c", "d")
head_width_mean<-c(1, 2, 1, 1.5)
sample_size<-c(5,2,3,7)
hw_2<-c(1.5)
ss_2<-c(4)
hw_2<-c(1.5, 1.5, 1.5, 0)
ss_2<-c(4,4,4,0)
SD<-c(1,3,2,4)
SD_2<-c(4,2,3,0)
df.3<-data.frame(spec_names, head_width_mean, sample_size, SD, SD_2)
df.3<-data.frame(spec_names, head_width_mean, sample_size, hw_2, ss_2, SD, SD_2)
grand.mean(df.1$head_width_mean, df.1$sample_size)

df.2<-df.1 %>% remove_rownames() %>% column_to_rownames(var="spec_names")

apply(a3, 1, grand.mean)
vector<-df.2[,1]
vector<-df.2$head_width_mean
a1<-c(df.2$head_width_mean,df.2$hw_2)
a2<-c(df.2$sample_size, df.2$ss_2)
a3<-cbind(a1,a2)

grand.mean(c(df.2[1,1],df.2[1,2]))
c(df.2[1,1],df.2[1,3])

for (i in df.2) {
  print(i)
}
for (row in 1:nrow(df.2)) {
  print(i)
}

apply(df.2, 1, grand.mean(row) {
  head_width_mean <- row["head_width_mean"]
  hw_2 <- row["hw_2"]
  #do something cool here
})

data.frame(df.2[1:2], lapply(df.2[3:4], grand.mean) )
#or
cbind(wifi[1:3], lapply(wifi[4:9], grand.mean) )


grand.mean(r$a[2:3], c(r$a[4:5]))



df.2[c(1,3)]


apply(c(df.2[,1], df.2[,3]), 1, mean)

c(split(df.1,df.1$spec_names),
  split(df.1$hw_2,df.1$spec_names))[df.1$spec_names]

r<-split(df.1,df.1$spec_names)
r$a
grand.mean(r$a[2:5]
           )

apply(df.2[c(1,3)], 1, mean)

cbind(df.2[c(1,3)], df.2[c(2,4)])

df.2[1,]
grand.mean(c(1,2), c(2,2))

c(df.2[1,c(1,3)])


apply(df.2[c(1,3)])



for (i in 1:nrow(df.2[1,3])) {
  g<-print(df.2[i,])
  # do more things with the data frame...
}
g
# Construct a 5x6 matrix
X <- matrix(rnorm(30), nrow=5, ncol=6)

# Sum the values of each column with `apply()`
apply(X, 2, sum)




grand.sd   <- function(S, M, N) {sqrt(weighted.mean(S^2 + M^2, N) -
                                        weighted.mean(M, N)^2)}
fj<-c()
for (i in 1:4) {
  print(grand.sd(c(df.3[i,6], df.3[i,7]), c(df.3[i,2], df.3[i,4]), c(df.3[i,3], df.3[i,5])))
}
for (i in 1:4) {
  fj<-append(fj,print(grand.sd(c(df.3[i,6], df.3[i,7]), c(df.3[i,2], df.3[i,4]), c(df.3[i,3], df.3[i,5]))))
}

squared <- function(x) {
  m <- c()
  for(i in 1:x) {
    y <- i * i
    m <- c(m, y)
  }
  return(m)
}
x<-c(1,2)
append(x,3)

v<-c(1,1,2,2)
sd(v)






