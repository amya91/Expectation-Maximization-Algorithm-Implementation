library(mvtnorm)
library(ggplot2)
library(scales)

#reading data
dat <- read.csv("Data/semeion.csv", header=FALSE)

#extracting labels from data
labels <- as.matrix(dat[,257:266])%*%as.matrix(0:9)

dat <- as.matrix(dat[,1:256])

n = nrow(dat) 
d = ncol(dat)
k = 10
q_vals = c(0,2,4,6)
n_iter = 40
AIC_qvals = vector()

#function to initialize the data to clusters using kmeans
initialize = function(mat,q) {
  
  set.seed(0)
  # initial cluster centers using kmeans
  initial_clust = kmeans(mat, k, iter.max = 20, nstart = 5)
  
  #Initailizing Gamma
  m = matrix(c(c(1:n),initial_clust$cluster),n,2)
  gamma = matrix(rep(0),1,n*k)
  gamma[,t(as.matrix((m[,1]-1)*k+m[,2]))] = 1
  gamma = t(matrix(gamma,k,n))
  sigma = array(dim = c(k,d,d))
  likelihood = rep(0, n_iter)
  assign("gamma", gamma, .GlobalEnv)
  assign("sigma", sigma, .GlobalEnv)
  assign("likelihood", likelihood, .GlobalEnv)
  m_step(mat,q)
  
}

#function for e step
e_step = function(mat){
  p = matrix(0,n,k)
  
  # calculating the Px matrix first using the dmvnorm function
  for(k in 1:k) {
    p[,k] = pi[k]*dmvnorm(mat, mu[k,], sigma[k,,], log = FALSE)
  }
  
  gamma = p/rowSums(p)
  assign("p",p, .GlobalEnv)
  assign("gamma", gamma, .GlobalEnv)
}

#function for m step
m_step = function(mat,q){
  
  mu = t(gamma)%*%mat/colSums(gamma)
  pi = colSums(gamma)/nrow(mat)
  assign("mu", mu, .GlobalEnv)
  assign("pi", pi, .GlobalEnv)
  
  # Estimating sigma_k matrix for each iteration
  for(k in 1:k) { # for each cluster
    
    # Calculate cov matrix for spectral decomposition
    mat_temp = t(mat)-mu[k,]
    cov_k = mat_temp%*%((gamma[,k])*t(mat_temp))/sum(gamma[,k])
    
    # calculating eigenvectors and eigenvalues  
    eigen_vecs = eigen(cov_k, symmetric=TRUE) # calculating eigenvectors of this matrix
    
    # Compute sigma_square
    sig_sq = (sum(eigen_vecs$values[q+1:d], na.rm = T)/(d-q)) # calculating sig_sq
    
    # checking to see if the number of principal components is zero
    if(q==0) {
      sigma[k,,] = sig_sq * diag(d)
    } 
    else {
      Vq = eigen_vecs$vectors[,1:q]
      diag_mat = diag(sqrt(eigen_vecs$values[1:q]-sig_sq))
      
      Wq = Vq %*% diag_mat
      sigma[k,,] = Wq %*% t(Wq) + (sig_sq * diag(d))
    }
    assign("sigma", sigma, .GlobalEnv)
  }
}

likelihood_fn = function() {
  
  likelihood[iter] = sum(log(rowSums(p)))
  assign("likelihood", likelihood, .GlobalEnv)
}


# declaring lists to store matrices and arrays for different Q values
gamma_qvals = vector("list", length(q_vals))
mu_qvals = vector("list", length(q_vals))
likelihood_qvals = vector("list", length(q_vals))
sigma_qvals = vector("list", length(q_vals))


#-----running EM algorithm for 40 iterations-----#
for (i in seq_along(q_vals)) {
  q = q_vals[i]
  print(c('Running EM Algorithm for q = ',q))
  initialize(dat,q)
  
  for (iter in 1:n_iter) {
    m_step(dat,q)
    e_step(dat)
    likelihood_fn()
    print(paste("Log-likelihood value for iteration ",iter,' is',round(likelihood[iter])))
  }
  likelihood_qvals[[i]] = likelihood
  gamma_qvals[[i]] = gamma
  mu_qvals[[i]] = mu
  sigma_qvals[[i]] = sigma
}

#------Plotting Log-likelihood vs Iteration plots for all q-------#

for (i in seq_along(q_vals)) {
  
  plot_data = as.data.frame(cbind(likelihood_qvals[[i]],c(1:n_iter)))
  colnames(plot_data) = c('log-likelihood','iteration-number')
  
  #Saving plots to disk
  png(paste0('Log-likelihood plot for ',q_vals[i],' principal components.png'), width = 900, height = 600, res = 120)
  print(ggplot(plot_data, aes(y=plot_data$`log-likelihood`, x=plot_data$`iteration-number`)) +
          geom_point(color="blue") +
          labs(title = paste("Log-likelihood v. Iteration for ",q_vals[i],' principal components'),
               x = "Iteration",
               y = "Log-likelihood") +
          scale_x_discrete(limits=c(1:n_iter)))
  dev.off()
}

#-------Calculating AIC and finding best q with least AIC--------#

for (i in seq_along(q_vals)) {
  #calculate AIC
  AIC = -2*tail(likelihood_qvals[[i]],1) + 2*(d*q_vals[i] + 1 - (q_vals[i]*(q_vals[i]-1)/2))
  AIC_qvals[i] = AIC
  print(paste('For',q,'Principal components, the value of AIC is',AIC))
}
best_q = q_vals[which(AIC_qvals==min(AIC_qvals))]
best_AIC = min(AIC_qvals)
print(paste('Principal components =',q_vals[which(AIC_qvals==min(AIC_qvals))],'minimizes AIC'))

#-----Visualization of Clusters-----#
dev.new(width=7,height=3.5)
par(mai=c(0.05,0.05,0.05,0.05),mfrow=c(10,6))

for(z in 1:10){
  image(t(matrix(mu_qvals[which(AIC_qvals==min(AIC_qvals))][[1]][z,], byrow=TRUE,16,16)[16:1,]),col=gray(0:1),axes=FALSE)
  box()
  for(j in 1:5){
    tempX = rmvnorm(1, mean <- mu_qvals[which(AIC_qvals==min(AIC_qvals))][[1]][z,], sigma_qvals[which(AIC_qvals==min(AIC_qvals))][[1]][z,,])
    image(t(matrix(tempX, byrow=TRUE,16,16)[16:1,]),col=gray(0:1),axes=FALSE)
    box()
  }
}

#-----Accuracy Assesment-----#

# Find most common digits for each cluster
final_label = rep(0,n)
for(i in 1:n) {
  final_label[i] = which.max(gamma_qvals[4][[1]][i,])
}

# Cluster Splittig
groups = split(labels, final_label)
perc = lapply(groups, function(group){
  return(sort(table(group), decreasing=TRUE)[1])
})

#Mis-classification Matrix
acc = matrix(0,4,k)

for(k in 1:k) {
  acc[1,k] = round(as.integer(names(perc[[k]])))#Most common digit
  acc[2,k] = round(as.integer(perc[[k]][[1]]))#No of occurences
  acc[3,k] = round(as.integer(length(groups[[k]])))#Total Observations
  acc[4,k] = (as.numeric(1 - (acc[2,k] / acc[3,k])))#Mis-classification rate
}
o_mis_rate = 1 - sum(acc[2,]) / sum(acc[3,])
acc[4,] = percent(acc[4,])
rownames(acc) = c('Most Common Digit','# of Occurances','Total Observations','Mis-classification rate')
print(acc)
print(paste('The overall mis-classification rate is',percent(o_mis_rate)))

