# Exact distribution of time T to observe T to observe "ATG" when positions 
# are independent and every letter is identically distributed

probA <- 0.1
probC <- 0.2
probG <- 0.3
probT <- 0.4

InitProb <- c(1, 0, 0, 0)

P <- matrix(0, nrow = 4, ncol = 4)
P[1,1:2] <- c(1-probA, probA)
P[2,1:3] <- c(probC+probG, probA, probT)
P[3,] <- c(probC+probT, probA, 0, probG)
P[4,4] <- 1

# let's go step by step
InitProb %*% P
InitProb %*% P %*% P
InitProb %*% P %*% P %*% P    
# Makes sense that after 3 steps we are able to reach 'ATG' state (see notes)
InitProb %*% P %*% P %*% P %*% P %*% P %*% P %*% P %*% P %*% P

result <- InitProb %*% P
probs <- c()
for (i in 1:10000) {
  probs <- append(probs, result[4])
  result <- result %*% P
}
result
plot(probs)
# This is the accumulative distribution so to make the distribution:
differences <- diff(probs)  # different between the sucesivos elements on probs (probs[2] - probs[1], probs[3]- probs[2],...)
sum(differences)

plot(differences)  # This is almost our distribution
points(dgeom(1:length(differences), prob = 8/1000), col = "red")  # We can see it is a geometric distribution

distribution <- c(0, differences) # it was missing a first 0
plot(distribution)

# Now to calculate the mean time for getting ATG
mean2 <- sum((1:length(distribution)*distribution))
mean2
mean3 <- sum(1-probs)
mean3
