import numpy as np
import matplotlib.pyplot as plt


size = 1000
biases = np.linspace(0,1,size)
priors = np.ones(size)/len(biases)

# generate random p
p = np.random.rand()
print("Coin bias, p(X=1) = {}".format(p))

# define number of iterations
num_iterations = 1000


# iterate through 
for curr_iter in range(num_iterations):

    plt.plot(biases, priors, color = "blue", alpha = 0.1) #alpha = (curr_iter/num_iterations))

    print("iter {}".format(curr_iter))
    c = np.random.choice([0, 1], p = [1-p, p]) # flip coin
    evidence = np.sum(np.multiply(biases, priors)) # estimated probaility of heads

    for i in range(len(biases)):
        # prosterior = prior * likelihood / evidence
        if c == 1: # prior * bias / evidence
            priors[i] = priors[i] * biases[i] / evidence
        else: # prior * 1-bias / evidence
            priors[i] = priors[i] * (1-biases[i]) / evidence

plt.plot(biases, priors, color = "red") #alpha = (curr_iter/num_iterations))
plt.title("Estimating Coin Bias of {}".format(p))
plt.xlabel("heads bias")
plt.ylabel("P(heads bias | flips)")

plt.show()

### APPLYING TO GENOMES
# Dependency is a bias of 1?
# Coin is flipped for each co-appearence of two genes
# 


