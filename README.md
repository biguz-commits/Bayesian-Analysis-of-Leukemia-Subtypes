
Section I :  Dataset Introduction

In the realm of medical research, particularly in the study of leukemia, understanding the molecular and genetic differences between various subtypes is crucial for developing targeted therapies and improving patient outcomes. Leukemia, a type of blood cancer, can be classified into different subtypes based on the French-American-British (FAB) classification system. This project aims to leverage Bayesian modeling techniques to analyze a leukemia dataset with the primary objective of comparing posterior expectations between different FAB groups.

The provided dataset includes various molecular markers and their corresponding levels in patients with different subtypes of leukemia. In the above we are just considering only 4 groups (M0, M1,M2 and M4), because of the presence of multiple Na’s in the others.

The main goals of this analysis are three: 
	Verify the posterior probability distribution of each parameter, i.e.  we are going to take the parameters Ω  , θ_0  ,θ_j  ,τ^2  and derive their posterior from the Gibbs output. The Gibbs is made necessary by the fact that the marginal distribution of the data is not available in a closed form, so we are going to sample directly from the full conditional of each parameter.
	Comparison of Posterior Expectation between groups, i.e.  we are going to take the posterior of θ_j [vectors of expectations] with mean for every p covariate, and from these vectors of expectations we are going to extract the mean value of a protein (a covariate) in a given group and compare it to the mean of the same protein of the other three groups.
	Shrinkage: what influence have the groups within each other regarding the values of the covariates, so how far is the the posterior expectation of θ_j from (y_j ) ̅.
 



Section II :  Statistical Introduction

Bayesian modelling has become an indispensable tool in modern statistics and data science due to its coherent framework for incorporating prior knowledge and updating beliefs with new data. Among the various Bayesian models, the Multivariate Normal Hierarchical Model (MVNHM) stands out for its flexibility and applicability to a wide range of complex, real-world problems. This project extends the traditional Multivariate Normal Model, by integrating it with the classical hierarchical Normal Model, enhancing its capability to handle intricate dependencies and structures inherent in high-dimensional data. 
