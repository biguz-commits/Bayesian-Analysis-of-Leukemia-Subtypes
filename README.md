
### **Section I: Dataset Introduction**

In medical research, particularly in the study of leukemia, understanding molecular and genetic differences between subtypes is critical for developing targeted therapies and improving patient outcomes. Leukemia, a type of blood cancer, is classified into subtypes based on the **French-American-British (FAB) classification system**. 

This project uses **Bayesian modeling techniques** to analyze a leukemia dataset, focusing on comparing **posterior expectations** between different FAB groups.

#### **Dataset Details**
- The dataset contains molecular markers and their levels in patients with different leukemia subtypes.
- For this analysis, we focus on **four groups**: **M0, M1, M2, and M4**, due to the high prevalence of missing values (NaNs) in other subtypes.

#### **Goals of the Analysis**
1. **Posterior Probability Distribution**  
   - Derive the posterior distributions of the parameters (**Ω, θ₀, θⱼ, τ²**) using Gibbs sampling.  
   - Gibbs sampling is necessary because the marginal distribution of the data is not available in closed form, so sampling from the full conditional distributions is used.

2. **Comparison of Posterior Expectations Between Groups**  
   - Analyze the posterior of **θⱼ** (vectors of expectations) for each covariate.  
   - Extract the mean values of specific proteins (covariates) within a given group and compare these means across the other three groups.

3. **Shrinkage Analysis**  
   - Investigate the influence of group membership on covariate values.  
   - Assess how far the posterior expectation of **θⱼ** is from the observed group mean (**ȳⱼ**).

---

### **Section II: Statistical Introduction**

**Bayesian modeling** has become a cornerstone of modern statistics and data science, offering a coherent framework to incorporate prior knowledge and update beliefs with new data. 

#### **Multivariate Normal Hierarchical Model (MVNHM)**
- The **MVNHM** is a powerful Bayesian model known for its flexibility and applicability in handling high-dimensional, complex data.
- This project extends the traditional **Multivariate Normal Model** by integrating it with the **classical hierarchical Normal Model**, enabling it to:
  - Handle intricate dependencies between variables.
  - Address the structure and complexity inherent in high-dimensional datasets.

By leveraging these advanced Bayesian methods, the project aims to uncover meaningful patterns and relationships within the leukemia dataset, facilitating deeper insights into its molecular and genetic landscape. 

