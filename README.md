# PCA_cluster
Clusters positions along a biological sequence (genome or proteins) according to the values of the observables provided in an input matrix. Applications: Chromatin states, protein domain decomposition, dynamical domains

Program PCA_cluster, author Ugo Bastolla ubastolla@cbm.csic.es

Please cite:
Sequeira-Mendes, J., Aragüez, I., Peiró, R., Mendez-Giraldez, R., Zhang, X., Jacobsen, S. E., Bastolla, U., & Gutierrez, C. (2014). The Functional Topography of the Arabidopsis Genome Is Organized in a Reduced Number of Linear Motifs of Chromatin States. The Plant cell, 26(6), 23512366. https://doi.org/10.1105/tpc.114.124578

DESCRIPTION

PCA_cluster assigns hidden states to positions along a biological sequence as a function of observables along the sequence. Possible applications are chromatin states (possible observables are histone modifications, DNA methylation, GC content), structural domains of a protein with known structure (observables are eigenvectors of the contact matrix), dynamical domains (possible observables are eigenvectors of dynamical couplings computed with the normal modes of the elastic network model), etc.

Given a matrix where each line represents a position ordered along a biological sequence and each column represent an observable, the program performs the following steps.

1) Principal Component Analysis of the O observables, which are substituted by n (input parameter) principal components. It is advisable to experiment different values of n. Smaller n accelerates the computation and reduces the risk of overfitting the HMM model, but it may eliminate useful information. In my experience, n=3 or n=4 (default) is sufficient for determining robust chromatin states.

Output file: <file>.pca (<file> is the name of the input file). It contains the pairwise correlations between the O input observables, the fraction of variance explained by each PC (individual and cumulative) and the contribution of each of the O observables (one per row) to the n PCs (one per column).

Input parameter: 
   -n <number of PCA used for clustering> (default 4)

2) Infer the HMM with s hidden states between s=1 and smax. For each s, the program determines the HMM parameters of each of the s states that depend on the n PCs that are considered in the model (number of parameters: s times n average values plus n(n+1)/2 elements of the covariance matrix, plus s weights and s(s-1) transition probabilities). The inference of parameters is done iteratively until convergence starting from two randomly drawn initial parameter sets and the similarity between the final parameters is evaluated, if it is above a threshold the computation is stopped, the parameters that yield maximum likelihood are selected.

After inferring the parameters, the program assigns each sample to the cluster that maximizes the posterior probability (out of s possibe ones), it computes the similarity between the profiles of each cluster (average values and standard errors of the O observables) and it joins two clusters if their similarity in terms is above a threshold. Therefore, the final number of clusters may be smaller than s. Finally, the program assesses the decomposition by computing its BIC and AIC scores (likelihoods penalized by the number of parameters).

Output files:
In the file names below, <n> indicates the number of PCs (default 4)

<file>_n<n>_HMM.log
log file with summary information

<file>_n<n>_HMM_bic.dat
For each value of s, it reports the final scores: Likelihood (the larger the better), BIC, AIC (the smaller the better) with respect to s=1, cluster score (the smaller the better), mean square difference between the average PCs of all pairs of clusters (the larger the better) and maximum cosine between any two different clusters (the smaller the better).

<file>_n<n>_HMM_clusters.txt
The most important file: for each position in the sequence it reports the most likely state (=cluster) index, ranked according to the average value of the first principal component so that cluster 0 has the largest average value of PC0, and the and the posterior probability to belong to the state computed with the HMM.
Each pair of columns represents a value of the number of clusters. They are printed in decreasing order from nclust_max to nclust_min. If the optimal number of clusters is different from the maximum one, the first pair of columns always represents the optimal number of clusters.

<file>_n<n>_HMM_cluster_profiles.dat
For each studied number of clusters (=states) and each cluster, it reports the average value and the standard error of the O observables and the n PCs and also the number and percentage of elements of each cluster.

<file>_n<n>_HMM_cluster_similarity.dat
For each studied number of clusters (=states) and each pair of clusters it reports the similarity between the profiles mentioned above, measured as cosine between the n dimensional vectors of the average PCs and cosine between the n dimensional vectors of the binarized average PCs. Pairs of clusters with similarity larger than a threshold are merged.

The above files contain results for each value of s from 2 until the maximal one that maximizes the user defined score or until the value smax.

Input parameters:

   -smin <Minimum number of clusters> default 1
   -smax <Maximum number of states> default 12
   -stop <stopping criterion> Allowed: aic, bic, sep, score
   -sim <Maximum similarity of two clusters>, default 0.900
   -diff <Minimum number of different variables in two clusters> (default 0)
   -noorder: Do not take into account transition probability (no HMM)

The stopping criteria are. minimum AIC, minimum BIC, minimum cluster score, maximum separation between clusters.
If -noorder is set, the transition probabilities between different clusters are not considered and the model is a Gaussian mixture model instead of a HMM

COMPILATION

Copy PCA_cluster.zip in a new directory, then run
>unzip PCA_cluster.zip
>make
>mv PCA_cluster <path>
where <path> is one of the directories of your local path
>rehash

USAGE

PCA_cluster <input_file> (mandatory)
Options:
   -h help
   -n <number of PCA used for clustering> (default 4)
   -smin <Minimum number of states> default 1
   -smax <Maximum number of states> default 12
   -stop <stopping criterion> Allowed: aic, bic, sep, score
   -noorder: Do not take into account positional correlations (HMM)
   -sim <Maximum similarity of two clusters>, default 0.900
   -diff <Minimum number of different variables in two clusters> (default 0)
   -split <file> file with N lines (as samples) to be splitted
   -print  Print n principal components of all elements
   -matrix Print matrix with selected observables used for PCA

FORMAT of Input file

Matrix with N rows (samples) and M columns (variables).

The first line must contain variable names: ### name1 name2 name3...

The second line selects variables used for PCA (1=selected, 0=omitted): ### 1 0 1 ...

In this way, the file can contain additional information such as, for chromatin states, chromosome and genomic coordinate of each position.

The github repository contains the example file Properties_genome_w250_norm.dat.gz
It was obtained dividing the genome in bins of width 500bp, averaging the Z scores of each of the M genomic and epigenomic observables in the bin, and smoothing in 5 neighboring bins (d=-2,-1,0,1,2) with exponentially decaying weights w(d)=exp(-|d|*A).
