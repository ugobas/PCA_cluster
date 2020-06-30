float Average_linkage_coord(float **Clust_sim, int *cluster, int *Nnew,
			    int N, float **x, int *Nele, int Nvar, float Thr);
float Cosine(float *v1, float *v2, int N);
int Silouhette(int *cluster, int *numclus, int nclust,
	       float **x_VarSam, int Nsam, int Nvar);
