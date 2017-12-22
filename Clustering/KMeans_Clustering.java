import java.util.*;  // Needed for Scanner class and for Vector class
import java.io.*;  // Needed for File class


public class KMeans_Clustering extends Clustering {

	/**************************************************************
     ********************** INSTANCE VARIABLES ********************
     **************************************************************/

	//max number of clusters you can have 
    private int numClusters;
    private String fileName;
    
	/***************************************************************

	Creates an initially empty KMeans_Clustering.

	A set of genes and experiments are determined from the specified String representing the name of a 
	file. Genes and experiments are read-in from the tab-delimited file. Initially, the constructed 
	KMeans_Clustering consists of the specified number of clusters, but each cluster has not yet been 
	assigned any genes.

	***************************************************************/

	public KMeans_Clustering(String fileName, int numClusters){

		super(fileName);
		this.numClusters = numClusters;
		this.fileName = fileName;

		//create the specified number of clusters, but don't assign any genes
		for (int i = 0; i < numClusters; i++){

			//create an empty cluster and add to clusters vector 
			Cluster initCluster = new Cluster();
			clusters.add(initCluster);

		}

		// System.out.println("number of empty clusters in clustering: " + getNumClusters());

	}

	//Initializes each cluster so that each cluster contains zero genes
	public void initializeAllClusters(){

		//loop through the clusters and initialize to 0 genes
		for (int i = 0; i < numClusters; i++){

			//initialize all clusters to 0 genes
			clusters.get(i).initialize();

		}
		
	}

	//Return a collection of the mean (average) expression vectors for all of the clusters.
	public Vector<Vector<Double>> getMeansOfAllClusters(){
		Vector<Vector<Double>> meansOfAllClusters = new Vector<Vector<Double>>();

		//loop through the collection of clusters
		for (int i = 0; i < numClusters; i++){

			//get the mean of the current cluster and add it to meansOfAllClusters
			Vector<Double> clusterMean = clusters.get(i).getClusterMean();
			meansOfAllClusters.add(clusterMean);

		}

		return meansOfAllClusters;
	}

	/*
	If any clusters are empty (contain zero genes), then genes are moved from clusters containing multiple genes.

	Assumes all genes have been assigned to a cluster. For any empty clusters, a cluster with mulitple genes is 
	found randomly and one gene is removed from the multiple gene cluster and added to the empty cluster. When 
	the method completes, no cluster contains zero genes.
	*/
	public void populateEmptyClusters(){

		//number of clusters never changes so can use the numClusters variable
		//loop through clusters
		for (int i = 0; i < numClusters; i++){

			//look at cluster
			Cluster possEmptyCluster = clusters.get(i);

			//check the size of the cluster. if 0, execute another loop
			if (possEmptyCluster.getSizeOfCluster() == 0){
				// System.out.println("the cluster is empty!");

				//set a boolean to find a random cluster with a size > 1
				boolean b = false;
				Random rand = new Random(); //to find random cluster
				while (b == false){

					//generate random index
					int n = rand.nextInt(numClusters);

					//make sure you don't get the index of the empty cluster
					if (n != i){

						Cluster currCluster = clusters.get(n);
						int clusterSize = currCluster.getSizeOfCluster();
						//if the cluster found has more than one element then take a gene and add to empty cluster
						if (clusterSize > 1){

							// System.out.println("relacing");


							Gene geneToSwitch = currCluster.getGene(clusterSize - 1); //grab last gene in vector so it doesn't have to reindex

							currCluster.removeGene(geneToSwitch);
							possEmptyCluster.addGene(geneToSwitch);
							b = true;

						}	
					}
				}
			}
		}

		//loop through to see find all empty clusters 
		//if so, then take a random find cluster with genes (since multiple it needs to be more than 1) while loop
		//then add gene to empty cluster 


	}

	/*
	Assigns each gene to a random cluster.

	Each gene is randomly assigned to one of the clusters. When the method completes, no cluster should contain 
	zero genes.
	*/
	public void randomlyAssignGenesToClusters(){

		Random rand = new Random();

		for (int i = 0; i < getNumGenes(); i++) {

			//find a random cluster then add gene to cluster
			int n = rand.nextInt(numClusters);

			clusters.get(n).addGene(genes.get(i));

		}

		//call to make sure no clusters are empty
		populateEmptyClusters();

		// System.out.println(clusters.toString());

		//loop through all the genes and assign to random cluster with random number
		//call populateEmptyClusters

	}

	/*
	Assigns each gene to the cluster whose mean expression vector is closest to the gene.

	The parameter means corresponds to a collection of Vectors (analagous to a 2D array). Each entry of 
	means corresponds to the set of average expression values for a particular cluster. The method returns 
	true if the cluster assignments are not improving, i.e., we have achieved a locally optimal clustering. 
	The method returns false if the cluster assignments are an improvement over previous cluster assignments, 
	i.e., better clusters than previously have been found. The measure used for comparing one clustering to 
	another clustering is the sum of the distance of each gene from its cluster's mean.
	*/
	public boolean assignGenesToClusters(Vector<Vector<Double>> means){

		//reinitialize the clusters
		initializeAllClusters();

		//loop through genes 
		for (int i = 0; i < getNumGenes(); i++) {

			//grab the gene and get the distance from it to all the means
			Gene gene = genes.get(i);

			// System.out.println(means.get(0));
			// System.out.println(gene.getExpressionVector());

			//set a base
			double minDistance = gene.distanceToExpressionVector(means.get(0));
			int minIndex = 0;
			//find cluster it best fits in
			for (int j = 1; j < means.size(); j++) {

				double distance = gene.distanceToExpressionVector(means.get(j));
				if (distance < minDistance){

					//if you find a new min replace old min and index
					minDistance = distance;
					minIndex = j;
				}
				
			}

			//assign the gene to its new cluster
			clusters.get(minIndex).addGene(gene);

			
		}

		//recalculate cluster mean
		Vector<Vector<Double>> newClusterMeans = getMeansOfAllClusters();

		if (newClusterMeans.equals(means)) return true;

		return false;

	}

	/*
	Performs k-means clustering.
	Initially, genes are randomly assigned to clusters. Then, iteratively, the clustering is improved until a 
	local optima is reached (until clusterings are no longer improving). In each iteration, first the mean 
	(average) expression vector is calculated for each cluster. Second, each gene is assigned to the cluster whose 
	mean (average) expression vector is closest to the gene's expression vector.
	*/
	public void kMeans(){

		//assign to random clusters
		randomlyAssignGenesToClusters();
		Vector<Vector<Double>> means = new Vector<Vector<Double>>();

		boolean loopBool = false;

		while(loopBool == false){
			means = getMeansOfAllClusters();
			loopBool = assignGenesToClusters(means);

		}

	}

	//print the information out
	public void printClusters(){

		System.out.println("\n");

		for (int i = 0; i < numClusters; i++){

			//get the number of genes in the cluster
			Cluster clusterForInfo = clusters.get(i);
			int numGenesInCluster = clusterForInfo.getSizeOfCluster();

			System.out.println("Cluster # " + i + " containing " + numGenesInCluster + " genes.");
			System.out.println(clusterForInfo.toString());

		}

	}

	public static void main(String [] args) {

		KMeans_Clustering kClustering = new KMeans_Clustering(args[0], Integer.valueOf(args[1]));

		// kClustering.randomlyAssignGenesToClusters();

		// kClustering.assignGenesToClusters(kClustering.getMeansOfAllClusters());
		
		kClustering.kMeans();
		kClustering.printClusters();


	}
}
