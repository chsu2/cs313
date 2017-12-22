import java.util.*;  // Needed for Scanner class and for Vector class
import java.io.*;  // Needed for File class


public class CAST_Clustering extends Clustering {

	/**************************************************************
     ********************** INSTANCE VARIABLES ********************
     **************************************************************/

	//max number of clusters you can have 
    private double affinityThreshold;
    private String fileName;
    private int numberOfGenes;

    
	//Creates an initially empty CAST_Clustering.
		/*
		A set of genes and experiments are determined from the specified String representing 
		the name of a file. Genes and experiments are read-in from the tab-delimited file. Initially, 
		the constructed CAST_Clustering is empty, i.e., none of the genes are assigned to clusters.
	    */

	public CAST_Clustering(String fileName, double affinityThreshold){

		super(fileName);
		this.affinityThreshold = affinityThreshold;
		this.fileName = fileName;
		numberOfGenes = getNumGenes(); //keep track of the original number of genes in the file
		
	}


	//Adds to the specified Cluster all unassigned genes that are closer to the genes in the Cluster, 
	//on average, than the affinity threshold.
	//when returns 0, reach convergence 

	public int addGenesWithHighAffinity(Cluster current){

		//get the mean of the cluster inputted
		Vector<Double> clusterMean = current.getClusterMean();

		//initialize count 
		int count = 0;

		//keep a temp array so you can loop through properly
		Vector<Gene> temp = (Vector) genes.clone();


		//loop through the genes and see if their values are within the affinityThreshold
		for (int i = 0; i < getNumGenes() ; i++) {
			
			Gene gene = temp.get(i);

			//if the gene is less than the affinity threshold and not in the current cluster
			if (gene.distanceToExpressionVector(clusterMean) < affinityThreshold){

				//add to current cluster and remove from unused genes
				current.addGene(gene);
				genes.remove(gene);
				count++;
			}
		}

		return count;
	}


	//Removes from the specified Cluster any genes that are farther from the genes in the Cluster, 
	//on average, than the affinity threshold.
	//when returns 0, reach convergence 

	public int removeGenesWithLowAffinity(Cluster current){

		//get the mean of the cluster inputted
		Vector<Double> clusterMean = current.getClusterMean();

		//initialize count
		int count = 0;

		//loop through the genes in cluster and see if their values are within the affinityThreshold
		for (int i = 0; i < current.getSizeOfCluster(); i++) {
			
			Gene gene = current.getGene(i);

			//if the gene is less than the affinity threshold
			if (gene.distanceToExpressionVector(clusterMean) > affinityThreshold){

				//remove from current cluster and add back to unused genes
				current.removeGene(gene);
				genes.add(gene);

				count++;
			}
		}

		return count;
	}


	//Performs CAST clustering.
	/*	Repeats the following until all genes are assigned to some cluster. Chooses a gene not already 
		assigned to a cluster and assigns it to a new cluster. Then, iteratively (1) adds to the new 
		cluster any unassigned genes that are closer to the genes in the new cluster, on average, than 
		the affinity threshold and (2) removes from the new cluster any genes that are farther from the 
		genes in the new cluster, on average, than the affinity threshold.
	*/
	public void cast(){

		//keep track of the genes used
		int numberOfGenesUsed = 0;

		//grab a random unused gene
		Random rand = new Random();

		while (numberOfGenesUsed != numberOfGenes){

			//get a random unused gene in the genes vector and create a new cluster with it
			Cluster initCluster = new Cluster();
			int n = rand.nextInt(getNumGenes());
			initCluster.addGene(genes.get(n));

			//remove once used for initial cluster
			genes.removeElementAt(n);

			//base cases
			int addConvergence = addGenesWithHighAffinity(initCluster);
			int removeConvergence = removeGenesWithLowAffinity(initCluster);

			//loop till convergence
			while(addConvergence != 0 || removeConvergence != 0){
				addConvergence = addGenesWithHighAffinity(initCluster);
				removeConvergence = removeGenesWithLowAffinity(initCluster);
			}

			//add the clusters to the cluster
			clusters.add(initCluster);

			//increase the number of genes used 
			numberOfGenesUsed += initCluster.getSizeOfCluster();
		}
	}





	//print the information out
	public void printClusters(){

		System.out.println("\n");

		for (int i = 0; i < getNumClusters(); i++){

			//get the number of genes in the cluster
			Cluster clusterForInfo = clusters.get(i);
			int numGenesInCluster = clusterForInfo.getSizeOfCluster();

			//print the information
			System.out.println("Cluster # " + i + " containing " + numGenesInCluster + " genes.");
			System.out.println(clusterForInfo.toString());

		}

	}

	public static void main(String [] args) {

		CAST_Clustering cClustering = new CAST_Clustering(args[0], Integer.valueOf(args[1]));

		cClustering.cast();
		cClustering.printClusters();
	}
}
