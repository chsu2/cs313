import java.util.*;  // Needed for Scanner class and for Vector class
import java.io.*;  // Needed for File class


public class Hierarchical_Clustering extends Clustering {

	/**************************************************************
     ********************** INSTANCE VARIABLES ********************
     **************************************************************/

	//max number of clusters you can have 
    private int numClusters;
    private String fileName;
    
	//Creates an initially empty Hierarchical_Clustering.
	//from clustering:
		/** 
	     * Creates an initially empty Clustering.
	     * 
	     * A set of genes and experiments are determined from the specified String representing the name of a file.
	     * Genes and experiments are read-in from the tab-delimited file. Initially, the constructed
	     * Clustering is empty.
	     *
	     * @param   fileName   the name of a tab-delimited text file containing gene and experiment data
	     */

	public Hierarchical_Clustering(String fileName, int numClusters){

		super(fileName);
		this.numClusters = numClusters;
		this.fileName = fileName;

	}

	//Assigns each gene to its own unique cluster.
	public void initiallyAssignOneGeneToEachCluster(){

		//produce a vector with all the gene information
		for (int i = 0; i < getNumGenes(); i++){

			//create an empty cluster to add one gene to 
			Cluster initCluster = new Cluster();
			initCluster.addGene(genes.get(i));	

			//add cluster to the clusters vector
			clusters.add(initCluster);
		}
	}

	//Identifies and merges together the two closest clusters.
	public void mergeTwoClosestClusters(){

		//assign to first indices 
		int indexOfClosest1 = 0;
		int indexOfClosest2 = 1;

		//set minimum distance for a base case
		double minimumDistance = 100000000;

		int numberOfClusters = getNumClusters();

		//nested loop to compare all distances together and find the two closest clusters
		for (int i = 0; i < numberOfClusters - 1; i++){
			for (int j = i + 1; j < numberOfClusters ; j++) {

				double currDistance = clusters.get(i).getDistanceToCluster(clusters.get(j));

				if (currDistance < minimumDistance){
					//replace maximum
					minimumDistance = currDistance;
					indexOfClosest1 = i;
					indexOfClosest2 = j;
				} 
			}
		}

		//now merge the two closest clusters and remove the second cluster 
		clusters.get(indexOfClosest1).absorbCluster(clusters.get(indexOfClosest2));
		clusters.remove(indexOfClosest2);
	}

	//Performs centroid-linkage hierarchical clustering
	public void hierarchical(){

		//first assign each gene to one cluster
		initiallyAssignOneGeneToEachCluster();

		System.out.println("Merging clusters: ");
		while (getNumClusters() > numClusters){

			//for printing purposes
			int currentNumClusters = getNumClusters();
			if (currentNumClusters%10 == 0) System.out.print(currentNumClusters + " ");

			//merge clusters until numClusters is reached
			mergeTwoClosestClusters();
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


	//main method
	public static void main(String [] args) {

		Hierarchical_Clustering hClustering = new Hierarchical_Clustering(args[0], Integer.valueOf(args[1]));
		hClustering.hierarchical();
		hClustering.printClusters();
	}
}
