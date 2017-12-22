public class Test_KMeans {

    public static void main(String[] args) {

	if (args.length < 2) {
            System.out.println("\nTo execute this program, two command line arguments corresponding\nto a file name and the number of clusters are required. For example,\n");
            System.out.println("\tjava Test_KMeans data/yeast_10.txt 2\n\n");
            System.exit(0);
        }
	
	String fileName = args[0];
	int numClusters = Integer.parseInt(args[1]);
	KMeans_Clustering c = new KMeans_Clustering(fileName, numClusters);
	c.kMeans();
	System.out.println(c.toString());

    }

}

