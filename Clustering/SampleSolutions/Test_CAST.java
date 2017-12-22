public class Test_CAST {

    public static void main(String[] args) {

        if (args.length < 2) {
            System.out.println("\nTo execute this program, two command line arguments corresponding\nto a file name and an affinity threshold are required. For example,\n");
            System.out.println("\tjava Test_CAST data/yeast_10.txt 3.5\n\n");
            System.exit(0);
        }
        
        String fileName = args[0];
        double affinityThreshold = Double.parseDouble(args[1]);
        CAST_Clustering cc = new CAST_Clustering(fileName, affinityThreshold);
        cc.cast();
        System.out.println(cc.toString());

    }

}

