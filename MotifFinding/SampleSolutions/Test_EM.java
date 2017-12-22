public class Test_EM {
  
  public static void main(String[] args) {

    if (args.length < 3) {
      System.err.println("\nTo execute this program, three command line arguments are required corresponding to a FASTA file name, the length of the desired motif, and the number of EM iterations. For example,\n");
      System.err.println("\tjava Test_EM modA.txt 16 500\n\n");
      return;
    }

    EM_MotifSearch_Sols m = new EM_MotifSearch_Sols(args[0], Integer.parseInt(args[1]));
    m.run_EM_multiple_times(Integer.parseInt(args[2]));
    
    // Output results
    System.out.println("\n" + m.motifInstancesToString());
    System.out.println(m.matrixToString());
    System.out.println("Consensus sequence: " + m.getConsensusSequence() + "\n");
    System.out.println("Information content: " + m.getInformationContentOfMatrix() + "\n");
  } 
  
}

