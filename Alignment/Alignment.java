import java.io.File;
import java.io.FileNotFoundException;
import java.io.*;
import java.util.Scanner;
import java.io.PrintWriter;
import java.util.Hashtable;
import java.util.Vector;
import java.util.*;
import java.awt.Point;
import java.lang.*;

/*************************************************************
  * An instance of the Alignment class represents an 
  * optimal pairwise alignment for two genomic sequences.
  *************************************************************/
public class Alignment {
  
  /*********************************************************
    ****************** INSTANCE VARIABLES *******************
    *********************************************************/
  
  private String seq1;
  private String seq2;
  private int[][] alignmentTable;
  private Point[][] backtrackTable;
  private int optimalScore;
  private String alignment;
  private double pValue;
  private Point localOptimalLocation;
  private Hashtable<String, Integer> matrixScoring = new Hashtable<String, Integer>();
  private Vector<Integer> alignScores;
  
  // Parameters
  private int matchScore;
  private int mismatchScore;
  private int linearGapScore;
  private int alphaGapScore;
  private int betaGapScore;
  
  // Options
  private boolean useLocalAlignment;  // As opposed to global alignment, which is default
  private boolean useMatrixScoring;  // As opposed to fixed match/mismatch scoring, which is default
  private boolean useAffineGapScoring;  // As opposed to linear gap scoring, which is default
  private boolean useFastAlignment;  // As opposed to quadratic time alignment, which is default
  private int numGaps;  // When performing a "fast" alignment, consider at least this many gaps
  
  
  
  /*********************************************************
    ****************** CONSTRUCTORS *************************
    *********************************************************/
  
  /**
   * Creates an Alignment from the genomic sequences found 
   * in the two specified FASTA files.
   *
   * @param   file1   a <code>File</code> object referring to a FASTA file containing a genomic sequence
   * @param   file2   a <code>File</code> object referring to a FASTA file containing a genomic sequence
   */
  public Alignment(File file1, File file2) {
    this(Alignment.readSequenceFromFile(file1), Alignment.readSequenceFromFile(file2));
  }
  
  /**
   * Creates an Alignment from the two specified genomic sequences.
   *
   * @param   sequence1   a genomic sequence to be aligned
   * @param   sequence2   a genomic sequence to be aligned
   */
  public Alignment(String sequence1, String sequence2) {
    seq1 = sequence1;
    seq2 = sequence2;
    if (seq1.charAt(seq1.length()-1) == '*') seq1 = seq1.substring(0, seq1.length()-1);
    if (seq2.charAt(seq2.length()-1) == '*') seq2 = seq2.substring(0, seq2.length()-1);
    setGlobalAlignment();  // Global alignment is default
    setFixedScoring(5, -4);  // Fixed match/mismatch scoring is default
    setLinearGaps(-6);  // Linear gap scoring is default
    useFastAlignment = false;  // Quadratic alignment is default
    alignmentTable = new int[seq1.length() + 1][seq2.length() + 1];  // Create alignment table
    backtrackTable = new Point[seq1.length() + 1][seq2.length() + 1];  // Create backtrack table
  }
  
  
  
  /*********************************************************
    ****************** PUBLIC INSTANCE METHODS **************
    *********************************************************/
  
  /**
   * Returns the first of two genomic sequences in this <code>Alignment</code>.
   *
   * @return   a <code>String</code> corresponding to the first of two genomic sequences in the pairwise alignment
   */
  public String sequence1() {
    return seq1;
  }
  
  /**
   * Returns the second of two genomic sequences in this <code>Alignment</code>.
   *
   * @return   a <code>String</code> corresponding to the second of two genomic sequences in the pairwise alignment
   */
  public String sequence2() {
    return seq2;
  }
  /**
   * Computes the optimal pairwise alignment of two genomic sequences.
   * <p>
   * Either the optimal global or optimal local pairwise alignment is computed.
   * Computation of the optimal pairwise alignment includes determining the
   * optimal pairwise alignment score as well as the corresponding alignment.
   */
  public void computeAlignment() {
    if (!useLocalAlignment) {  // Global alignment
      if (useFastAlignment) globalAlignment_FAST();  // Linear time alignment, with limited gaps considered
      else globalAlignment();  // Quadratic time alignment, with any number of gaps considered
      alignment = globalBacktrack();
    } else {  // Local alignment
      localAlignment();
      alignment = localBacktrack();
    }
  }
  
  /**
   * Returns the optimal pairwise alignment score for this <code>Alignment</code>.
   *
   * @return   the optimal pairwise alignment score
   */
  public int getAlignmentScore() {
    return optimalScore;
  }
  
  /**
   * Returns a <code>String</code> representation of the alignment table generated during
   * alignment computation. This method is used primarily for debugging and is only useful
   * for <em>small</em> alignment tables, i.e., when aligning short sequences.
   *
   * @return   a <code>String</code> representation of the alignment table
   */
  public String alignmentTableToString() {
    return tableToString(alignmentTable);
  }
  
  /**
   * Returns a <code>String</code> representation of the backtrack table generated during
   * alignment computation. This method is used primarily for debugging and is only useful
   * for <em>small</em> backtrack tables, i.e., when aligning short sequences.
   *
   * @return   a <code>String</code> representation of the backtrack table
   */
  public String backtrackTableToString() {
    return tableToString(backtrackTable);
  }
  
  /**
   * Returns the optimal pairwise alignment.
   *
   * @return   a <code>String</code> representing the optimal pairwise alignment
   */
  public String getAlignment() {
    return alignment;
  }
  
  /**
   * Returns the p-value of this <code>Alignment</code>.
   * <p>
   * The p-value of this <code>Alignment</code> is the probability (between 0.0 and 1.0) that the 
   * optimal pairwise alignment score of two random sequences is greater than or equal to the
   * optimal pairwise alignment score for this <code>Alignment</code>. A p-value close to 1.0 suggests
   * that an alignment was likely to have occurred merely by chance. A p-value close to 0.0 suggests
   * that an alignment was unlikely to have occurred by chance. In this case (especially if the p-value
   * is less than about 0.05), the alignment is <em>significant</em> and the sequences are deemed
   * similar.
   *
   * @return   the p-value of the optimal pairwise alignment
   */
  public double getPValue() {
    return pValue;
  }
  
  /**
   * Indicate that a <em>global</em> pairwise alignment should be performed.
   */
  public void setGlobalAlignment() {
    this.useLocalAlignment = false;
  }
  
  /**
   * Indicate that a <em>local</em> pairwise alignment should be performed.
   */
  public void setLocalAlignment() {
    this.useLocalAlignment = true;
  }
  
  /**
   * Indicate that a fast linear-time pairwise alignment should be performed.
   * <p>
   * For two sequences of length n, rather than use an O(n^2) algorithm that identifies 
   * the optimal alignment with any number of gaps, a <em>FAST</em> alignment runs in 
   * O(numGaps*n) time where <em>numGaps</em> is the number of gaps considered. 
   * This option is only available for <em>global</em> alignments.
   * 
   * @param   numGaps   at least this many gaps are considered when computing the optimal alignment
   */
  public void setFastAlignment(int numGaps) {
    this.useFastAlignment = true;
    this.numGaps = numGaps;
  }
  
  /**
   * Score gaps in this <code>Alignment</code> using a <em>linear</em> model.
   * <p>
   * With a linear model for scoring gaps, every gap is penalized the
   * same amount as specified by the <code>linearGapScore</code> parameter.
   *
   * @param   linearGapScore   the (negative) contribution to the alignment score of each gap
   */
  public void setLinearGaps(int linearGapScore) {
    this.linearGapScore = linearGapScore;
    this.useAffineGapScoring = false;
  }
  
  /**
   * Score gaps in this <code>Alignment</code> using an <em>affine</em> model.
   * <p>
   * With an affine model for scoring gaps, the first gap in a sequence of 
   * consecutive gaps is penalized by the <code>alphaGapScore</code> parameter 
   * whereas subsequent gaps in a sequence of consecutive gaps are penalized
   * by the <code>betaGapScore</code> parameter.
   * <p>
   * Affine gap scoring is meant to model empirical biological evidence that the <em>existence</em>
   * of a gap is more significant than the <em>length</em> of the gap. It is expensive, biologically,
   * to add to or splice from a genomic sequence, but the length of the addition or deletion is less important.
   * 
   * @param   alphaGapScore   the (negative) contribution to the alignment score of initiating each sequence of gaps
   * @param   betaGapScore   the (negative) contribution to the alignment score of extending each sequence of gaps
   */
  public void setAffineGaps(int alphaGapScore, int betaGapScore) {
    this.alphaGapScore = alphaGapScore;
    this.betaGapScore = betaGapScore;
    this.useAffineGapScoring = true;
  }
  
  /**
   * When aligning two characters, one from each genomic sequence, if the two characters
   * are identical then the alignment score of the two characters should be the
   * <code>match</code> score. If the two characters differ then the alignment score of
   * the two characters should be the <code>mismatch</code> score.
   * 
   * @param   matchScore   the (positive) contribution to the alignment score of aligning two identical characters
   * @param   mismatchScore   the (negative) contribution to the alignment score of aligning two different characters
   */
  public void setFixedScoring(int matchScore, int mismatchScore) {
    this.matchScore = matchScore;
    this.mismatchScore = mismatchScore;
    this.useMatrixScoring = false;
  }
  
  /**
   * When aligning two characters, one from each genomic sequence, the alignment score of the
   * two characters should be determined from a matrix of scores found in a file with the 
   * specified name.
   * <p>
   * In <em>fixed</em> scoring, all pairs of identical characters (e.g., G|G, C|C, T|T) are scored 
   * the same and all pairs of different characters (e.g., G|C, G|T, C|T) are scored the same. 
   * However, with genomic sequences, not all pairs of characters are equally similar or dissimilar. 
   * In <em>matrix</em> scoring, different pairs of identical characters (e.g., G|G, C|C, T|T) may
   * be scored differently and different pairs of mismatching characters (e.g., G|C, G|T, C|T) may
   * be scored differently. For example, adenine (A) and guanine (G) are both purines whereas 
   * cytosine (C) and thymine (T) are both pyrimidines. Since adenine is more similar to guanine
   * than to thymine, an adenine aligned with a guanine (A|G) should not penalize an alignment as much
   * as an adenine aligned with a thymine (A|T). Analogously for protein sequences, two different
   * hydrophobic amino acids aligned together might not penalize an alignment as much as a hydrophobic
   * amino acid aligned with a hydrophilic amino acid.
   * <p>
   * In <em>matrix</em> scoring,
   * the alignment score of every possible pair of characters is specified in a matrix that must
   * be read in from a file. For DNA sequences, since there are 4 characters in the DNA alphabet,
   * there are 16 possible pairs of characters and the matrix contains 16 entries. For protein
   * sequences, since there are 20 characters in the protein alphabet, there are 400 possible pairs
   * of characters and the matrix contains 400 entries. 
   *
   * @param   fileName   the name of a file containing a matrix of alignment scores for all pairs of characters
   */
  public void setMatrixScoring(String fileName) {
    
    try{
      BufferedReader br = new BufferedReader(new FileReader(fileName));
      // read the first two lines separately since it's file name and blank line
      br.readLine();
      br.readLine();
      //store the order of the characters
      String line = br.readLine();
      line = line.replaceAll("\\s","");
      String[] chars = line.split(""); //where chars[0] is a space
      // read the rest of the file into one string
      // for (int i = 0; i < chars.length; i++) System.out.println(chars[i]);

      while ((line = br.readLine()) != null){ 
        line = line.replaceAll("   "," "); //replace triple space with single space
        line = line.replaceAll("  "," "); //replace double space with single space
        int val=0;
        String key;
        String[] split_line = line.split(" ");
        int j = 0; //index of row character
        // System.out.println("length of split_line: " + split_line.length);

        for (int i = 1; i<split_line.length; i++){ //i is index of row character

          // System.out.println("i: " + i + "\tj: " + j);
          key = split_line[0] + chars[j];
          val = Integer.parseInt(split_line[i]);
          matrixScoring.put(key, val);
          j++;
        }
      }
      br.close();
    }
    catch (FileNotFoundException e){
      e.printStackTrace();
    }
    catch (IOException e){
      e.printStackTrace();
    }
    
    this.useMatrixScoring = true;
    
  }
  
  /**
   * Estimates the p-value of this <code>Alignment</code>.
   * <p>
   * The p-value of this <code>Alignment</code> is the probability (between 0.0 and 1.0) that the 
   * optimal pairwise alignment score of two random sequences is greater than or equal to the
   * optimal pairwise alignment score for this <code>Alignment</code>. A p-value close to 1.0 suggests
   * that an alignment was likely to have occurred merely by chance. A p-value close to 0.0 suggests
   * that an alignment was unlikely to have occurred by chance. In this case (especially if the p-value
   * is less than about 0.05), the alignment is <em>significant</em> and the sequences are deemed
   * similar.
   * <p>
   * A p-value for an alignment of two sequences can be estimated as follows. Randomly generate 1000
   * pairs of sequences by randomly permuting the original
   * two sequences. For each of the 1000 pairs of random sequences, determine the optimal
   * pairwise alignment score. These 1000 scores approximate an extreme value distribution.
   * Calculate the <em>mean</em> and <em>standard_deviation</em> of the 1000 scores. 
   * The <em>mean</em> and <em>standard_deviation</em>
   * can be used to calculate two parameters, <em>mu</em> and <em>beta</em>, representing the extreme 
   * value distribution. <em>mu</em> can be calculated as <em>mean</em> - (0.5772*<em>beta</em>).
   * <em>beta</em> can be calculated as <em>standard_deviation</em> * &radic;6 / &pi;.
   * Finally, the p-value can be calculated as 
   * 1.0 - <em>e</em>^(-<em>e</em>^(-(<em>x</em>-<em>mu</em>)/<em>beta</em>))
   * where <em>x</em> is the optimal pairwise alignment score of the original pair of sequences.
   */
  public void calculatePValue() {
    if (useFastAlignment && (Math.abs(seq1.length()-seq2.length()) > numGaps)) this.pValue=Double.NaN; 
    else{ //actually calculate it
      alignScores = new Vector<Integer>();
      //randomly permute two sequences 1000 times
      for (int i=0; i<1000; i++){
        String seq1p = SequenceOps.randomPermutation(seq1);
        String seq2p = SequenceOps.randomPermutation(seq2);
        Alignment temp = new Alignment(seq1p, seq2p);
        //calculate optimal pairwise alignment for each pair and store:
        temp.computeAlignment();
        alignScores.add(temp.getAlignmentScore());
      }
      //calculate mu and beta with mean and std
      double beta = getVectorStandardDeviation(alignScores) * Math.sqrt(6) / Math.PI;
      double mu = getVectorMean(alignScores) - (0.5772*beta);
      double currAlignmentScore = getAlignmentScore();
      double inside = (-1.0)*((currAlignmentScore-mu)/beta);
      this.pValue = 1.0 - Math.exp((-1.0)*Math.exp(inside));
    }
  }
  
  /**
   * Returns a <code>String</code> representation of this <code>Alignment</code>.
   *
   * @return   a <code>String</code> representation of this <code>Alignment</code>
   */
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("\n***********************************************\n");
    if (!useLocalAlignment && !useFastAlignment) sb.append("Alignment:\t" + "Global" + "\n");
    else if (!useLocalAlignment && useFastAlignment) sb.append("Alignment:\t" + "Global FAST" + "\n");
    else sb.append("Alignment:\t" + "Local" + "\n");
    if (!useMatrixScoring) sb.append("Scoring:\t" + "Fixed (match=" + matchScore + ", mismatch=" + mismatchScore + ")" + "\n");
    else sb.append("Scoring:\t" + "Matrix" + "\n");
    if (!useAffineGapScoring) sb.append("Gap scoring:\t" + "Linear (score=" + linearGapScore + ")" + "\n");
    else sb.append("Gap scoring:\t" + "Affine (alpha=" + alphaGapScore + ", beta=" + betaGapScore + ")" + "\n");
    sb.append("***********************************************\n");
    
    sb.append("\nSequence1: " + seq1 + "\n" + "Sequence2: " + seq2 + "\n");
    sb.append("\n" + getAlignment());
    if (!useLocalAlignment)  // Global alignment
      sb.append("Optimal global alignment score: " + getAlignmentScore() + "\n");
    else  // Local alignment
      sb.append("Optimal local alignment score: " + getAlignmentScore() + "\n");
    sb.append("p-value: " + getPValue() + "\n\n");
    return sb.toString();
  }
  
  
  
  /*********************************************************
    ****************** PUBLIC CLASS METHODS *****************
    *********************************************************/
  
  /**
   * Outputs a histogram of optimal pairwise alignment scores to a file.
   * <p>
   * A histogram with 101 bins is created from the set of alignment scores
   * stored in the Vector <code>v</code>. The histogram indicates the number
   * of alignment scores corresponding to each of the 101 bins. The histogram
   * is output to the specified file <code>fileName</code>. The first column
   * of the output file represents the x-axis of the histogram, i.e., the
   * 101 bins corresponding to 101 possible alignment scores. The second
   * column of the output file indicates the number of alignment scores corresponding
   * to each bin. The third column is a normalized version of the second
   * column, i.e., each entry in the third column is the corresponding
   * entry in the second column divided by the total number of alignment
   * scores. The fourth column represents a mathematical function - an
   * extreme value distribution - that approximates the third column. The
   * extreme value distribution is determined from the mean and standard
   * deviation of the set of alignment scores in <code>v</code>.
   *
   * @param   fileName   the name of a file to which the histogram will be output
   * @param   v          a <code>Vector</code> of optimal pairwise alignment scores
   */
  public static void outputHistogramOfRandomAlignmentScores(String fileName, Vector<Integer> v) {
    // Generate histogram
    int binsInHistogram = 100;
    double mean = getVectorMean(v);
    double std_dev = getVectorStandardDeviation(v);
    int min = (int)(mean) - binsInHistogram/2;
    int max = (int)(mean) + binsInHistogram/2;
    Vector<Integer> histogram = new Vector<Integer>();
    for (int i=0; i<max-min+1; i++) histogram.add(0);  // Initialize histogram entries
    double total = 0.0;
    for (int i=0; i<v.size(); i++) {  // Populate histogram
      if ((v.get(i) >= min) && (v.get(i) <= max)) {
        int index = v.get(i) - min;
        histogram.set(index, histogram.get(index) + 1);
        total += 1.0;
      }
    }
    
    // Output histogram to file
    try {
      PrintWriter writer = new PrintWriter(new File(fileName));
      writer.println("Score" + "\t" + "Actual # of alignments with score" + "\t" + "Actual % of alignments with score" + "\t" + "Estimated % of alignments with score");
      double beta = std_dev * Math.sqrt(6.0) / Math.PI;
      double mu = mean - (0.5772 * beta);
      for (int i=0; i<histogram.size(); i++) {
        double value = (1.0/beta) * (Math.pow(Math.E, -((min+i - mu)/beta))) * (Math.pow(Math.E, -Math.pow(Math.E, -((min+i - mu)/beta))));
        writer.println((min+i) + "\t" + histogram.get(i) + "\t" + (double)(histogram.get(i))/(double)(total) + "\t" + value);
      }
      writer.close();
    } catch (FileNotFoundException e) {
      System.err.println("Error - the file " + fileName + " could not be opened and written to.");
    }
  }
  
  
  
  /*********************************************************
    ****************** PRIVATE INSTANCE METHODS *************
    *********************************************************/
  
  /**
   * Computes the optimal <em>global</em> pairwise alignment score for two genomic sequences.
   * <p>
   * An alignment table (i.e., 2D integer array) is populated as part of the computation of
   * the optimal alignment score.
   * A backtracking table (i.e., 2D <code>Point</code> array) is also populated so that the
   * optimal alignment (and not just the score) can be re-constructed using the method
   * <code>globalBacktrack</code>.
   */
  private void globalAlignment() {
    // Initialize first entry (0, 0) in tables
    alignmentTable[0][0] = 0;  //  Alignment table
    backtrackTable[0][0] = new Point(-1, -1);  // Backtrack table
    
    // Initialize first column in tables
    for (int i=1; i<alignmentTable.length; i++) {
      alignmentTable[i][0] = i*getGapScore(i, 0, true);  // Alignment table
      backtrackTable[i][0] = new Point(i-1, 0);  // Backtrack table
    }
    
    // Initialize first row in tables
    for (int i=1; i<alignmentTable[0].length; i++) {
      alignmentTable[0][i] = i*getGapScore(0, i, false);  // Alignment table
      backtrackTable[0][i] = new Point(0, i-1);  // Backtrack table
    }
    
    // Fill in the remaining entries in the two tables
    // In a global alignment, each entry is the maximum of 3 possible values,
    // which are computed based on 3 previously determined table entries.
    for (int i=1; i<alignmentTable.length; i++) {
      for (int j=1; j<alignmentTable[0].length; j++) {
        int gapAbove = alignmentTable[i-1][j] + getGapScore(i, j, true);
        int gapLeft = alignmentTable[i][j-1] + getGapScore(i, j, false);
        int match = alignmentTable[i-1][j-1] + getScore(seq1.charAt(i-1), seq2.charAt(j-1));
        int maximum = max(gapAbove, gapLeft, match);
        
        // Assign value to Alignment table
        alignmentTable[i][j] = maximum;
        
        // Assign value to Backtrack table
        if (match == maximum) backtrackTable[i][j] = new Point(i-1, j-1);
        else if (gapAbove == maximum) backtrackTable[i][j] = new Point(i-1, j);
        else backtrackTable[i][j] = new Point(i, j-1);
      }
    }
    // Keep track of optimal global alignment score
    optimalScore = alignmentTable[alignmentTable.length - 1][alignmentTable[0].length - 1];
  }
  
  /**
   * Computes the optimal <em>local</em> pairwise alignment score for two genomic sequences.
   * <p>
   * An alignment table (i.e., 2D integer array) is populated as part of the computation of
   * the optimal alignment score.
   * A backtracking table (i.e., 2D <code>Point</code> array) is also populated so that the
   * optimal alignment (and not just the score) can be re-constructed using the method
   * <code>localBacktrack</code>.
   */
  private void localAlignment() {
    
    // Initialize first column in tables with zeroes
    for (int i=1; i < alignmentTable.length; i++) {
      alignmentTable[i][0] = 0;  // Alignment table
      backtrackTable[i][0] = new Point(i-1, 0);  // Backtrack table
    }
    
    // Initialize first row in tables with zeroes
    for (int i=1; i < alignmentTable[0].length; i++) {
      alignmentTable[0][i] = 0;  // Alignment table
      backtrackTable[0][i] = new Point(0, i-1);  // Backtrack table
    }
    
    // Step 2: fill in alignmentTable and backtrackTable
    int localOptimal = 0;
    localOptimalLocation = new Point(0, 0);
    
    for (int i = 1; i < alignmentTable.length; i++) {
      for (int j = 1; j < alignmentTable[0].length; j++) {
        int gapAbove = alignmentTable[i-1][j] + getGapScore(i, j, true);
        int gapLeft = alignmentTable[i][j-1] + getGapScore(i, j, false);
        int match = alignmentTable[i-1][j-1] + getScore(seq1.charAt(i-1), seq2.charAt(j-1));
        int maximum = max(gapAbove, gapLeft, match);
        if (maximum < 0){
         maximum = 0;
        }
//        
//        System.out.println("comparing " + seq2.charAt(j -1) + " and " + seq1.charAt(i - 1));
//        System.out.println("maximum: " + maximum + " index: " + i + ", " + j);

        // Assign value to Alignment table
        alignmentTable[i][j] = maximum;
        
//        System.out.println(alignmentTable[i][j]);
        
        // Assign value to Backtrack table
        if (match == maximum || maximum == 0) backtrackTable[i][j] = new Point(i-1, j-1);
        else if (gapAbove == maximum) backtrackTable[i][j] = new Point(i-1, j);
        else backtrackTable[i][j] = new Point(i, j-1);
//        System.out.println("backtrack (x,y): " + backtrackTable[i][j]  + "\n");
        if (localOptimal < maximum){ 
          localOptimal = maximum; //reset maximum
          localOptimalLocation.setLocation(i, j);
          //System.out.println(localOptimal);
        }
      }
    }
    
    this.optimalScore = localOptimal;
    
  }
  
  /**
   * Computes the optimal <em>global</em> pairwise alignment score for two genomic sequences when considering
   * at least a small fixed number of gaps.
   * <p>
   * An alignment table (i.e., 2D integer array) is partially populated as part of the computation of
   * the optimal alignment score.
   * For k gaps, runs in O(k*n) time.
   * A backtracking table (i.e., 2D <code>Point</code> array) is also partially populated so that the
   * optimal alignment (and not just the score) can be re-constructed using the method
   * <code>globalBacktrack</code>.
   */
  private void globalAlignment_FAST() {
    
    // Initialize first entry (0, 0) in tables
    alignmentTable[0][0] = 0;  //  Alignment table
    backtrackTable[0][0] = new Point(-1, -1);  // Backtrack table
    
    // Initialize first column in tables
    for (int i=1; i<= numGaps; i++) {
      alignmentTable[i][0] = i*getGapScore(i, 0, true);  // Alignment table
      backtrackTable[i][0] = new Point(i-1, 0);  // Backtrack table
    }
 
    // Initialize first row in tables
    for (int i=1; i<=numGaps; i++) {
      alignmentTable[0][i] = i*getGapScore(0, i, false);  // Alignment table
      backtrackTable[0][i] = new Point(0, i-1);  // Backtrack table
    }
    
    // Fill in the remaining entries in the two tables
    // In a global alignment, each entry is the maximum of 3 possible values,
    // which are computed based on 3 previously determined table entries.
    for (int i=1; i<alignmentTable.length; i++) { //loops through every row in the table
      for (int j=i-numGaps; j<=i+numGaps; j++) { //loop through the shorter row rather than the whole row
        if (j<1) j = 1;
        if (j == alignmentTable[0].length) break;
        
        int gapAbove = alignmentTable[i-1][j] + getGapScore(i, j, true);
        int gapLeft = alignmentTable[i][j-1] + getGapScore(i, j, false);
        int match = alignmentTable[i-1][j-1] + getScore(seq1.charAt(i-1), seq2.charAt(j-1));
        int maximum;
        if (j == i+numGaps){ //don't compare above
          maximum = Math.max(gapLeft, match);
        } 
        else if (j == i-numGaps){ //don't compare left
          maximum = Math.max(gapAbove, match);
        }
        else{
          maximum = max(gapAbove, gapLeft, match);
        }
        
        // Assign value to Alignment table
        alignmentTable[i][j] = maximum;
        
        // Assign value to Backtrack table
        if (match == maximum) backtrackTable[i][j] = new Point(i-1, j-1);
        else if (gapAbove == maximum) backtrackTable[i][j] = new Point(i-1, j);
        else backtrackTable[i][j] = new Point(i, j-1);
      }
    }
    // Keep track of optimal global alignment score
    this.optimalScore = alignmentTable[alignmentTable.length - 1][alignmentTable[0].length - 1];
  }
  
  /**
   * Once an alignment <em>score</em> has been computed, the backtracking table is searched to
   * determine the optimal <em>alignment</em>.
   * <p>
   * The optimal <em>global</em> alignment is re-contructed by searching through the
   * backtracking table. The optimal alignment is returned as a <code>String</code>.
   *
   * @return   a <code>String</code> representing the optimal global alignment
   */
  private String globalBacktrack() {
    
    if (useFastAlignment && (Math.abs(seq1.length()-seq2.length()) > numGaps)) return "Too many gaps - no alignment possible.\n";
    
    StringBuilder s1_alignment = new StringBuilder();  // First of two sequences
    StringBuilder alignment = new StringBuilder();     // Alignment of two sequences
    StringBuilder s2_alignment = new StringBuilder();  // Second of two sequences
    
    // Begin backtracking search at optimal global alignment score, i.e., at last entry in table
    int i = alignmentTable.length - 1;  // Last row in table
    int j = alignmentTable[0].length - 1;  // Last column in table
    
    // Backtrack through table to beginning of alignment at position (0,0)
    while ((i != 0) || (j != 0)) {
      Point p = backtrackTable[i][j];  // Current position in table
      if ((p.x == i-1) && (p.y == j-1)) { // match or mismatch
        s1_alignment.insert(0, seq1.charAt(i-1));
        alignment.insert(0, getMatchMismatchChar(seq1.charAt(i-1), seq2.charAt(j-1)));
        s2_alignment.insert(0, seq2.charAt(j-1));
        i--;
        j--;
      } else if ((p.x == i-1) && (p.y == j)) { // gap above
        s1_alignment.insert(0, seq1.charAt(i-1));
        alignment.insert(0, " ");
        s2_alignment.insert(0, "-");
        i--;
      } else { // gap left
        s1_alignment.insert(0, "-");
        alignment.insert(0, " ");
        s2_alignment.insert(0, seq2.charAt(j-1));
        j--;
      }
    }
    
    // When displaying alignment, each line should be no more than some fixed length (60)
    StringBuilder formattedAlignment = new StringBuilder();
    int lineLength = 60;
    int k = 0;
    while (k < s1_alignment.length()) {
      if (k+lineLength > s1_alignment.length()) {  // Last line of alignment
        formattedAlignment.append("\t" + s1_alignment.substring(k) + "\n");
        formattedAlignment.append("\t" + alignment.substring(k) + "\n");
        formattedAlignment.append("\t" + s2_alignment.substring(k) + "\n\n");
      } else {
        formattedAlignment.append("\t" + s1_alignment.substring(k, k+lineLength) + "\n");
        formattedAlignment.append("\t" + alignment.substring(k, k+lineLength) + "\n");
        formattedAlignment.append("\t" + s2_alignment.substring(k, k+lineLength) + "\n\n");
      }
      k += lineLength;
    }
    return formattedAlignment.toString();
  }
  
  /**
   * Once an alignment <em>score</em> has been computed, the backtracking table is searched to
   * determine the optimal <em>alignment</em>.
   * <p>
   * The optimal <em>local</em> alignment is re-contructed by searching through the
   * backtracking table. The optimal alignment is returned as a <code>String</code>.
   *
   * @return   a <code>String</code> representing the optimal local alignment
   */
  private String localBacktrack() {
    
   if (useFastAlignment && (Math.abs(seq1.length()-seq2.length()) > numGaps)) return "Too many gaps - no alignment possible.\n";
    
    StringBuilder s1_alignment = new StringBuilder();  // First of two sequences
    StringBuilder alignment = new StringBuilder();     // Alignment of two sequences
    StringBuilder s2_alignment = new StringBuilder();  // Second of two sequences
    
    // Begin backtracking search at optimal local alignment score, i.e., at localOptimalLocation
    int i = localOptimalLocation.x; 
    int j = localOptimalLocation.y;
    
    // Backtrack through table to beginning of alignment at position (0,0)
    while ((i != 0) || (j != 0)) {
      Point p = backtrackTable[i][j];  // Current position in table
      if(alignmentTable[i][j] == 0) break; // we reached a zero, end alignment
      if ((p.x == i-1) && (p.y == j-1)) { // match or mismatch
        s1_alignment.insert(0, seq1.charAt(i-1));
        alignment.insert(0, getMatchMismatchChar(seq1.charAt(i-1), seq2.charAt(j-1)));
        s2_alignment.insert(0, seq2.charAt(j-1));
        i--;
        j--;
      } else if ((p.x == i-1) && (p.y == j)) { // gap above
        s1_alignment.insert(0, seq1.charAt(i-1));
        alignment.insert(0, " ");
        s2_alignment.insert(0, "-");
        i--;
      } else { // gap left
        s1_alignment.insert(0, "-");
        alignment.insert(0, " ");
        s2_alignment.insert(0, seq2.charAt(j-1));
        j--;
      }
    }
    
    // When displaying alignment, each line should be no more than some fixed length (60)
    StringBuilder formattedAlignment = new StringBuilder();
    int lineLength = 60;
    int k = 0;
    while (k < s1_alignment.length()) {
      if (k+lineLength > s1_alignment.length()) {  // Last line of alignment
        formattedAlignment.append("\t" + s1_alignment.substring(k) + "\n");
        formattedAlignment.append("\t" + alignment.substring(k) + "\n");
        formattedAlignment.append("\t" + s2_alignment.substring(k) + "\n\n");
      } else {
        formattedAlignment.append("\t" + s1_alignment.substring(k, k+lineLength) + "\n");
        formattedAlignment.append("\t" + alignment.substring(k, k+lineLength) + "\n");
        formattedAlignment.append("\t" + s2_alignment.substring(k, k+lineLength) + "\n\n");
      }
      k += lineLength;
    }
    return formattedAlignment.toString();
    
  }
  
  /**
   * Returns the gap score for a particular position in a pairwise alignment.
   * <p>
   *
   *
   * @param   i            the index for sequence 1, i.e., the row of the alignment/backtrack tables
   * @param   j            the index for sequence 2, i.e., the column of the alignment/backtrack tables
   * @param   isGapAbove   for affine gap scoring, indicates an extension of a gap from above (as opposed to a gap from the left)
   * @return               the gap score associated with the indices <code>i</code> and <code>j</code>
   */
  private int getGapScore(int i, int j, boolean isGapAbove) {
    if (!useAffineGapScoring)
      return linearGapScore;
    else {  // Affine gap scoring
      if (isGapAbove){
        //check if cell above had a gap:
        if (backtrackTable[i-1][j].equals(new Point(i-2, j)) || backtrackTable[i-1][j].equals(new Point(i-1, j-1))){
          return betaGapScore;
        }
        return alphaGapScore;
      }
      else{
      //check if cell to left had a gap:
        if (backtrackTable[i][j-1].equals(new Point(i-1, j-1)) || backtrackTable[i][j-1].equals(new Point(i, j-2))){
          return betaGapScore;
        }
        return alphaGapScore;
      }
    }
  }
  
  /**
   * Returns the alignment score for two characters.
   *
   * @param   c1   the first of two aligned characters
   * @param   c2   the second of two aligned characters
   * @return       the score associated with aligning the two characters <code>c1</code> and <code>c2</code>
   */
  private int getScore(char c1, char c2) {
    if (!useMatrixScoring) {  // Fixed match/mismatch scoring
      if (c1 == c2) return matchScore;
      return mismatchScore;
    } else {  // A matrix of match/mismatch scores
      String key = Character.toString(c1)+Character.toString(c2);
      return matrixScoring.get(key);
      
    }
  }
  
  /**
   * Determines how aligned characters are displayed.
   * <p>
   * If two aligned characters are identical, a "|" is shown in the alignment.
   * If two aligned characters are similar but not identical (i.e., match score
   * greater than zero), a ":" is shown in the alignment. If two aligned
   * characters are dissimilar, a " " is shown in the alignment.
   *
   * @param   c1   the first of two aligned characters
   * @param   c2   the second of two aligned characters
   * @return       a character indicating whether <code>c1</code> and <code>c2</code> are identical or similar or dissimalar
   */
  private char getMatchMismatchChar(char c1, char c2) {
    if (c1 == c2) return '|';  // Characters are identical
    else if (getScore(c1, c2) > 0) return ':';  // Characters are similar
    return ' ';  // Characters are dissimilar
  }
  
  /**
   * Returns a <code>String</code> representation of a table (i.e., 2D array) of integers.
   *
   * @param   t   a 2-dimenstional array of integers
   * @return   a <code>String</code> representation of a 2D integer array
   */
  private String tableToString(int[][] t) {
    StringBuilder sb = new StringBuilder("\t");
    for (int i=0; i<seq2.length(); i++)
      sb.append("\t" + seq2.charAt(i));
    sb.append("\n");
    for (int i=0; i<t.length; i++) {
      if (i != 0) sb.append(seq1.charAt(i-1));
      
      for (int j=0; j<t[i].length; j++) {
        sb.append("\t" + t[i][j]);
      }
      sb.append("\n");
    }
    return sb.toString();
  }
  
  /**
   * Returns a <code>String</code> representation of a table (i.e., 2D array) of <code>Points</code>.
   *
   * @param   t   a 2-dimenstional array of <code>Points</code>
   * @return   a <code>String</code> representation of a 2D <code>Point</code> array
   */
  private String tableToString(Point[][] t) {
    StringBuilder sb = new StringBuilder("\t");
    for (int i=0; i<seq2.length(); i++)
      sb.append("\t" + seq2.charAt(i));
    sb.append("\n");
    for (int i=0; i<t.length; i++) {
      if (i != 0) sb.append(seq1.charAt(i-1));
      
      for (int j=0; j<t[i].length; j++) {
        if (t[i][j] == null) sb.append("\t" + "-1" + ":" + "-1");
        else sb.append("\t" + t[i][j].x + ":" + t[i][j].y);
      }
      sb.append("\n");
    }
    return sb.toString();
  }
  
  
  
  /*********************************************************
    ****************** PRIVATE CLASS METHODS ****************
    *********************************************************/
  
  /**
   * Reads in and returns a genomic sequence from the specified FASTA file.
   *
   * @param   f   a <code>File</code> object referring to a FASTA file containing a genomic sequence
   * @return      the genomic sequence read in from a FASTA file
   */
  private static String readSequenceFromFile(File f) {
    StringBuilder sequence = new StringBuilder();
    try {
      Scanner reader = new Scanner(f);
      String header = reader.nextLine();  // Header line of FASTA file
      if ((header.length() == 0) || (header.charAt(0) != '>')) {
        System.err.println("Error - first line of file " + f + " is not in FASTA format.");
        return sequence.toString();
      }
      while (reader.hasNext()) {  // continue until we reach end of file
        sequence.append(reader.nextLine());
      }
      reader.close();
    } catch (FileNotFoundException e) {
      System.err.println("Error - the file " + f + " could not be found and opened.");
      return sequence.toString();
    }
    return sequence.toString();
  }
  
  /**
   * Returns the maximum of three integers.
   *
   * @param   a   the first of three integers
   * @param   b   the second of three integers
   * @param   c   the third of three integers
   * @return      the maximum of three specified integers
   */
  private static int max(int a, int b, int c) {
    return Math.max(a, Math.max(b, c));
  }
  
  /**
   * Returns the mean (i.e., the average value) of a <code>Vector</code> of numbers.
   *
   * @param   v   a <code>Vector</code> of numbers
   * @return      the mean of the numbers in a <code>Vector</code>
   */
  private static double getVectorMean(Vector<Integer> v) {
    double sum = 0.0;
    double count = 0.0;
    Iterator<Integer> itr = v.iterator();
    while(itr.hasNext()){
      sum += itr.next();
      count++;
    }
    return sum/count;
    
  }
  
  /**
   * Returns the sample standard deviation of a <code>Vector</code> of numbers.
   *
   * @param   v   a <code>Vector</code> of numbers
   * @return      the standard deviation of the numbers in a <code>Vector</code>
   */
  private static double getVectorStandardDeviation(Vector<Integer> v) {
    double mean = getVectorMean(v);
    double sigma = 0.0;
    double count = 0.0;
    Iterator<Integer> itr = v.iterator();
    while(itr.hasNext()){
      double diff = Math.pow((itr.next() - mean), 2);
      sigma += diff;
      count++;
    }
    return Math.sqrt(sigma/(count-1));
    
  }
  
  public Vector<Integer> getAlignScores(){
   return alignScores; 
  }
  
  
  /*********************************************************
    ****************** MAIN METHOD **************************
    *********************************************************/
  
  /**
   * The <code>main</code> method creates an optimal pairwise alignment for two genomic sequences.
   * <p>
   * The <code>main</code> method expects the <code>Alignment</code> program is
   * executed with exactly two command line arguments: the names of two FASTA files each
   * containing a genomic sequence. The <code>main</code> method creates an 
   * <code>Alignment</code> object based on the two genomic sequences and computes
   * the optimal pairwise alignment of the sequences.
   * 
   * @param   args   an array of Strings representing any command line arguments
   */
  public static void main(String [] args) {
    
    if (args.length < 2) {
      System.err.println("When executing this program, please enter the name of two files,");
      System.err.println("each containing a sequence. The program will align the two sequences.");
      System.err.println("\tjava Alignment file1.txt file2.txt\n");
      //System.exit(0);
    } else {
      
      Alignment a = new Alignment(new File(args[0]), new File(args[1]));
      
      // Global alignment with fixed scoring and linear gap scoring
      a.setGlobalAlignment();
      a.setFixedScoring(5, -4);
      a.setLinearGaps(-6);
      a.computeAlignment();
      a.calculatePValue();
      System.out.println(a);

       // Local alignment with fixed scoring and linear gap scoring
       a.setLocalAlignment();
       a.setFixedScoring(5, -4);
       a.setLinearGaps(-6);
       a.computeAlignment();
       a.calculatePValue();
       System.out.println(a);

       // Global alignment with matrix scoring and linear gap scoring
       a.setGlobalAlignment();
       a.setMatrixScoring("data/Blosum62_Matrix.txt");
       a.setLinearGaps(-6);
       a.computeAlignment();
       a.calculatePValue();
       System.out.println(a);
      
       // Global alignment with fixed scoring and affine gap scoring
       a.setGlobalAlignment();
       a.setFixedScoring(5, -4);
       a.setAffineGaps(-11, -2);
       a.computeAlignment();
       a.calculatePValue();
       System.out.println(a);
       
       // Global alignment *FAST* considering at least 2 gaps
       a.setGlobalAlignment();
       a.setFixedScoring(5, -4);
       a.setLinearGaps(-6);
       a.setFastAlignment(2);
       a.computeAlignment();
       a.calculatePValue();
       System.out.println(a);
       
       //calculate the p-value for two longer sequences
       Alignment a2 = new Alignment("data/protein_seq1.txt", "data/protein_seq2.txt");
       a2.setGlobalAlignment();
       a2.setFixedScoring(5, -4);
       a2.setLinearGaps(-6);
       a2.computeAlignment();
       a2.calculatePValue();
       outputHistogramOfRandomAlignmentScores("outputHistogramData.tsv", a2.getAlignScores());
       
    }
    
  }
  
}

