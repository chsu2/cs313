/* Caroline Hsu
 * November 16, 2017
 */

import java.util.*;
import java.io.*;

public class EM_MotifSearch extends MotifSearch{
  
  
  public EM_MotifSearch(String fileName, int motifLength) {

    //call super from the parent function
    super(fileName, motifLength);
  }
  
  //The initial random seed step in the EM algorithm.
  //For each sequence, randomly determine the start index of a motif instance in the sequence.
  
  public void setRandomLocationsForMotifInstances() {

    //create a random int to find random locations for the motif instances 
    Random rand = new Random();

    //want a max index so we don't go out of range
    int maxIndex = 0;

    //loops through the sequences and finds a random index
    for (int i = 0; i < numSequences; i++){

      maxIndex = sequences.get(i).length() - motifLength + 1; //be careful of indexing here
      int intg = rand.nextInt(maxIndex);

      //save the location in the instanceLocations vector
      instanceLocations.set(i,intg);
    }
    
  }
  
  //The Maximization step in the EM algorithm.
  //Based on the motif instances in the sequences, creates a matrix motif model. The 
  //resulting matrix model should be updated with pseudocounts so that no entries in the 
  //matrix correspond to 0.0.
  
  public void determineMatrixModel(){
  
    //a vector of strings to store the motifs    
    Vector<String> motifs = new Vector<String>();

    //clear so that you can find new motifs
    // motifs.clear();

    //loop through the sequences 
    for (int i = 0; i < numSequences; i++){

      //get the instance location 
      int startLoc = instanceLocations.get(i);

      //add the motifs at the given instance locations to the vector 
      motifs.add(sequences.get(i).substring(startLoc, startLoc + motifLength)); // be careful of indexing here

    }
    
    //make a loop as long as the motif length
    for (int i = 0; i < motifLength; i++){

      double countA = 0;
      double countC = 0;
      double countT = 0;
      double countG = 0;
      
      // System.out.println("loop number: " + i);

      // System.out.println("matrix elements before replaced: ");
      // System.out.println("a: " + matrix[0][i]);
      // System.out.println("c: " + matrix[1][i]);
      // System.out.println("g: " + matrix[2][i]);
      // System.out.println("t: " + matrix[3][i]);
      
      
      //then loop through the jth character in the motifs and count up all the characters for that index
      for (int j = 0; j < numSequences; j++){
        
        char current = motifs.get(j).charAt(i);
        
        if (current == 'A') {

          countA++;

        } else if (current == 'C') {

          countC++;

        } else if (current == 'T') {

          countT++;

        } else if (current == 'G') {

          countG++;

        } 
      }
      
      //add the appropriate values to the motif matrix
      matrix[0][i] = countA / ((double) numSequences); // first row is A
      matrix[1][i] = countC / ((double) numSequences); // (they're in alaphetical order)
      matrix[2][i] = countG / ((double) numSequences);
      matrix[3][i] = countT / ((double) numSequences);

      // System.out.println("matrix elements after replaced: ");
      // System.out.println("a: " + matrix[0][i]);
      // System.out.println("c: " + matrix[1][i]);
      // System.out.println("g: " + matrix[2][i]);
      // System.out.println("t: " + matrix[3][i]);

    }


    //call this function so that no 0 is actually a 0
    addPseudocountsToMatrix();
    
  }

  //Given a candidate motif instance, returns the score (probability) of that instance based on the matrix model.
  //The length of the motif instance specified by String s must be the same as the number of columns in the matrix model.
  /*Parameters:
   s - a String corresponding to a motif instance
   Returns:
   the score (probability) of the motif instance matching the matrix model
   */
  
  public double getScoreForMotifInstance(String s){
    
    double scoreSoFar = 1; 
    //omg...start at 1 because anything multiplied by 0 is 0!!! fooool!!!
    
    //loop through and calculate score from matrix model
    for (int k = 0; k < s.length(); k++) {

          char current = s.charAt(k);

          if (current == 'A') {

            scoreSoFar *= matrix[0][k];

          } else if (current == 'C') {

            scoreSoFar *= matrix[1][k];

          } else if (current == 'G') {

            scoreSoFar *= matrix[2][k];

          } else if (current == 'T') {

            scoreSoFar *= matrix[3][k];

          }
          
        }
    
    //return the score
    return scoreSoFar;
  }
  
  
  //The Expectation step in the EM algorithm.
  //Based on the matrix model, identifies motif instances in the sequences. One motif instance 
  //is identified in each sequence. For each sequence, the motif instance that best matches 
  //the model is chosen.
  //ugh
  public void determineMotifInstances(){
    
    //loop through all the sequences
    for (int i = 0; i < numSequences; i++) {

      //start some counts
      int bestStart = -5;
      double bestScore = -7;
      
      //so our array works
      int maxIndex = sequences.get(i).length() - motifLength + 1;
      
      //loop through the motif at i and sees which substring has the best score from the matrix

      for (int j = 0; j < maxIndex; j++){

        //grab the motif 
        String substring = sequences.get(i).substring( j, j + motifLength);

        //grab the score
        double currentScore = getScoreForMotifInstance(substring);

        //compare 
        if (currentScore > bestScore) {
          bestScore = currentScore;
          bestStart = j;
        }
      }

      //set the new best motif based on current matrix model
      instanceLocations.set(i,bestStart);

    }
  }
  
  
//  Returns the information content associated with the matrix model.
//The information content of a matrix model is described in Task 2 of 
  //Exercise 6. The method getNucleotideContent may be useful in determining 
  //the background frequency of different nucleotides.

//Returns:
//the information content of the matrix model
  public double getInformationContentOfMatrix(){
    
    //to store the information content
    double informationContent = 0;
    
    //see how many times the nucleotide appears in the sequences 
    double bFreqA = getNucleotideContent('A');
    double bFreqC = getNucleotideContent('C');
    double bFreqG = getNucleotideContent('G');
    double bFreqT = getNucleotideContent('T');
    
    //find the information content for each idex of the motif
    for (int i = 0; i < motifLength; i ++){

      double aVal = matrix[0][i];
      double cVal = matrix[1][i];
      double gVal = matrix[2][i];
      double tVal = matrix[3][i];
     
      double aPart = aVal * (Math.log(aVal/bFreqA)/Math.log(2));
      double cPart = cVal * (Math.log(cVal/bFreqC)/Math.log(2));
      double gPart = gVal * (Math.log(gVal/bFreqG)/Math.log(2));
      double tPart = tVal * (Math.log(tVal/bFreqT)/Math.log(2));
      
      informationContent += aPart + cPart + gPart + tPart;
      
    }

    //return the content
    return informationContent;
  }
  
//  Executes the EM (Expectation Maximization) algorithm.
//  Initially, the EM algorithm is randomly seeded, i.e., one motif instance is randomly chosen in each sequence. Then, the Maximization and Expectation steps are alternately repeated until convergence. The EM algorithm converges when the information content of the matrix model no longer improves.
  public void EM(){

    //first set randomm locations for the motifs
    setRandomLocationsForMotifInstances();

    //for testing purposes 
    // int count = 0;

    //for comparisons 
    double infoContent = getInformationContentOfMatrix();
    double lastInfoContent = 0;

    //run till the info content doesn't improve
    while (lastInfoContent != infoContent) {

      // System.out.println("ENTERING LOOP NUMBER: " + count);

      //determines a motif matrix based on the motif instances 
      determineMatrixModel();

      //EM step: finds best motifs in each sequence based on the matrix model 
      determineMotifInstances();

      //set the last info content to what it was previously 
      lastInfoContent = infoContent;

      //replace the infoContent with the new info content of the current matrix
      infoContent = getInformationContentOfMatrix();

      //for testing purposes
      // count++;
    }
  }
  
  
//  Executes the EM (Expectation Maximization) algorithm multiple times.
//  The number of times that the algorithm is executed is specified by the integer parameter. The best motif, as determined by information content, is identified over all executions of the algorithm. Upon completion of this method, this EM_MotifSearch should correspond to the best motif (including matrix model and motif instances) identified over all executions of the algorithm.
  
  //can just save the instances and then compute the matrix and ifoContent from that
  public void run_EM_multiple_times(int iterations){

    //keep track of the best motif matrix, best scores, and best instances
    double bestScore = -10000000;
    double currentScore = -1;
    Vector<Integer> bestInstances = new Vector<Integer>();

    //loop through the number of times you are told to
    for (int i = 0; i < iterations; i++) {

      // System.out.println("ENTERING ITERATION NUMBER: " + i);

      //run EM
      EM();

      //get the info content of the matrix
      currentScore = getInformationContentOfMatrix();
      // System.out.println("CURRENT SCORE: " + currentScore);

      //if it's better than the previous best
      if (currentScore > bestScore) {

        bestInstances.clear(); //omg clear OUTSIDE OF THE LOOOOOOP
        
        //copy over the best instances 
        bestInstances = getInstanceLocations();

        //replace the best score for reference 
        bestScore = currentScore;
      }
    }

    //replace the instance variables
  
    instanceLocations.clear();

    for (int l = 0; l < numSequences; l++) {
      instanceLocations.add(bestInstances.get(l));
    }

    //recompute matrix
    determineMatrixModel();

  }
  
  public static void main(String[] args){
      if (args.length < 3) {
      System.err.println("\nTo execute this program, three command line arguments are required corresponding to a FASTA file name, the length of the desired motif, and the number of EM iterations. For example,\n");
      System.err.println("\tjava Test_EM modA.txt 16 500\n\n");
      return;
    }

    EM_MotifSearch m = new EM_MotifSearch(args[0], Integer.parseInt(args[1]));
    m.run_EM_multiple_times(Integer.parseInt(args[2]));

    // m.EM();
    
    // Output results
    System.out.println("\n" + m.motifInstancesToString());
    System.out.println(m.matrixToString());
    System.out.println("Consensus sequence: " + m.getConsensusSequence() + "\n");
    System.out.println("Information content: " + m.getInformationContentOfMatrix() + "\n");
  } 
  
  
  
  
}