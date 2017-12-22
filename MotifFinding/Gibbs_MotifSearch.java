/* Jess Abramson and Caroline Hsu
 * November 16, 2017
 */


import java.util.*;
import java.io.*;

public class Gibbs_MotifSearch extends EM_MotifSearch{



	//constructor!! 
	public Gibbs_MotifSearch(String fileName, int motifLength){

		super(fileName, motifLength);
	}


	/*
	Returns the index of a randomly sampled value in a Vector.
	One value from the Vector is chosen at random and the value's index (not the 
	value itself) is returned. The value is not chosen uniformly at random, but rather 
	via sampling, i.e., higher values are more likely to be chosen and lower values 
	are less likely to be chosen.

	One approach for randomly sampling a collection of values proceeds as follows:

		  - Normalize the values in the collection so that they sum to 1.0. The values 
			now represent a probability distribution.
		  - Convert the values from a probability distribution to a cumulative distribution. 
		  	In a cumulative distribution, the value at index i represents the sum of all values 
		  	at indices less than or equal to i in the probability distribution. The final value 
		  	in a cumulative distribution should be 1.0 since the sum of all values in a 
		  	probability distribution is 1.0.
		  - Generate a number uniformly at random between 0.0 and 1.0. Return the index of 
			the smallest value in the cumulative distribution that is at least as big as 
			the random number.

	Parameters:
		values - a Vector of decimal numbers to be sampled

	Returns:
		the index of a randomly sampled value from the Vector
	*/
	public int getIndexViaSampling(Vector<Double> values){

		//instead of choosing biggest, choose random 
		//but higher scores are most likely to be chosen 
		int randomIndex = -1;
		double sumOfVectorVals = 0;

		//find the sum of all the values in the vector inputted
		for (int i = 0; i < values.size(); i++) {
			sumOfVectorVals += values.get(i);
		}


		double valsOfPrevious = 0;

		//get random number
		Random rand = new Random();
		double randomDouble = rand.nextDouble();


		//normalize vector and turn into cumulative distribution
		for (int i =0; i < values.size(); i++) {

			double normVal = values.get(i)/sumOfVectorVals;
			valsOfPrevious += normVal;

			values.set(i, valsOfPrevious);

			if (valsOfPrevious > randomDouble){
				randomIndex = i;
				break;

			}
		}

		return randomIndex;
	}


	/*
	The Expectation step in the EM algorithm.
	Based on the matrix model, identifies motif instances in the sequences. One motif 
	instance is identified in each sequence. For each sequence, a motif instance is 
	chosen by sampling the scores of each possible motif instance in that sequence. 
	The score of each possible motif instance is based on the matrix model.

	Overrides:
		determineMotifInstances in class EM_MotifSearch
	*/
	public void determineMotifInstances(){

		//keep track of all the scores in the sequence in like vector 
		//then pass vector into index sampling  

		for (int i = 0; i < numSequences; i++) {
	      
	      //so our array works
	      int maxIndex = sequences.get(i).length() - motifLength + 1;

	      //create a vector of doubles to store score
	      Vector<Double> scoresOfMotifs = new Vector<Double>();
	      
	      //loop through the motif at i and store the score for each in the vector of scores
	      
	      for (int j = 0; j < maxIndex; j++){

	        //grab the motif 
	        String substring = sequences.get(i).substring(j,j+motifLength);

	        //grab the score
	        double currentScore = getScoreForMotifInstance(substring);

	        //add score to vector of scores
	        scoresOfMotifs.add(currentScore);
	        }

		    //find random score with helper function 
		    int randInstanceIndex = getIndexViaSampling(scoresOfMotifs);

		    //set the new best motif based on current matrix model
		    instanceLocations.set(i, randInstanceIndex);
	      }

	    }
	
	


	public static void main(java.lang.String[] args){

	if (args.length < 3) {
      System.err.println("\nTo execute this program, three command line arguments are required corresponding to a FASTA file name, the length of the desired motif, and the number of Gibbs sampling iterations. For example,\n");
      System.err.println("\tjava Test_Gibbs modA.txt 16 500\n\n");
      return;
    }

    Gibbs_MotifSearch m = new Gibbs_MotifSearch(args[0], Integer.parseInt(args[1]));
    m.run_EM_multiple_times(Integer.parseInt(args[2]));
    
    // Output results
    System.out.println("\n" + m.motifInstancesToString());
    System.out.println(m.matrixToString());
    System.out.println("Consensus sequence: " + m.getConsensusSequence() + "\n");
    System.out.println("Information content: " + m.getInformationContentOfMatrix() + "\n");
	}
	
}
	  