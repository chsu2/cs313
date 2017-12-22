
import java.io.*;
import java.util.*;
import java.lang.*;

public class SequenceOps{

	//create hashtable
	private static Hashtable <String, String> codons = new Hashtable <String, String>();

	/*Returns a String corresponding to the genomic sequence found 
	in the specified FASTA file*/

	public static String sequenceFromFastaFile(String fileName){

		//create bufferedreader and filereader to read in file
		BufferedReader br = null;
		FileReader fr = null;

		//create stringbuilder for the genome sequence and first line
		StringBuilder genomeSequence = new StringBuilder();
		StringBuilder firstLine = new StringBuilder();

		//string for the current line of the file
		String currentLine;

		try {

			//create bufferedreader and filereader to read in file
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);

			//read in lines in file with a loop
			while ((currentLine = br.readLine()) != null){

				//don't append first line
				if (currentLine.contains(">")){

					firstLine.append(currentLine);

				}

				else {
					//add to genomeSequence if not first line
					genomeSequence.append(currentLine);
				}
			}

		// catch exceptions 
		} catch (IOException e){

			e.printStackTrace();
		
		//when done, close readers
		} finally {

			try {

				if (br != null) br.close();

				if (fr != null) fr.close();

			} catch (IOException ex){

				ex.printStackTrace();
			}
		}
		
		// System.out.println(firstLine);
		// System.out.println(genomeSequence);

		return genomeSequence.toString();

	}

	/*Create a helper function to read in the translation text for the codons and 
	populate a global hashtable that will be accessed later by a different function
	

	wanted to call this at the top of the class so i could make the hashtable a populated
	global variable but couldn't figure it out*/

	private static void fileToHashtable(String fileName){

		//create bufferedreader and filereader to read in file
		BufferedReader br = null;
		FileReader fr = null;

		//create stringbuilder for the genome sequence and first line
		StringBuilder firstLine = new StringBuilder();

		//string for the current line of the file
		String currentLine;

		try {

			//create bufferedreader and filereader to read in file
			fr = new FileReader(fileName);
			br = new BufferedReader(fr);

			//do not need the first line of the document in this case
			br.readLine();

			//read in lines in file with a loop
			//split into array with values separated by tabs 
			//then add to hashtable 
			while ((currentLine = br.readLine()) != null){

				String [] temp = currentLine.split("\t");
				// System.out.println("the codon: " + temp[0] +"\nthe key: " + temp[1]);
				codons.put(temp[0], temp[1]);

			}

		// catch exceptions 
		} catch (IOException e){

			e.printStackTrace();
		
		//when done, close readers
		} finally {

			try {

				if (br != null) br.close();

				if (fr != null) fr.close();

			} catch (IOException ex){

				ex.printStackTrace();
			}
		}

	}

	/*Takes the specified genomic sequence s and returns a String 
	representing the sequence in FASTA format. Each line in the 
	returned String should be 60 characters in length, except 
	possibly the last line. The returned String need not include 
	the FASTA header line beginning with the character '>'.*/

	public static String sequenceToFastaFormat(String s){
		
		//create an empty stringbuilder
		StringBuilder fastaFormat = new StringBuilder();

		//keep an index for the parsing
		int index = 0;

		//period for parsing
		int period = 60;

		//loop through the string
		while (index < s.length()){

			//append 60 characters at a time to the stringbuilder. use Math.min for the upper bound
			fastaFormat.append(s.substring(index, Math.min(index + period, s.length())));

			//append new line character
			fastaFormat.append("\n");

			//increase the index
			index += period;
		}

		// System.out.println("fasta format: \n" + fastaFormat.toString());

		return fastaFormat.toString();

	}

	/*Returns a String corresponding to the reversed version of 
	the specified sequence s. For example, if s represents the 
	sequence ACGGACTGC then the method should return the string 
	CGTCAGGCA.*/

	public static String reverse(String s){

		//convert string to stringBuilder
		StringBuilder reverse = new StringBuilder(s);

		//use stringBuilder's method to reverse string
		reverse.reverse();

		return reverse.toString();

	}

	/*Returns a String corresponding to the complement of the 
	specified nucleotide sequence s. For example, if s represents 
	the sequence ACGGACTGC then the method should return the string 
	TGCCTGACG.*/

	public static String complement(String s){

		//in case the file is miswritted
		s.toUpperCase();

		//create stringbuilder to add to
		StringBuilder complement = new StringBuilder();

		/*loop through and switch:
			A --> T
			T --> A
			C --> G	
			G --> C 
			*/

		for (int i = 0; i < s.length(); i++) {
			
			char currentChar = s.charAt(i);

			if (currentChar == 'A') complement.append('T');

			if (currentChar == 'T') complement.append('A');

			if (currentChar == 'C') complement.append('G');

			if (currentChar == 'G') complement.append('C');


		}

		return complement.toString();

	}

	/*Returns a String corresponding to the reverse complement of 
	the specified nucleotide sequence s. For example, if s 
	represents the sequence ACGGACTGC then the method should return 
	the string GCAGTCCGT.*/

	public static String reverseComplement(String s){

		//use the functions we wrote earlier to get the reverse complement easily

		return complement(reverse(s));


	}

	/*Returns the GC content of the specified nucleotide sequence 
	s. For example, if s represents the sequence ACGGACTGC then the 
	method should return 0.6666666666666666.

	GC content of:

	- Escherichia coli:		0.5047480343799055 
	- Human Chromosome 22: 	0.4791780185638117

	*/

	public static double GC_content(String s){

		//start a count
		double gcCount = 0;

		for (int i = 0; i < s.length(); i++) {
			
			char currentChar = s.charAt(i);

			//if G or C add to the count
			if (currentChar == 'C' || currentChar == 'G') gcCount++;

		}

		double content = gcCount/s.length();

		return content;

	}

	/*Returns a random permutation of the specified sequence s. The 
	returned String should have all the same characters as the input 
	String s, only the order of the characters should be determined 
	randomly. Each invocation of the method on a particular String 
	s should result in a (almost certainly) different permutation of 
	s.*/

	public static String randomPermutation(String s){

		/* Using the Fisher-Yates shuffle. should have runtime O(n) 
		Fisher-Yates essentially puts all the elements into a hat and then */

		//create a stringbuilder to act as an array

		StringBuilder randomPerm = new StringBuilder(s);

		Random r = new Random();

		//permutates through the string backwards
		for (int i = randomPerm.length() - 1; i > 0 ; i --) {
			
			//grab a random number 0 through i
			int index = r.nextInt(i);

			//swap 
			char temp = randomPerm.charAt(index);

			//set the character at index index to the character at position i in s
			randomPerm.setCharAt(index, randomPerm.charAt(i));
			randomPerm.setCharAt(i, temp);

		}

		return randomPerm.toString();

	}

	/*Returns a random genomic sequence of the specified length with 
	the expected specified GC content. Each nucleotide in the random 
	sequence can be generated independently. The probability that 
	each nucleotide is a G or C (as opposed to A or T) is specified 
	by the GC_content parameter. Each nucleotide has an equal 
	probability of being A or T. Each nucleotide has an equal 
	probability of being G or C. Note: each invocation of the 
	randomSequence method may not generate a random sequence with a 
	GC content equal to that specified by the GC_content parameter. 
	However, the expected GC content will be equal to the GC_content 
	parameter, i.e., if you invoked the method infinitely many times 
	then the GC content of all randomly generated sequences should 
	equal the GC_content parameter.*/

	public static String randomSequence(int length, double GC_content){


		//create two stringbuilders with respective nucleotides
		StringBuilder at = new StringBuilder("AT");
		StringBuilder cg = new StringBuilder("CG");

		//create a random number for assigning 0 or 1
		Random r = new Random();

		StringBuilder notRandomString = new StringBuilder();

		//tells you how many numbers should be GC. may always be an underapproximation though
		//since we round down
		int numOfGC = (int) (length * GC_content); 

		//GC loop
		for (int i = 0; i < numOfGC ; i++) {

			int index = r.nextInt(2);
			notRandomString.append(cg.charAt(index));
			
		}

		//AT loop
		for (int i = 0; i < length - numOfGC ; i++) {

			int index = r.nextInt(2);

			notRandomString.append(at.charAt(index));
			
		}

		//use randomPermutation to randomize the string

		return randomPermutation(notRandomString.toString());

	}

	/*Returns a random genomic sequence with the same length as s 
	and the same expected GC content as s. Note: each invocation of 
	the randomSampling method may not generate a random sequence 
	with a GC content equal to that of s. However, the expected GC 
	content will be equal to that of s, i.e., if you invoked the 
	method infinitely many times then the GC content of all randomly 
	generated sequences should equal that of s.*/

	public static String randomSampling(String s){

		//figure out the GC content of the input string 
		double gcContent = GC_content(s);
		// System.out.println("the gc content in this randomSampling is " + String.valueOf(gcContent));
		//then use randomSequence and return a random sampling


		return randomSequence(s.length(), gcContent);

	}

	/*Assuming the input sequence s is a genomic sequence of exactly 
	3 nucleotides, the method translates the codon, i.e., it returns 
	a character representing the amino acid corresponding to the 
	codon. For example, if the input sequence s is CUA then the 
	method should return the character L corresponding to the amino 
	acid Leucine. A table of codon translations is available in HTML 
	format and in tab-delimited text format. Note: the expected 
	running time of this method should not be linear in the size of 
	the table, i.e., the method should not compare s to each of the 
	64 codons. Instead, the expected running time of the method 
	should be constant. Consider what data structures would be 
	appropriate to store the table in order to support a constant 
	expected running time of this method.*/

	//use a hashtable. look how to store an HTML file in a hashtable 
	public static char translateCodon(String s){

		//create hashtable to reference 
		fileToHashtable("translation.txt");

		//return value at that key	`
		return codons.get(s).charAt(0);

	}

	/*Assuming the length of the input sequence s is divisible by 
	3, the method returns a translated version of the sequence, 
	i.e., every 3 nucleotides in s are translated into an amino acid, 
	and the sequence of amino acids is returned. For example, if s 
	represents the sequence AUGGCCUUUCGAUAG then the method should 
	return the string MAFR*.*/

	public static String translateORF(String s){

		fileToHashtable("translation.txt");

		//create a stringbuilder to append the translation onto
		StringBuilder translation = new StringBuilder();

		//number of amino acids there should be 
		int numOfSubStrings = s.length()/3;

		//loop through and translate string
		for (int i = 0; i < numOfSubStrings ; i++ ) {

			//translate the codon using the hashtable and then add to the string builder
			translation.append(codons.get(s.substring(i * 3, (i + 1) * 3)));
			
		}

		//get rid of all white space from string when returning it
		return translation.toString().replaceAll("\\s+","");

	}

	public static void main(String[] args) {

		/*TESTING*/
		
		// SequenceOps test = new SequenceOps();

		// String yeastFile = "yeastGenome.txt";
		// String hemoglobinFile = "hemoglobin.txt";
		// String humanC22 = "humanChromosome22.txt";
		// String ecoli = "ecoli.txt";

		// String hemoglobin = "ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA";

		// System.out.println(test.sequenceFromFastaFile(hemoglobinFile));
		// System.out.println("\n");

		// System.out.println(test.sequenceToFastaFormat(hemoglobin));

		// System.out.println("The complement of hemoglobin: " + test.complement(hemoglobin));
		// System.out.println("The reverse complement of hemoglobin: " + test.reverseComplement(hemoglobin));

		// System.out.println("the GC content of ACGGACTGC is: " + test.GC_content("ACGGACTGC"));

		// System.out.println(test.GC_content(test.sequenceFromFastaFile(hemoglobinFile)));

		// System.out.println(test.randomPermutation("i love katherine so much i just hope that this permutes the way i want it to"));

		System.out.println(translateORF("AUGGCCUUUCGAUAGGG"));


	}
	

}