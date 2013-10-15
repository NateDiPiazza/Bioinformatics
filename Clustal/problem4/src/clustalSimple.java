import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

/**
 * Problem 4 HW#2 CS 576
 * ClustalW like Multiple Sequence Alignment program
 * @author Nate DiPiazza
 *
 */
public class clustalSimple {
  /**
   * @param args
   */
  public static void main(String[] args) throws IOException {
    int numSeq = 0;
    String finalMSA = "";
    /*
     * array to index the children of all internal nodes
     *-1 = root or uninitialized, else ArrayList<Sequence> sequences
     *used as id in the following arrays
    */
    int[] internal;
    //child arrays contain indexed 'pointers' to parent
    int[] childI;
    int[] childJ;
    //matrix used to store the distance values
    //Note: removed or invalid entries are designated to infinity
    double[][] distanceMatrix;
    ArrayList<Sequence> sequences = new ArrayList<Sequence>();
    //String filename = "sequence_set.txt"; //"fasta1.txt"; //
    if(args.length != 1){
      System.out.println("Usage: java clustalSimple <FASTA filename>");
      System.exit(-1);
    }
    String filename = args[0];
      FileReader fr1 = new FileReader(filename);
      BufferedReader br = new BufferedReader(fr1);
      String line = "";
      //File parser block/////////////////////////////////////////////////////
      while ((line = br.readLine()) != null) {
        //FASTA flag
        if(line.charAt(0) == '>') {
          String seqName = line.substring(1);
          Sequence sequence = new Sequence(seqName);
          numSeq++;
          line = br.readLine(); //get actual sequence string
          sequence.setAlignment(line);
          sequences.add(sequence);
        } else {
          System.out.println("Error: file is not in FASTA format");
          System.exit(-1);
        } 
      }
      ////////////////////////////////////////////////////////////////////////
      int numLeaves = 2*numSeq-1; // 2n-1
      distanceMatrix = new double[numLeaves][numLeaves];
      //initialize all distances to infinity
      for (int i = 0; i < distanceMatrix.length; i++) {
        for (int j = 0; j < distanceMatrix.length; j++) {
            distanceMatrix[i][j] = Double.POSITIVE_INFINITY;
        }
      }
      internal = new int[numLeaves];
      childI = new int[numLeaves];
      childJ = new int[numLeaves];
      double min = Double.POSITIVE_INFINITY;
      int minI = -1;
      int minJ = -1;
      //set unused flag as -1
      for (int i = 0; i < internal.length; i++) {
        internal[i] = -1;
        childI[i] = -1;
        childJ[i] = -1;
      }
      for (int i = 0; i < numSeq; i++) {
        for (int j = 0; j < numSeq; j++) {
            if(i != j){
              String[] stringAlign = getGlobalAlignment(sequences.get(i).getAlignment(), sequences.get(j).getAlignment());
              //String[] stringAlign = getGlobalAlignment(sequenceList.get(i), sequenceList.get(j));
              distanceMatrix[i][j] = getDistance(stringAlign[0], stringAlign[1]);
              //if min update running totals
              if(distanceMatrix[i][j] < min){
                min = distanceMatrix[i][j];
                minI = i;
                minJ = j;
              }
            }
        }
      }
      /*Now that initial distances have been determined, loop through
       *UPGMA algorithm until the guide-tree has been build.
       *The sequence object and parent,childI, and childJ can be used to MSA it.*/
      while(numSeq < numLeaves){
          //internal node printout format; ex: s1-s2 s1
          String newName = sequences.get(minI).getName()+ "-" + sequences.get(minJ).getName();
          String result = printMSA(sequences.get(minI).getAlignment(), sequences.get(minJ).getAlignment());
          finalMSA = result;
          Sequence newSeq = new Sequence(newName);
          newSeq.setAlignment(result);
          sequences.add(newSeq);
          //remove winner from matrix
          distanceMatrix[minI][minJ] = Double.POSITIVE_INFINITY;
          distanceMatrix[minJ][minI] = Double.POSITIVE_INFINITY;
          //print output to console
          /*System.out.printf("%s %s\n", sequences.get(numSeq).getName(), sequences.get(minI).getName());
          System.out.printf("%s %s\n", sequences.get(numSeq).getName(), sequences.get(minJ).getName());
          System.out.printf("%s", result);*/
          internal[minI] = numSeq; //numSeq == id for new node
          internal[minJ] = numSeq;
          childI[numSeq] = minI;
          childJ[numSeq] = minJ;
          //update distanceMatrix
          for (int i = 0; i < numSeq; i++) {
            if(internal[i] == -1){
                //System.out.println("need to update " + distanceMatrix[minJ][i] + " " + distanceMatrix[minI][i]);
                int jElem = getNumElements(numSeq, childI);
                int iElem = getNumElements(numSeq, childJ);
                //distanceMatrix[i][numSeq] = (iElem*distanceMatrix[i][minI] + jElem*distanceMatrix[i][minJ])/(iElem+jElem);
                distanceMatrix[i][numSeq] = (iElem*distanceMatrix[minI][i] + jElem*distanceMatrix[minJ][i])/(iElem+jElem);
                //remove components of winning internal node from distance matrix
                distanceMatrix[i][minI] = Double.POSITIVE_INFINITY;
                distanceMatrix[i][minJ] = Double.POSITIVE_INFINITY;
                distanceMatrix[minI][i] = Double.POSITIVE_INFINITY;
                distanceMatrix[minJ][i] = Double.POSITIVE_INFINITY;
            }
          }
          if(anyNodes(distanceMatrix) == false){
            break;
          }
          //now find another clustering candidate
          numSeq++;
          min = Double.POSITIVE_INFINITY;
          for (int i = 0; i < numSeq; i++) {
            for (int j = 0; j < numSeq; j++) {
                if(distanceMatrix[i][j] < min){
                    min = distanceMatrix[i][j];
                    minI = i;
                    minJ = j;
                  }
            }
          }
      } //end while
      
      //guide tree built. use it to create a MSA
      /*String test1 = "AGC\nACC";
      String test2 = "AAAC\n-GAC";
      String result = printMSA(test1, test2);*/
      finalMSA = finalMSA.substring(0, finalMSA.length()-2);
      System.out.println(finalMSA);
  } //end main
  
  public static double getDistance(String sequence1, String sequence2) {
    double score = 0.0;
    char[] seq1 = sequence1.toCharArray();
    char[] seq2 = sequence2.toCharArray();
    for (int i = 0; i < seq2.length; i++) {
      if(seq1[i] == seq2[i])
        score++;
    }
    double similarity = (score/(double) seq1.length);
    return 1-similarity;
  }
  
  public static boolean anyNodes(double[][] distMatrix){
    for (int i = 0; i < distMatrix.length; i++) {
      for (int j = 0; j < distMatrix[0].length; j++) {
        double dist = distMatrix[i][j];
        if(dist != Double.POSITIVE_INFINITY)
          return true;
      }
    }
    return false;
  }
  
  public static int getNumElements(int index, int[] child){
      int numElem = 0;
      index = child[index];
      while(index != -1){
          index = child[index];
          numElem++;
      }
      return numElem;
  }
  
  /**
   * Compute the pairwise global sequence 
   * alignment using n.w. dynamic programming.
   * The method then builds the alignment strings
   * by back tracing path pointers.
   * The result is printed to the console.
   * @param seq1
   * @param seq2
   * @return 
   * @return void
   */
  private static String[] getGlobalAlignment(String seq1, String seq2){
    int xSize = seq1.length() + 1;
    int ySize = seq2.length() + 1;
    SequenceMatrix[][] seqMatrix = new SequenceMatrix[ySize][xSize];
    //Initialize first row and column of matrix to zero
    for (int i = 0; i < seqMatrix.length; i++) {
      for (int j = 0; j < seqMatrix[0].length; j++) {
        SequenceMatrix entry = new SequenceMatrix(0);
        seqMatrix[i][j] = entry;
      }
    }
    int max = 0;
    int maxI = 0;
    int maxJ = 0;
    //Fill in rest of matrix from top to bottom, left to right
    for (int i = 1; i < seqMatrix.length; i++) {
      for (int j = 1; j < seqMatrix[0].length; j++) {
        /**
         * Recurrence:
         * F(j,i) = max{ 
         *   F(j-1,i-1)+s(Xi,Yj), 
         *   F(j-1,i)-d, 
         *   F(j,i-1)-d }
         */
        int mVal = 0;
        int gap2 = 0; 
        int gap1 = 0;
        //System.out.printf("current [%d,%d], ", i, j);
        mVal = seq1.charAt(j-1) == seq2.charAt(i-1) ? seqMatrix[i-1][j-1].score + 2 : seqMatrix[i-1][j-1].score + 1;
        gap1 = seqMatrix[i-1][j].score;
        gap2 = seqMatrix[i][j-1].score;
        //tie goes to mVal
        if(mVal >= gap1 && mVal >= gap2){
          //System.out.printf("%c and %c\n", seq1.charAt(j-1), seq2.charAt(i-1));
          max = mVal;
          maxI = i-1;
          maxJ = j-1;
        }
          else if (gap1 >= gap2) {
            max = gap1;
          //maxI = i;
          //maxJ = j-1;
            maxI = i-1;
          maxJ = j;
        }
        else {
          max = gap2;
          //maxI = i-1;
          //maxJ = j;
          maxI = i;
          maxJ = j-1;
        }
        seqMatrix[i][j].score = max;
        seqMatrix[i][j].setBackPointer(maxI, maxJ);
        
      }
      //System.out.println("");
    }
    //dumpMatrix(seqMatrix);
    //System.out.println();
    return buildAlignment(seqMatrix, seq1, seq2);
  }
  
  public static String[] buildAlignment(SequenceMatrix[][] seqMatrix, String seq1, String seq2){
    char[] seq1Char = seq1.toCharArray();
    char[] seq2Char = seq2.toCharArray();
    String tempSeq1 = "";
    String tempSeq2 = "";
    int i = seqMatrix.length-1;
    int j = seqMatrix[0].length-1;
    //Trace back from F(m, n) to F(0, 0) to recover alignment
    while(true){
      //System.out.printf("Pointer: [%d,%d], Index: %d,%d\n", seqMatrix[i][j].getBackPointer()[0], seqMatrix[i][j].getBackPointer()[1], i,j);
      if(seqMatrix[i][j].getBackPointer()[0] < 0 || seqMatrix[i][j].getBackPointer()[1] < 0)
        break;
      //Gap1
      if( (i > seqMatrix[i][j].getBackPointer()[0]) && (j == seqMatrix[i][j].getBackPointer()[1]) ){
        tempSeq1 = tempSeq1 + "-";
        tempSeq2 = tempSeq2 + seq2Char[i-1];
        i = seqMatrix[i][j].getBackPointer()[0];  
      }
      //Gap2
      else if( (i == seqMatrix[i][j].getBackPointer()[0]) && (j > seqMatrix[i][j].getBackPointer()[1]) ){
        tempSeq2 = tempSeq2 + "-";
        tempSeq1 = tempSeq1 + seq1Char[j-1];
        //System.out.printf("Gap2: %c\n", seq1Char[j-1]);
        j = seqMatrix[i][j].getBackPointer()[1];
      }
      else {
        //System.out.printf("Match: %c\n", seq1Char[j-1]);
        tempSeq2 = tempSeq2 + seq2Char[i-1];
        tempSeq1 = tempSeq1 + seq1Char[j-1];
        int temp = seqMatrix[i][j].getBackPointer()[0];
        j = seqMatrix[i][j].getBackPointer()[1];
        i = temp;
      }
      if(seqMatrix[i][j].getBackPointer()[0] == null || seqMatrix[i][j].getBackPointer()[1] == null){
        //if not at [0,0] fill in the last gap
        //if not at [0,0] fill in the last gap
        if ((j == 0) && (i > 0)) {
          tempSeq1 = tempSeq1 + "-";
          tempSeq2 = tempSeq2 + seq2Char[i-1];
        }
        if ((i == 0) && (j > 0)) {
          tempSeq2 = tempSeq2 + "-";
          tempSeq1 = tempSeq1 + seq1Char[j-1];
        }
        break;
      }
        
    }
    // strings are in reverse order call method to reverse and print out to the console
    //System.out.printf("%2s\n%2s", reverseSeq(tempSeq1), reverseSeq(tempSeq2));
    //return reverseSeq(tempSeq1) + "\n" + reverseSeq(tempSeq2);
    String[] align = {reverseSeq(tempSeq1), reverseSeq(tempSeq2)};
    return align;
  }
  /*
   * helper method used to
   * reverse sequence
   * @param seq
   * @return reversed string
   */
  public static String reverseSeq(String seq){
    char[] seqArray = seq.toCharArray();
    char[] revArray = new char[seq.toCharArray().length];
    int j = seqArray.length-1;
    for (int i = 0; i < seqArray.length; i++) {
      char c = seqArray[i];
      revArray[j--] = c;
      
    }
    return String.valueOf(revArray);
  }
  
  public static String printMSA(String[] seq1, String[] seq2){
    Profile profile1 = new Profile();
    Profile profile2 = new Profile();
    for (int i = 0; i < seq1.length; i++) {
      profile1.add(seq1[i]);
    }
    for (int i = 0; i < seq2.length; i++) {
      profile2.add(seq2[i]);
    }
    String msa = getMultipleAlignment(profile1, profile2);
    return msa;
  }
  //single alignment case
  public static String printMSA(String seq1, String seq2){
    Profile profile1 = new Profile();
    Profile profile2 = new Profile();
    String[] sSeq1 = seq1.split("\n");
    String[] sSeq2 = seq2.split("\n");
    for (int i = 0; i < sSeq1.length; i++) {
      profile1.add(sSeq1[i]);
    }
    for (int i = 0; i < sSeq2.length; i++) {
      profile2.add(sSeq2[i]);
    }
    String msa = getMultipleAlignment(profile1, profile2);
    return msa;
    }
  
  private static String getMultipleAlignment(Profile profile1, Profile profile2) {
    int xSize = profile1.length() + 1;
    int ySize = profile2.length() + 1;
    MultipleSequence[][] seqMatrix = new MultipleSequence[ySize][xSize];
    // Initialize first row and column of matrix to zero
    for (int i = 0; i < seqMatrix.length; i++) {
      for (int j = 0; j < seqMatrix[0].length; j++) {
        MultipleSequence entry = new MultipleSequence(0);
        seqMatrix[i][j] = entry;
      }
    }
    double max = 0;
    int maxI = 0;
    int maxJ = 0;
    // Fill in rest of matrix from top to bottom, left to right
    for (int i = 1; i < seqMatrix.length; i++) {
      for (int j = 1; j < seqMatrix[0].length; j++) {
        double mVal = 0;
        double gap2 = 0;
        double gap1 = 0;
        mVal = profile1.getCol(j - 1).equals(profile2.getCol(i - 1)) ? seqMatrix[i - 1][j - 1].score
            + seqMatrix[i - 1][j - 1].getColScore(profile1.getCol(j - 1),
                profile2.getCol(i - 1))
            : seqMatrix[i - 1][j - 1].score
                + seqMatrix[i - 1][j - 1].getColScore(profile1.getCol(j - 1),
                    profile2.getCol(i - 1));
        // gap cost = 0
        gap1 = seqMatrix[i - 1][j].score;
        gap2 = seqMatrix[i][j - 1].score;
        // tie goes to mVal
        if (mVal >= gap1 && mVal >= gap2) {
          max = mVal;
          maxI = i - 1;
          maxJ = j - 1;
        } else if (gap1 >= gap2) {
          max = gap1;
          maxI = i - 1;
          maxJ = j;
        } else {
          max = gap2;
          maxI = i;
          maxJ = j - 1;
        }
        seqMatrix[i][j].score = max;
        seqMatrix[i][j].setBackPointer(maxI, maxJ);

      }
    }
    String result = printMultipleAlignment(seqMatrix, profile1, profile2);
    return result;
  }
  
  public static String printMultipleAlignment(MultipleSequence[][] seqMatrix, Profile profile1,
      Profile profile2) {
    String result = "";
    Profile alignPro1 = new Profile(colLen(profile1));
    Profile alignPro2 = new Profile(colLen(profile2));
    int i = seqMatrix.length - 1;
    int j = seqMatrix[0].length - 1;
    // Trace back from F(m, n) to F(0, 0) to recover alignment
    while (true) {
      if (seqMatrix[i][j].getBackPointer()[0] < 0
          || seqMatrix[i][j].getBackPointer()[1] < 0)
        break;
      // Gap1
      if ((i > seqMatrix[i][j].getBackPointer()[0])
          && (j == seqMatrix[i][j].getBackPointer()[1])) {
        String gap = getGap(colLen(profile1));
        alignPro1.updateAlign(gap);
        alignPro2.updateAlign(profile2.getCol(i - 1));
        i = seqMatrix[i][j].getBackPointer()[0];
      }
      // Gap2
      else if ((i == seqMatrix[i][j].getBackPointer()[0])
          && (j > seqMatrix[i][j].getBackPointer()[1])) {
        String gap = getGap(colLen(profile2));
        alignPro2.updateAlign(gap);
        alignPro1.updateAlign(profile1.getCol(j - 1));
        j = seqMatrix[i][j].getBackPointer()[1];
      } else {
        alignPro2.updateAlign(profile2.getCol(i - 1));
        alignPro1.updateAlign(profile1.getCol(j - 1));
        int temp = seqMatrix[i][j].getBackPointer()[0];
        j = seqMatrix[i][j].getBackPointer()[1];
        i = temp;
      }
      if (seqMatrix[i][j].getBackPointer()[0] == null
          || seqMatrix[i][j].getBackPointer()[1] == null) {
        // if not at [0,0] fill in the last gap
        if ((j == 0) && (i > 0)) {
          String gap = getGap(colLen(profile1));
          alignPro1.updateAlign(gap);
          alignPro2.updateAlign(profile2.getCol(i - 1));
        }
        if ((i == 0) && (j > 0)) {
          String gap = getGap(colLen(profile2));
          alignPro2.updateAlign(gap);
          alignPro1.updateAlign(profile1.getCol(j - 1));
        }
        //leave loop
        break;
      }

    } //end infinite loop 
    for (i = 0; i < alignPro1.newAlign.length; i++) {
      //System.out.println(reverseSeq(alignPro1.newAlign[i]));
      result += reverseSeq(alignPro1.newAlign[i]) + "\n";
    }
    for (i = 0; i < alignPro2.newAlign.length; i++) {
      //System.out.println(reverseSeq(alignPro2.newAlign[i]));
      result += reverseSeq(alignPro2.newAlign[i]) + "\n";
    }
    return result;
  } // end print alignment
  
  /**
   * build gap string which is the
   * number of the profile's column
   * @param colLen
   * @return gap
   */
  public static String getGap(int colLen) {
    String gap = "";
    for (int i = 0; i < colLen; i++) {
      gap = gap + "-";
    }
    return gap;
  }
  //assumption is columns are equal
  public static int colLen(Profile p) {
    return p.getCol(0).length();
  }
  
  public static double maxScore(MultipleSequence[][] seqMatrix){
    return seqMatrix[seqMatrix.length-1][seqMatrix[0].length-1].score;
  }


}
