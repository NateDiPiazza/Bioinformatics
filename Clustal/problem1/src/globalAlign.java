import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

/**
 * <h2>Needleman-Wunsch Global Alignment program</h2>
 * <br/>
 * Scoring: <li>match=2</li> <li>mismatch=1</li> <li>gap=0</li>
 * <p>
 * <br/>
 * Tiebreaker ranking:  match/mismatch, gap1, gap2
 * Arbitrarily chosen tiebreakers.
 * </p><br/>
 * CS 576 UW
 * Instructor: Sushmita Roy
 * @param args
 * @throws IOException
 * @author Nate DiPiazza
 */
public class globalAlign {

  public static void main(String[] args) throws IOException {
    //String fileName1 = "seq1.txt";
    //String fileName2 = "seq2.txt";
    if( (args.length < 2) || (args.length > 2) ){
      System.out.println("Usage: java globalAlign <sequence 1 filename> <sequence 2 filename>");
      System.exit(-1);
    }
    FileReader fr1 = new FileReader(fileName1);
    BufferedReader br = new BufferedReader(fr1);
    //should only be one line, so no need for while loop
    String seq1 = br.readLine();
    FileReader fr2 = new FileReader(fileName2);
    br = new BufferedReader(fr2);
    String seq2 = br.readLine();
    getGlobalAlignment(seq1, seq2);
    br.close();
  }
  
  /**
   * Compute the pairwise global sequence 
   * alignment using n.w. dynamic programming.
   * The method then builds the alignment strings
   * by back tracing path pointers.
   * The result is printed to the console.
   * @param seq1
   * @param seq2
   * @return void
   */
  private static void getGlobalAlignment(String seq1, String seq2){
    int xSize = seq1.length() + 1;
    int ySize = seq2.length() + 1;
    Sequence[][] seqMatrix = new Sequence[ySize][xSize];
    //Initialize first row and column of matrix to zero
    for (int i = 0; i < seqMatrix.length; i++) {
      for (int j = 0; j < seqMatrix[0].length; j++) {
        Sequence entry = new Sequence(0);
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
         *   F(j,i) = max{ 
         *   F(j-1,i-1)+s(Xi,Yj), 
         *   F(j-1,i)-d, 
         *   F(j,i-1)-d }
         */
        int mVal = 0;
        int gap2 = 0; 
        int gap1 = 0;
        mVal = seq1.charAt(j-1) == seq2.charAt(i-1) ? seqMatrix[i-1][j-1].score + 2 : seqMatrix[i-1][j-1].score + 1;
        gap1 = seqMatrix[i-1][j].score;
        gap2 = seqMatrix[i][j-1].score;
        //tie goes to mVal
        if(mVal >= gap1 && mVal >= gap2){
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
    }
    //dumpMatrix(seqMatrix); //debug
    //System.out.println();
    printAlignment(seqMatrix, seq1, seq2);
  }
  /**
   * Print out the scoring matrix
   * @param seqMatrix
   */
  public static void dumpMatrix(Sequence[][] seqMatrix){
    for (int i = 0; i < seqMatrix.length; i++) {
      for (int j = 0; j < seqMatrix[0].length; j++) {
        System.out.printf("%4d ", seqMatrix[i][j].score);
      }
      System.out.println();
    }
  }
  /**
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
  
  public static void printAlignment(Sequence[][] seqMatrix, String seq1, String seq2){
    char[] seq1Char = seq1.toCharArray();
    char[] seq2Char = seq2.toCharArray();
    String tempSeq1 = "";
    String tempSeq2 = "";
    int i = seqMatrix.length-1;
    int j = seqMatrix[0].length-1;
    //Trace back from F(m, n) to F(0, 0) to recover alignment
    while(true){
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
        j = seqMatrix[i][j].getBackPointer()[1];
      }
      else {
        tempSeq2 = tempSeq2 + seq2Char[i-1];
        tempSeq1 = tempSeq1 + seq1Char[j-1];
        int temp = seqMatrix[i][j].getBackPointer()[0];
        j = seqMatrix[i][j].getBackPointer()[1];
        i = temp;
      }
      if(seqMatrix[i][j].getBackPointer()[0] == null || seqMatrix[i][j].getBackPointer()[1] == null){
        //if not at [0,0] fill in the last gap
        if(j == 0){
          tempSeq1 = tempSeq1 + "-";
          tempSeq2 = tempSeq2 + seq2Char[i-1];
        }
        if(i == 0){
          tempSeq2 = tempSeq2 + "-";
          tempSeq1 = tempSeq1 + seq1Char[j-1];
        }
        break;
      }
        
    }
    // strings are in reverse order call method to reverse and print out to the console
    System.out.printf("%2s\n%2s", reverseSeq(tempSeq1), reverseSeq(tempSeq2));
  }

}
