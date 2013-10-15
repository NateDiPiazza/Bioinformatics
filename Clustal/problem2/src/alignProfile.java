import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class alignProfile {
  /**
   * Compute the sum of scores for two alignments
   * @param args
   * @throws IOException
   */
  public static void main(String[] args) throws IOException {
    //String fileName1 = "pro4.txt";
    if((args.length < 2) || (args.length > 2)){
      System.out.println("Usage: java alignProfile <profile1 filename> <profile2 filename>");
      System.exit(-1);
    }
    String fileName1 = args[0];
    String fileName2 = args[1];
    FileReader fr1 = new FileReader(fileName1);
    BufferedReader br = new BufferedReader(fr1);
    Profile profile1 = new Profile();
    String sequence = "";
    while ((sequence = br.readLine()) != null) {
      profile1.add(sequence);
    }
    //String fileName2 = "pro5.txt";
    FileReader fr2 = new FileReader(fileName2);
    br = new BufferedReader(fr2);
    Profile profile2 = new Profile();
    while ((sequence = br.readLine()) != null) {
      profile2.add(sequence);
    }
    br.close();
    getGlobalAlignment(profile1, profile2);
  }

  private static void getGlobalAlignment(Profile profile1, Profile profile2) {
    int xSize = profile1.length() + 1;
    int ySize = profile2.length() + 1;
    Sequence[][] seqMatrix = new Sequence[ySize][xSize];
    // Initialize first row and column of matrix to zero
    for (int i = 0; i < seqMatrix.length; i++) {
      for (int j = 0; j < seqMatrix[0].length; j++) {
        Sequence entry = new Sequence(0);
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
    printAlignment(seqMatrix, profile1, profile2);
  }

  public static void printAlignment(Sequence[][] seqMatrix, Profile profile1,
      Profile profile2) {
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
      System.out.println(reverseSeq(alignPro1.newAlign[i]));
    }
    for (i = 0; i < alignPro2.newAlign.length; i++) {
      System.out.println(reverseSeq(alignPro2.newAlign[i]));
    }
    
  } // end print alignment
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
  /**
   * Print out the scoring matrix
   * @param seqMatrix
   */
  public static void dumpMatrix(Sequence[][] seqMatrix) {
    for (int i = 0; i < seqMatrix.length; i++) {
      for (int j = 0; j < seqMatrix[0].length; j++) {
        System.out.printf("%6.3f ", seqMatrix[i][j].score);
      }
      System.out.println();
    }
  }
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
  
  public static double maxScore(Sequence[][] seqMatrix){
    return seqMatrix[seqMatrix.length-1][seqMatrix[0].length-1].score;
  }

}
