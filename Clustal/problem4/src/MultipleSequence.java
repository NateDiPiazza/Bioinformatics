
public class MultipleSequence {
  double score = 0;
  Integer[] backPointer = new Integer[2];
  //constructor
  public MultipleSequence(double score) {
  this.score = score;
  }
  
  public void setBackPointer(int i, int j){
  backPointer[0] = i;
  backPointer[1] = j;
  }

  public Integer[] getBackPointer() {
  return backPointer;
  }
  /**
   * Compute the score average of 
   * all pairwise alignments given
   * a pair of sequence columns
   * @param col1
   * @param col2
   * @return averaged score
   */
  public double getColScore(String col1, String col2){
    double score = 0;
    int total = 0;
    char[] col1Arr = col1.toCharArray();
    char[] col2Arr = col2.toCharArray();
    for (int i = 0; i < col1Arr.length; i++) {
      char letter1 = col1Arr[i];
      for (int j = 0; j < col2Arr.length; j++) {
        char letter2 = col2Arr[j];
        total++;
        //gap penalty = zero
        if( !((letter1 == '-') || (letter2 == '-')) ){
          if(letter1 == letter2){
            score += 2; //Match score
          } else {
            score++; //Mismatch score
          }
        }
      }
    }
    
    return (score/(double) total);
  }
  
  
}
