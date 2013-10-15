//used to create the DP alignment matrix
public class SequenceMatrix {
  
  String seq1;
  String seq2;
  int score = 0;
  Integer[] backPointer = new Integer[2];
  //constructor
  public SequenceMatrix(int score) {
    //this.seq1 = seq1;
  //this.seq2 = seq2;
  this.score = score;
  }
  
  public void setBackPointer(int i, int j){
  backPointer[0] = i;
  backPointer[1] = j;
  }

  public Integer[] getBackPointer() {
  return backPointer;
  }
  
  
}
