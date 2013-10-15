import java.util.ArrayList;
import java.util.Iterator;


/**
 * @author Owner
 *
 */
public class Profile extends ArrayList<String> {
  
  /**
   * 
   */
  private static final long serialVersionUID = 1L;
  String[] newAlign;
  
  public Profile(){}
  
  public Profile(int initCap) {
    newAlign = new String[initCap];
    for (int i = 0; i < newAlign.length; i++) {
      newAlign[i] = "";
    }
  }
  /**
   * @param profile column index [0..n]
   * @return transposed column as a String
   */
  public String getCol(int colNum){
     String col = "";
     Iterator<String> proIter = super.iterator();
     while (proIter.hasNext()) {
      String sequence = proIter.next();
      //concatenate each column entry
      col = col + sequence.substring(colNum, colNum+1);
     }
    return col;
  }
  //override 
  public int length(){
   int len = super.get(0).length();
   return len;
  }
  
  public void updateAlign(String entry){
    char[] splitEntry = entry.toCharArray();
    for (int i = 0; i < splitEntry.length; i++) {
      newAlign[i] = newAlign[i] + splitEntry[i];
    }
  }
  
  
  
}
