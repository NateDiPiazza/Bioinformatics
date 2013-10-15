import java.util.ArrayList;


public class Sequence {
  private String name;
  private String alignment;
  public boolean isLeaf;
 public ArrayList<String> msAlign = new ArrayList<String>();
  
  public Sequence(String name) {
    this.name = name;
  }
  
  public String getAlignment() {
    return alignment;
  }

  public void setAlignment(String alignment) {
    this.alignment = alignment;
  }

  public String getName() {
    return name;
  }

  public void setName(String name) {
    this.name = name;
  }
  
}
