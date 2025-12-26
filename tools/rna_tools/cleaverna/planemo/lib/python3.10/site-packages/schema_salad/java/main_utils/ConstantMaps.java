package ${package}.utils;

import java.util.HashMap;

public class ConstantMaps {
  // declare as HashMap for clone().
  public static final HashMap<String, String> vocab = new HashMap();
  public static final HashMap<String, String> rvocab = new HashMap();

  static {
${vocab}

${rvocab}
  }
}
