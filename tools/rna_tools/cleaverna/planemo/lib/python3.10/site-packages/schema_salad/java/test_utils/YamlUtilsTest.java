package ${package}.utils;

import java.util.Map;
import org.junit.Assert;
import org.junit.Test;

public class YamlUtilsTest {
  @Test
  public void testSimpleLoad() {
    final String yamlStr = "moo: cow\nbark: dog\n";
    final Map<String, Object> loaded = YamlUtils.mapFromString(yamlStr);
    Assert.assertEquals(loaded.get("moo"), "cow");
    Assert.assertEquals(loaded.get("bark"), "dog");
  }
}
