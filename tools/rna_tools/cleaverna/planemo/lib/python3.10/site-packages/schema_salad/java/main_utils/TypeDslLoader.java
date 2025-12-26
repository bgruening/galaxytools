package ${package}.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class TypeDslLoader<T> implements Loader<T> {
  private final Loader<T> innerLoader;
  private final Integer refScope;
  private static final Pattern TYPE_DSL_REGEX = Pattern.compile("^([^\\[?]+)(\\[\\])?(\\?)?$");

  public TypeDslLoader(final Loader<T> innerLoader, final Integer refScope) {
    this.innerLoader = innerLoader;
    this.refScope = refScope;
  }

  private Object resolve(
      final String doc_, final String baseUri, final LoadingOptions loadingOptions) {
    final Matcher m = TYPE_DSL_REGEX.matcher(doc_);
    if (m.matches()) {
      final String first =
          loadingOptions.expandUrl(m.group(1), baseUri, false, true, this.refScope);
      Object second = null;
      Object third = null;
      if (m.group(2) != null && m.group(2).length() > 0) {
        HashMap<String, Object> resolvedMap = new HashMap<String, Object>();
        resolvedMap.put("type", "array");
        resolvedMap.put("items", first);
        second = resolvedMap;
      }
      if (m.group(3) != null && m.group(3).length() > 0) {
        third = Arrays.asList("null", second != null ? second : first);
      }
      if (third != null) {
        return third;
      } else {
        return second != null ? second : first;
      }
    } else {
      return doc_;
    }
  }

  public T load(
      final Object doc_,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    Object doc = doc_;
    if (doc instanceof List) {
      final List<Object> docList = (List<Object>) doc;
      final List<Object> r = new ArrayList<Object>();
      for (final Object d : docList) {
        if (d instanceof String) {
          Object resolved = this.resolve((String) d, baseUri, loadingOptions);
          if (resolved instanceof List) {
            List<Object> resolvedList = (List<Object>) resolved;
            for (Object i : resolvedList) {
              if (!r.contains(i)) {
                r.add(i);
              }
            }
          } else {
            if (!r.contains(resolved)) {
              r.add(resolved);
            }
          }
        } else {
          r.add(d);
        }
      }
      doc = docList;
    } else if (doc instanceof String) {
      doc = this.resolve((String) doc, baseUri, loadingOptions);
    }
    return this.innerLoader.load(doc, baseUri, loadingOptions);
  }
}
