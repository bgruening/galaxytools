package ${package}.utils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

public class IdMapLoader<T> implements Loader<T> {
  private final Loader<T> innerLoader;
  private final String mapSubject;
  private final String mapPredicate;

  public IdMapLoader(
      final Loader<T> innerLoader, final String mapSubject, final String mapPredicate) {
    this.innerLoader = innerLoader;
    this.mapSubject = mapSubject;
    this.mapPredicate = mapPredicate;
  }

  public T load(
      final Object doc_,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    Object doc = doc_;
    if (doc instanceof Map) {
      final Map<String, Object> docMap = (Map<String, Object>) doc;
      final List<Object> asList = new ArrayList();
      for (final String key : docMap.keySet()) {
        final Object el = docMap.get(key);
        if (el instanceof Map) {
          final Map<String, Object> v2 = new HashMap<String, Object>((Map<String, Object>) el);
          v2.put(this.mapSubject, key);
          asList.add(v2);
        } else {
          if (this.mapPredicate != null) {
            final Map<String, Object> v3 = new HashMap<String, Object>();
            v3.put(this.mapPredicate, el);
            v3.put(this.mapSubject, key);
            asList.add(v3);
          } else {
            throw new ValidationException("No mapPredicate");
          }
        }
      }
      doc = asList;
    }
    return this.innerLoader.load(doc, baseUri, loadingOptions);
  }
}
