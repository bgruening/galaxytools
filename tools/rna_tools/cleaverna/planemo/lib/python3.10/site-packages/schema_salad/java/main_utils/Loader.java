package ${package}.utils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public interface Loader<T> {

  T load(
      final Object doc,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot);

  default T load(final Object doc, final String baseUri, final LoadingOptions loadingOptions) {
    return load(doc, baseUri, loadingOptions, null);
  }

  default T documentLoad(
      final String doc, final String baseUri, final LoadingOptions loadingOptions) {
    return load(doc, baseUri, loadingOptions);
  }

  default T documentLoad(
      final Map<String, Object> doc_, final String baseUri_, final LoadingOptions loadingOptions_) {
    Map<String, Object> doc = doc_;
    LoadingOptions loadingOptions = loadingOptions_;
    if (doc.containsKey("$namespaces")) {
      final Map<String, String> namespaces = (Map<String, String>) doc.get("$namespaces");
      loadingOptions =
          new LoadingOptionsBuilder().copiedFrom(loadingOptions).setNamespaces(namespaces).build();
      doc = copyWithoutKey(doc, "$namespaces");
    }
    String baseUri = baseUri_;
    if (doc.containsKey("$base")) {
      baseUri = (String) doc.get("$base");
    }
    if (doc.containsKey("$graph")) {
      return load(doc.get("$graph"), baseUri, loadingOptions);
    } else {
      return load(doc, baseUri, loadingOptions, baseUri);
    }
  }

  default T documentLoad(
      final List<Object> doc, final String baseUri, final LoadingOptions loadingOptions) {
    return load(doc, baseUri, loadingOptions);
  }

  default T documentLoadByUrl(final String url, final LoadingOptions loadingOptions) {
    if (loadingOptions.idx.containsKey(url)) {
      Object result = loadingOptions.idx.get(url);
      if (result instanceof String) {
        return documentLoad((String) result, url, loadingOptions);
      } else if (result instanceof Map) {
        return documentLoad((Map<String, Object>) result, url, loadingOptions);
      }
      return load(result, url, loadingOptions);
    }

    final String text = loadingOptions.fetcher.fetchText(url);
    try {
      Map<String, Object> resultMap = YamlUtils.mapFromString(text);
      loadingOptions.idx.put(url, resultMap);
      final LoadingOptionsBuilder urlLoadingOptions =
          new LoadingOptionsBuilder().copiedFrom(loadingOptions).setFileUri(url);
      return documentLoad(resultMap, url, urlLoadingOptions.build());
    } catch (ClassCastException e) {
      List<Object> resultList = YamlUtils.listFromString(text);
      loadingOptions.idx.put(url, resultList);
      final LoadingOptionsBuilder urlLoadingOptions =
          new LoadingOptionsBuilder().copiedFrom(loadingOptions).setFileUri(url);
      return documentLoad(resultList, url, urlLoadingOptions.build());
    }
  }

  default T loadField(
      final Object val_, final String baseUri, final LoadingOptions loadingOptions) {
    Object val = val_;
    if (val instanceof Map) {
      Map<String, Object> valMap = (Map<String, Object>) val;
      if (valMap.containsKey("$import")) {
        if (loadingOptions.fileUri == null) {
          throw new ValidationException("Cannot load $import without fileuri");
        }
        return documentLoadByUrl(
            loadingOptions.fetcher.urlJoin(loadingOptions.fileUri, (String) valMap.get("$import")),
            loadingOptions);
      } else if (valMap.containsKey("$include")) {
        if (loadingOptions.fileUri == null) {
          throw new ValidationException("Cannot load $import without fileuri");
        }
        val =
            loadingOptions.fetcher.fetchText(
                loadingOptions.fetcher.urlJoin(
                    loadingOptions.fileUri, (String) valMap.get("$include")));
      }
    }
    return load(val, baseUri, loadingOptions);
  }

  default Map<String, Object> copyWithoutKey(final Map<String, Object> doc, final String key) {
    final Map<String, Object> result = new HashMap();
    for (final Map.Entry<String, Object> entry : doc.entrySet()) {
      if (!entry.getKey().equals(key)) {
        result.put(entry.getKey(), entry.getValue());
      }
    }
    return result;
  }

  static <T> T validateOfJavaType(final Class<T> clazz, final Object doc) {
    if (!clazz.isInstance(doc)) {
      String className = "null";
      if (doc != null) {
        className = doc.getClass().getName();
      }
      final String message =
          String.format(
              "Expected object with Java type of %s but got %s", clazz.getName(), className);
      throw new ValidationException(message);
    }
    return (T) doc;
  }
}
