package ${package}.utils;

public class NullLoader implements Loader<Object> {

  public Object load(
      final Object doc,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    if (doc != null) {
      throw new ValidationException("Expected null");
    }
    return doc;
  }
}
