package ${package}.utils;

public class ExpressionLoader implements Loader<String> {

  public ExpressionLoader() {
  }

  public String load(
      final Object doc_,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    if (doc_ instanceof String) {
	    return (String) doc_;
    } else {
        throw new ValidationException("Expected a string.");
    }
  }
}
