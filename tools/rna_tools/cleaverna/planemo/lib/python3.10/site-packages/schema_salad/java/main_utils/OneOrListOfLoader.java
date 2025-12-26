package ${package}.utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Optional;


public class OneOrListOfLoader<T> implements Loader<OneOrListOf<T>> {
  private final Loader<T> oneLoader;
  private final Loader<List<T>> listLoader;

  public OneOrListOfLoader(Loader<T> oneLoader, Loader<List<T>> listLoader) {
      this.oneLoader = oneLoader;
      this.listLoader = listLoader;
  }

  public OneOrListOf<T> load(
      final Object doc,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    final List<ValidationException> errors = new ArrayList();
    try {
      return OneOrListOf.oneOf(this.oneLoader.load(doc, baseUri, loadingOptions, docRoot));
    } catch (ValidationException e) {
      errors.add(e);
    }
    try {
      return OneOrListOf.listOf(this.listLoader.load(doc, baseUri, loadingOptions, docRoot));
    } catch (ValidationException e) {
      errors.add(e);
    }
    throw new ValidationException("Failed to one or list of of type", errors);
  }
}
