package ${package}.utils;

import java.util.Optional;


public class OptionalLoader<T> implements Loader<Optional<T>> {
  private final Loader<T> itemLoader;

  public OptionalLoader(Loader<T> itemLoader) {
    this.itemLoader = itemLoader;
  }

  public Optional<T> load(
      final Object doc,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    if(doc == null) {
      return Optional.empty();
    }
    return Optional.of(itemLoader.load(doc, baseUri, loadingOptions, docRoot));
  }
}
