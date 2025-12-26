package ${package}.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class UnionLoader implements Loader<Object> {
  private final ArrayList<Loader> alternates;

  public UnionLoader(List<Loader> alternates) {
    this.alternates = new ArrayList<Loader>(alternates);
  }

  public UnionLoader(Loader[] alternates) {
    this(Arrays.asList(alternates));
  }

  public void addLoaders(List<Loader> loaders) {
    this.alternates.addAll(loaders);
  }

  public void addLoaders(Loader[] loaders) {
    this.addLoaders(Arrays.asList(loaders));
  }

  public Object load(
      final Object doc,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    final List<ValidationException> errors = new ArrayList();
    for (final Loader loader : this.alternates) {
      try {
        return loader.load(doc, baseUri, loadingOptions, docRoot);
      } catch (ValidationException e) {
        errors.add(e);
      }
    }
    throw new ValidationException("Failed to match union type", errors);
  }
}
