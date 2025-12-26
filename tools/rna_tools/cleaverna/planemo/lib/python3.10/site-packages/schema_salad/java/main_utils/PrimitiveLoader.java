package ${package}.utils;

public class PrimitiveLoader<T> implements Loader<T> {
  private Class<T> clazz;

  public PrimitiveLoader(Class<T> clazz) {
    this.clazz = clazz;
  }

  public T load(
      final Object doc,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    return Loader.validateOfJavaType(this.clazz, doc);
  }
}
