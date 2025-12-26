package ${package}.utils;

import java.lang.reflect.Method;
import java.lang.ReflectiveOperationException;
import java.util.Arrays;
import java.util.List;

public class EnumLoader<T extends Enum> implements Loader<T>{
  private final Class<T> symbolEnumClass;

  public EnumLoader(final Class<T> symbolEnumClass) {
    this.symbolEnumClass = symbolEnumClass;
  }

  public T load(
      final Object doc,
      final String baseUri,
      final LoadingOptions loadingOptions,
      final String docRoot) {
    final String docString = Loader.validateOfJavaType(String.class, doc);
    try {
      final Method m = symbolEnumClass.getMethod("fromDocumentVal", String.class);
      final T val = (T) m.invoke(null, docString);
      return val;
    } catch (final ReflectiveOperationException e) {
      final Throwable cause = e.getCause();
      if (cause instanceof RuntimeException) {
        throw (RuntimeException) cause;
      }
      throw new RuntimeException(e);
    }
  }
}
