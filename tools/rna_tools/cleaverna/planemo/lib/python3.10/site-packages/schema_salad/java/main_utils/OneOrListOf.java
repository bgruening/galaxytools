package ${package}.utils;

import java.util.List;
import java.util.Optional;

public class OneOrListOf<T> {
    private Optional<T> object;
    private Optional<List<T>> objects;

    private OneOrListOf(final T object, final List<T> objects) {
        this.object = Optional.ofNullable(object);
        this.objects = Optional.ofNullable(objects);
    }

    public static <T> OneOrListOf<T> oneOf(T object) {
        return new OneOrListOf(object, null);
    }

    public static <T> OneOrListOf<T> listOf(List<T> objects) {
        assert objects != null;
        return new OneOrListOf(null, objects);
    }

    public boolean isOne() {
        return this.getOneOptional().isPresent();
    }

    public boolean isList() {
        return this.getListOptional().isPresent();
    }

    public Optional<T> getOneOptional() {
        return this.object;
    }
    
    public Optional<List<T>> getListOptional() {
        return this.objects;
    }

    public T getOne() {
        return this.getOneOptional().get();
    }

    public List<T> getList() {
        return this.getListOptional().get();
    }

}
