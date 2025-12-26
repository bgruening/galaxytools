namespace ${project_name};
internal class PrimitiveLoader<T> : ILoader<T>
{
    public T Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        if (doc == null)
        {
            throw new ValidationException("Expected non null");
        }

        if (doc is not T)
        {
            throw new ValidationException($"Expected object with type of {typeof(T)} but got {doc.GetType()}");
        }

        return (T)doc;
    }

    object ILoader.Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot)
    {
        return Load(doc,
                    baseuri,
                    loadingOptions,
                    docRoot)!;
    }
}

