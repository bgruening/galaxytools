namespace ${project_name};

internal class EnumLoader<T> : ILoader<T> where T : IEnumClass<T>
{
    public T Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        if (doc is not string)
        {
            throw new ValidationException("Expected a string");
        }

        if (T.Contains((string)doc))
        {
            T val = T.Parse((string)doc);
            return val;
        }
        else
        {
            throw new ValidationException(
                $"Symbol not contained in {typeof(T).Name} Enum, expected one of {string.Join(", ", T.Symbols())}");
        }
    }

    object ILoader.Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot)
    {
        return Load(doc,
                    baseuri,
                    loadingOptions,
                    docRoot)!;
    }
}

