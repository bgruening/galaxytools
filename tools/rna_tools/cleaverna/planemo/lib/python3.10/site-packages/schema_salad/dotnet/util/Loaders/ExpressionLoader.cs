namespace ${project_name};

internal class ExpressionLoader : ILoader<string>
{
    public string Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        if (doc is string docString)
        {
            return docString;
        }
        else
        {
            throw new ValidationException("Expected a string");
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

