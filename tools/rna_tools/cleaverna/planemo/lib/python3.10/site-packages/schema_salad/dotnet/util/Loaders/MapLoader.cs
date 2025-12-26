using System.Collections;

namespace ${project_name};

internal class MapLoader<T> : ILoader<Dictionary<string, T>>
{
    private readonly ILoader valueLoader;
    private readonly string? container;
    private readonly bool? noLinkCheck;

    public MapLoader(in ILoader valueLoader, in string? container, in bool? noLinkCheck)
    {
        this.valueLoader = valueLoader;
        this.container = container;
        this.noLinkCheck = noLinkCheck;
    }

    public Dictionary<string, T> Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        if (doc == null)
        {
            throw new ValidationException("Expected non null");
        }

        if (doc is not IDictionary)
        {
            throw new ValidationException("Expected list");
        }

        LoadingOptions innerLoadingOptions = loadingOptions;
        if (this.container != null || this.noLinkCheck != null)
        {
            innerLoadingOptions = new LoadingOptions(copyFrom: loadingOptions, container: this.container, noLinkCheck: this.noLinkCheck);
        }

        IDictionary docDictionary = (IDictionary)doc;
        Dictionary<string, T> returnValue = new();
        List<ILoader> loaders = new()
        {
            this,
            valueLoader
        };
        ILoader<object> unionLoader = new UnionLoader(loaders);
        List<ValidationException> errors = new();

        foreach (KeyValuePair<string, T> item in docDictionary)
        {
            try
            {
                dynamic loadedField = unionLoader.LoadField(item.Value, baseuri, innerLoadingOptions);
                returnValue[item.Key] = loadedField;
            }
            catch (ValidationException e)
            {
                errors.Add(e);
            }
        }

        if (errors.Count > 0)
        {
            throw new ValidationException("", errors);
        }

        return returnValue;
    }

    object ILoader.Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot)
    {
        return Load(doc, baseuri, loadingOptions, docRoot);
    }
}
