using System.Collections;

namespace ${project_name};

internal class ArrayLoader<T> : ILoader<List<T>>
{
    private readonly ILoader itemLoader;

    public ArrayLoader(in ILoader itemLoader)
    {
        this.itemLoader = itemLoader;
    }

    public List<T> Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        if (doc == null)
        {
            throw new ValidationException("Expected non null");
        }

        if (doc is not IList)
        {
            throw new ValidationException("Expected list");
        }

        IList docList = (IList)doc;
        List<T> returnValue = new();
        List<ILoader> loaders = new()
        {
            this,
            itemLoader
        };
        ILoader<object> unionLoader = new UnionLoader(loaders);
        List<ValidationException> errors = new();

        foreach (object? e1 in docList)
        {
            try
            {
                dynamic loadedField = unionLoader.LoadField(e1, baseuri, loadingOptions);
                bool flatten = !String.Equals("@list", loadingOptions.container);
                if (flatten && loadedField is IList)
                {
                    returnValue.AddRange((List<T>)loadedField);
                }
                else
                {
                    returnValue.Add(loadedField);
                }
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
