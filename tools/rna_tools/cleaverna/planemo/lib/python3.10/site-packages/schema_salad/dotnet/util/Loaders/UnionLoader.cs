namespace ${project_name};

internal class UnionLoader : ILoader<object>
{
    private readonly List<ILoader> alternatives;

    public UnionLoader(in List<ILoader> alternatives)
    {
        this.alternatives = alternatives;
    }

    public void addLoaders(IEnumerable<ILoader> loaders)
    {
        this.alternatives.AddRange(loaders);
    }

    public object Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        List<ValidationException> errors = new();

        foreach (ILoader loader in this.alternatives)
        {
            try
            {
                return loader.Load(doc, baseuri, loadingOptions, docRoot);
            }
            catch (ValidationException e)
            {
                errors.Add(e);
            }
        }

        throw new ValidationException("Failed to match union type:", errors);
    }
}
