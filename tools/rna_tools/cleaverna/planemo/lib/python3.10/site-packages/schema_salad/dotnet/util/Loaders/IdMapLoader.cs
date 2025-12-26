using System.Collections;

namespace ${project_name};

internal class IdMapLoader : ILoader<object>
{
    readonly ILoader inner;
    readonly string mapSubject;
    readonly string? mapPredicate;

    public IdMapLoader(ILoader inner, string mapSubject, string? mapPredicate = null)
    {
        this.inner = inner;
        this.mapSubject = mapSubject;
        this.mapPredicate = mapPredicate;
    }

    public object Load(in object doc_, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        object doc = doc_;
        if (doc_ is IDictionary docDict)
        {
            List<object> r = new();
            Dictionary<string, object> d = docDict.Cast<dynamic>().ToDictionary(entry => (string)entry.Key, entry => entry.Value);
            foreach (string? k in d.Keys)
            {
                object val = d[k];
                if (val is IDictionary dictionary)
                {
                    Dictionary<string, object> v2 = new(
                        dictionary.Cast<dynamic>().ToDictionary(entry => (string)entry.Key, entry => entry.Value))
                    {
                        [mapSubject] = k
                    };
                    r.Add(v2);
                }
                else
                {
                    if (mapPredicate != null)
                    {
                        Dictionary<string, object> v3 = new()
                        {
                            [mapPredicate] = val,
                            [mapSubject] = k
                        };
                        r.Add(v3);
                    }
                    else
                    {
                        throw new ValidationException("No mapPredicate was specified");
                    }
                }
            }

            doc = r;
        }

        return inner.Load(doc, baseuri, loadingOptions);
    }

    object ILoader.Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot)
    {
        return Load(doc,
                    baseuri,
                    loadingOptions,
                    docRoot)!;
    }
}
