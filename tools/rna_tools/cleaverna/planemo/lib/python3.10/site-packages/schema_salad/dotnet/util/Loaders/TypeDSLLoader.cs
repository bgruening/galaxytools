using System.Collections;
using System.Text.RegularExpressions;

namespace ${project_name};

internal class TypeDSLLoader : ILoader<object>
{
    readonly ILoader inner;
    readonly int refScope;

    private static readonly Regex typeDSLRegex = new(@"^([^\\[?]+)(\\[\\])?(\\?)?$");

    public TypeDSLLoader(in ILoader inner, in int refScope)
    {
        this.inner = inner;
        this.refScope = refScope;
    }

    private object Resolve(in string doc_, in string baseuri, in LoadingOptions loadingOptions)
    {
        Match m = typeDSLRegex.Match(doc_);
        if (m.Success)
        {
            string first = loadingOptions.ExpandUrl(m.Groups[1].ToString(), baseuri, false, true, this.refScope);
            object? second = null;
            object? third = null;
            if (m.Groups.Count >= 3 && m.Groups[2].Length > 0)
            {
                Dictionary<string, object> resolveMap = new()
                {
                    { "type", "array" },
                    { "items", first }
                };
                second = resolveMap;
            }

            if (m.Groups.Count >= 4 && m.Groups[3].Length > 0)
            {
                third = new List<object> { "null", second ?? first };
            }

            return third ?? second ?? first;
        }

        return doc_;
    }

    public object Load(in object doc_, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        object doc = doc_;
        if (doc is IList)
        {
            List<object> docList = (List<object>)doc;
            List<object> r = new();
            foreach (object d in docList)
            {
                if (d is string dString)
                {
                    object resolved = Resolve(dString, baseuri, loadingOptions);
                    if (resolved is IList)
                    {
                        List<object> resolvedList = (List<object>)resolved;
                        foreach (object i in resolvedList)
                        {
                            if (!r.Contains(i))
                            {
                                r.Add(i);
                            }
                        }
                    }
                    else
                    {
                        if (!r.Contains(resolved))
                        {
                            r.Add(resolved);
                        }
                    }
                }
                else
                {
                    r.Add(d);
                }
            }

            doc = docList;
        }
        else if (doc is string docString)
        {
            doc = Resolve(docString, baseuri, loadingOptions);
        }

        return (object)inner.Load(doc, baseuri, loadingOptions);
    }

    object ILoader.Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot)
    {
        return Load(doc,
                    baseuri,
                    loadingOptions,
                    docRoot)!;
    }
}
