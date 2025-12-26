using System.Collections;

namespace ${project_name};

internal class UriLoader : ILoader<object>
{
    readonly ILoader inner;
    readonly bool scopedID;
    readonly bool vocabTerm;
    readonly int? scopedRef;
    readonly bool? noLinkCheck;

    public UriLoader(in ILoader inner, in bool scopedID, in bool vocabTerm, in int? scopedRef, in bool? noLinkCheck)
    {
        this.inner = inner;
        this.scopedID = scopedID;
        this.vocabTerm = vocabTerm;
        this.scopedRef = scopedRef;
        this.noLinkCheck = noLinkCheck;
    }

    public object Load(in object doc_, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot = null)
    {
        LoadingOptions innerLoadingOptions = loadingOptions;
        if (this.noLinkCheck != null)
        {
            innerLoadingOptions = new LoadingOptions(copyFrom: loadingOptions, noLinkCheck: this.noLinkCheck);
        }
        object doc = doc_;
        if (doc is IList)
        {
            List<object> docList = (List<object>)doc_;
            List<object> docWithExpansion = new();
            foreach (object val in docList)
            {
                if (val is string valString)
                {
                    docWithExpansion.Add(innerLoadingOptions.ExpandUrl(valString, baseuri, scopedID, vocabTerm, scopedRef));
                }
                else
                {
                    docWithExpansion.Add(val);
                }
            }

            doc = docWithExpansion;
        }
        else if (doc is string docString)
        {
            doc = innerLoadingOptions.ExpandUrl(docString, baseuri, scopedID, vocabTerm, scopedRef);
        }

        return (object)inner.Load(doc, baseuri, innerLoadingOptions);
    }

    object ILoader.Load(in object doc, in string baseuri, in LoadingOptions loadingOptions, in string? docRoot)
    {
        return Load(doc,
                    baseuri,
                    loadingOptions,
                    docRoot)!;
    }
}
