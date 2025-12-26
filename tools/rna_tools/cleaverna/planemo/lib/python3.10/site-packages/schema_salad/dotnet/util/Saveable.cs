namespace ${project_name};
using System.Collections;
using OneOf;
using OneOf.Types;

public interface ISaveable
{
    public static abstract ISaveable FromDoc(object doc, string baseUri, LoadingOptions loadingOptions, string? docRoot = null);
    public abstract Dictionary<object, object> Save(bool top, string baseUrl, bool relativeUris);

    public static object Save(object? val_, bool top = true, string baseurl = "", bool relativeUris = true)
    {
        object? val = val_;
        if (val is IOneOf oneOfVal)
        {
            return Save(oneOfVal.Value, top, baseurl, relativeUris);
        }

        if (val is None)
        {
            return null!;
        }

        if (val_ is null)
        {
            return null!;
        }

        if (val is IEnumClass)
        {
            val = val.ToString()!;
        }

        if (val is ISaveable valSaveable)
        {
            return valSaveable.Save(top, baseurl, relativeUris);
        }

        if (val is IList valList)
        {
            List<object> r = new();
            foreach (object v in valList)
            {
                r.Add(Save(v, false, baseurl, relativeUris));
            }

            return r;
        }

        if (val is IDictionary valDict)
        {
            Dictionary<object, object> newDict = new();
            foreach (DictionaryEntry entry in valDict)
            {
                newDict[entry.Key] = Save(entry.Value, false, baseurl, relativeUris);
            }

            return newDict;
        }

        return val!;
    }

    public static object SaveRelativeUri(object? uri_, bool scopedId, bool relativeUris, int? refScope, string baseUrl = "")
    {
        object? uri = uri_;

        if (uri is IOneOf oneOfVal)
        {
            return SaveRelativeUri(oneOfVal.Value, scopedId, relativeUris, refScope, baseUrl);
        }

        if (uri is None)
        {
            return null!;
        }

        if (uri is null)
        {
            return null!;
        }

        if (uri is IEnumClass)
        {
            uri = uri.ToString()!;
        }

        if (relativeUris == false || (uri is string @string && @string == baseUrl))
        {
            return uri;
        }

        if (uri is IList uriList)
        {
            List<object> r = new();
            foreach (object v in uriList)
            {
                r.Add(SaveRelativeUri(v, scopedId, relativeUris, refScope, baseUrl));
            }

            return r;
        }
        else if (uri is string uriString)
        {
            Uri uriSplit = new(uriString, UriKind.RelativeOrAbsolute);
            Uri baseSplit = new(baseUrl, UriKind.RelativeOrAbsolute);
            if (((uriSplit.IsAbsoluteUri && uriSplit.AbsolutePath.Length < 1)
                && (baseSplit.IsAbsoluteUri && baseSplit.AbsolutePath.Length < 1)))
            {
                throw new ValidationException("Uri or baseurl need to contain a path");
            }

            if (uriSplit.IsAbsoluteUri && baseSplit.IsAbsoluteUri && uriSplit.Scheme == baseSplit.Scheme && uriSplit.Host == baseSplit.Host)
            {
                if (uriSplit.AbsolutePath != baseSplit.AbsolutePath)
                {
                    string p = Path.GetRelativePath(Path.GetDirectoryName(baseSplit.AbsolutePath)!, uriSplit.AbsolutePath);
                    if (uriSplit.Fragment.Length > 0)
                    {
                        p = p + "#" + uriSplit.FragmentWithoutFragmentation();
                    }

                    return p;
                }

                string baseFrag = baseSplit.FragmentWithoutFragmentation() + "/";
                if (refScope != null)
                {
                    List<string> sp = baseFrag.Split('/').ToList();
                    int i = 0;
                    while (i < refScope)
                    {
                        sp.RemoveAt(sp.Count - 1);
                        i += 1;
                    }

                    baseFrag = string.Join('/', sp);
                }

                if (uriSplit.FragmentWithoutFragmentation().StartsWith(baseFrag))
                {
                    return uriSplit.FragmentWithoutFragmentation().Substring(baseFrag.Length);
                }
                else
                {
                    return uriSplit.FragmentWithoutFragmentation();
                }
            }
            else
            {
                return Save(uri, false, baseUrl);
            }

        }

        throw new ValidationException("uri needs to be of type List or String");
    }
}
