namespace ${project_name};

internal static class UriExtensions
{
    /// Because DotNets URI leaves in the # when parsing Fragments, we need a function that removes the # character
    public static string FragmentWithoutFragmentation(this Uri uri)
    {
        return uri.Fragment.Length > 0 ? uri.Fragment.Substring(1) : "";
    }

    public static string FragmentWithoutFragmentation(this UriBuilder uri)
    {
        return uri.Fragment.Length > 0 ? uri.Fragment.Substring(1) : "";
    }
}
