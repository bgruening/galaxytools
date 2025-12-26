// Code implemented after https://github.com/common-workflow-language/schema_salad/blob/main/schema_salad/fetcher.py
namespace ${project_name};

public interface IFetcher
{
    string FetchText(string uri);
    bool CheckExists(string uri);
    string Urljoin(string baseUrl, string url);
    protected static readonly string[] Schemes = new string[] { "file", "http", "https", "mailto" };
}

public class DefaultFetcher : IFetcher
{

    private readonly HttpClient client;

    public DefaultFetcher()
    {
        this.client = new HttpClient();
    }

    public DefaultFetcher(HttpClient client)
    {
        this.client = client;
    }

    public bool CheckExists(string uri)
    {
        throw new NotImplementedException();
    }

    public string FetchText(string uri)
    {
        UriBuilder split = Utilities.Split(uri);
        string scheme = split.Scheme;
        if (IFetcher.Schemes.Contains(scheme))
        {
            if ((new[] { "http", "https" }).Contains(scheme))
            {
                try
                {
                    HttpResponseMessage response = client.GetAsync(uri).Result;
                    response.EnsureSuccessStatusCode();
                    return response.Content.ReadAsStringAsync().Result;
                }
                catch (Exception e)
                {
                    throw new ValidationException($"Error fetching {uri}: {e.Message} ");
                }
            }
            else if (scheme == "file")
            {
                try
                {
                    string fileContent = System.IO.File.ReadAllText(Path.GetFullPath(split.Path));
                    return fileContent;
                }
                catch (Exception e)
                {
                    throw new ValidationException($"Error reading file {uri}: {e.Message} ");
                }
            }
        }

        throw new ValidationException($"Unsupported scheme {scheme} in url: {uri}");
    }

    public string Urljoin(string baseUrl, string url)
    {
        if (url.StartsWith("_:"))
        {
            return url;
        }

        Uri baseUri = new(baseUrl);
        Uri uri = new(url, UriKind.Relative);

        return new Uri(baseUri, uri).ToString();
    }
}
