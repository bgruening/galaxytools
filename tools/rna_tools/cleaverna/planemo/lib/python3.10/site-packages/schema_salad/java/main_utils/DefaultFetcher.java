package ${package}.utils;

import java.io.IOException;
import java.net.URI;
import java.util.Arrays;
import java.util.Scanner;

public class DefaultFetcher implements Fetcher {

  public String urlJoin(final String baseUrl, final String url) {
    if (url.startsWith("_:")) {
      return url;
    }

    final URI baseUri = Uris.toUri(baseUrl);
    final URI uri = Uris.toUri(url);
    if (baseUri.getScheme() != null
        && !baseUri.getScheme().equals("file")
        && "file".equals(uri.getScheme())) {
      throw new ValidationException(
          String.format(
              "Not resolving potential remote exploit %s from base %s".format(url, baseUrl)));
    }
    String result = baseUri.resolve(uri).toString();
    if (result.startsWith("file:")) {
      // Well this is gross - needed for http as well?
      result = "file://" + result.substring("file:".length());
    }
    return result;
  }

  public String fetchText(final String url) {
    final URI uri = Uris.toUri(url);
    final String scheme = uri.getScheme();
    if (Arrays.asList("http", "https", "file").contains(scheme)) {
      Scanner scanner;
      try {
        scanner = new Scanner(uri.toURL().openStream(), "UTF-8").useDelimiter("\\A");
      } catch (IOException e) {
        throw new ValidationException("Error fetching %s: %s.".format(url, e));
      }
      String result = scanner.next();
      scanner.close();
      return result;
    }
    throw new ValidationException("Unsupported scheme in URL: %s".format(url));
  }
}
