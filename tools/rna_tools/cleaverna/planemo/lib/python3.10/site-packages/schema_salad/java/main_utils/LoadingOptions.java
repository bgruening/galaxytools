package ${package}.utils;

import java.net.URI;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

public class LoadingOptions {
  Fetcher fetcher;
  String fileUri;
  Map<String, String> namespaces;
  List<String> schemas;
  Boolean noLinkCheck;
  String container;
  Map<String, Object> idx;
  Map<String, String> vocab;
  Map<String, String> rvocab;

  LoadingOptions(
      final Fetcher fetcher,
      final String fileUri,
      final Map<String, String> namespaces,
      final List<String> schemas,
      final Boolean noLinkCheck,
      final String container,
      final Map<String, Object> idx) {
    this.fetcher = fetcher;
    this.fileUri = fileUri;
    this.namespaces = namespaces;
    this.schemas = schemas;
    this.noLinkCheck = noLinkCheck;
    this.container = container;
    this.idx = idx;

    if (namespaces != null) {
      this.vocab = (Map<String, String>) ConstantMaps.vocab.clone();
      this.rvocab = (Map<String, String>) ConstantMaps.rvocab.clone();
      for (Map.Entry<String, String> namespaceEntry : namespaces.entrySet()) {
        this.vocab.put(namespaceEntry.getKey(), namespaceEntry.getValue());
        this.rvocab.put(namespaceEntry.getValue(), namespaceEntry.getKey());
      }
    } else {
      this.vocab = (Map<String, String>) ConstantMaps.vocab;
      this.rvocab = (Map<String, String>) ConstantMaps.rvocab;
    }
  }

  public String expandUrl(
      String url_,
      final String baseUrl,
      final boolean scopedId,
      final boolean vocabTerm,
      final Integer scopedRef) {
    // NOT CONVERTING this - doesn't match type declaration
    // if not isinstance(url, str):
    //    return url
    String url = url_;
    if (url.equals("@id") || url.equals("@type")) {
      return url;
    }

    if (vocabTerm && this.vocab.containsKey(url)) {
      return url;
    }

    if (!this.vocab.isEmpty() && url.contains(":")) {
      String prefix = url.split(":")[0];
      if (this.vocab.containsKey(prefix)) {
        url = this.vocab.get(prefix) + url.substring(prefix.length() + 1);
      }
    }

    Uris.UriSplit split = Uris.split(url);
    final String scheme = split.scheme;
    final boolean hasFragment = stringHasContent(split.fragment);
    if (scheme != null
        && ((scheme.length() > 0
                && (scheme.equals("http") || scheme.equals("https") || scheme.equals("file")))
            || url.startsWith("$(")
            || url.startsWith("${"))) {
      // pass
    } else if (scopedId && !hasFragment) {
      final Uris.UriSplit splitbase = Uris.split(baseUrl);
      final String frg;
      if (stringHasContent(splitbase.fragment)) {
        frg = splitbase.fragment + "/" + split.path;
      } else {
        frg = split.path;
      }
      String pt;
      if (!splitbase.path.equals("")) {
        pt = splitbase.path;
      } else {
        pt = "/";
      }
      url = Uris.unsplit(splitbase.scheme, splitbase.netloc, pt, splitbase.query, frg);
    } else if (scopedRef != null && !hasFragment) {
      final Uris.UriSplit splitbase = Uris.split(baseUrl);
      final ArrayList<String> sp = new ArrayList(Arrays.asList(splitbase.fragment.split("/")));
      int n = scopedRef;
      while (n > 0 && sp.size() > 0) {
        sp.remove(sp.size()-1);
        n -= 1;
      }
      sp.add(url);
      final String fragment = String.join("/", sp);
      url = Uris.unsplit(splitbase.scheme, splitbase.netloc, splitbase.path, splitbase.query, fragment);
    } else {
      url = this.fetcher.urlJoin(baseUrl, url);
    }

    if (vocabTerm) {
      split = Uris.split(url);
      if (stringHasContent(split.scheme)) {
        if (this.rvocab.containsKey(url)) {
          return this.rvocab.get(url);
        }
      } else {
        throw new ValidationException("Term '{}' not in vocabulary".format(url));
      }
    }
    return url;
  }

  static boolean stringHasContent(final String s) {
    return s != null && s.length() > 0;
  }
}
