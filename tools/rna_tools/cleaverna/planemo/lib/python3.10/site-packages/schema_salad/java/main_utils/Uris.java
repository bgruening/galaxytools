package ${package}.utils;

import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.charset.StandardCharsets;

public class Uris {

	// Emulate Python's urlsplit.
	public static class UriSplit {
		String scheme;
		String netloc;
		String path;
		String query;
		String fragment;

		public UriSplit(String scheme, String netloc, String path, String query, String fragment) {
			this.scheme = scheme;
			this.netloc = netloc;
			this.path = path;
			this.query = query;
			this.fragment = fragment;
		}

		public String toString() {
			return String.format("UriSplit[%s,%s,%s,%s,%s]", this.scheme, this.netloc, this.path, this.query,
					this.fragment);
		}

	}

	public static String fileUri(final String path) {
		return fileUri(path, false);
	}

	public static String fileUri(final String path, final boolean splitFrag) {
		if (path.equals("file://")) {
			return path;
		}
		String frag;
		String urlPath;
		if (splitFrag) {
			final String[] pathsp = path.split("#", 2);
			// is quoting this?
			urlPath = Uris.quote(pathsp[0]);
			if (pathsp.length == 2) {
				frag = "#" + Uris.quote(pathsp[1]);
			} else {
				frag = "";
				urlPath = Uris.quote(path);
			}
		} else {
			urlPath = Uris.quote(path);
			frag = "";
		}
		if (urlPath.startsWith("//")) {
			return "file:" + urlPath + frag;
		} else {
			return "file://" + urlPath + frag;
		}
	}

	public static UriSplit split(final String uriString) {
		try {
			final URI uri = new URI(uriString);
			return new Uris.UriSplit(uri.getScheme(), uri.getAuthority(), uri.getPath(), uri.getQuery(),
					uri.getFragment());
		} catch (URISyntaxException e) {
			return new Uris.UriSplit(null, null, uriString, null, null);
		}
	}

	public static String unsplit(final String scheme, final String netloc, final String path, final String query,
			final String fragment) {
		try {
			return new URI(scheme, netloc, path, query, fragment).toString();
		} catch (URISyntaxException e) {
			if (scheme == null && path.startsWith("_:")) {
				String uri = path;
				if (fragment != null && fragment.length() > 0) {
					uri += "#" + fragment;
				}
				return fragment;
			}
			throw new RuntimeException(e);
		}
	}

	public static URI toUri(final String url) {
		try {
			return new URI(url);
		} catch (URISyntaxException e) {
			throw new RuntimeException(e);
		}
	}

	public static String quote(final String uri) {
		try {
			return java.net.URLDecoder.decode(uri, StandardCharsets.UTF_8.name());
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException(e);
		}
	}

	public static String unquote(final String uri) {
		try {
			return java.net.URLEncoder.encode(uri, StandardCharsets.UTF_8.name());
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException(e);
		}
	}

	public static String shortname(final String input_id) {
		try {
			final URI uri = new URI(input_id);
			final String fragment = uri.getFragment();
			if (fragment != null) {
				String[] fragment_elements = fragment.split("/");
				return fragment_elements[fragment_elements.length - 1];
			} else {
				String[] path_elements = uri.getPath().split("/");
				return path_elements[path_elements.length - 1];
			}
		} catch (URISyntaxException e) {
			return input_id;
		}
	}
}
