import { LoadingOptions, documentLoadByUrl, TypeGuards, ValidationException } from '../Internal'
import * as URI from 'uri-js'

export interface Loader {
  load: (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string) => Promise<any>
}

export async function loadField (val: any, fieldType: Loader, baseuri: string, loadingOptions: LoadingOptions): Promise<any> {
  if (TypeGuards.isDictionary(val)) {
    if ('$import' in val) {
      if (loadingOptions.fileUri == null) {
        throw Error('Cannot load $import without fileuri')
      }
      return await documentLoadByUrl(fieldType, loadingOptions.fetcher.urljoin(loadingOptions.fileUri, val.$import), loadingOptions)
    } else if ('$include' in val) {
      if (loadingOptions.fileUri == null) {
        throw Error('Cannot load $import without fileuri')
      }
      val = await loadingOptions.fetcher.fetchText(loadingOptions.fetcher.urljoin(loadingOptions.fileUri, val.$include))
    }
  }
  return await fieldType.load(val, baseuri, loadingOptions)
}

export function expandUrl (url: string, baseUrl: string, loadingOptions: LoadingOptions, scopedId = false, vocabTerm = false, scopedRef?: number): string {
  if (['@id', '@type'].includes(url)) {
    return url
  }

  if (vocabTerm && url in loadingOptions.vocab) {
    return url
  }

  if (loadingOptions.vocab != null && url.includes(':')) {
    const prefix = url.split(':')[0]
    if (prefix in loadingOptions.vocab) {
      url = loadingOptions.vocab[prefix] + url.slice(prefix.length + 1)
    }
  }

  const split = URI.parse(url)
  if ((split.scheme != null && ['http', 'https', 'file'].includes(split.scheme)) || url.startsWith('$(') || url.startsWith('${')) {
  } else if (scopedId && split.fragment === undefined) {
    const splitbase = URI.parse(baseUrl)
    let frg = ''
    if (splitbase.fragment != null) {
      frg = splitbase.fragment + '/' + (split.path ?? '')
    } else {
      frg = split.path ?? ''
    }
    const pt = splitbase.path ?? '/'
    const parts = {
      scheme: splitbase.scheme,
      userinfo: undefined,
      host: splitbase.host,
      port: undefined,
      path: pt,
      query: splitbase.query,
      fragment: frg,
      reference: undefined,
      error: undefined
    }

    url = URI.serialize(parts)
  } else if (scopedRef != null && split.fragment === undefined) {
    const splitbase = URI.parse(baseUrl)
    const sp = splitbase.fragment?.split('/') ?? []
    let n = scopedRef
    while (n > 0 && sp?.length > 0) {
      sp.pop()
      n -= 1
    }
    sp.push(url)
    const parts = {
      scheme: splitbase.scheme,
      userinfo: undefined,
      host: splitbase.host,
      port: undefined,
      path: splitbase.path,
      query: splitbase.query,
      fragment: sp.join('/'),
      reference: undefined,
      error: undefined
    }
    url = URI.serialize(parts)
  } else {
    url = loadingOptions.fetcher.urljoin(baseUrl, url)
  }

  if (vocabTerm) {
    const split = URI.parse(url)
    if (split.scheme !== undefined) {
      if (url in loadingOptions.rvocab) {
        return loadingOptions.rvocab[url]
      } 
    } else {
      throw new ValidationException(`Term '${url}' not in vocabulary`)
    }
  }
  return url
}
