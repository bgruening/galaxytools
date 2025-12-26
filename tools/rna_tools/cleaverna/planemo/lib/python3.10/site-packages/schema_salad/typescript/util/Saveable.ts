import { LoadingOptions, Dictionary, TypeGuards } from './Internal'
import * as URI from 'uri-js'
import path from 'path'

// eslint-disable-next-line @typescript-eslint/no-extraneous-class
export abstract class Saveable {
  loadingOptions: LoadingOptions
  constructor(loadingOptions?: LoadingOptions) {
    this.loadingOptions = loadingOptions ?? new LoadingOptions({})
  }
  static async fromDoc (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<Saveable> {
    throw new Error('Not Implemented')
  }
  abstract save (top: boolean, baseUrl: string, relativeUris: boolean): Dictionary<any>
}

export function save (val: any, top: boolean = true, baseUrl: string = '', relativeUris: boolean = true): any {
  if (val instanceof Saveable) {
    return val.save(top, baseUrl, relativeUris)
  }

  if (Array.isArray(val)) {
    const r = []
    for (const v of val) {
      r.push(save(v, false, baseUrl, relativeUris))
    }
    return r
  }

  if (TypeGuards.isDictionary(val)) {
    const newDict: Dictionary<any> = {}
    for (const key in val) {
      newDict[key] = save(val[key], false, baseUrl, relativeUris)
    }
    return newDict
  }

  return val
}

export function saveRelativeUri (uri: any, baseUrl: string='', scopedId: boolean, relativeUris: boolean, refScope?: number): any {
  if (relativeUris === false || uri === baseUrl) {
    return uri
  }
  if (Array.isArray(uri)) {
    const r = []
    for (const v of uri) {
      r.push(saveRelativeUri(v, baseUrl, scopedId, relativeUris, refScope))
    }
    return r
  } else if (typeof uri === 'string') {
    const uriSplit = URI.parse(uri)
    const baseSplit = URI.parse(baseUrl)
    if (uriSplit.path == null || baseSplit.path == null) {
      throw new Error('uri or baseurl need to contain a path.')
    }
    if (uriSplit.scheme === baseSplit.scheme && uriSplit.host === baseSplit.host) {
      if (uriSplit.path !== baseSplit.path) {
        let p = path.relative(path.dirname(baseSplit.path), uriSplit.path)
        if (uriSplit.fragment != null) {
          p = p + '#' + uriSplit.fragment
        }
        return p
      }

      if (baseSplit.fragment == null) {
        baseSplit.fragment = ''
      }
      let basefrag = baseSplit.fragment + '/'
      if (refScope != null) {
        const sp = basefrag.split('/')
        let i = 0
        while (i < refScope) {
          sp.pop()
          i += 1
        }
        basefrag = sp.join('/')
      }
      if (uriSplit.fragment == null) {
        uriSplit.fragment = ''
      }
      if (uriSplit.fragment.startsWith(basefrag)) {
        return uriSplit.fragment.slice(basefrag.length)
      } else {
        return uriSplit.fragment
      }
    } else {
      return save(uri, false, baseUrl)
    }
  }
}

export function prefixUrl (url: string, namespaces: Dictionary<string>): string {
  for (const k in namespaces) {
    if (url.startsWith(namespaces.k)) {
      return k + ':' + url.slice(namespaces.k.length)
    }
  }
  return url
}

/**
 * Compute the shortname of a fully qualified identifier.
 * See https://w3id.org/cwl/v1.2/SchemaSalad.html#Short_names. 
 *
 */
export function shortname (inputId: string): string {
  const parsedId = URI.parse(inputId)
  if (parsedId.fragment != null) {
    const fragmentSplit = parsedId.fragment.split('/')
    return fragmentSplit[fragmentSplit.length - 1]
  } else if (parsedId.path != null) {
    const pathSplit = parsedId.path.split('/')
    return pathSplit[pathSplit.length - 1]
  } else {
    return inputId
  }
}
