import { LoadingOptions, Loader, TypeGuards, LoaderInstances } from '../Internal'
import * as Internal from  '../Internal'
import * as yaml from 'js-yaml'
import * as URL from 'url'

export async function documentLoad (loader: Loader, doc: unknown, baseuri: string, loadingOptions: LoadingOptions): Promise<any> {
  if (typeof doc === 'string') {
    return await documentLoadByUrl(loader, loadingOptions.fetcher.urljoin(baseuri, doc), loadingOptions)
  }

  if (Array.isArray(doc)) {
    return await loader.load(doc, baseuri, loadingOptions)
  }

  if (TypeGuards.isDictionary(doc)) {
    if (doc != null) {
      if ('$namespaces' in doc || '$schemas' in doc) {
        loadingOptions = new LoadingOptions({ copyFrom: loadingOptions, namespaces: doc.$namespaces ?? undefined, schemas: doc.$schemas ?? undefined })
        delete doc.$schemas
        delete doc.$namespaces
      }

      if ('$base' in doc) {
        baseuri = doc.$base
      }

      if ('$graph' in doc) {
        return await loader.load(doc.$graph, baseuri, loadingOptions)
      } else {
        return await loader.load(doc, baseuri, loadingOptions, baseuri)
      }
    }
  }

  throw new Error('Reached unexpected path')
}

export async function documentLoadByUrl (loader: Loader, url: string, loadingOptions: LoadingOptions): Promise<void> {
  if (url in loadingOptions.idx) {
    return await documentLoad(loader, loadingOptions.idx[url], url, loadingOptions)
  }
  const text = await loadingOptions.fetcher.fetchText(url)
  const result = yaml.load(text)
  loadingOptions.idx[url] = result
  loadingOptions = new LoadingOptions({ copyFrom: loadingOptions, fileUri: url })
  return await documentLoad(loader, result, url, loadingOptions)
}

export async function loadDocument (doc: any, baseuri?: string, loadingOptions?: LoadingOptions): Promise<${root_loader_type}> {
  if (baseuri == null) {
    baseuri = URL.pathToFileURL(process.cwd() + '/').toString()
  }

  if (loadingOptions == null) {
    loadingOptions = new LoadingOptions({})
  }

  return await documentLoad(LoaderInstances.${root_loader}, doc, baseuri, loadingOptions)
}

export async function loadDocumentByString (doc: string, uri: string, loadingOptions?: LoadingOptions): Promise<${root_loader_type}> {
  const result = yaml.load(doc)

  if (loadingOptions == null) {
    loadingOptions = new LoadingOptions({ fileUri: uri })
  }
  loadingOptions.idx[uri] = result

  return await documentLoad(LoaderInstances.${root_loader}, result, uri, loadingOptions)
}
