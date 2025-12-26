import { expandUrl, Loader, LoadingOptions } from '../Internal'

export class _URILoader implements Loader {
  inner: Loader
  scopedID: boolean
  vocabTerm: boolean
  scopedRef?: number
  noLinkCheck?: boolean

  constructor (inner: Loader, scopedID: boolean, vocabTerm: boolean, scopedRef?: number, noLinkCheck?: boolean) {
    this.inner = inner
    this.scopedID = scopedID
    this.vocabTerm = vocabTerm
    this.scopedRef = scopedRef
    this.noLinkCheck = noLinkCheck
  }

  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<any> {
    if (this.noLinkCheck !== undefined) {
      loadingOptions = new LoadingOptions({ copyFrom: loadingOptions, noLinkCheck: this.noLinkCheck })
    }
    if (Array.isArray(doc)) {
      const newDoc: any[] = []
      for (const val of doc) {
        if (typeof val === 'string') {
          newDoc.push(expandUrl(val, baseuri, loadingOptions, this.scopedID, this.vocabTerm, this.scopedRef))
        } else {
          newDoc.push(val)
        }
      }
      doc = newDoc
    } else if (typeof doc === 'string') {
      doc = expandUrl(doc, baseuri, loadingOptions, this.scopedID, this.vocabTerm, this.scopedRef)
    }
    return await this.inner.load(doc, baseuri, loadingOptions)
  }
}
