import { Dictionary, expandUrl, Loader, LoadingOptions } from '../Internal'

export class _TypeDSLLoader implements Loader {
  typeDSLRegex = /^([^[?]+)(\[\])?(\?)?$/
  inner: Loader
  refScope?: number

  constructor (inner: Loader, refScope?: number) {
    this.inner = inner
    this.refScope = refScope
  }

  resolve (doc: string, baseuri: string, loadingOptions: LoadingOptions): Array<Dictionary<string> | string> | Dictionary<string> | string {
    const m = this.typeDSLRegex.exec(doc)
    if (m != null) {
      const group1 = m[1]
      if (group1 == null) {
        throw Error()
      }

      const first = expandUrl(group1, baseuri, loadingOptions, false, true, this.refScope)
      var second
      var third
      if (m[2] != null) {
        second = { type: 'array', items: first }
      }
      if (m[3] != null) {
        third = ['null', second ?? first]
      }
      return third ?? second ?? first
    }
    return doc
  }

  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<any> {
    if (Array.isArray(doc)) {
      const r: any[] = []
      for (const d of doc) {
        if (typeof d === 'string') {
          const resolved = this.resolve(d, baseuri, loadingOptions)
          if (Array.isArray(resolved)) {
            for (const i of resolved) {
              if (!r.includes(i)) {
                r.push(i)
              }
            }
          } else {
            if (!r.includes(resolved)) {
              r.push(resolved)
            }
          }
        } else {
          r.push(d)
        }
      }
      doc = r
    } else if (typeof doc === 'string') {
      doc = this.resolve(doc, baseuri, loadingOptions)
    }
    return await this.inner.load(doc, baseuri, loadingOptions)
  }
}
