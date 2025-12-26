import { Dictionary, TypeGuards, Loader, LoadingOptions, ValidationException } from '../Internal'

export class _IdMapLoader implements Loader {
  inner: Loader
  mapSubject: string
  mapPredicate?: string

  constructor (inner: Loader, mapSubject: string, mapPredicate?: string) {
    this.inner = inner
    this.mapSubject = mapSubject
    this.mapPredicate = mapPredicate
  }

  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<any> {
    if (TypeGuards.isDictionary(doc)) {
      const r: any[] = []
      for (var k of Object.keys(doc)) {
        const val = doc[k]
        if (TypeGuards.isDictionary(val)) {
          const v2 = Object.assign({}, val)
          v2[this.mapSubject] = k
          r.push(v2)
        } else {
          if (this.mapPredicate != null) {
            const v3: Dictionary<any> = {}
            v3[this.mapPredicate] = val
            v3[this.mapSubject] = k
            r.push(v3)
          } else {
            throw new ValidationException('No mapPredicate')
          }
        }
      }
      doc = r
    }
    return await this.inner.load(doc, baseuri, loadingOptions)
  }
}
