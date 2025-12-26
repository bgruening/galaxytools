import { Loader, LoadingOptions, Dictionary, TypeGuards, ValidationException } from '../Internal'

export class _SecondaryDSLLoader implements Loader {
  inner: Loader

  constructor (inner: Loader) {
    this.inner = inner
  }

  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<any> {
    const r: Array<Dictionary<any>> = []
    if (Array.isArray(doc)) {
      for (const d of doc) {
        if (typeof d === 'string') {
          if (d.endsWith('?')) {
            r.push({ pattern: d.slice(0, -1), required: false })
          } else {
            r.push({ pattern: d })
          }
        } else if (TypeGuards.isDictionary(d)) {
          const newDict: Dictionary<any> = {}
          if ('pattern' in d) {
            newDict.pattern = d.pattern
            delete d.pattern
          } else {
            throw new ValidationException(`Missing pattern in secondaryFiles specification entry: ${JSON.stringify(d)}`)
          }
          if ('required' in d) {
            newDict.required = d.required
            delete d.required
          }
          if (Object.keys(d).length > 0) {
            throw new ValidationException(`Unallowed values in secondaryFiles specification entry: ${JSON.stringify(d)}`)
          }
          r.push(newDict)
        } else {
          throw new ValidationException('Expected a string or sequence of (strings or mappings).')
        }
      }
    } else if (TypeGuards.isDictionary(doc)) {
      const newDict: Dictionary<any> = {}
      if ('pattern' in doc) {
        newDict.pattern = doc.pattern
        delete doc.pattern
      } else {
        throw new ValidationException(`Missing pattern in secondaryFiles specification entry: ${JSON.stringify(doc)}`)
      }
      if ('required' in doc) {
        newDict.required = doc.required
        delete doc.required
      }
      if (Object.keys(doc).length > 0) {
        throw new ValidationException(`Unallowed values in secondaryFiles specification entry: ${JSON.stringify(doc)}`)
      }
      r.push(newDict)
    } else if (typeof doc === 'string') {
      if (doc.endsWith('?')) {
        r.push({ pattern: doc.slice(0, -1), required: false })
      } else {
        r.push({ pattern: doc })
      }
    } else {
      throw new ValidationException('Expected str or sequence of str')
    }
    return await this.inner.load(r, baseuri, loadingOptions, docRoot)
  }
}
