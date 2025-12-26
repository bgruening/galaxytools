import { Loader, LoadingOptions, ValidationException } from '../Internal'

export class _PrimitiveLoader implements Loader {
  typeGuard: (val: any) => boolean

  constructor (tp: (val: any) => boolean) {
    this.typeGuard = tp
  }

  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<any> {
    if (!this.typeGuard(doc)) {
      throw new ValidationException(`Expected a ${this.typeGuard.name} but got ${(typeof doc)}`)
    }
    return doc
  }
}
