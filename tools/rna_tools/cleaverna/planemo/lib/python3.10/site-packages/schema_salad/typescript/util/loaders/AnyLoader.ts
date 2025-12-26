import { Loader, LoadingOptions, ValidationException } from '../Internal'

export class _AnyLoader implements Loader {
  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string | undefined): Promise<any> {
    if (doc != null) {
      return doc
    }
    throw new ValidationException('Expected non-null')
  }
}
