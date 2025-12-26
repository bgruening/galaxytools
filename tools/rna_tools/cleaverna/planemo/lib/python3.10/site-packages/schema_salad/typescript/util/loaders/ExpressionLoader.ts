import { Loader, LoadingOptions, ValidationException } from '../Internal'

export class _ExpressionLoader implements Loader {

  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<any> {
    if (typeof doc !== 'string') {
      throw new ValidationException('Expected a str')
    }
    return doc
  }
}
