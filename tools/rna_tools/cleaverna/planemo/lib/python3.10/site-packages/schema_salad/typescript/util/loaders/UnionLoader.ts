import { Loader, LoadingOptions, Saveable, ValidationException } from '../Internal'

export class _UnionLoader implements Loader {
  alternates: Loader[]

  constructor (alternates: Loader[]) {
    this.alternates = alternates
  }

  addLoaders(loaders: Loader[]) {
    this.alternates = this.alternates.concat(loaders);
  }

  async load (doc: any, baseuri: string, loadingOptions: LoadingOptions, docRoot?: string): Promise<Saveable> {
    const errors: ValidationException[] = []
    for (const t of this.alternates) {
      try {
        return await t.load(doc, baseuri, loadingOptions, docRoot)
      } catch (e) {
        if (e instanceof ValidationException) {
          errors.push(new ValidationException(`tried ${t.constructor.name} but`, [e]))
        } else {
          throw e
        }
      }
    }
    throw new ValidationException('', errors).withBullet('-')
  }
}
