export class ValidationException extends Error {
  static indentPerLevel = 2

  children: ValidationException[]
  bullet: string = ''

  constructor (message: string);
  constructor (message: string, _children: ValidationException[]);
  constructor (message: string, _children = new Array<ValidationException>()) {
    super(message)
    Object.setPrototypeOf(this, new.target.prototype)
    const children = new Array<ValidationException>()
    for (const child of _children) {
      children.push(...child.simplify())
    }
    this.children = children
  }

  withBullet (bullet: string): ValidationException {
    this.bullet = bullet
    return this
  }

  simplify (): ValidationException[] {
    if (this.toString().length > 0) {
      return new Array(this)
    }
    return this.children
  }

  summary (level: number): string {
    const spaces = new Array(level * ValidationException.indentPerLevel).join(' ')
    return spaces + this.bullet + this.message
  }

  prettyStr (level = 0): string {
    const parts = new Array<string>()
    let nextLevel = level

    if (this.message != null && this.message.length > 0) {
      parts.push(this.summary(level))
      nextLevel++
    }

    for (const child of this.children) {
      parts.push(child.prettyStr(nextLevel))
    }

    const ret = parts.join('\n')
    return ret
  }

  override toString (): string {
    return this.prettyStr()
  }
}
