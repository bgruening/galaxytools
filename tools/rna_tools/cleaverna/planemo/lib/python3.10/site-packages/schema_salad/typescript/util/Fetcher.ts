import { ValidationException } from './Internal'
import fetch from 'node-fetch'
import * as fs from 'fs'
import * as URI from 'uri-js'

// Code implemented after https://github.com/common-workflow-language/schema_salad/blob/main/schema_salad/fetcher.py
export abstract class Fetcher {
  abstract fetchText (url: string, contentTypes?: string[]): Promise<string>
  abstract checkExists (url: string): boolean
  abstract urljoin (baseUrl: string, url: string): string

  static schemes = ['file', 'http', 'https', 'mailto']
}

export class DefaultFetcher extends Fetcher {
  async fetchText (urlString: string): Promise<string> {
    // TODO: cache
    const split = URI.parse(urlString)
    const scheme = split.scheme ?? ''
    if (Fetcher.schemes.includes(scheme)) {
      if (['http', 'https'].includes(scheme)) {
        try {
          // TODO: content types
          const result = await fetch(new URL(urlString))
          if (!result.ok) {
            throw Error(`HTTP Error Response: ${result.status} ${result.statusText}`)
          }
          return await result.text()
        } catch (e) {
          if (e instanceof Error) {
            throw new ValidationException(`Error fetching ${urlString}: ${e.message}`)
          } else {
            throw e
          }
        }
      } else if (scheme === 'file') {
        try {
          return fs.readFileSync(split.path ?? '', { encoding: 'utf8' })
        } catch (e) {
          if (e instanceof Error) {
            throw new ValidationException(`Error reading file ${urlString}: ${e.message}`)
          } else {
            throw e
          }
        }
      }
    }
    throw new ValidationException(`Unsupported scheme ${scheme} in url: ${urlString}`)
  }

  checkExists (url: string): boolean {
    throw new Error('Not implemented.')
  }

  urljoin (baseUrlString: string, urlString: string): string {
    if (urlString.startsWith('_:')) {
      return urlString
    }
    const baseUrl = URI.parse(baseUrlString)
    const url = URI.parse(urlString)
    if (baseUrl.scheme != null && baseUrl.scheme !== 'file' && url.scheme === 'file') {
      throw new ValidationException(`Not resolving potential remote exploit ${urlString} from base ${baseUrlString}`)
    }
    // TODO: Windows specific join?
    return new URL(urlString, baseUrlString).toString()
  }
}
