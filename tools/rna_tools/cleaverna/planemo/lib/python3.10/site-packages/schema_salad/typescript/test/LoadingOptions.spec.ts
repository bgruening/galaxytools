import { assert } from 'chai'
import { Fetcher, LoadingOptions } from '../util/Internal'

class TestFetcher implements Fetcher {
  checkExists(url: string): boolean {
    return true
  }
  async fetchText(url: string, contentTypes?: string[] | undefined): Promise<string> {
    return "TestFetcher"
  }
  urljoin(baseUrl: string, url: string): string {
    return `${baseUrl}/${url}`
  }
}

describe('Test LoadingOptions', () => {
  describe('copyFrom', () => {
    const original = new LoadingOptions({fetcher: new TestFetcher()})
    it('should have the same Fetcher as the original', async () => {
      const copy = new LoadingOptions({copyFrom:original})
      assert.equal(copy.fetcher,original.fetcher)
    })
    it('fetcher should take precedence over copyFrom', async () => {
      const fetcher = new TestFetcher()
      const copy = new LoadingOptions({fetcher,copyFrom:original})
      assert.equal(copy.fetcher,fetcher)
    })
  })
})
