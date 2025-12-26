import { assert } from 'chai'
import { DefaultFetcher, ValidationException } from '../util/Internal'
import sinon from 'sinon'
import * as fetchModule from 'node-fetch'
import { Response } from 'node-fetch'
import path from 'path'
import URL from 'url'

describe('Test Fetcher', () => {
  const fet = new DefaultFetcher()
  describe('Test fetchText()', () => {
    afterEach(function () {
      sinon.restore()
    })

    it('Should fetch text from http urls', async () => {
      sinon.stub(fetchModule, 'default').returns(new Promise((resolve) => resolve(new Response('test', { status: 200 }))))
      assert.equal(await fet.fetchText('http://www.example.com'), 'test')
    })

    it('Should fetch text from https urls', async () => {
      sinon.stub(fetchModule, 'default').returns(new Promise((resolve) => resolve(new Response('test', { status: 200 }))))
      assert.equal(await fet.fetchText('https://www.example.com'), 'test')
    })

    it('Should fetch text from files', async () => {
      const filepath = URL.pathToFileURL(path.resolve('./src/test/data/test.txt')).toString()
      assert.equal(await fet.fetchText(filepath), 'test\n')
    })

    it('Throw a 404 exception', async () => {
      sinon.stub(fetchModule, 'default').returns(new Promise((resolve) => resolve(new Response('test', { status: 404 }))))
      let err
      try {
        await fet.fetchText('https://www.example.com')
      } catch (e) {
        err = e
      }
      assert.exists(err)
      assert.isTrue(err instanceof ValidationException)
      assert.equal((err as ValidationException).message, 'Error fetching https://www.example.com: HTTP Error Response: 404 Not Found')
    })

    it('Throw an invalid schema exception', async () => {
      let err
      try {
        await fet.fetchText('invalidscheme://www.example.com')
      } catch (e) {
        err = e
      }
      assert.exists(err)
      assert.isTrue(err instanceof ValidationException)
      assert.equal((err as ValidationException).message, 'Unsupported scheme invalidscheme in url: invalidscheme://www.example.com')
    })
  })

  describe('Test urlJoin()', () => {
    it('Should correctly join urls', async () => {
      assert.equal(fet.urljoin('http://example.com/base', 'one'), 'http://example.com/one')
      assert.equal(fet.urljoin('http://example.com/base', 'two'), 'http://example.com/two')
      assert.equal(fet.urljoin('http://example.com/base', '#three'), 'http://example.com/base#three')
      assert.equal(fet.urljoin('http://example.com/base', 'four#five'), 'http://example.com/four#five')
      assert.equal(fet.urljoin('http://example.com/base', '_:five'), '_:five')
    })

    it('Should throw a remote exploit exception', async () => {
      let err
      try {
        fet.urljoin('http://example.com/base', 'file:///test/test.txt')
      } catch (e) {
        err = e
      }
      assert.exists(err)
      assert.isTrue(err instanceof ValidationException)
      assert.equal((err as ValidationException).message, 'Not resolving potential remote exploit file:///test/test.txt from base http://example.com/base')
    })
  })
})
