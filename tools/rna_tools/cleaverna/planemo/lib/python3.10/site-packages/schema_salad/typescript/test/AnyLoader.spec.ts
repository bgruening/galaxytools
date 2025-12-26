import chai from 'chai'
import chaiAsPromised from 'chai-as-promised'

import { _AnyLoader, LoadingOptions } from '../util/Internal'

chai.use(chaiAsPromised)
const assert = chai.assert

const loader = new _AnyLoader()

describe('Test AnyLoader', () => {
  it('Should load the documents', async () => {
    assert.deepEqual(await loader.load({}, '', new LoadingOptions({})), {})
    assert.equal(await loader.load(2, '', new LoadingOptions({})), 2)
    assert.equal(await loader.load('2', '', new LoadingOptions({})), '2')
    assert.deepEqual(await loader.load([2], '', new LoadingOptions({})), [2])
  })

  it('Should throw an exception on null-like input', async () => {
    await assert.isRejected(loader.load(undefined, '', new LoadingOptions({})), 'Expected non-null')
    await assert.isRejected(loader.load(null, '', new LoadingOptions({})), 'Expected non-null')
  })
})
