import { assert } from 'chai'
import { shortname } from '../util/Internal'

describe('Test Utilities', () => {
  describe('shortname', () => {
    it('test1', async () => {
      assert.equal(
        shortname('file:/Users/jdidion/projects/cwlScala/target/test-classes/CommandLineTools/conformance/#anon_enum_inside_array_inside_schemadef.cwl/first/user_type_2/species/homo_sapiens'),
        'homo_sapiens'
      )
    }),
    it('test2', async () => {
      assert.equal(
        shortname('file:///home/michael/cwljava/src/test/resources/org/w3id/cwl/cwl1_2/utils/valid_anon_enum_inside_array_inside_schemadef.cwl#vcf2maf_params/ncbi_build/GRCh37'),
        'GRCh37'
      )
    }),
    it('test3', async () => {
      assert.equal(
        shortname('http://example.com/foo'), 'foo'
      )
    }),
    it('test4', async () => {
      assert.equal(
        shortname('http://example.com/#bar'), 'bar'
      )
    }),
    it('test5', async () => {
      assert.equal(
        shortname('http://example.com/foo/bar'), 'bar'
      )
    }),
    it('test6', async () => {
      assert.equal(
        shortname('http://example.com/foo#bar'), 'bar'
      )
    }),
    it('test7', async () => {
      assert.equal(
        shortname('http://example.com/#foo/bar'), 'bar'
      )
    }),
    it('test8', async () => {
      assert.equal(
        shortname('http://example.com/foo#bar/baz'), 'baz'
      )
    })
  })
})
