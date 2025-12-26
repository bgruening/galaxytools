import { assert } from 'chai'
import { TypeGuards } from '../util/Internal'

describe('Test Typeguards', () => {
  describe('Int', () => {
    it('Should return true', () => {
      assert.equal(TypeGuards.Int(2), true)
      assert.equal(TypeGuards.Int(0), true)
    })

    it('Should return false', () => {
      assert.equal(TypeGuards.Int(2.2), false)
      assert.equal(TypeGuards.Int('2.2'), false)
      assert.equal(TypeGuards.Int([2]), false)
      assert.equal(TypeGuards.Int({}), false)
      assert.equal(TypeGuards.Int(null), false)
      assert.equal(TypeGuards.Int(undefined), false)
    })
  })
  describe('Float', () => {
    it('Should return true', () => {
      assert.equal(TypeGuards.Float(2.0), true)
      assert.equal(TypeGuards.Float(2), true)
      assert.equal(TypeGuards.Float(0), true)
    })

    it('Should return false', () => {
      assert.equal(TypeGuards.Float([2]), false)
      assert.equal(TypeGuards.Float('2.2'), false)
      assert.equal(TypeGuards.Float({}), false)
      assert.equal(TypeGuards.Float(null), false)
      assert.equal(TypeGuards.Float(undefined), false)
    })
  })

  describe('Bool', () => {
    it('Should return true', () => {
      assert.equal(TypeGuards.Bool(true), true)
      assert.equal(TypeGuards.Bool(false), true)
    })

    it('Should return false', () => {
      assert.equal(TypeGuards.Bool([1]), false)
      assert.equal(TypeGuards.Bool('1'), false)
      assert.equal(TypeGuards.Bool(1), false)
      assert.equal(TypeGuards.Bool({}), false)
      assert.equal(TypeGuards.Bool(null), false)
      assert.equal(TypeGuards.Bool(undefined), false)
    })
  })

  describe('String', () => {
    it('Should return true', () => {
      assert.equal(TypeGuards.String('2.2'), true)
      assert.equal(TypeGuards.String(''), true)
      assert.equal(TypeGuards.String('test'), true)
    })

    it('Should return false', () => {
      assert.equal(TypeGuards.String([2]), false)
      assert.equal(TypeGuards.String(2), false)
      assert.equal(TypeGuards.String({}), false)
      assert.equal(TypeGuards.String(null), false)
      assert.equal(TypeGuards.String(undefined), false)
    })
  })

  describe('Undefined', () => {
    it('Should return true', () => {
      assert.equal(TypeGuards.Undefined(undefined), true)
      assert.equal(TypeGuards.Undefined(null), true)
    })

    it('Should return false', () => {
      assert.equal(TypeGuards.Undefined([1]), false)
      assert.equal(TypeGuards.Undefined('1'), false)
      assert.equal(TypeGuards.Undefined(1), false)
      assert.equal(TypeGuards.Undefined(1.1), false)
      assert.equal(TypeGuards.Undefined({}), false)
    })
  })

  describe('Dictionary', () => {
    it('Should return true', () => {
      assert.equal(TypeGuards.isDictionary({}), true)
      assert.equal(TypeGuards.isDictionary({ test: 'test' }), true)
    })

    it('Should return false', () => {
      assert.equal(TypeGuards.isDictionary([]), false)
      assert.equal(TypeGuards.isDictionary('1'), false)
      assert.equal(TypeGuards.isDictionary(1), false)
      assert.equal(TypeGuards.isDictionary(1.1), false)
      assert.equal(TypeGuards.isDictionary(undefined), false)
      assert.equal(TypeGuards.isDictionary(null), false)
    })
  })
})
