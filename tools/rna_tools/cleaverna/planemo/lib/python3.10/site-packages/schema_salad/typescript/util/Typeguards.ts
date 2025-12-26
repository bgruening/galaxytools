import { Dictionary } from './Internal'

export function Int (doc: any): boolean {
  return typeof doc === 'number' && Number.isInteger(doc)
}

export function Float (doc: any): boolean {
  return typeof doc === 'number'
}

export function Bool (doc: any): boolean {
  return typeof doc === 'boolean'
}

export function String (doc: any): boolean {
  return typeof doc === 'string'
}

export function Undefined (doc: any): boolean {
  return doc == null
}

export function isDictionary (doc: any): doc is Dictionary {
  return (doc != null && typeof doc === 'object' && !Array.isArray(doc))
}
