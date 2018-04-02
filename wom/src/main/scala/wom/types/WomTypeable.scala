package wom.types

@typeclass
trait WomTypeable[A] {
  def womType(a: A): WomType
}

