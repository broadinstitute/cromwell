package wdl.transforms.base.ast2wdlom

import simulacrum.typeclass

@typeclass
trait GenericAstNodeMaker[A] {
  def generify(a: A): GenericAstNode
}

@typeclass
trait GenericAstMaker[A] {
  def generify(a: A): GenericAst
}

@typeclass
trait GenericAstListMaker[A] {
  def generify(a: A): GenericAstList
}

@typeclass
trait GenericTerminalMaker[A] {
  def generify(a: A): GenericTerminal
}
