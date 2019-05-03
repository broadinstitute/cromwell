package wdl.transforms.base.ast2wdlom

import cats.instances.vector._
import cats.syntax.either._
import cats.syntax.traverse._
import common.validation.Checked._
import common.Checked
import common.transforms.CheckedAtoB

trait GenericAstNode {
  def astListAsVector: Checked[Vector[GenericAstNode]] = this match {
    case list: GenericAstList => list.astNodeList.toVector.validNelCheck
    // Parsing for a list of terminals returns a single terminal element if there was only one
    case t: GenericTerminal => Vector(t).validNelCheck
    case _ => s"Invalid target for astListAsVector: ${getClass.getSimpleName}".invalidNelCheck
  }

  def allTerminals: Set[GenericTerminal] = this match {
    case t: GenericTerminal => Set(t)
    case ast: GenericAst => ast.getAttributes.values.toSet[GenericAstNode] flatMap { _.allTerminals }
    case astList: GenericAstList => astList.astNodeList.toSet[GenericAstNode] flatMap { _.allTerminals }
  }

  def firstTerminal: Option[GenericTerminal] = {
    def foldFunction(acc: Option[GenericTerminal], next: GenericTerminal): Option[GenericTerminal] = acc match {
      case Some(t) if t.getLine > next.getLine || (t.getLine == next.getLine && t.getColumn > next.getColumn) => Some(next)
      case None => Some(next)
      case _ => acc
    }

    allTerminals.foldLeft[Option[GenericTerminal]](None)(foldFunction)
  }

  // We would actually like to get the extent that is covered by this AST. A representation
  // like:  (startLine, startColumn, endLine, endColumn) would be efficient, while conveying
  // all the information needed for downstream analysis phases. However, it turns out that
  // getting accurate information out of Hermes is not that simple. For now, we just
  // get the initial source line, which is -more or less- accurate.
  def getSourceLine: Option[Int] = {
    firstTerminal map {t => t.getLine }
  }

  def lineAndColumnString = firstTerminal map { t => s" at line ${t.getLine} column ${t.getColumn}"} getOrElse("")
}

trait GenericAst extends GenericAstNode {
  def getAttribute(attr: String): GenericAstNode
  def getAttributes: Map[String, GenericAstNode]
  def getName: String

  private def getAttributeAsAstNodeVector(attr: String, optional: Boolean): Checked[Vector[GenericAstNode]] = {
    Option(getAttribute(attr)) match {
      case Some(attributeNode) => attributeNode.astListAsVector
      case None if optional => Vector.empty.validNelCheck
      case None => s"No expected attribute '$attr' found".invalidNelCheck
    }
  }

  /**
    * Will get an attribute on this Ast as an AstNode and then convert that into a single element of
    * the required type.
    */
  def getAttributeAs[A](attr: String)(implicit toA: CheckedAtoB[GenericAstNode, A]): Checked[A] = {
    val attribute = Option(getAttribute(attr))
    // Note: if you see one of these in the wild and can recreate it, you might want to try switching in this more complete message:
    //  "No attribute '$attr' found on Ast '$getName'. Did you mean: ${getAttributes.keys.mkString(", ")}"
    attribute.map(toA.run).getOrElse(s"No expected attribute '$attr' found".invalidNelCheck)
  }

  /**
    * Will get an attribute on this Ast as an AstList and then convert that into a vector of Ast
    * @param attr The attribute to read from this Ast
    */
  def getAttributeAsVector[A](attr: String, optional: Boolean = false)(implicit toA: CheckedAtoB[GenericAstNode, A]): Checked[Vector[A]] = {
    for {
      asVector <- getAttributeAsAstNodeVector(attr, optional)
      // This toValidated/toEither dance is necessary to
      // (1) collect all errors from the traverse as an ErrorOr, then
      // (2) convert back into a Checked for the flatMap
      result <- asVector.traverse(item => toA.run(item).toValidated).toEither
    } yield result
  }

  /**
    * Will get an attribute on this Ast as an AstList and then convert that into a vector of Ast
    * @param attr The attribute to read from this Ast
    */
  def getAttributeAsVectorF[A](attr: String, optional: Boolean = false)(toA: GenericAstNode => Checked[A]): Checked[Vector[A]] = {
    for {
      asVector <- getAttributeAsAstNodeVector(attr, optional)
      // This toValidated/toEither dance is necessary to
      // (1) collect all errors from the traverse as an ErrorOr, then
      // (2) convert back into a Checked for the flatMap
      result <- asVector.traverse(item => toA(item).toValidated).toEither
    } yield result
  }

  /**
    * Gets an attribute on this Ast as an Optional Ast, returns an empty Option if the attribute is empty.
    */
  def getAttributeAsOptional[A](attr: String)(implicit toA: CheckedAtoB[GenericAstNode, A]): Checked[Option[A]] =  {
    Option(getAttribute(attr)) match {
      case None => None.validNelCheck
      case Some(attribute) => toA.run(attribute).map(Option.apply)
    }
  }
}

trait GenericAstList extends GenericAstNode {
  def astNodeList: Seq[GenericAstNode]
}

trait GenericTerminal extends GenericAstNode {
  def getSourceString: String
  def getTerminalStr: String
  def getLine: Int
  def getColumn: Int
}
