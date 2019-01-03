package cromwell.engine.workflow.lifecycle.execution.ejea

private[ejea] sealed trait ExpectOne[+A] {
  def checkIt(block: A => Any): Unit = throw new IllegalStateException("An ExpectOne must have exactly one element for checkIt to work")
  def hasExactlyOne: Boolean
  def foundOne[B >: A](theFoundOne: B) = this match {
    case NothingYet => GotOne(theFoundOne)
    case GotOne(theOriginalOne) => GotTooMany(List(theOriginalOne, theFoundOne))
    case GotTooMany(theOnes) => GotTooMany(theOnes :+ theFoundOne)
  }
}

private[ejea] case object NothingYet extends ExpectOne[scala.Nothing] {
  override def hasExactlyOne = false
}

private[ejea] case class GotOne[+A](theOne: A) extends ExpectOne[A] {
  override def checkIt(block: A => Any): Unit = { block(theOne); () }
  override def hasExactlyOne = true
}

private[ejea] case class GotTooMany[+A](theOnes: List[A]) extends ExpectOne[A] {
  override def hasExactlyOne = false
}