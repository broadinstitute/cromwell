package cromwell.engine.workflow.lifecycle.execution.ejea

private[ejea] sealed trait ExpectOne[+A] {
  def checkIt(block: A => Any): Unit = throw new IllegalStateException("This ExpectOne must have exactly one element for checkIt to work")
  def checkLatest(block: A => Any): Unit = throw new IllegalStateException("This ExpectOne must have at least one element for checkLatest to work")
  def hasAtLeastOne: Boolean
  def hasExactlyOne: Boolean
  def foundOne[B >: A](theFoundOne: B) = this match {
    case NothingYet => GotOne(theFoundOne)
    case GotOne(theOriginalOne) => GotTooMany(List(theOriginalOne, theFoundOne))
    case GotTooMany(theOnes) => GotTooMany(theOnes :+ theFoundOne)
  }
}

private[ejea] case object NothingYet extends ExpectOne[scala.Nothing] {
  override def hasAtLeastOne = false
  override def hasExactlyOne = false
}

private[ejea] case class GotOne[+A](theOne: A) extends ExpectOne[A] {
  override def checkIt(block: A => Any): Unit = { block(theOne); () }
  override def checkLatest(block: A => Any): Unit = checkIt(block)
  override def hasAtLeastOne = true
  override def hasExactlyOne = true
}

private[ejea] case class GotTooMany[+A](theOnes: List[A]) extends ExpectOne[A] {
  override def checkLatest(block: A => Any): Unit = { block(theOnes.last); ()  }
  override def hasAtLeastOne = true
  override def hasExactlyOne = false
}