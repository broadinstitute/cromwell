package cromwell.backend.validation

import wom.types._

object ContinueOnReturnCode {
  val validWdlTypes = Set[WomType](WomArrayType(WomIntegerType), WomBooleanType, WomIntegerType)
}

/**
  * Decides if a call/job continues upon a specific return code.
  */
sealed trait ContinueOnReturnCode extends ReturnCode {

  /**
    * Returns true if the call is a success based on the return code.
    *
    * @param returnCode Return code from the process / script.
    * @return True if the call is a success.
    */
  final override def continueFor(returnCode: Int): Boolean =
    this match {
      case ContinueOnReturnCodeFlag(continue) => continue || returnCode == 0
      case _ => super.continueFor(returnCode)
    }
}

/**
  * Continues based on a generic true / false flag, that when false, only zero return codes continue.
  * @param continue If true, all return codes are valid for continuing.
  */
case class ContinueOnReturnCodeFlag(continue: Boolean) extends ContinueOnReturnCode {
  override def toString = continue.toString
}
