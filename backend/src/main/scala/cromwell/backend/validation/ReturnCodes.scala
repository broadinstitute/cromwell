package cromwell.backend.validation

import wom.types._

object ReturnCodes {
  val validWdlTypes = Set[WomType](WomArrayType(WomIntegerType), WomStringType, WomIntegerType)
}

/**
 * Decides if a call/job continues upon a specific return code.
 */
sealed trait ReturnCodes extends ReturnCode {

  /**
   * Returns true if the call is a success based on the return code.
   *
   * @param returnCode Return code from the process / script.
   * @return True if the call is a success.
   */
  final override def continueFor(returnCode: Int): Boolean =
    this match {
      case ReturnCodesString(continue) => continue.equals("*") || returnCode == 0
      case _ => super.continueFor(returnCode)
    }
}

/**
 * Continues based on a string, if "*" all return codes continue.
 * @param returnCode If "*", all return codes are valid for continuing.
 */
case class ReturnCodesString(returnCode: String) extends ReturnCodes {
  override def toString = returnCode
}
