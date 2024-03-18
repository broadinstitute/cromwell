package cromwell.backend.validation

import wom.types._

object ReturnCodes {
  val validWdlTypes = Set[WomType](WomArrayType(WomIntegerType), WomStringType, WomIntegerType)
}

/**
 * Decides if a call/job continues upon a specific return code.
 */
sealed trait ReturnCodes {

  /**
   * Returns true if the call is a success based on the return code.
   *
   * @param returnCode Return code from the process / script.
   * @return True if the call is a success.
   */
  final def continueFor(returnCode: Int): Boolean =
    this match {
      case ReturnCodesString(continue) => continue.equals("*") || returnCode == 0
      case ReturnCodesSet(returnCodes) => returnCodes.contains(returnCode)
    }
}

/**
 * Continues based on a string, if "*" all return codes continue.
 * @param returnCode If "*", all return codes are valid for continuing.
 */
case class ReturnCodesString(returnCode: String) extends ReturnCodes {
  override def toString = returnCode
}

/**
 * Continues only if the call/job return code is found in returnCodes.
 * @param returnCodes Inclusive set of return codes that specify a job success.
 */
case class ReturnCodesSet(returnCodes: Set[Int]) extends ReturnCodes {
  override def toString = returnCodes match {
    case single if single.size == 1 => returnCodes.head.toString
    case multiple => s"[${multiple.mkString(",")}]"
  }
}
