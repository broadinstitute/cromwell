package cromwell.backend.validation

import wdl4s.types.{WdlArrayType, WdlBooleanType, WdlIntegerType, WdlType}

object ContinueOnReturnCode {
  val validWdlTypes = Set[WdlType](WdlArrayType(WdlIntegerType), WdlBooleanType, WdlIntegerType)
}

/**
  * Decides if a call/job continues upon a specific return code.
  */
sealed trait ContinueOnReturnCode {
  /**
    * Returns true if the call is a success based on the return code.
    *
    * @param returnCode Return code from the process / script.
    * @return True if the call is a success.
    */
  final def continueFor(returnCode: Int): Boolean = {
    this match {
      case ContinueOnReturnCodeFlag(continue) => continue || returnCode == 0
      case ContinueOnReturnCodeSet(returnCodes) => returnCodes.contains(returnCode)
    }
  }
}

/**
  * Continues based on a generic true / false flag, that when false, only zero return codes continue.
  * @param continue If true, all return codes are valid for continuing.
  */
case class ContinueOnReturnCodeFlag(continue: Boolean) extends ContinueOnReturnCode {
  override def toString = continue.toString
}

/**
  * Continues only if the call/job return code is found in returnCodes.
  * @param returnCodes Inclusive set of return codes that specify a job success.
  */
case class ContinueOnReturnCodeSet(returnCodes: Set[Int]) extends ContinueOnReturnCode {
  override def toString = returnCodes match {
    case single if single.size == 1 => returnCodes.head.toString
    case multiple => s"[${multiple.mkString(",")}]"
  }
}
