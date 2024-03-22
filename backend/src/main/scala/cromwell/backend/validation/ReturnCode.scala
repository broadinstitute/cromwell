package cromwell.backend.validation

trait ReturnCode {

  /**
   * Returns true if the call is a success based on the return code.
   *
   * @param returnCode Return code from the process / script.
   * @return True if the call is a success.
   */
  def continueFor(returnCode: Int): Boolean =
    this match {
      case ReturnCodeSet(returnCodes) => returnCodes.contains(returnCode)
    }
}

/**
 * Continues only if the call/job return code is found in returnCodes.
 * @param returnCodes Inclusive set of return codes that specify a job success.
 */
case class ReturnCodeSet(returnCodes: Set[Int]) extends ReturnCode {
  override def toString = returnCodes match {
    case single if single.size == 1 => returnCodes.head.toString
    case multiple => s"[${multiple.mkString(",")}]"
  }

  /**
   * Returns true if the call is a success based on the return code.
   *
   * @param returnCode Return code from the process / script.
   * @return True if the call is a success.
   */
  override def continueFor(returnCode: Int): Boolean = super.continueFor(returnCode)
}
