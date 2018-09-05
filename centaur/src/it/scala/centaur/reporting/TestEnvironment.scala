package centaur.reporting

/**
  * Information about a test.
  *
  * @param name The test name.
  * @param retries The total number of retries.
  * @param attempt The zero based attempt.
  */
case class TestEnvironment(name: String, retries: Int, attempt: Int)
