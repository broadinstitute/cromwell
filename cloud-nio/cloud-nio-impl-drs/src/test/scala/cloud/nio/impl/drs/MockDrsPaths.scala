package cloud.nio.impl.drs

import scala.concurrent.duration._

object MockDrsPaths {
  val marthaUrl = "http://mock.martha"

  val mockDrsConfig: DrsConfig = DrsConfig(MockDrsPaths.marthaUrl, 1, 1.seconds, 1.seconds, 1d, 0d)

  val mockToken = "mock.token"

  private val drsPathPrefix = "drs://drs-host"

  val drsRelativePath = "drs-host/4d427aa3-5640-4f00-81ae-c33443f84acf/f3b148ac-1802-4acc-a0b9-610ea266fb61"

  val drsReplacedChar = "drs-host/4d427aa3_5640_4f00_81ae_c33443f84acf/f3b148ac-1802-4acc-a0b9-610ea266fb61"

  val gcsRelativePathWithFileName = "drs-host/d7c75399-bcd3-4762-90e9-434de005679b/file.txt"

  val drsPathResolvingGcsPath = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-c33443f84acf"

  val drsPathWithNonPathChars = s"$drsPathPrefix/4d427aa3_5640_4f00_81ae_c33443f84acf"

  val drsPathResolvingWithFileName = s"$drsPathPrefix/d7c75399-bcd3-4762-90e9-434de005679b"

  val drsPathResolvingToNoGcsPath = s"$drsPathPrefix/226686cf-22c9-4472-9f79-7a0b0044f253"

  val drsPathNotExistingInMartha = s"$drsPathPrefix/5e21b8c3-8eda-48d5-9a04-2b32e1571765"
}
