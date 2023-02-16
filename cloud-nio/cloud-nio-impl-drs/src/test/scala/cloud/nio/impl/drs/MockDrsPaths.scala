package cloud.nio.impl.drs

import scala.concurrent.duration._

object MockDrsPaths {
  val drsResolverUrl = "http://mock.drshub"

  val mockDrsConfig: DrsConfig = DrsConfig(MockDrsPaths.drsResolverUrl, 1, 1.seconds, 1.seconds, 1d, 0d)

  val mockToken = "mock.token"

  val DrsLocalizationPathsContainer = "drs_localization_paths"

  private val drsPathPrefix = "drs://drs-host"

  val drsRelativePath = "drs-host/4d427aa3-5640-4f00-81ae-c33443f84acf/f3b148ac-1802-4acc-a0b9-610ea266fb61"

  val drsReplacedChar = "drs-host/4d427aa3_5640_4f00_81ae_c33443f84acf/f3b148ac-1802-4acc-a0b9-610ea266fb61"

  val gcsRelativePathWithFileName = "drs-host/d7c75399-bcd3-4762-90e9-434de005679b/file.txt"

  val gcsRelativePathWithFileNameFromLocalizationPath = s"$DrsLocalizationPathsContainer/dir/subdir/file.txt"

  val gcsRelativePathWithFileNameFromAllThePaths = s"$DrsLocalizationPathsContainer/dir/subdir/file.txt"

  val drsPathResolvingGcsPath = s"$drsPathPrefix/4d427aa3-5640-4f00-81ae-c33443f84acf"

  val drsPathWithNonPathChars = s"$drsPathPrefix/4d427aa3_5640_4f00_81ae_c33443f84acf"

  val drsPathResolvingWithFileName = s"$drsPathPrefix/d7c75399-bcd3-4762-90e9-434de005679b"

  val drsPathResolvingWithLocalizationPath = s"$drsPathPrefix/1e7ecfa6-2a77-41d7-a251-38a2f4919842"

  val drsPathResolvingWithAllThePaths = s"$drsPathPrefix/0524678a-365e-42f3-a1e7-e4c6ac499b35"

  val drsPathResolvingToNoGcsPath = s"$drsPathPrefix/226686cf-22c9-4472-9f79-7a0b0044f253"

  val drsPathNotExistingInDrsResolver = s"$drsPathPrefix/5e21b8c3-8eda-48d5-9a04-2b32e1571765"
}
