package cromwell.backend.google.pipelines.v2alpha1

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common._
import cromwell.backend.standard.StandardAsyncExecutionActorParams
import cromwell.core.path.DefaultPathBuilder
import wom.core.FullyQualifiedName
import wom.values.{GlobFunctions, WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomUnlistedDirectory}

import scala.language.postfixOps

class PipelinesApiAsyncBackendJobExecutionActor(standardParams: StandardAsyncExecutionActorParams) extends cromwell.backend .google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor(standardParams) {
  
  // The original implementation assumes the WomFiles are all WomMaybePopulatedFiles and wraps everything in a PipelinesApiFileInput
  // In v2 we can differentiate files from directories 
  override protected def pipelinesApiInputsFromWomFiles(jesNamePrefix: String,
                                                        remotePathArray: Seq[WomFile],
                                                        localPathArray: Seq[WomFile],
                                                        jobDescriptor: BackendJobDescriptor): Iterable[PipelinesApiInput] = {
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath: WomMaybeListedDirectory, localPath), index) =>
        maybeListedDirectoryToPipelinesParameters(jesNamePrefix, remotePath, index, localPath.valueString)
      case ((remotePath: WomUnlistedDirectory, localPath), index) =>
        Seq(PipelinesApiDirectoryInput(s"$jesNamePrefix-$index", getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
      case ((remotePath: WomMaybePopulatedFile, localPath), index) =>
        maybePopulatedFileToPipelinesParameters(jesNamePrefix, remotePath, index, localPath.valueString)
      case ((remotePath, localPath), index) =>
        Seq(PipelinesApiFileInput(s"$jesNamePrefix-$index", getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
    }
  }

  // The original implementation recursively finds all non directory files, in V2 we can keep directory as is
  override protected def callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.fullyQualifiedInputs mapValues {
    womFile =>
      womFile collectAsSeq {
        case womFile: WomFile => womFile
      }
  } map identity // <-- unlazy the mapValues

  // Simply create a PipelinesApiDirectoryOutput in v2 instead of globbing
  override protected def generateUnlistedDirectoryOutputs(unlistedDirectory: WomUnlistedDirectory, optional: Boolean): List[PipelinesApiOutput] = {
    val destination = callRootPath.resolve(unlistedDirectory.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(unlistedDirectory.value, runtimeAttributes.disks)
    val directoryOutput = PipelinesApiDirectoryOutput(makeSafeJesReferenceName(unlistedDirectory.value), destination, relpath, disk, optional)
    List(directoryOutput)
  }

  // Delocalize the glob directory as a PipelinesApiDirectoryOutput instead of using * pattern match
  override def generateJesGlobFileOutputs(womFile: WomGlobFile): List[PipelinesApiOutput] = {
    val globName = GlobFunctions.globName(womFile.value)
    val globDirectory = globName + "/"
    val globListFile = globName + ".list"
    val gcsGlobDirectoryDestinationPath = callRootPath.resolve(globDirectory)
    val gcsGlobListFileDestinationPath = callRootPath.resolve(globListFile)

    val (_, globDirectoryDisk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)

    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      PipelinesApiDirectoryOutput(makeSafeJesReferenceName(globDirectory), gcsGlobDirectoryDestinationPath, DefaultPathBuilder.get(globDirectory), globDirectoryDisk, optional = false),
      // The glob list file:
      PipelinesApiFileOutput(makeSafeJesReferenceName(globListFile), gcsGlobListFileDestinationPath, DefaultPathBuilder.get(globListFile), globDirectoryDisk, optional = false)
    )
  }
  
  private def maybePopulatedFileToPipelinesParameters(jesNamePrefix: String, maybePopulatedFile: WomMaybePopulatedFile, index: Int, localPath: String) = {
    val secondaryFiles = maybePopulatedFile.secondaryFiles.flatMap({ secondaryFile =>
      pipelinesApiInputsFromWomFiles(secondaryFile.valueString, List(secondaryFile), List(relativeLocalizationPath(secondaryFile)), jobDescriptor)
    })

    Seq(PipelinesApiFileInput(s"$jesNamePrefix-$index", getPath(maybePopulatedFile.valueString).get, DefaultPathBuilder.get(localPath), workingDisk)) ++ secondaryFiles
  }
  
  private def maybeListedDirectoryToPipelinesParameters(jesNamePrefix: String, womMaybeListedDirectory: WomMaybeListedDirectory, index: Int, localPath: String) = womMaybeListedDirectory match {
     // If there is a path, simply localize as a directory
    case WomMaybeListedDirectory(Some(path), _, _) =>
      List(PipelinesApiDirectoryInput(s"$jesNamePrefix-$index", getPath(path).get, DefaultPathBuilder.get(localPath), workingDisk))

    // If there is a listing, recurse and call pipelinesApiInputsFromWomFiles on all the listed files
    case WomMaybeListedDirectory(_, Some(listing), _) if listing.nonEmpty =>
      listing.flatMap({
        case womFile: WomFile if isAdHocFile(womFile) =>
          pipelinesApiInputsFromWomFiles(makeSafeJesReferenceName(womFile.valueString), List(womFile), List(fileName(womFile)), jobDescriptor)
        case womFile: WomFile =>
          pipelinesApiInputsFromWomFiles(makeSafeJesReferenceName(womFile.valueString), List(womFile), List(relativeLocalizationPath(womFile)), jobDescriptor)
      })
    case _ => List.empty
  }
}
