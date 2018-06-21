package cromwell.backend.google.pipelines.v2alpha1

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common._
import cromwell.backend.google.pipelines.common.io.PipelinesApiWorkingDisk
import cromwell.backend.standard.StandardAsyncExecutionActorParams
import cromwell.core.path.DefaultPathBuilder
import wom.core.FullyQualifiedName
import wom.expression.FileEvaluation
import wom.values.{GlobFunctions, WomFile, WomGlobFile, WomMaybeListedDirectory, WomMaybePopulatedFile, WomSingleFile, WomUnlistedDirectory}

class PipelinesApiAsyncBackendJobExecutionActor(standardParams: StandardAsyncExecutionActorParams) extends cromwell.backend .google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor(standardParams) {
  
  // The original implementation assumes the WomFiles are all WomMaybePopulatedFiles and wraps everything in a PipelinesApiFileInput
  // In v2 we can differentiate files from directories 
  override protected def pipelinesApiInputsFromWomFiles(inputName: String,
                                                        remotePathArray: Seq[WomFile],
                                                        localPathArray: Seq[WomFile],
                                                        jobDescriptor: BackendJobDescriptor): Iterable[PipelinesApiInput] = {
    (remotePathArray zip localPathArray) flatMap {
      case (remotePath: WomMaybeListedDirectory, localPath) =>
        maybeListedDirectoryToPipelinesParameters(inputName, remotePath, localPath.valueString)
      case (remotePath: WomUnlistedDirectory, localPath) =>
        Seq(PipelinesApiDirectoryInput(inputName, getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
      case (remotePath: WomMaybePopulatedFile, localPath) =>
        maybePopulatedFileToPipelinesParameters(inputName, remotePath, localPath.valueString)
      case (remotePath, localPath) =>
        Seq(PipelinesApiFileInput(inputName, getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
    }
  }

  // The original implementation recursively finds all non directory files, in V2 we can keep directory as is
  override protected def callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.localInputs map {
    case (key, womFile) =>
      key -> womFile.collectAsSeq({
        case womFile: WomFile if !inputsToNotLocalize.contains(womFile) => womFile
      })
  }

  // Simply create a PipelinesApiDirectoryOutput in v2 instead of globbing
  override protected def generateUnlistedDirectoryOutputs(unlistedDirectory: WomUnlistedDirectory, fileEvaluation: FileEvaluation): List[PipelinesApiOutput] = {
    val destination = callRootPath.resolve(unlistedDirectory.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(unlistedDirectory.value, runtimeAttributes.disks)
    val directoryOutput = PipelinesApiDirectoryOutput(makeSafeReferenceName(unlistedDirectory.value), destination, relpath, disk, fileEvaluation.optional, fileEvaluation.secondary)
    List(directoryOutput)
  }

  // De-localize the glob directory as a PipelinesApiDirectoryOutput instead of using * pattern match
  override def generateGlobFileOutputs(womFile: WomGlobFile): List[PipelinesApiOutput] = {
    val globName = GlobFunctions.globName(womFile.value)
    val globDirectory = globName + "/"
    val globListFile = globName + ".list"
    val gcsGlobDirectoryDestinationPath = callRootPath.resolve(globDirectory)
    val gcsGlobListFileDestinationPath = callRootPath.resolve(globListFile)

    val (_, globDirectoryDisk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)

    // We need both the glob directory and the glob list:
    List(
      // The glob directory:
      PipelinesApiDirectoryOutput(makeSafeReferenceName(globDirectory), gcsGlobDirectoryDestinationPath, DefaultPathBuilder.get(globDirectory), globDirectoryDisk, optional = false, secondary = false),
      // The glob list file:
      PipelinesApiFileOutput(makeSafeReferenceName(globListFile), gcsGlobListFileDestinationPath, DefaultPathBuilder.get(globListFile), globDirectoryDisk, optional = false, secondary = false)
    )
  }
  
  private def maybePopulatedFileToPipelinesParameters(inputName: String, maybePopulatedFile: WomMaybePopulatedFile, localPath: String) = {
    val secondaryFiles = maybePopulatedFile.secondaryFiles.flatMap({ secondaryFile =>
      pipelinesApiInputsFromWomFiles(secondaryFile.valueString, List(secondaryFile), List(relativeLocalizationPath(secondaryFile)), jobDescriptor)
    })

    Seq(PipelinesApiFileInput(inputName, getPath(maybePopulatedFile.valueString).get, DefaultPathBuilder.get(localPath), workingDisk)) ++ secondaryFiles
  }
  
  private def maybeListedDirectoryToPipelinesParameters(inputName: String, womMaybeListedDirectory: WomMaybeListedDirectory, localPath: String) = womMaybeListedDirectory match {
     // If there is a path, simply localize as a directory
    case WomMaybeListedDirectory(Some(path), _, _) =>
      List(PipelinesApiDirectoryInput(inputName, getPath(path).get, DefaultPathBuilder.get(localPath), workingDisk))

    // If there is a listing, recurse and call pipelinesApiInputsFromWomFiles on all the listed files
    case WomMaybeListedDirectory(_, Some(listing), _) if listing.nonEmpty =>
      listing.flatMap({
        case womFile: WomFile if isAdHocFile(womFile) =>
          pipelinesApiInputsFromWomFiles(makeSafeReferenceName(womFile.valueString), List(womFile), List(fileName(womFile)), jobDescriptor)
        case womFile: WomFile =>
          pipelinesApiInputsFromWomFiles(makeSafeReferenceName(womFile.valueString), List(womFile), List(relativeLocalizationPath(womFile)), jobDescriptor)
      })
    case _ => List.empty
  }

  override def generateSingleFileOutputs(womFile: WomSingleFile, fileEvaluation: FileEvaluation) = {
    val (relpath, disk) = relativePathAndAttachedDisk(womFile.value, runtimeAttributes.disks)
    // If the file is on a custom mount point, resolve it so that the full mount path will show up in the cloud path
    // For the default one (cromwell_root), the expectation is that it does not appear
    val mountedPath = if (!disk.mountPoint.isSamePathAs(PipelinesApiWorkingDisk.Default.mountPoint)) disk.mountPoint.resolve(relpath) else relpath
    // Normalize the local path (to get rid of ".." and "."). Also strip any potential leading / so that it gets appended to the call root
    val normalizedPath = mountedPath.normalize().pathAsString.stripPrefix("/")
    val destination = callRootPath.resolve(normalizedPath)
    val jesFileOutput = PipelinesApiFileOutput(makeSafeReferenceName(womFile.value), destination, relpath, disk, fileEvaluation.optional, fileEvaluation.secondary)
    List(jesFileOutput)
  }
}
