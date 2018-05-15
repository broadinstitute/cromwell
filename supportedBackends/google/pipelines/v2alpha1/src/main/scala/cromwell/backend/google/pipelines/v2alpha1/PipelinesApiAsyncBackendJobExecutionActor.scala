package cromwell.backend.google.pipelines.v2alpha1

import cromwell.backend.BackendJobDescriptor
import cromwell.backend.google.pipelines.common._
import cromwell.backend.standard.StandardAsyncExecutionActorParams
import cromwell.core.path.DefaultPathBuilder
import wom.core.FullyQualifiedName
import wom.values.{GlobFunctions, WomFile, WomGlobFile, WomMaybeListedDirectory, WomUnlistedDirectory}

import scala.language.postfixOps

class PipelinesApiAsyncBackendJobExecutionActor(standardParams: StandardAsyncExecutionActorParams) extends cromwell.backend .google.pipelines.common.PipelinesApiAsyncBackendJobExecutionActor(standardParams) {
  override protected def jesInputsFromWomFiles(jesNamePrefix: String,
                                               remotePathArray: Seq[WomFile],
                                               localPathArray: Seq[WomFile],
                                               jobDescriptor: BackendJobDescriptor): Iterable[PipelinesApiInput] = {
    def maybeListedDirectory(womMaybeListedDirectory: WomMaybeListedDirectory, index: Int, localPath: String) = womMaybeListedDirectory match {
      case WomMaybeListedDirectory(Some(path), _, _) if isAdHocFile(womMaybeListedDirectory) =>
        List(PipelinesApiDirectoryInput(s"$jesNamePrefix-$index", getPath(path).get, DefaultPathBuilder.get(localPath), workingDisk))
      case WomMaybeListedDirectory(_, Some(listing), _) if listing.nonEmpty =>
        listing.flatMap({
          case womFile: WomFile if isAdHocFile(womFile) =>
            jesInputsFromWomFiles(makeSafeJesReferenceName(womMaybeListedDirectory.valueString), List(womFile), List(fileName(womFile)), jobDescriptor) 
          case womFile: WomFile =>
            jesInputsFromWomFiles(makeSafeJesReferenceName(womMaybeListedDirectory.valueString), List(womFile), List(relativeLocalizationPath(womFile)), jobDescriptor) 
        })
      case WomMaybeListedDirectory(Some(path), _, _) =>
        List(PipelinesApiDirectoryInput(s"$jesNamePrefix-$index", getPath(path).get, DefaultPathBuilder.get(localPath), workingDisk))
      case _ => List.empty
    }
    
    (remotePathArray zip localPathArray zipWithIndex) flatMap {
      case ((remotePath: WomMaybeListedDirectory, localPath), index) =>
        maybeListedDirectory(remotePath, index, localPath.valueString)
      case ((remotePath: WomUnlistedDirectory, localPath), index) =>
        Seq(PipelinesApiDirectoryInput(s"$jesNamePrefix-$index", getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
      case ((remotePath, localPath), index) =>
        Seq(PipelinesApiFileInput(s"$jesNamePrefix-$index", getPath(remotePath.valueString).get, DefaultPathBuilder.get(localPath.valueString), workingDisk))
    }
  }

  override protected def callInputFiles: Map[FullyQualifiedName, Seq[WomFile]] = jobDescriptor.fullyQualifiedInputs mapValues {
    womFile =>
      womFile collectAsSeq {
        case womFile: WomFile => womFile
      }
  } map identity // <-- unlazy the mapValues
  
 override protected def generateUnlistedDirectoryOutputs(unlistedDirectory: WomUnlistedDirectory, optional: Boolean): List[PipelinesApiOutput] = {
    val destination = callRootPath.resolve(unlistedDirectory.value.stripPrefix("/"))
    val (relpath, disk) = relativePathAndAttachedDisk(unlistedDirectory.value, runtimeAttributes.disks)
    val directoryOutput = PipelinesApiDirectoryOutput(makeSafeJesReferenceName(unlistedDirectory.value), destination, relpath, disk, optional)
    List(directoryOutput)
  }
  
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
}
