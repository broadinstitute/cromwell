package cromwell.engine.finalcall

import java.nio.file.{FileAlreadyExistsException, Files}

import cromwell.engine.backend.{BackendCallJobDescriptor, WorkflowDescriptor}
import cromwell.logging.WorkflowLogger
import cromwell.util.{SimpleExponentialBackoff, TryUtil}
import wdl4s._
import wdl4s.types.{WdlArrayType, WdlStringType}
import wdl4s.values.{WdlArray, WdlString, WdlValue}

import scala.concurrent.duration._
import scala.util.Try

/**
  * HACK: Copies files within this process, until we can figure out a way to generate an command line that can copy
  * to/from a source/destination via a generated external process.
  */
object FinalCallCopyHack {
  private val CopiesParam = "pairs_of_src_dst_to_copy"

  /** Generate an inert task, since maybeCopyFiles will do the actual copying _before_ this task runs. */
  def wdlSource(finalCallName: String) =
    s"""
       |task $finalCallName {
       |  Array[Array[String]] $CopiesParam
       |  command {
       |    echo running final call $finalCallName
       |    echo TODO: copy these files out of process
       |    cat $${write_tsv($CopiesParam)}
       |  }
       |}
     """.stripMargin

  /** Takes a sequence of copies and converts it into a CallInputs, the inputs later compatible with maybeCopyFiles. */
  def toCallInputs(copies: Seq[FinalCallCopy]): CallInputs = {
    Map(
      CopiesParam -> WdlArray(WdlArrayType(WdlArrayType(WdlStringType)),
        copies map { copy =>
          WdlArray(WdlArrayType(WdlStringType),
            Seq(WdlString(copy.source.toString), WdlString(copy.destination.toString)))
        }
      )
    )
  }

  /**
    * HACK: Possibly copy the files now, in process, since the copying files via an external process is not implemented,
    * yet.
    */
  def maybeCopyFiles(logger: WorkflowLogger, jobDescriptor: BackendCallJobDescriptor): Unit = {
    import FinalCall._
    if (jobDescriptor.call.unqualifiedName.isFinalCall) {
      val callFqn = jobDescriptor.call.fullyQualifiedName
      val copies = toCopies(jobDescriptor.workflowDescriptor, jobDescriptor.locallyQualifiedInputs)
      logger.info(s"Copying files (via in process APIs) for $callFqn")
      TryUtil.sequence(copies map tryCopyFile(logger), s"Unable to copy files for $callFqn: ")
      // ignore the result of the try on this final call
    }
  }

  /** Convert our hacked call inputs _back_ into a sequence of final call copies. */
  private def toCopies(workflowDescriptor: WorkflowDescriptor, callInputs: CallInputs): Seq[FinalCallCopy] = {
    val copyPairs = callInputs.get(CopiesParam) collect {
      case WdlArray(WdlArrayType(WdlArrayType(WdlStringType)), arrCopies) => arrCopies collect {
        case WdlArray(WdlArrayType(WdlStringType), arrSrcDst) => arrSrcDst
      }
    } getOrElse Seq.empty

    copyPairs map toCopy(workflowDescriptor)
  }

  /**
    * Use the filesystems embedded within the workflow descriptor to create a Path that knows how to copy to/from our
    * source/destination.
    */
  private def toCopy(workflowDescriptor: WorkflowDescriptor)(pairSrcDst: Seq[WdlValue]): FinalCallCopy = {
    pairSrcDst match {
      case Seq(source, destination) =>
        val srcPath = workflowDescriptor.toFileSystemPath(source.valueString)
        val dstPath = workflowDescriptor.toFileSystemPath(destination.valueString)
        FinalCallCopy(srcPath, dstPath)
    }
  }

  /** Try to copy a file, repeatedly. NOTE: Blocks the current thread while trying. */
  private def tryCopyFile(logger: WorkflowLogger)(copy: FinalCallCopy): Try[_] = {
    val src = copy.source
    val dest = copy.destination

    def tryCopy(): Unit = {
      logger.info(s"Trying to copy file $src to $dest")
      Files.createDirectories(dest.getParent)
      Files.copy(src, dest)
    }

    TryUtil.retryBlock(
      fn = (retries: Option[Unit]) => tryCopy(),
      retryLimit = Option(5),
      backoff = SimpleExponentialBackoff(5.seconds, 10.seconds, 1.1D),
      logger = logger,
      failMessage = Option(s"Failed to copy file $src to $dest"),
      /*
        TODO: Fix partial copies after we switch to CloudStorageFileSystemProvider

        https://github.com/GoogleCloudPlatform/gcloud-java/blob/95c98afa22f2eadfb61d5beff89bf0ecd65f03e3/gcloud-java-contrib/gcloud-java-nio/src/main/java/com/google/gcloud/storage/contrib/nio/CloudStorageFileSystemProvider.java

        If GCS NIO is released before this hack is fixed, replace checking for FAEE with
        StandardCopyOption.REPLACE_EXISTING. Right now skipping FAEE may lead to incomplete copies if the first attempt
        fails.
         */
      isFatal = (t: Throwable) => t.isInstanceOf[FileAlreadyExistsException]
    ) recover {
      case _: FileAlreadyExistsException =>
        logger.info(s"Tried to copy the same file multiple times. Skipping subsequent copies for $src")
    }
  }
}
