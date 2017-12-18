package cromwell.backend.impl.jes

import com.google.api.services.genomics.model.{LocalCopy, PipelineParameter}
import com.typesafe.config.ConfigFactory
import cromwell.backend.impl.jes.io.JesAttachedDisk
import cromwell.core.path.Path
import net.ceedubs.ficus.Ficus._


sealed trait JesParameter {
  def name: String
  def toGooglePipelineParameter: PipelineParameter
  def toGoogleRunParameter: String
}

object JesInput {
  // This is used for CWL conformance testing on PAPI. The current assumption is that the input files have already
  // been staged up to this path.
  lazy val defaultGcsPrefix = {
    val rawPrefix = ConfigFactory.load().as[Option[String]]("papi.default-input-gcs-prefix")
    rawPrefix map { p => if (p.endsWith("/")) p else p + "/" }
  }

  private[jes] def prefixIfNecessary(url: String): String = defaultGcsPrefix match {
    case Some(prefix) if !url.startsWith("gs://") => prefix + url
    case _ => url
  }
}

sealed trait JesInput extends JesParameter

final case class JesFileInput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesInput {
  def toGooglePipelineParameter = {
    new PipelineParameter().setName(name).setLocalCopy(
      new LocalCopy().setDisk(mount.name).setPath(local.pathAsString)
    )
  }

  val toGoogleRunParameter: String = JesInput.prefixIfNecessary(gcs)

  def containerPath: Path = mount.mountPoint.resolve(local)
}

final case class JesLiteralInput(name: String, value: String) extends JesInput {
  def toGooglePipelineParameter = new PipelineParameter().setName(name)
  val toGoogleRunParameter: String = value
}

final case class JesFileOutput(name: String, gcs: String, local: Path, mount: JesAttachedDisk) extends JesParameter {
  def toGooglePipelineParameter = {
    new PipelineParameter().setName(name).setLocalCopy(
      new LocalCopy().setDisk(mount.name).setPath(local.pathAsString)
    )
  }
  val toGoogleRunParameter: String = gcs
}
