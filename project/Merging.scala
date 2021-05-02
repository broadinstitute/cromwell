import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtassembly.{MergeStrategy, PathList}

object Merging {
  //noinspection SameParameterValue
  // Based on https://stackoverflow.com/questions/24363363/how-can-a-duplicate-class-be-excluded-from-sbt-assembly#57759013
  private def excludeFromJar(jarName: String): sbtassembly.MergeStrategy = new sbtassembly.MergeStrategy {
    override def name: String = "excludeFromJar"

    override def apply(tempDir: File, path: String, files: Seq[File]): Either[String, Seq[(File, String)]] = {
      val filteredFiles = files flatMap { file =>
        val (source, _, _, isFromJar) = sbtassembly.AssemblyUtils.sourceOfFileForMerge(tempDir, file)
        if (isFromJar && source.getName != jarName) Option(file -> path) else None
      }
      Right(filteredFiles)
    }
  }

  val customMergeStrategy: Def.Initialize[String => MergeStrategy] = Def.setting {
    case PathList(ps@_*) if ps.last == "project.properties" =>
      // Merge/Filter project.properties files from Google jars that otherwise collide at merge time.
      MergeStrategy.filterDistinctLines
    case PathList(ps@_*) if ps.last == "logback.xml" =>
      MergeStrategy.first
    // AWS SDK v2 configuration files - can be discarded
    case PathList(ps@_*) if Set("codegen.config" , "service-2.json" , "waiters-2.json" , "customization.config" , "examples-1.json" , "paginators-1.json").contains(ps.last) =>
      MergeStrategy.discard
    case x@PathList("META-INF", path@_*) =>
      path map {
        _.toLowerCase
      } match {
        case "spring.tooling" :: _ =>
          MergeStrategy.discard
        case "io.netty.versions.properties" :: Nil =>
          MergeStrategy.first
        case "maven" :: "com.google.guava" :: _ =>
          MergeStrategy.first
        case _ =>
          val oldStrategy = (assembly / assemblyMergeStrategy).value
          oldStrategy(x)
      }
    case x@PathList("OSGI-INF", path@_*) =>
      path map {
        _.toLowerCase
      } match {
        case "l10n" :: "bundle.properties" :: Nil =>
          MergeStrategy.concat
        case _ =>
          val oldStrategy = (assembly / assemblyMergeStrategy).value
          oldStrategy(x)
      }
    case "asm-license.txt" | "module-info.class" | "overview.html" | "cobertura.properties" =>
      MergeStrategy.discard
    case PathList("mime.types") =>
      MergeStrategy.last
    case PathList("scala", "annotation", "nowarn.class" | "nowarn$.class") =>
      // scala-collection-compat shading bug: https://github.com/scala/scala-collection-compat/issues/426
      excludeFromJar("scala-collection-compat_2.12-2.4.2.jar")
    case x =>
      val oldStrategy = (assembly / assemblyMergeStrategy).value
      oldStrategy(x)
  }
}
