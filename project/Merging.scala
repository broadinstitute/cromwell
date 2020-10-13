import sbt._
import sbtassembly.AssemblyPlugin.autoImport._
import sbtassembly.{MergeStrategy, PathList}

object Merging {
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
        case "spring.tooling" :: xs =>
          MergeStrategy.discard
        case "io.netty.versions.properties" :: Nil =>
          MergeStrategy.first
        case "maven" :: "com.google.guava" :: xs =>
          MergeStrategy.first
        case _ =>
          val oldStrategy = (assemblyMergeStrategy in assembly).value
          oldStrategy(x)
      }
    case x@PathList("OSGI-INF", path@_*) =>
      path map {
        _.toLowerCase
      } match {
        case "l10n" :: "bundle.properties" :: Nil =>
          MergeStrategy.concat
        case _ =>
          val oldStrategy = (assemblyMergeStrategy in assembly).value
          oldStrategy(x)
      }
    case "asm-license.txt" | "module-info.class" | "overview.html" | "cobertura.properties" =>
      MergeStrategy.discard
    case PathList("mime.types") =>
      MergeStrategy.last
    case x =>
      val oldStrategy = (assemblyMergeStrategy in assembly).value
      oldStrategy(x)
  }
}
