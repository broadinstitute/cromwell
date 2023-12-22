/*
 * Copyright 2018 Amazon.com, Inc. or its affiliates.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the copyright holder nor the names of its
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
 *  BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 *  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *  THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 *  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 *  IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */
package cromwell.backend.impl.aws.io


import wom.values._
import cromwell.backend.io._
import cromwell.backend.standard._
import scala.concurrent.Future
import java.nio.file.Paths
import cromwell.core.path.DefaultPathBuilder


trait AwsBatchGlobFunctions extends GlobFunctions {
    
  def standardParams: StandardExpressionFunctionsParams

  

  /**
    * Returns a list of path from the glob.
    *
    * The paths are read from a list file based on the pattern.
    *
    * @param pattern The pattern of the glob. This is the same "glob" passed to globPath().
    * @return The paths that match the pattern.
    */
  override def glob(pattern: String): Future[Seq[String]] = {
    // get access to globName()
    import GlobFunctions._
    
    // GOAL : 
    //  - get config (backend / runtime / ...) here to evaluate if efsMntPoint is set & if efs delocalization is set. 
    //  - according to those values : write the pattern as s3:// or as local path. 
    //  - get the wf id from the config settings.

    // this function reads in the globfile and locates globbed files : "local" or NIO access is needed to the files. 

    // for now : hard coded as local at mount point /mnt/efs.
    val wfid_regex = ".{8}-.{4}-.{4}-.{4}-.{12}".r
    val wfid = callContext.root.toString.split("/").toList.filter(element => wfid_regex.pattern.matcher(element).matches()).lastOption.getOrElse("")
    val globPatternName = globName(s"${pattern}-${wfid}")
    val globbedDir = Paths.get(pattern).getParent.toString
    val listFilePath = if (pattern.startsWith("/mnt/efs/")) {
        DefaultPathBuilder.get(globbedDir + "/." + globPatternName + ".list")
    } else {
        callContext.root.resolve(s"${globbedDir}/.${globPatternName}.list".stripPrefix("/"))
    }
    asyncIo.readLinesAsync(listFilePath.toRealPath()) map { lines =>
      lines.toList map { fileName =>
        // again : this should be config based...
        if (pattern.startsWith("/mnt/efs/")) {
            s"${globbedDir}/.${globPatternName}/${fileName}"
        } else {
            callContext.root.resolve(s"${globbedDir}/.${globPatternName}/${fileName}".stripPrefix("/")).pathAsString
        }
      }
    }
  }
  
}