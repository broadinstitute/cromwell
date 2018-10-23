cwlVersion: "v1.0"
class: "CommandLineTool"
label: "Example trivial wrapper for Java 7 compiler"
hints: 
  - dockerPull: "java:7-jdk"
    class: "DockerRequirement"
baseCommand: "javac"
arguments: 
  - "-d"
  - "$(runtime.outdir)"
inputs: 
  - type: "File"
    inputBinding: 
      position: 1
    id: "file:///home/dan/wdl4s/arguments.cwl#src"
outputs: 
  - type: "File"
    outputBinding: 
      glob: "*.class"
    id: "file:///home/dan/wdl4s/arguments.cwl#classfile"
id: "file:///home/dan/wdl4s/arguments.cwl"
name: "file:///home/dan/wdl4s/arguments.cwl"
