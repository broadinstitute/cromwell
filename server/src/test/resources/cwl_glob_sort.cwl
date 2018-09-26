cwlVersion: v1.0
class: CommandLineTool
requirements:
  - class: InlineJavascriptRequirement
hints:
  DockerRequirement:
    dockerPull: "debian:stretch-slim"
inputs: []
baseCommand: [touch, z, y, x, w, c, b, a]
outputs:
  letters:
    type: string
    outputBinding:
      glob: '?'
      outputEval: |
        ${ return self.sort(function(a,b) { return a.location > b.location ? 1 : (a.location < b.location ? -1 : 0) }).map(f => f.basename).join(" ") }
