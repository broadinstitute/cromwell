cwlVersion: v1.0
$graph:
- id: one
  class: CommandLineTool

  inputs:
    r:
      type: float

  outputs:
    rSquared:
      type: float
      outputBinding:
        glob: stdout.txt
        loadContents: true
        outputEval: $(parseFloat(self[0].contents))
    rCopy:
      type: float
      outputBinding:
        outputEval: $(inputs.r)

  arguments:
  - valueFrom: "echo $(inputs.r * inputs.r)"
    shellQuote: false

  stdout: stdout.txt

- id: two
  class: CommandLineTool

  inputs:
    rSquared: float
    pi: float

  outputs:
    area:
      outputBinding:
        glob: stdout.txt
        loadContents: true
        outputEval: $(parseInt(self[0].contents))
      type: int
    rSquaredCopy:
      outputBinding:
        outputEval: $(inputs.rSquared)
      type: float

  arguments:
  - valueFrom: "echo $(inputs.rSquared * inputs.pi)"
    shellQuote: false

  stdout: stdout.txt

# Throwing in some ExpressionTools for fun. Currently these not only don't cache but aren't even listed in the metadata
# except in the original workflow blob (see #3499). The current thinking is that these should appear in metadata with a
# call-like presentation but that they don't need to be eligible for caching since ExpressionTools are supposed to be
# very lightweight and not worth the expense of hashing and potential cache hit copying.
- id: three
  class: ExpressionTool

  inputs:
    rSquared: float
    pi: float

  outputs:
    area: int
    rSquaredCopy: float

  expression: |
    ${
    return {"area": parseInt(inputs.pi * inputs.rSquared),
            "rSquaredCopy": inputs.rSquared };
    }


- id: cwl-cache-within-workflow
  class: Workflow

  requirements:
  - class: DockerRequirement
    dockerPull: "ubuntu:latest"
  - class: ShellCommandRequirement

  inputs:
    radius: float
    pi: float

  outputs:
    area:
      type: int
      outputSource: re-bar/area
    area-expression:
      type: int
      outputSource: re-baz/area

  steps:
    foo:
      run: "#one"
      in:
        r: radius
        pi: pi
      out: [rSquared, rCopy]

    bar:
      run: "#two"
      in:
        rSquared: "foo/rSquared"
        pi: pi
      out: [area, rSquaredCopy]

    re-bar:
      run: "#two"
      in:
        rSquared: "bar/rSquaredCopy"
        pi: pi
      out: [area, rSquaredCopy]

    baz:
      run: "#three"
      in:
        rSquared: "foo/rSquared"
        pi: pi
      out: [area, rSquaredCopy]

    re-baz:
      run: "#three"
      in:
        rSquared: "baz/rSquaredCopy"
        pi: pi
      out: [area, rSquaredCopy]
