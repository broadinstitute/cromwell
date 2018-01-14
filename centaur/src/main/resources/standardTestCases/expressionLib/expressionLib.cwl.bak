class: CommandLineTool
cwlVersion: v1.0
requirements:
  InlineJavascriptRequirement:
    expressionLib:
      - "function foo() { return 2; }"
inputs: []
outputs:
  out: stdout
arguments: [echo, $(foo())]
stdout: whatever.txt
