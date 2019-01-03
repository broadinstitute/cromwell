cwlVersion: v1.0
$graph:
- id: cwl_input_json
  class: Workflow
  inputs: []
  outputs:
    final_output:
      type: File
      outputSource: round/output_file
  steps:
    make:
      run: "#makefile"
      in: []
      out: [fileoutput]
    round:
      run: "#roundtrip"
      in:
        input_record: 
          source: "make/fileoutput"
      out: [output_file]
          
- id: makefile
  class: CommandLineTool
  requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: DockerRequirement
    dockerPull: "ubuntu:latest"
  inputs: []
  outputs:
    fileoutput:
      type:
        fields:
        - name: input_file
          type: File
        name: input_record
        type: record
  arguments:
    - valueFrom: >
         echo foo > foo && echo '{ "fileoutput": { "input_file": {"path": "$(runtime.outdir)/foo", "class": "File"} } }' > cwl.output.json
      shellQuote: false

- id: roundtrip
  class: CommandLineTool
  hints:
  - class: DockerRequirement
    dockerPull: "stedolan/jq:latest"
  inputs:
  - id: input_record
    type:
      fields:
      - name: input_file
        type: File
      name: input_record
      type: record
  outputs:
  - id: output_file
    type: File
  requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
  - class: InitialWorkDirRequirement
    listing:
    - entry: $(JSON.stringify(inputs))
      entryname: cwl.inputs.json
  arguments:
  # Round-trips the file referenced in cwl.input.json to cwl.output.json. Also ls it in the command to make sure it's there.
    - valueFrom: >
         INPUT_FILE=\$(cat cwl.inputs.json | jq -r '.. | .path? // empty') && ls $INPUT_FILE && echo "{\"output_file\": {\"path\": \"\$INPUT_FILE\", \"class\": \"File\"} }" > cwl.output.json
      shellQuote: false
