cwlVersion: v1.0
$namespaces:
  dx: https://www.dnanexus.com/cwl#
$graph:
- id: inputdir_zero_doesnt_localize
  class: Workflow
  hints:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: google/cloud-sdk:slim
  - class: dx:InputResourceRequirement
    indirMin: 0
  inputs: []
  outputs:
  - id: errors
    type: File
    outputSource: findFile/errors
  steps:
  - id: echo
    in: []
    out:
    - id: echoOut
    - id: echoArrayOut
    - id: echoRecordOut
    - id: echoMaybeOut
    run:
      class: CommandLineTool
      inputs: []
      outputs:
      - id: echoOut
        outputBinding:
          glob: lonely_file
        type: File
        secondaryFiles: [.also]
      - id: echoArrayOut
        outputBinding:
          glob: '*.txt*'
          outputEval: $([self[0]])  
        type:
          items: File
          type: array
      - id: echoRecordOut
        secondaryFiles: [.also]
        type:
          type: record
          name: foo
          fields:
            - name: a
              type: File
              outputBinding:
                glob: file_in_object
      - id: echoMaybeOut
        type: 
        - File
        - 'null'
        secondaryFiles: [.also]
        outputBinding:
          glob: maybe_file
      arguments:
       - valueFrom: |
           echo "lonely_file" > lonely_file
           echo "lonely_file.also" > lonely_file.also
           echo "file_in_array.txt" > file_in_array.txt
           echo "file_in_array.txt.also" > file_in_array.txt.also
           echo "file_in_object" > file_in_object
           echo "file_in_object.also" > file_in_object.also
           echo "maybe_file" > maybe_file
           echo "maybe_file.also" > maybe_file.also
         shellQuote: false
  - id: findFile
    in:
    - id: f
      source: "#inputdir_zero_doesnt_localize/echo/echoOut"
    - id: g
      source: "#inputdir_zero_doesnt_localize/echo/echoArrayOut"
    - id: h
      source: "#inputdir_zero_doesnt_localize/echo/echoRecordOut"
    - id: i
      source: "#inputdir_zero_doesnt_localize/echo/echoMaybeOut"
    out:
    - id: errors
    run:
      inputs:
      - id: f
        type: File
        secondaryFiles: [.also]
      - id: i
        type: File?
        secondaryFiles: [.also]
      - id: g
        type:
          type: array
          items: File
        secondaryFiles: [.also]
      - id: h
        secondaryFiles: [.also]
        type:
          type: record
          name: foo
          fields:
            - name: a
              type: File
      outputs:
      - id: errors
        type: string
        outputBinding:
          glob: errors.txt
          loadContents: true
          outputEval: $(self[0].contents)
      class: CommandLineTool
      requirements:
      - class: ShellCommandRequirement
      - class: InlineJavascriptRequirement
      arguments:
      # Unfortunately we can't use the same trick as in draft3_nio_file_papi1 because in papi 2 the instance metadata doesn't seem to contain
      # the PAPI operation ID. So instead simply check that the input files are not there. This is not as neat as the PAPI1 version
      # but should check the right thing as long as the files are localized (in general, but not here) somewhere in the same directory as where the command
      # runs. If that changes then the directory where "find" searches should be updated 
        - valueFrom: |
            touch errors.txt
            find . -name lonely_file >> errors.txt
            find . -name lonely_file.also >> errors.txt
            find . -name file_in_array >> errors.txt
            find . -name file_in_array.also >> errors.txt
            find . -name file_in_object >> errors.txt
            find . -name file_in_object.also >> errors.txt
            find . -name maybe_file >> errors.txt
            find . -name maybe_file.also >> errors.txt
          shellQuote: false
