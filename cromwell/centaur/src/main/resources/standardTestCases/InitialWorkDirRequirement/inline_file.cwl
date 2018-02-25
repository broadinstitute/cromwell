# Little tool that uses a prime sieve to generate primes up to 100

cwlVersion: v1.0
$graph:
- id: iwdr_inline_file
  class: CommandLineTool
  baseCommand: ['python3', 'prime_sieve.py', '100']
  requirements:
      - class: DockerRequirement
        dockerPull: "python:3.5.0"
      - class: InitialWorkDirRequirement
        listing:
            - entryname: 'prime_sieve.py'
              entry: |
                  import sys
                  import math

                  limit = int(sys.argv[1])
                  sieve = [True for i in range(limit)]
                  for i in range(2, math.floor(limit / 2)):
                      if sieve[i]:
                          for j in range(i * 2, limit, i):
                              sieve[j] = False

                  result = "["
                  for i in range(2, limit):
                      if sieve[i]:
                          if result != "[":
                              result += ", "
                          result += str(i)
                  result += "]"

                  print(result)

  stdout: "primes"
  inputs: []
  outputs:
      - id: prime_list
        type: string
        outputBinding:
          glob: primes
          loadContents: true
          outputEval: $(self[0].contents.trim())
