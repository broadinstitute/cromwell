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

package main

import (
        "bufio"
        "fmt"
        "os"
        "os/exec"
)

// spawnoninput will spawn a new command every time new input is detected
// it operates on stdin. It is useful with grep --line-buffered to spawn
// on specific strings in the inbound stream. The actual line of input
// on stdin will be sent in as the last argument to the process, so
// "spawnoninput my_snazzy_program -vv -q -l -a" will result in the arguments
// -vv, -q, -l, -a, <single line from stdin>" being sent to the program
func main() {
        cmdName := os.Args[1]
        cmdArgs := os.Args[2:]

        reader := bufio.NewReaderSize(os.Stdin, 80)
        scanner := bufio.NewScanner(reader)
        for scanner.Scan() {
                execute(cmdName, cmdArgs, scanner.Text())
        }
        if err := scanner.Err(); err != nil {
                fmt.Fprintln(os.Stderr, "error reading stdin: ", err)
        }
}

func execute(cmdName string, cmdArgs []string, line string) {
        argsWithLine := append(cmdArgs, line)
        cmd := exec.Command(cmdName, argsWithLine...)
        cmdOut, err := cmd.Output()
        if err != nil {
                fmt.Fprintln(os.Stderr, "error running command: ", err)
        }
        fmt.Println(string(cmdOut))
}
