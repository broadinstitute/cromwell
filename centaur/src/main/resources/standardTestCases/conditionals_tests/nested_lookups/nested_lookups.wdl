workflow nested_lookups {
  Int i = 27
  Int a = 82
  if(true) {
    call mirror as throwaway1 { input: i = i } # Make sure this 'i' OGIN doesn't duplicate the nested m1's 'i' OGIN
    call mirror as throwaway2 { input: i = i } # Make sure this 'i' OGIN doesn't duplicate throwaway1's, or the nested m1's 'i' OGIN
    if(true) {
      if(true) {
        call mirror as m1 { input: i = i }
        call mirror as needs_wf_input
        Int b = a

        Int f = 100000
      }
    }
    if(true) {
      if(false) {
        Int? f1 = f
      }
    }
  }

  Int c = select_first([b, i])

  if(true) {
    Int? throwaway3 = m1.out # Make sure this 'm1.out' OGIN doesn't duplicate the nested m2's 'm1.out' OGIN
    Int? throwaway4 = m1.out # Make sure this 'm1.out' OGIN doesn't duplicate throwaway3's, or the nested m2's 'm1.out' OGIN
    Int? throwaway5 = b # Make sure this 'b' OGIN doesn't duplicate the nested e's 'b' OGIN
    Int? throwaway6 = b # Make sure this 'b' OGIN doesn't duplicate throwaway5's, or the nested e's 'b' OGIN

    if(true) {
      if(true) {
        call mirror as m2 { input: i = select_first([m1.out, 5]) + 1 }
        Int d = c
        Int? e = b
        Int? f2 = f1
      }
    }
  }

  output {
    Int? m1_out = m1.out
    Int? m2_out = m2.out
    Int? needs_wf_input_out = needs_wf_input.out

    Int? b_out = b
    Int c_out = c
    Int? d_out = d
    Int? e_out = e

    Int? f1_out = f1
    Int? f2_out = f2
  }
}

task mirror {
  Int i

  command {
    echo ${i}
  }
  output {
    Int out = read_int(stdout())
  }
  runtime {
    docker: "ubuntu:latest"
  }
}
