task t {
 String s
 command {...}
 output {
   String o = "yahh"
   Int o2 = 5
 }
}

workflow w2 {
 String inp
 call t { input: s = inp }
 output {
   String result = t.o
 }
}