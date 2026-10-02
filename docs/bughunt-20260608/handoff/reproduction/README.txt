GPU-free parser check (proves Bug #4). Build & run:
  g++ -O2 -std=c++17 parse_repro2.cpp -o parse_repro2
  ./parse_repro2 ../test_case_NaX/force_field.def
Expect: NEW(buggy) Nmixrule=0 ; NEW+FIX Nmixrule=11
