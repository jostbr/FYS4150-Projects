
-------------------------------------------------------------------------
In order to obtain the files
- result10.txt
- results100.txt
- result1000.txt

run tridiag.cpp with command line arguments 10 100 1000.
-------------------------------------------------------------------------

-------------------------------------------------------------------------
In order to obtain timing_general.txt and timing_special.txt, make calls 
time_algorithm("general", 23, 5) and time_algorithm("special", 23, 5) while having
changed the base number (for N), within \texttt{time_algorithm()} from 10 to 2).
-------------------------------------------------------------------------

-------------------------------------------------------------------------
In order to obtain runtime_summary.txt, call time_algorithm("general", 3, 100),
time_algorithm("special", 3, 100) (while using base 10 for N in time_algorithm())
and lu_solver().
-------------------------------------------------------------------------
