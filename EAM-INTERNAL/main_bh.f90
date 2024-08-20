PROGRAM BH

  USE BH_MOD
  USE FUNCTION_PAR

  implicit none
  Integer :: n_iter
  Integer :: n_err = 0

  Call read_bh
  Call drive_bh(n_err,n_iter)


END PROGRAM BH
