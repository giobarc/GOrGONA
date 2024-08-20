MODULE HISTO
  
  Integer, Parameter :: nbarsmax=1000
  Integer :: nbars(2)
  Integer :: nisto
  Integer :: order_parameter(2)
  Real*8 :: parameter_value(2)
  Real*8 :: legamimisti
  Real*8 :: hist_step(2)
  Real*8 :: start_histo(2)
  Real*8 :: weight
  Real*8 :: perc,perc_old
  Real*8 :: histogram(nbarsmax+1,nbarsmax+1),perc_histo(nbarsmax+1,nbarsmax+1)
  Real*8 :: histo_sig555(nbarsmax+1),perc_sig555(nbarsmax+1)
  Real*8 :: histo_misti(nbarsmax+1),perc_misti(nbarsmax+1)
  Real*8 :: hi_l(2),hi_s(2)
  Real*8 :: nsig555
  Real*8 :: sig555,sig322,sig433,sig321,sig331,sig422,sig421,sig333,sig554,sig666
  Real*8 :: sig533,sig211,sig544,sig200,sig311,sig300,sig555_old
  Real*8 :: legamimisti_old,frac_contact
  Integer :: renew_histo

END MODULE HISTO
