module atm_pars
implicit none

include 'atmpars.fi'

real(8),allocatable::  xa(:), ya(:)  !atmospheric lon-grid and lat-grid

integer, allocatable:: atm_mask(:,:)

endmodule atm_pars