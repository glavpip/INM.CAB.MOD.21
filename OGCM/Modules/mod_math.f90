module math_tools
use constants
implicit none

#ifndef __INTEL_COMPILER

contains

function cosd(x)
    real(4):: x 
    real(4):: cosd

    cosd = cos((x/180.0)*pi)
end function

function sind(x)
    real(4):: x 
    real(4):: sind

    sind = sin((x/180.0)*pi)
end function

function dcosd(x)
    real(8):: x 
    real(8):: dcosd

    dcosd = dcos((x/180.0d0)*dpi)
end function

function dsind(x)
    real(8):: x 
    real(8):: dsind

    dsind = dsin((x/180.0d0)*dpi)
end function

function dasind(x)
    real(8):: x 
    real(8):: dasind

    dasind = dasin(x)/dpi*180.0d0
end function

function dacosd(x)
    real(8):: x 
    real(8):: dacosd

    dacosd = dacos(x)/dpi*180.0d0
end function

function dtand(x)
    real(8):: x 
    real(8):: dtand

    dtand = dtan((x/180.0d0)*dpi)
end function

#endif

end module math_tools
