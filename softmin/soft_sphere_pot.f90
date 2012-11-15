! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!                            
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

subroutine soft_sphere_pot(dimen,npart,x,diams,energy,force) ! Subroutine that calculates the force and potential.

  implicit none

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! npart = num particles, x = coords, diams = particle diameters, !
  ! dimen = dimension (2 or 3), force = particle force array,      !
  ! energy = potential energy.                                     !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer, intent(in) :: npart,dimen
  double precision, intent(in) :: x(dimen*npart)
  double precision, intent(out) :: force(dimen*npart), energy
  double precision, intent(in) :: diams(npart)
  double precision :: rij,Odij,Orij,xij,yij,zij,fr,ptemp
  integer :: i,j,i2,i2m1,j2,j2m1,i3,i3m2,i3m1,j3,j3m2,j3m1
  
  force = 0.0D0

  energy = 0.0D0

  if (dimen == 2) then

    do i=1,npart-1
      
      i2 = 2*i
      i2m1 = i2-1
      
      do j=i+1,npart

        j2 = 2*j
        j2m1 = j2-1

        xij = x(i2m1)-x(j2m1)
        xij = xij-dnint(xij) ! PBC
        yij = x(i2)-x(j2)
        yij = yij-dnint(yij) ! PBC
        rij = dsqrt(xij*xij+yij*yij)
        Orij = 1.0D0/rij
        Odij = 2.0D0/(diams(i)+diams(j))
        
        if (Orij > Odij) then
          
          ptemp = (1.0D0-rij*Odij)
          energy = energy + ptemp*ptemp
          fr=Odij*Orij*ptemp
          force(i2m1) = force(i2m1) - fr*xij
          force(j2m1) = force(j2m1) + fr*xij
          force(i2) = force(i2) - fr*yij
          force(j2) = force(j2) + fr*yij
          
        end if
      end do
    end do
    
  else if (dimen == 3) then

    do i=1,npart-1
      
      i3 = 3*i
      i3m1 = i3-1
      i3m2 = i3-2
      
      do j=i+1,npart

        j3 = 3*j
        j3m1 = j3-1
        j3m2 = j3-2

        xij = x(i3m2)-x(j3m2)
        xij = xij-dnint(xij) ! PBC
        yij = x(i3m1)-x(j3m1)
        yij = yij-dnint(yij) ! PBC
        zij = x(i3)-x(j3)
        zij = zij-dnint(zij) ! PBC
        rij = dsqrt(xij*xij+yij*yij+zij*zij)
        Orij = 1.0D0/rij
        Odij = 2.0D0/(diams(i)+diams(j))
        
        if (Orij > Odij) then
          
          ptemp = (1.0D0-rij*Odij)
          energy = energy + ptemp*ptemp
          fr=Odij*Orij*ptemp
          force(i3m2) = force(i3m2) - fr*xij
          force(j3m2) = force(j3m2) + fr*xij
          force(i3m1) = force(i3m1) - fr*yij
          force(j3m1) = force(j3m1) + fr*yij
          force(i3) = force(i3) - fr*zij
          force(j3) = force(j3) + fr*zij
          
        end if
      end do
    end do
    
  else
    
    print *,'Wrong dimension in soft_sphere_pot subroutine.'
    stop
    
  end if
  
  energy = energy/2.0D0
  
end subroutine soft_sphere_pot
