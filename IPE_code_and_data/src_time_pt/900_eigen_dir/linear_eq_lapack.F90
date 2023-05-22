subroutine lapack_linear_eq_typ1(N,NRHS,A,b,x)
implicit none
integer,parameter :: wp = 8
integer,intent(in) :: N
integer,intent(in) :: NRHS
real(wp),intent(in) :: A(N,N)
real(wp),intent(in) :: b(N,NRHS)
real(wp),intent(out) :: x(N,NRHS)


!.. Local Scalars ..
real(wp) :: RCOND
integer ::           I, IFAIL, INFO, J 
character ::        EQUED
!.. Local Arrays ..
real(wp) ::  AF(N,N),  &
                 BERR(NRHS), C(N), FERR(NRHS), R(N), &
                 WORK(4*N)
integer ::       IPIV(N), IWORK(N)

integer :: LDA,LDAF,LDB,LDX,NRHSMX

   LDA=N
   LDAF=N
   LDB=N
   LDX=N
   NRHSMX=NRHS
         CALL DGESVX('Equilibration','No transpose',N,NRHS,A,LDA,AF, &
                    LDAF,IPIV,EQUED,R,C,B,LDB,X,LDX,RCOND,FERR,BERR, &
                    WORK,IWORK,INFO)

end subroutine lapack_linear_eq_typ1

