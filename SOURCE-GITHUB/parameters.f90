!*******************************************************************************
! Data Ultima Modifica: 12 novembre 2001                                       *
!                                                                              *
! Autore: ARNALDO RAPALLO                                                      *
!                                                                              *
! Consiglio Nazionale delle Ricerche (CNR):                                    *
! Istituto di Studi Chimico Fisici di Macromolecole Sintetiche e Naturali(IMAG)*
! Via de Marini 6                                                              *
! CAP 16149 Genova (GE)                                                        *
! tel.     +390106475868                                                       *
! Fax.     +390106475880                                                       *
! e-mail   rapallo@imag.ge.cnr.it                                              *
!*******************************************************************************

MODULE PARAMETERS
Implicit None

Real(8), PARAMETER :: pigreco=dacos(-1.d0)

! numero massimo di atomi
Integer, Parameter :: maxpar=5600

! numero che viene considerato il massimo Real(8)
Real(8), Parameter :: maxreal=1.E+300_8

! numero che viene considerato il massimo Integer
Integer, Parameter :: maxint=2000000000

! massimo consentito per i rifiuti di mossa
Integer, Parameter :: maxfail=100

! parametro di tolleranza sulla determinazione di un minimo locale
Real(8), Parameter :: ftoll=1.e-7_8

! passo temporale per il quenching
Real(8), Parameter :: tstep=2.e-15_8

! attrito per la parte iniziale del quenching
Real(8), Parameter :: eta_andersen=1.e15_8

Real(8), PARAMETER :: pi=3.141592653589793
Real(8), PARAMETER :: cbol=8.62d-05
Real(8), PARAMETER :: evsujoule=1.60219d-19
Real(8), PARAMETER :: uasukg=1./6.023d26
Real(8), PARAMETER :: angsum=1.d-10
Integer, Parameter :: nmax0=maxpar,nvmax0=maxpar
Real(8), Parameter :: rac2=1.41421356237309504880_8
Real(8), Parameter :: rac3=1.73205080756887729352_8
Real(8), Parameter :: rac5=2.23606797749978969640_8
Real(8), Parameter :: rac8=2._8*rac2
Real(8), Parameter :: epsi=1.e-6_8

CONTAINS
Real*8 FUNCTION rndd()
      	REAL*8 zz
	call RANMAR(zz)
	rndd = zz
END FUNCTION rndd

END MODULE PARAMETERS
