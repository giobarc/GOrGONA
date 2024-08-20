!*******************************************************************************
! Data Ultima Modifica: 3 dicembre 2001                                        *
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

! MODULO PER LA DEFINIZIONE DI GRANDEZZE E PARAMETRI PER LA FUNZIONE
!===================================================================
MODULE FUNCTION_PAR
USE PARAMETERS
USE BH_MOD
Implicit None

SAVE

Real(8) :: ecoh(2),p(3),q(3),a(3),qsi(3),arete(3)
Real(8) :: a5(3),a4(3),a3(3),x5(3),x4(3),x3(3)
Real(8) :: dist(3)
Real(8) :: cutoff_start(3),cutoff_end(3)
Real(8) :: rat(2)
Real(8) :: mass(2)
Real(8) :: cutz(2,2)
Real(8) :: Rij(nmax0*(nmax0-1)/2)
Real(8) :: e_meno_pr(nmax0*(nmax0-1)/2)
Real(8) :: e_meno_dueqr(nmax0*(nmax0-1)/2)
Real(8) :: fx(nmax0),fy(nmax0),fz(nmax0)
Real(8) :: dfx(nmax0),dfy(nmax0),dfz(nmax0)
Real(8) :: dispx(nvmax0,nmax0),dispy(nvmax0,nmax0),dispz(nvmax0,nmax0)
Real(8) :: stress_atom(maxpar,3),pres_atom(maxpar)
Real(8) :: emetal
Real(8) :: sigma_lj,rmet_lj,epsi_lj,epsi_lj_k
! potential energy
Real(8) :: ener,ecinet,etot,ener_old,esub_old,deltaen,boltz_fact,ener_debug,ener_att
Integer :: nvois(nmax0),ivois(nvmax0,nmax0)
Integer :: imet(2)
Integer :: znmax(2)
Integer :: itype(nmax0)
Real(8) :: energ(nmax0)

Real*8 :: ddd(2,3,3,3)
Real*8 :: first_layer(2), amgo(2)
Real*8:: nv1(maxpar)
Real*8 :: e_met_mgo
Logical :: is_nan

Character(2) :: cutoff

END MODULE FUNCTION_PAR
