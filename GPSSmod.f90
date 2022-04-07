module GPSS
 implicit none
 integer, parameter:: max_asm    =  300
 integer, parameter:: max_bin    =   50
 integer, parameter:: max_drn    =   30
 integer, parameter:: max_fac    =   50
 integer, parameter:: max_fam    =  300
 integer, parameter:: max_ga1    =  300
 integer, parameter:: max_ga2    =   10
 integer, parameter:: max_mfac   =   20
 integer, parameter:: max_se     =   20
 integer, parameter:: max_pla    =   20
 integer, parameter:: max_polvec =  300
 integer, parameter:: max_pol    =   20
 integer, parameter:: max_sm     = 1024
 integer, parameter:: max_sbv    =   40
 integer, parameter:: max_src    =   20
 integer, parameter:: max_sto    =   40
 integer, parameter:: max_str    =   40
 integer, parameter:: max_uchf1  =  300
 integer, parameter:: max_uchf2  =   10
 integer, parameter:: max_ucht   =   10
 integer, parameter:: max_tx     =  300    ! к-во tx
 integer, parameter:: max_tx2    =   30    ! вторая размерность TX
 integer, parameter:: max_el     =  300
 integer, parameter:: max_state  = 9200
 integer, parameter:: max_tab    =  100

 integer  asm(max_asm,10)
 integer  bin(max_bin,8)
 real     binsta(max_bin,2)
 real*8   drn(max_drn),DFACT(max_drn),DMODUL,DCONST(max_drn)
 integer  fac(max_fac,3)
 integer  fam(max_fam,3)
 integer  gathf(max_ga1,10)
 integer  gatht(max_ga2)
 integer  LSE,MFAC(max_mfac,2),MBV(max_mfac),SE(max_se,3)
 integer  PLAMA(max_pla,2)
 integer  POLVEC(max_polvec),POLC,POL(max_pol,2)
 integer  LSM , SBV(max_sbv) , SM(max_sm,2)
 integer  SRC(max_src,2) , NTXC
 integer  sto(max_sto,2)
 integer  STRAMA(max_str,2)
 integer  UCHF(max_uchf1,max_uchf2,2)
 integer  UCHT(max_ucht,2)
 !un
 integer tx(max_tx,30),al(max_tx,2),lal,ltx,el(max_el,2)
 integer lel,lev,lfam,naddr,state(max_state),ok
 integer n,t,it,rt





 ! par
 integer el1,tx1,tx2,fam1,src1,asm1,fac1,mfac1,se1,sto1
 integer sm1,ngate1,ngate2,gatht1,gathf1,ucht1,uchf1,kend,STATE1
 integer efac,emfac,esto,egate1,egate2,egathf,egatht,euchf,bmfac
 integer bsto,bgate2,pol1,bin1,tab1,drn1,ind,outd,savi,savo
 ! end par

! parg
!      common /par/ el1,tx1,tx2,fam1,src1,asm1,fac1,mfac1,se1,sto1,
!     +sm1,Agate1,Agate2,gatht1,gathf1,ucht1,uchf1,kend,state1,
!     +efac,emfac,esto,egate1,egate2,egathf,egatht,euchf,bmfac,
!     +bsto,bgate2,pol1,bin1,tab1,drn1,ind,outd,savi,savo

contains
      subroutine activ1(*)
!     ***
!     ***          call activ1(&1006)
!     ***
!     *** function :  search out the next event or the next
!     ***             scheduled tx activation
!     *** parameter:  exit1 = exit to final analysis section
!     ***
      implicit integer (a-z)
!
!     clear the line pointers
!     =======================
      lev = 0
      ltx = 0
!
!     set the simulator clock
!     =======================
      t1 = t - rt
      t  = n + 1
!
!     search out the next event
!     =========================
      do 100 j = 1 , lel
      if(el(j,2).eq.0.or.el(j,2).ge.t) goto 100
      t = el(j,2)
      lev = j
      naddr = el(lev,1)
100   continue
!
!     search out the next tx activation
!     =================================
      do 200 j = 1 , lal
      if(al(j,2).le.0.or.al(j,2).ge.t) goto 200
      t = al(j,2)
      ltx = j
      naddr = al(ltx,1)
      lev = 0
200   continue
!
!     check for simulation halt
!     =========================
      if(t.gt.n) return 1
!
!     set the user clock
!     ==================
      rt = t - t1
!
!     family membership
!     =================
      lfam = 0
      if(ltx.eq.0) goto 400
      if(tx(ltx,2).eq.0) return
      do 300 j = 1 , fam1
      if(tx(ltx,1).eq.fam(j,1)) goto 350
300   continue
      return
350   lfam = j
      return
!
!     clear the event's entry
!     =======================
400   el(lev,1) = 0
      el(lev,2) = 0
!
!     reset the list-end pointer, lel
!     ===============================
450   if(el(lel,1).ne.0.or.lel.eq.1) return
      lel = lel-1
      goto 450
      end subroutine activ1


      subroutine activ2(*,*)
!     ***
!     ***          call activ2(&1001,&1003)
!     ***
!     *** function  : search for conditioned activations
!     *** parameters: exit1 = sheduled activation exit
!     ***             exit2 = conditioned activation exit
!     ***
      implicit integer ( a - z )
!
!     clear the event pointer
!     =======================
      lev = 0
!
!     let txs try the type-2 gates
!     ============================
      if(it.eq.t) goto 200
      it = t
      do 100 i = bgate2 , egate2
100   state(i) = 1
!
!     search for a blocked tx
!     =======================
200   do 201 ja = 1 , lal
      if(al(ja,2).ge.0) goto 201
      k = - al(ja,2)
      if(k.gt.kend) goto 201
      if(state(k).NE.0) goto 250
201   continue
      it = 0
      return 1
!
!     assemble
!     ========
250   polc = 1
      polvec(polc)=ja
      ja = ja + 1
      if(ja.gt.lal) goto 350
      do 300 jb = ja , lal
      if(al(jb,2).NE.-k) goto 300
      polc = polc + 1
      polvec(polc) = jb
300   continue
!
!     policy selection
!     ================
350   call policy(k)
!
!     conclusion
!     ==========
      if(k.le.egate2) ok = 1
      if(k.lt.bgate2.or.k.gt.egate2) it = 0
!
!     activate the tx
!     ===============
      naddr = al(ltx,1)
      lfam=0
      if(tx(ltx,2).eq.0) return 2
      do 400 j = 1 , fam1
      if(tx(ltx,1).eq.fam(j,1)) goto 450
400   continue
      return 2
450   lfam = j
      return 2
      end subroutine activ2



      subroutine advanc(AT,idn,*,iprint)
!     ***
!     ***
!     ***          call advanc(at,idn,&1005,iprint)
!     ***
!     *** function  : detain a tx
!     *** parameters: at  = detention time
!     ***             idn = target
!     ***
      implicit integer (a-z)
!
!     detain
!     ======
      al(ltx,1)=idn
      al(ltx,2)=t+at
      if(iprint.eq.0) return 1
      write(outd,3000) t,tx(ltx,1),tx(ltx,2),al(ltx,2)
3000  format(8h ADVANC:,3x,2hT=,i7,2x,2hTX,i5,1h,,i3,2x,&
      26h WILL BE DETAINED UNTIL T=,i7)
      return 1
      end subroutine advanc


      subroutine alloc(nst,ne,mark,id,line,lock,*,*,iprint)
!     ***
!     ***          call alloc(nst,ne,mark,id,line,lock,&1005,
!     ***                     &1006,iprint)
!     ***
!     *** function   :  acquisition of space in an addressible st
!     *** parameters :  nst  = storage number
!     ***               ne   = number of elements requested
!     ***               mark = segment code
!     ***               id   = alloc call's statement number
!     ***               line = storage address
!     ***               lock = lock flag
!     ***                      = 0: newly arriving txs try to
!     ***                           acquire space
!     ***                      = 1: newly arriving txs are locked
!     ***                           immediately
!     ***
      implicit integer (a - z)
!
!     determine the station number
!     ============================
      k = emfac + nst
!
!     watchdog
!     ========
      if(ok.eq.0) goto 50
      ok = 0
      if(ne.le.0) return
!
!     storage-space assignment by the strategy
!     ========================================
      call strata(nst,ne,*70)
      if(LSM.eq.0) goto 60
!
!     set the segment matrix
!     ======================
      if(sm(lsm,1).ne.0) goto 30
      do 100 i = 1,sm1
      j = lsm - i
      if(sm(j,1).ne.0) goto 20
100   continue
20    sm(lsm,1) = sm(j,1) - i
      sm(j,1) = i
30    if(sm(lsm,1).eq.ne) goto 40
      j = lsm + ne
      sm(j,1) = sm(lsm,1) - ne
      sm(j,2) = - 1
40    sm(lsm,1) = ne
      sm(lsm,2) = mark
      sto(nst,1)= sto(nst,1) + ne
      tx(ltx,8) = 0
      line = lsm - sbv(nst) + 1
      if(iprint.eq.0) return
      write(outd,3000) t,tx(ltx,1),tx(ltx,2),ne,nst,line
3000  format(8h ALLOC :,3x,2hT=,I7,2x,2hTX,i5,1h,,i3,2x,&
      9h ACQUIRES,i5,16h ELEMENTS IN STO,i3,12h STARTING AT,&
      8h ADDRESS,i5)
      return
!
!     first acquisition attempt
!     =========================
50    if(ne.le.0) return
      if(iprint.eq.0) goto 5001
      write(outd,3001) t,tx(ltx,1),tx(ltx,2),ne,nst
3001  format(8h ALLOC :,3x,2hT=,i7,2x,2hTX,i5,1h,,i3,2x,&
      9h REQUIRES,i5,16h ELEMENTS IN STO,i3)
5001  if(lock.gt.0) goto 60
      al(ltx,1) = id
      al(ltx,2) = -k
      tx(ltx,8) = t
      return 1
!
!     lock
!     ====

60    al(ltx,1) = ID
      al(ltx,2) = - kend - k
      if(tx(ltx,8).eq.0) tx(ltx,8) = t
      if(iprint.eq.0) return 1
      write(outd,3002) t,tx(ltx,1),tx(ltx,2),nst
3002  format(8h ALLOC :,3x,2hT=,i7,2x,2hTX,i5,1h,,i3,2x,&
      17h IS LOCKED AT STO,i3)
      return 1
!
!     storage error
!     =============
70    return 2
      end subroutine alloc

      subroutine arrive(nbn,ne,*,iprint)
!     ***
!     ***          call arrive(nbn, ne, &1006,iprint)
!     ***
!     *** function   :   deposit tokens in a bin
!     *** parameters :   nbn = bin number
!     ***                ne  = number of arriving tokens
!     ***
      implicit integer (a - z)
!
!     error control
!     =============
      do 100 i = 9 , 13
      if(tx(ltx,i).eq.0) goto 200
100   continue
      write(outd,3000) t,tx(ltx,1),tx(ltx,2),nbn
3000  format(1h0,24(1h+),26h ERROR: SUBROUTINE ARRIVE ,30(1h+)/1&
      24(1h+),3h T=,i7,2x,2hTX,i5,1h,,i3,13h ALREADY HAS ,&
      20(1h+),/1x,24(1h+),30h TOKENS IN 5 BINS, DEPOSIT IN ,&
      22(1h+)/1x,24(1h+),4h BIN,i3,14h NOT POSSIBLE ,35(1h+)/)
      return 1
!
!     set the tx matrix
!     =================
200   tx(ltx,i) = nbn
      tx(ltx,i+5) = t
!
!     set the bin matrix
!     ==================
      bin(nbn,7) = bin(nbn,7) + bin(nbn,1) * (t - bin(nbn,8))
      bin(nbn,8) = t
      bin(nbn,1) = bin(nbn,1) + ne
      if(bin(nbn,2).lt.bin(nbn,1)) bin(nbn,2) = bin(nbn,1)
      bin(nbn,3) = bin(nbn,3) + ne
      if(iprint.eq.0)  return
      write(outd,3001) t,tx(ltx,1),tx(ltx,2),ne,nbn
3001  format(8h ARRIVE:,3x,2hT=,i7,2x,2hTX,i5,1h,,i3,2x,&
      9h DEPOSITS,i5,17h TOKENS IN BIN NO,i3)
      return
      end subroutine arrive




      subroutine assemb(numass,nassem,*,iprint)
!     ***
!     ***          call assemb(numass,nassem,&1005,iprint)
!     ***
!     *** function  : reunite family members
!     *** parameters: numass = assembly station number
!     ***             nassem = number of txs to be assembled
!     ***
      implicit integer (a - z)
      if(iprint.eq.0) goto 5000
      write(outd,3000) t,tx(ltx,1),tx(ltx,2),numass
3000  format(8h ASSEMB:,3x,2hT=,i7,2x,2hTX,i5,1h,,i3,2x,&
      32h HAS REACHED ASSEMBLY STATION NO,i3)
!
!     test for family membership
!     ==========================
5000  if(lfam.eq.0) return
!
!     initiate an assembly
!     ====================
      if(asm(lfam,numass).eq.0) asm(lfam,numass)=nassem
!
!     decrement the count
!     ===================
      asm(lfam,numass) = asm(lfam,numass) - 1
!
!     fuse
!     ====
      if(asm(lfam,numass).gt.0) goto 200
      if(fam(lfam,2).gt.1) goto 150
      if(iprint.eq.0) return
      write(outd,3001)
3001  format(36x,44hAND LEAVES STATION AS THE LAST SURVIVING   ,&
      8h KINSMAN)
      return
150   continue
      if(iprint.eq.0) return
      write(outd,3002)
3002  format(36x,36hAND LEAVES THE STATION AFTER REUNION)
      return
!
!     annihilate
!     ==========
200   call termin(*250,iprint)
250   return 1
      end subroutine assemb



      subroutine bfit(nst , ne)
!     ***
!     ***          call bfit(nst,ne)
!     ***
!     *** function  : find a free space using the best-fit
!     ***             strategy
!     *** parameters: nst = storage number
!     ***             ne  = number of elements to be acquired
!     ***
      implicit integer (a-z)
!
!     initialization the search
!     =========================
      lsm = 0
      d = sto(nst,2)
      i = sbv(nst)
      ie = i + sto(nst,2) - 1
!
!     find a segment
!     ==============
10    if(sm(i,2).ne.-1) goto 20
      if(sm(i,1).lt.ne) goto 20
      if(sm(i,1)-ne.ge.d) goto 20
      d = sm(i,1) - ne
      lsm = i
20    i = i +sm(i,1)
      if(i.le.ie) goto 10
      return
      end subroutine bfit


      subroutine boxexp(min,max1,max2,ratio,rnum,random)
!     ***
!     ***          call boxexp(mi,max1,max2,ratio,rnum,random)
!     ***
!     *** function  : generate a random sequence uniform on one
!     ***             interval and exponential on an adjoining
!     ***             interval to the right of the first
!     *** parameters: min    = lower bound for the uniform
!     ***                      distribution
!     ***             max1   = upper bound for the uniform
!     ***                      distribution
!     ***             max2   = upper bound for the exponential
!     ***                      distribution
!     ***             ratio  = ration of terms to be drawn from
!     ***                      the interval with uniform
!     ***                      distribution
!     ***             rnum   = generator's identifier
!     ***             random = resulting random number
!     ***
      implicit none 
      integer rnum
      real min ,max1 , max2 , ratio, random
	  real z
!
!     generate a random number in the uniform interval
!     ================================================
      call unifrm(0.,1.,rnum,z)
      if(z.ge.ratio) goto 100
      random = min + ((max1 - min) / ratio) * z
      return
!
!     generate random number in the exponential interval
!     ==================================================
100   random = max1 - (((1. - ratio) * (max1 - min)) / ratio)&
      *alog((1. - z) / (1. - ratio))
      if(random.le.max2) return
      call unifrm(ratio,1.,rnum,z)
      goto 100
      end subroutine boxexp



      subroutine buffer(idn,*,iprint)
!     ***
!     ***          call buffer(idn,&1005,iprint)
!     ***
!     *** function  : interrupt a txn activation
!     *** parameter : idn=target
!     ***
      implicit integer (a-z)
!
!     deactivate
!     ==========
      al(ltx,1)=idn
      al(ltx,2)=t
      if(iprint.eq.0) return 1
      write(outd,3000) t,tx(ltx,1),tx(ltx,2)
3000  format(8h BUFFER:,3x,2hT=,i7,2x,2hTX,i5,1h,,i3,2x,&
      15h IS DEACTIVATED)
      return 1
      end subroutine buffer


      subroutine clear(nfa,*,*,iprint)
!     ***
!     ***          call clear(nfa,exit1,&1006,iprint)
!     ***
!     *** function  : free a facility
!     *** parameters: nfa   = facility number
!     ***             exit1 = preemption exit
!     ***
      implicit integer (a - z)
!
!     error check
!     ===========
      if(iabs(fac(nfa,1)).eq.ltx) goto 100
      write(outd,3000) t,tx(ltx,1),tx(ltx,2),nfa
3000  format(1h0,24(1h+),25h ERROR: SUBROUTINE CLEAR ,31(1h+)/1x,&
      24(1h+),3h T=,i7,2x,2hTX,i5,1h,,i3,19h DOESN'T OCCUPY FAC,&
      1x,i4,10(1h+)/)
      return 2
!
!     free
!     ====
100   do 150 i = 1 , 3
150   fac(nfa,i) = 0
      state(nfa) = 1
!
!     block a preempted tx
!     ====================
      if(tx(ltx,6).eq.0) goto 200
      al(ltx,1) = tx(ltx,5)
      al(ltx,2) = - nfa
      tx(ltx,8) = t
200   tx(ltx,5) = 0
      if(iprint.eq.0) goto 5001
      write(outd,3001) t,tx(ltx,1),tx(ltx,2),nfa
3001  format(8h CLEAR :,3x,2hT=,i7,2x,2hTX,i5,1h,,i3,2x,&
      11h LEAVES FAC,i3)
!
!     return
!     ======
5001  if(tx(ltx,6).eq.0) return
      return 1
      end subroutine clear


      subroutine cont
!     ***
!     ***          call cont
!     ***
!     *** function  : restore the system's state from a file with
!     ***             logical device number 'savi'
!     ***
      implicit integer (a-z)
      rewind savi
      read(savi,1) lal,lel
1     format(10i10)
      k=5
      do 2 i=1 , lel , 5
      if(k.gt.lel) k = lel
      read(savi,1) ((el(j,m),m=1,2),j=i,k)
2     k = k +5
      k = 5
      do 3 i = 1 , lal , 5
      if(k.gt.lal) k = lal
      read(savi,1) ((al(j,m),m = 1 , 2 ) , j = i , k )
3     k = k + 5
      do 4 i=1,lal
      k = 10
      do 4 j = 1,tx2,10
      if(k.gt.tx2) k = tx2
      read(savi,1) (tx(i,m),m=j,k)
4     k = k + 10
      k = 10
      do 5 i = 1 , state1 , 10
      if(k.gt.state1) k = state1
      read(savi,1) (state(m),m=i,k)
5     k = k + 10
      read(savi,1) t
      read(savi,1) rt
      do 6 i = 1 , fam1
      k = 10
      do 6 j = 1 , asm1 , 10
      if(k.gt.asm1) k = asm1
      read(savi,1) (asm(i,m),m=j,k)
6     k = k + 10
      k = 3
      do 7 i = 1 , fam1 , 3
      if(k.gt.fam1) k = fam1
      read(savi,1) ((fam(j,m),m=1,3),j=i,k)
7     k = k + 3
      do 8 i = 1 , fam1
      k = 10
      do 8 j = 1 , gathf1 , 10
      if(k.gt.gathf1) k = gathf1
      read(savi,1) (gathf(i,m),m=j,k)
8     k = k + 10
      do 9 i = 1 , fam1
      k = 5
      do 9 j = 1 , uchf1 , 5
      if(k.gt.uchf1) k = uchf1
      read(savi,1) ((uchf(i,m,l),l=1,2),m=j,k)
9     k = k + 5
      k = 3
      do 10 i = 1 , fac1 , 3
      if(k.gt.fac1) k = fac1
      read(savi,1) ((fac(j,m),m=1,3),j=i,k)
10    k = k + 3
      k = 5
      do 11  i = 1 , mfac1 , 5
      if(k.gt.mfac1) k = mfac1
      read(savi,1) ((mfac(j,m),m=1,2),j=i,k)
11    k = k + 5
      k = 10
      do 12  i = 1 , mfac1 , 10
      if(k.gt.mfac1) k = mfac1
      read(savi,1) (mbv(m),m=i,k)
12    k = k + 10
      k = 5
      do 13  i = 1 , mfac1 , 5
      if(k.gt.mfac1) k = mfac1
      read(savi,1) ((plama(j,m),m=1,2),j=i,k)
13    k = k + 5
      k = 3
      do 14  i = 1 , se1 , 3
      if(k.gt.se1) k = se1
      read(savi,1) ((se(j,m),m=1,3),j=i,k)
14    k = k + 3
      k = 5
      do 15  i = 1 , sto1 , 5
      if(k.gt.sto1) k = sto1
      read(savi,1) ((sto(j,m),m=1,2),j=i,k)
15    k = k + 5
      k = 10
      do 16 i = 1 , sto1 , 10
      if(k.gt.sto1) k = sto1
      read(savi,1) (sbv(m),m=i,k)
16    k = k + 10
      k = 5
      do 17  i = 1 , sm1 , 5
      if(k.gt.sm1) k = sm1
      read(savi,1) ((sm(j,m),m=1,2),j=i,k)
17    k = k + 5
      k = 5
      do 18  i = 1 , sto1 , 5
      if(k.gt.sto1) k = sto1
      read(savi,1) ((strama(j,m),m=1,2),j=i,k)
18    k = k + 5
      k = 5
      do 19  i = 1 , pol1 , 5
      if(k.gt.pol1) k = pol1
      read(savi,1) ((pol(j,m),m=1,2),j=i,k)
19    k = k + 5
      k = 10
      do 20  i = 1 , gatht1 , 10
      if(k.gt.gatht1) k = gatht1
      read(savi,1) (gatht(m),m=i,k)
20    k = k + 10
      k = 5
      do 21  i = 1 , ucht1 , 5
      if(k.gt.ucht1) k = ucht1
      read(savi,1) ((ucht(j,m),m=1,2),j=i,k)
21    k = k + 10
      read(savi,1) ntxc
      k = 5
      do 22  i = 1 , src1 , 5
      if(k.gt.src1) k = src1
      read(savi,1) ((src(j,m),m=1,2),j=i,k)
22    k = k + 5
      do 23  i = 1 , bin1
23    read(savi,1) (bin(i,m),m=1,8)
      k = 3
      do 25  i = 1 , bin1 , 3
      if(k.gt.bin1) k = bin1
      read(savi,24) ((binsta(j,m),m=1,2),j=i,k)
24    format(6d15.8)
25    k = k + 3
      k = 6
      do 26  i = 1 , drn1 , 6
      if(k.gt.drn1) k = drn1
      read(savi,24) (drn(m),m=i,k)
26    k = k + 6
      rewind savi
      return
      end subroutine cont



      subroutine depart(nbn,ne,*,iprint)
!     ***
!     ***          call depart(nbn, ne, &1006, iprint)
!     ***
!     *** function  :  withdraw tokens from a bin
!     *** parameters:  nbn = bin number
!     ***              ne = number of departing tokens
!     ***
      implicit integer (a - z)
!
!     error control
!     =============
      do 100  i = 9 , 13
      if(tx(ltx,i).eq.nbn) goto 200
100   continue
      write(outd,3000) t,tx(ltx,1),tx(ltx,2),nbn
3000  format(1h0,24(1h+),26h ERROR: SUBROUTINE DEPART ,30(1h+)/1x,&
      24(1h+),3h T=,i7,2x,2hTX,i5,1h,,i3,14h DOESN'T HAVE ,&
      19(1h+)/1x,24(1h+),25h ANY TOKENS IN BIN NUMBER,i3,1x,27(1h+)/)
      return 1
!
!     clear the bin from the caller's entry
!     =====================================
200   if(bin(nbn,1).lt.ne) goto 500
      inbin = tx(ltx,i+5)
      tx(ltx,i) = 0
      tx(ltx,i+5) = 0
!
!     set the bin matrix
!     ==================
      bin(nbn,7) = bin(nbn,7) + bin(nbn,1) * (t - bin(nbn,8))
      bin(nbn,8) = t
      bin(nbn,1) = bin(nbn,1) - ne
      bin(nbn,4) = bin(nbn,4) + ne
      if(inbin.eq.t) goto 300
      bin(nbn,6) = bin(nbn,6) + (t - inbin) * ne
      goto 400
300   bin(nbn,5) = bin(nbn,5) + ne
!
!     compute the mean values
!     =======================
400   if(t.eq.1) goto 450
      binsta(nbn,1) = bin(nbn,7) / ((bin(nbn,3) + bin(nbn,4)) / 2)
      binsta(nbn,2) = bin(nbn,7) / (bin(nbn,8) - 1.)
450   if(iprint.eq.0) return
      write(outd,3001) t,tx(ltx,1),tx(ltx,2),ne,nbn
3001  format(8h DEPART:,3x,2hT=,I7,2X,2HTX,i5,1h,,i3,2x,&
      10h WITHDRAWS,i5,19h TOKENS FROM BIN NO,i3)
      return
500   write(outd,3002) t,tx(ltx,1),tx(ltx,2),ne,nbn
3002  format(1h0,24(1h+),26h ERROR: SUBROUTINE DEPART ,30(1h+)/1x,&
      24(1h+),3h T=,i7,2x,2hTX,i5,1h,,i3,17h NUMBER OF TOKENS,i4,&
      11(1h+)/1x,24(1h+),32h EXCEEDS AVAILABLE SUPPLY IN BIN,i3,&
      20(1h+)/)
      return 1
      end subroutine depart


      INTEGER FUNCTION DYNPR(LTX1)
!     ***
!     *** FUNCTION  : RESET A TX'S PRIORITY
!     *** PARAMETER : LTX1 = LINE NUMBER OF THE TX WHOSE PRIORITY
!     ***             TO BE CHANGED
!     ***
      IMPLICIT INTEGER (A-Z)
      DYNPR = TX(LTX1,4)
      RETURN
      END FUNCTION DYNPR

      subroutine dynval(k,pcount,iprint)
!     ***
!     ***          call dynval(k,pcount,iprint)
!     ***
!     *** function  : reset the priorities of all txs waiting a
!     ***             station whose station number is k
!     *** parameters: k      = station number
!     ***             pcount = count of reassigned priorities
!     ***
!     ***
      implicit integer (a - z)
!
!     reassign the priorities
!     =======================
      PCOUNT = 0
      DO 100 I = 1,LAL
      IF(AL(I,2).NE.-K.AND.AL(I,2).NE.-K-KEND) GOTO 100
      PR = TX(I,4)
      TX(I,4) = DYNPR(I)
      PCOUNT = PCOUNT +1
      IF(IPRINT.EQ.0) GOTO 100
      WRITE(OUTD,3000) T,TX(I,1),TX(I,2),PR,TX(I,4)
3000  FORMAT(8H DYNVAL:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      13H OLD PRIORITY,I7,14H IS CHANGED TO,I7)
100   CONTINUE
      RETURN
      END subroutine dynval



      SUBROUTINE ENDBIN
!     ***
!     ***          CALL ENDBIN
!     ***
!     *** FUNCTION  : BRING A BIN'S DATA UP TO DATE
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     COMPLETE THE BIN'S DATA
!     =======================
      DO 100 NBN = 1 , BIN1
      IF(BIN(NBN,3).EQ.0) GOTO 100
      BIN(NBN,7) = BIN(NBN,7) + BIN(NBN,1) * (T - BIN(NBN,8))
      BIN(NBN,8) = T
!
!     COMPUTE THE MEAN VALUES
!     =======================
      IF(T.EQ.1) GOTO 100
      BINSTA(NBN,1) = BIN(NBN,7)/((BIN(NBN,3)+BIN(NBN,4))/2.)
      BINSTA(NBN,2) = BIN(NBN,7)/(BIN(NBN,8)-1.)
100   CONTINUE
      RETURN
      END SUBROUTINE ENDBIN


      SUBROUTINE ENTER(NST,NE,ID,LOCK,*,IPRINT)
!     ***
!     ***          CALL ENTER(NST,NE,ID,LOCK,&1005,IPRINT)
!     ***
!     *** FUNCTION  :  ACQUIRE SPACE IN A NON - ADDRESSIBLE STORAGE
!     *** PARAMETERS:  NST  = STORAGE NUMBER
!     ***              NE   = NUMBER OF ELEMENTS REQUESTED
!     ***              ID   = ENTER CALL'S STATEMENT NUMBER
!     ***              LOCK = LOCK FLAG
!     ***                     = 0: NEWLY ARRIVING TXS TRY TO ACQUIRE
!     ***                          SPACE
!     ***                     = 1: NEWLY ARRIVING TXS ARE LOCKED
!     ***                          IMMEDIATELY
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EMFAC + NST
!
!     WATCHDOG
!     ========
      IF(OK.EQ.0) GOTO 10
      OK = 0
!
!     TEST THE SPACE REQUEST
!     ======================
      IF(STO(NST,1)+NE.GT.STO(NST,2)) GOTO 20
!
!     ACQUIRE
!     =======
      STO(NST,1) = STO(NST,1) + NE
      TX(LTX,8) = 0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NE,NST
3000  FORMAT(8H ENTER :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      9H ACQUIRES,I5,16H ELEMENTS IN STO,I3)
      RETURN
!
!     FIRST ACQUISITION ATTEMPT
!     =========================
10    CONTINUE
      IF(IPRINT.EQ.0) GOTO 5001
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NE
3001  FORMAT(8H ENTER :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      9H REQUIRES,I5,13H STO ELEMENTS)
5001  IF(LOCK.GT.0) GOTO 20
      AL(LTX,1) = ID
      AL(LTX,2) = -K
      TX(LTX,8) = T
      RETURN 1
!
!     LOCK
!     ====
20    AL(LTX,1) = ID
      AL(LTX,2) = -KEND-K
      IF(TX(LTX,8).EQ.0)  TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3002) T,TX(LTX,1),TX(LTX,2),NST
3002  FORMAT(8H ENTER :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      17H IS LOCKED AT STO,I3)
      RETURN 1
      END SUBROUTINE ENTER



      SUBROUTINE ERLANG(MEAN,K,MIN,MAX,RNUM,RANDOM,*)
!     ***
!     ***          CALL ERLANG(MEAN,K,MIN,MAX,RNUM,RANDOM,&1006)
!     ***
!     *** FUNCTION  : GENERATE AN ERLANG-DISTRIBUTED RANDOM
!     ***             SEQUENCE
!     *** PARAMETERS: MEAN   = MEAN VALUE
!     ***             K      = DEGREE
!     ***             MIN    = INTERVAL'S LOWER BOUND
!     ***             MAX    = INTERVAL'S UPPER BOUND
!     ***             RNUM   = GENERATOR'S IDENTIFIER
!     ***             RANDOM = RESULTING RANDOM NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
      REAL MAX,MIN,MEAN,RANDOM,ALPHA,ALOG,R
   !   REAL ALPHA,R
!
!     ERROR CONTROL
!     =============
      IF(K-1.GE.0) GOTO 50
      WRITE(OUTD,3000) K
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE ERLANG ,30(1H+),/1X,&
      24(1H+),10H DEGREE K=,I3,16H SMALLER THEN 1 ,27(1H+)/)
      RETURN 1
!
!     GENERATE AN ERLANG-DISTRIBUTED RANDOM NUMBER
!     ============================================
50    ALPHA = FLOAT(K)/MEAN
100   R = 1.0
      DO 150  I = 1 , K
150   R = R * RN(RNUM)
      IF(R.EQ.0) GOTO 100
      RANDOM = -1.0/ALPHA*ALOG(R)
!
!     INTERVAL CHECK
!     ==============
      IF(RANDOM.LT.MIN.OR.RANDOM.GT.MAX) GOTO 100
      RETURN
      END SUBROUTINE ERLANG



      SUBROUTINE EVALUE(NG,TAB,NTAB,IPRINT)
!     ***
!     ***          CALL EVALUE(NG,TAB,NTAB,IPRINT)
!     ***
!     *** FUNCTION  : EVALUATE AND PRINT A FREQUENCY TABLE
!     *** PARAMETERS: NG     = NUMBER OF INTERVALS
!     ***             TAB    = FREQUENCY TABLE
!     ***             NTAB   = TABLE IDENTIFIER
!     ***             IPRINT = PRINT CONTROL
!     ***
      INTEGER ng,ntab,iprint
	  real TAB(MAX_TAB,4),TABX(MAX_TAB)

	  real XQU,VAR,SDV,SUMF,SUMX,SUMXQ,GBR,GM,F
	  integer i,j
!
!     INITIALIZE THE VARIABLES
!     ========================
      XQU = 0.
      VAR = 0.
      SDV = 0.
      SUMF = 0.
      SUMX = 0.
      SUMXQ = 0.
      GBR = TAB(2,1) - TAB(1,1)
      GM = TAB(1,1) - GBR / 2.
!
!     FIRST PASS THROUGH THE TABLE
!     ============================
      DO 500  I = 1 , NG
      F = TAB(I,2)
      SUMF = SUMF + F
      SUMX = SUMX + F * GM
      SUMXQ = SUMXQ + F * GM * GM
500   GM = GM + GBR
!
!     SECOND PASS THROUGH THE TABLE
!     =============================
      IF(SUMF.EQ.0) GOTO 700
      DO 600  I = 1 , NG
      TAB(I,3) = TAB(I,2) / SUMF
      IF(TAB(I,2).NE.0) TABX(I) = TAB(I,4) / TAB(I,2)
600   CONTINUE
!
!     COMPUTER MEAN VALUE, VARIANCE AND STANDARD DEVIATION
!     ====================================================
      XQU = SUMX / SUMF
      IF(SUMF.EQ.1.) VAR = 0.
      IF(SUMF.NE.1) VAR = (SUMXQ - SUMX * SUMX /SUMF) / (SUMF - 1.)
      SDV = SQRT(VAR)
!
!     PRINT THE RESULTS
!     =================
700   IF(IPRINT.GT.0) WRITE(OUTD,3000)
3000  FORMAT(1H1)
      WRITE(OUTD,3001) NTAB
3001  FORMAT(1H0,15X,5HTABLE,I3/)
      IF(IPRINT.EQ.0) GOTO 900
      WRITE(OUTD,3002)
3002  FORMAT(3X,1HI,12X,1HX,13X,4HF(X),10X,3HF/N,&
      11X,4HCUMY,11X,4HE(Y)/)
      DO 800  I = 1 , NG
      IF(TAB(I,2).EQ.0) GOTO 800
      WRITE(OUTD,3003) I,(TAB(I,J),J=1,4),TABX(I)
3003  FORMAT(1X,I3,3(5X,F10.3),5X,E10.3,5X,F10.3)
800   CONTINUE
900   WRITE(OUTD,3004) SUMF,XQU,VAR,SDV
3004  FORMAT(/2X,6H SUMF=,F10.1,2X,5H XQU=,F10.2,2X,5H VAR=,F10.1,&
      5H SDV=,F10.2//)
      RETURN
      END SUBROUTINE EVALUE

      subroutine EVENT(EVT,EVAD,NS,*,IPRINT)
!     ***
!     ***          CALL EVENT(EVT,EVAD,NS,&1006,IPRINT)
!     ***
!     *** FUNCTION  : SCHEDULE AN EVENT
!     *** PARAMETERS: EVT  = EVENT TIME
!     ***             EVPD = TARGET
!     ***             NS   = SOURCE NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     RESET ACTIVATION TIME
!     =====================
      IF(NS.EQ.0) GOTO 50
      IF(NS.LE.SRC1) GOTO 20
      WRITE(OUTD,3000) T,NS
3000  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE EVENT ,31(1H+)/1X,&
      24(1H+),3H T=,I7,14H SOURCE NUMBER,I3,1X,28(1H+)/1X,&
      24(1H+),37H EXCEEDS AVAILABLE NUMBER OF SOURCES ,19(1H+)/)
      RETURN 1
20    IF(EVT.EQ.0) GOTO 200
      IF(SRC(NS,1).EQ.0) GOTO 50
      I=SRC(NS,1)
      GOTO 150
!
!     SHEDULE
!     =======
50    DO 100  I = 1 , EL1
      IF(EL(I,1).EQ.0) GOTO 150
100   CONTINUE
      WRITE(OUTD,3001) T,EVT,EVAD
3001  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE EVENT ,31(1H+)/1X,&
      24(1H+),3H T=,I7,21H EVENT LIST OVERFLOW ,25(1H+)/&
      1X,24(1H+),9H FOR EVT=,I7,6H EVAD=,I3,1X,30(1H+)/)
      RETURN 1
150   EL(I,1) = EVAD
      EL(I,2) = EVT
      IF(LEL.LT.I) LEL = I
!
!     SCHEDULE A TX GENERATION
!     ========================
      IF(NS.GT.0) SRC(NS,1) = I
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3002) T,EVT,EL(I,1)
3002  FORMAT(8H EVENT :,3X,2HT=,I7,15X,&
      23H EVENT SCHEDULED FOR T=,I7,8H TARGET=,I3)
      RETURN
!
!     SHUT DOWN A SOURCE
!     ==================
200   I = SRC(NS,1)
      EL(I,1) = 0
      EL(I,2) = 0
      IF(IPRINT.EQ.0) GOTO 250
      WRITE(OUTD,3003) T,NS
3003  FORMAT(8H EVENT :,3X,2HT=,I7,15X,7H SOURCE,I3,&
      13H IS SHUT DOWN)
!
!     RESET THE LIST-END POINTER, LEL
!     ===============================
250   IF(EL(LEL,1).NE.0.OR.LEL.EQ.1) RETURN
      LEL = LEL - 1
      GOTO 250
      END subroutine EVENT

      SUBROUTINE FFIT(NST,NE)
!     ***
!     ***          CALL FFIT(NST ,NE)
!     ***
!     *** FUNCTION  : FIND A FREE SPACE USING THE FIRST-FIT
!     ***             STRATEGY
!     *** PARAMETERS: NST = STORAGE NUMBER
!     ***             NE  = NUMBER OF ELEMENTS TO BE ACQUIRED
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     INITIALIZE THE SEARCH
!     =====================
      LSM = 0
      I = SBV(NST)
      IE = I + STO(NST,2) - 1
!
!     FIND A SEGMENT
!     ==============
10    IF(SM(I,2).NE.-1) GOTO 20
      IF(SM(I,1).GE.NE) GOTO 30
20    I = I + SM(I,1)
      IF(I.LE.IE) GOTO 10
      RETURN
30    LSM = I
      RETURN 
      END SUBROUTINE FFIT


      SUBROUTINE FIFO
!     ***
!     ***          CALL FIFO
!     ***
!     *** FUNCTION : SELECT A TX ON THE BASE OF FIFO
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     INITIALIZE THE SEARCH
!     =====================
      LTX = POLVEC(1)
      IF(POLC.EQ.1) RETURN
!
!     COMPARE
!     =======
      DO 150  I = 2 , POLC
      LTX1 = POLVEC(I)
100   IF(TX(LTX,8)-TX(LTX1,8)) 150 , 110 , 130
110   IF(TX(LTX,1)-TX(LTX1,1)) 150 , 120 , 130
120   IF(TX(LTX,2)-TX(LTX1,2)) 150 , 150 , 130
130   LTX = LTX1
150   CONTINUE
      RETURN
      END SUBROUTINE FIFO


      SUBROUTINE FREE(NST,NE,KEY,LINE,*,*,IPRINT)
!     ***
!     ***          CALL FREE(NST,NE,KEY,LINE,EXIT1,&1006,IPRINT)
!     ***
!     *** FUNCTION  : FREE AN AREA IN ADDRESSIBLE STORAGE
!     *** PARAMETERS: NST   = STORAGE NUMBER
!     ***             NE    = NUMBER OF ELEMENTS TO BE FREED
!     ***             KEY   = FREEING KEY
!     ***             LINE  = REMAINDER ADDRESS
!     ***             EXIT1 = ABORTED-FREEING EXIT
!     ***
      IMPLICIT INTEGER (A - Z)
      IF(STRAMA(NST,1).EQ.0) GOTO 50
!
!     DETERMINE LSM
!     =============
      LSM = KEY + SBV(NST) - 1
      IF(STRAMA(NST,2).NE.0) CALL STRATF(NST,KEY,*1)
1     IF(LSM.EQ.0) GOTO 40
!
!     MARK REMAINDER
!     ==============
      LINE = 0
      IF(SM(LSM,1).LT.NE) GOTO 40
      IF(SM(LSM,1).EQ.NE) GOTO 10
      J = LSM + NE
      SM(J,1) = SM(LSM,1) - NE
      SM(J,2) = SM(LSM,2)
      LINE = J - SBV(NST) + 1
      SM(LSM,1) = NE
!
!     FREE THE STORAGE AREA
!     =====================
10    SM(LSM,2) = -1
      STO(NST,1) = STO(NST,1) - NE
      IF(IPRINT.EQ.0) GOTO 5000
      IADDR = LSM - SBV(NST) + 1
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),SM(LSM,1),NST,IADDR
3000  FORMAT(8H FREE  :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      6H FREES,I5,16H ELEMENTS IN STO,I3,21H, STARTING AT ADDRESS)
!
!     LUMB CONTIGUOUS FREE SEGMENTS INTO ONE SEGMENT
!     ==============================================
5000  J = LSM +SM(LSM,1)
      IF(SM(J,2).NE.-1.OR.J.GE.SBV(NST)+STO(NST,2)) GOTO 20
      SM(LSM,1) = SM(LSM,1)+SM(J,1)
      SM(J,1) = 0
      SM(J,2) = 0
20    IF(LSM.EQ.SBV(NST)) RETURN
      DO 100 I = 1, SM1
      J=LSM-I
      IF(SM(J,1).NE.0) GOTO 30
100   CONTINUE
30    IF(SM(J,2).NE.-1) RETURN
      SM(J,1) = SM(LSM,1) +I
      SM(LSM,1)=0
      SM(LSM,2)=0
      RETURN
!
!     SPECIAL EXITS
!     =============
40    WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NST
3001  FORMAT(1H0,24(1H+),26H WARNING: SUBROUTINE FREE ,30(1H+)/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,10H TOO MUCH ,&
      23(1H+),/1X,24(1H+),25H SPACE BEING FREED IN STO,&
      I3,1X,27(1H+)/)
      RETURN 1
50    WRITE(OUTD,3002) T,TX(LTX,1),TX(LTX,2),NST
3002  FORMAT(1H0,24(1H+),24H ERROR: SUBROUTINE FREE ,32(1H+),1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,17H NON-ADDRESSIBLE ,&
      16(1H+),/1X,24(1H+),8H STORAGE,I3,17H CANNOT BE FREED ,&
      28(1H+),1X,24(1H+),23H USING THIS SUBROUTINE ,33(1H+)/)
      RETURN 2
      END SUBROUTINE FREE



      subroutine GATE1(LOGEXP,GLOBAL,NGATE,ID,LOCK,*,IPRINT)
!     ***
!     ***          CALL GATE1(LOGEXP,GLOBAL,NGATE,ID,LOCK,&1006,IPRINT)
!     ***
!     *** FUNCTION  : LOCK TXS OR LET THEM PASS, DEPENDING ON THE
!     ***             VALUE OF A LOGICAL EXPRESSION
!     *** PARAMETERS: LOGEXP = WAIT CONDITION
!     ***             GLOBAL = PARAMETER CODE
!     ***                      = 0: THE LOGICAL EXPRESSION CONT
!     ***                           LOCAL PARAMETERS
!     ***                      = 1: THE LOGICAL EXPRESSION CONT
!     ***                           ONLY GLOBAL PARAMETERS
!     ***             NGATE  = TYPE-1 GATE NUMBER
!     ***             ID     = GATE1 CALL'S STATEMENT NUMBER
!     ***             LOCK   = LOCK FLAG
!     ***                      = 0: NEWLY ARRIVING TXS TRY THE
!     ***                      = 1: NEWLY ARRIVING TXS ARE LOCK
!     ***                           IMMEDIATELY
!     ***
      IMPLICIT INTEGER(A - Z)
      LOGICAL LOGEXP
      IF(IPRINT.EQ.0) GOTO 5000
	  if (LTX>MAX_TX .or. LTX<1) THEN
         write(*,*) "**ERROR** GATE1. T=",T," LTX=",LTX," MAX_TX=",MAX_TX," ngate=",NGATE
		 stop
	  end if
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NGATE
3000  FORMAT(8H GATE1 :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      20H ARRIVES AT GATE1 NO,I3)
!
!     DETERMINE THE STATION NUMBER
!     ============================
5000  K = ESTO + NGATE
!
!     WATCHDOG
!     ========
      IF(OK.EQ.0) GOTO 10
      OK = 0
!
!     LOCKING DECISION
!     ================
      IF(.NOT.LOGEXP) GOTO 20
      TX(LTX,8) = 0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3001)
3001  FORMAT(36X,32HAND PROCEEDS TO THE NEXT STATION)
      RETURN
!
!     NEW ARRIVAL
!     ===========
10    IF(LOCK.GT.0) GOTO 30
      AL(LTX,1) = ID
      AL(LTX,2) = - K
      TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3002)
3002  FORMAT(36X,14HAND IS BLOCKED)
      RETURN 1
!
!     TEST THE PARAMETER CODE
!     =======================
20    IF(GLOBAL.GT.0) GOTO 40
!     LOCK THE ACTIVE TX
!     ==================
30    AL(LTX,1) = ID
      AL(LTX,2) = - KEND - K
      IF(TX(LTX,8).EQ.0) TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3003)
3003  FORMAT(36X,13HAND IS LOCKED)
      RETURN 1
!
!     LOCK ALL STARTED TXS
!     ====================
40    AL(LTX,1) = ID
      AL(LTX,2) = - KEND - K
      IF(TX(LTX,8).EQ.0) TX(LTX,8) = T
      DO 100  I = 1 , LAL
      IF(AL(I,2).EQ.-K) AL(I,2) = - KEND - K
100   CONTINUE
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3003)
      RETURN 1
      END subroutine GATE1



      SUBROUTINE GATE2(LOGEXP,NGATE,ID,*,IPRINT)
!     ***
!     ***          CALL GATE2(LOGEXP,NGATE,ID,&1005,IPRINT)
!     ***
!     *** FUNCTION  : LOCK TXS OR LET THEM PASS, DEPENDING ON T
!     ***             OF A LOGICAL EXPRESSION
!     *** PARAMETERS: LOGEXP = WAIT CONDITION
!     ***             NGATE  = TYPE-2 GATE NUMBER
!     ***             ID     = GATE2 CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A - Z)
      LOGICAL LOGEXP
      IF(IPRINT.EQ.0) GOTO 5000
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NGATE
3000  FORMAT(8H GATE2 :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      20H ARRIVES AT GATE2 NO,I3)
!
!     DETERMINE THE STATION NUMBER
!     ============================
5000  K = EGATE1 + NGATE
!
!     BLOCKING DECISION
!     =================
      IF(OK.EQ.0) GOTO 250
      OK = 0
      IF(.NOT.LOGEXP) GOTO 200
100   TX(LTX,8) = 0
      IT = 0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3001)
3001  FORMAT(36X,32HAND PROCEEDS TO THE NEXT STATION)
      RETURN
!
!     MARK THE GATE INACCESSIBLE
!     ==========================
200   STATE(K) = 0
!
!     BLOCK
!     =====
250   AL(LTX,1) = ID
      AL(LTX,2) = - K
      IF(TX(LTX,8).EQ.0) TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3002)
3002  FORMAT(36X,14HAND IS BLOCKED)
      RETURN 1
      END SUBROUTINE GATE2


      SUBROUTINE GATHR1(NG,NGATH,ID,*,IPRINT)
!     ***
!     ***          CALL GATHR1(NG,NGATH,ID,&1005,IPRINT)
!     ***
!     *** FUNCTION  : GATHER A GROUP OF TXS WITH REGARD TO FAMILY
!     ***             MEMBERSHIP
!     *** PARAMETERS: NG    = NUMBER OF TXS TO COLLECT
!     ***             NGATH = TYPE-1 GATHER STATION NUMBER
!     ***             ID    = GATE1 CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     TEST FOR FAMILY MEMBERSHIP
!     ==========================
      IF(LFAM.EQ.0) RETURN
!
!     CLEAR THE BLOCKING TIME
!     =======================
      IF(TX(LTX,8).EQ.0) GOTO 100
      TX(LTX,8) = 0
      RETURN
!
!     DETERMINE THE STATION NUMBER
!     ============================
100   K = EGATE2 + FAM1 * (NGATH - 1) + LFAM
      IF(IPRINT.EQ.0) GOTO 5000
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NGATH
3000  FORMAT(8H GATHR1:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      29H ARRIVES AT GATHR1 STATION NO,I3)
!
!     LOCK
!     ====
5000  AL(LTX,1) = ID
      AL(LTX,2) = - KEND - K
      TX(LTX,8) = T
      GATHF(LFAM,NGATH) = GATHF(LFAM,NGATH) + 1
!
!     TEST THE COUNT
!     ==============
      IF(GATHF(LFAM,NGATH).LT.NG) GOTO 200
!
!     START THE GATHERED TXS
!     ======================
      IF(IPRINT.EQ.0) GOTO 5001
      WRITE(OUTD,3001)
3001  FORMAT(36X,29HAND COMPLETES THIS COLLECTION)
5001  CALL UNLOCK(K,IPRINT)
      GATHF(LFAM,NGATH) = 0
      RETURN 1
200   CONTINUE
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3002)
3002  FORMAT(36X,13HAND IS LOCKED)
      RETURN 1
      END SUBROUTINE GATHR1


      SUBROUTINE GATHR2 (NG,NGATH,ID,*,IPRINT)
!     ***
!     ***          CALL GATHR2(NG,NGATH,ID,&1005,IPRINT)
!     ***
!     *** FUNCTION  : GATHER A GROUP OF TXS WITHOUT REGARD
!     ***             TO FAMILY MEMBERSHIP
!     *** PARAMETERS: NG    = NUMBER OF TXS TO COLLECT
!     ***             NGATH = TYPE-2 GATHER STATION NUMBER
!     ***             ID    = GATE2 CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     CLEAR THE BLOCKING TIME
!     =======================
      IF(TX(LTX,8).EQ.0) GOTO 100
      TX(LTX,8) = 0
      RETURN
!
!     COMPUTE THE STATION NUMBER
!     ==========================
100   K = EGATHF + NGATH
      IF(IPRINT.EQ.0) GOTO 5000
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NGATH
3000  FORMAT(8H GATHR2:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      29H ARRIVES AT GATHR2 STATION NO,I3)
!
!     LOCK
!     ====
5000  AL(LTX,1) = ID
      AL(LTX,2) = - KEND - K
      TX(LTX,8) = T
      GATHT(NGATH) = GATHT(NGATH) + 1
!
!     TEST THE COUNT
!     ==============
      IF(GATHT(NGATH).LT.NG) GOTO 200
!
!     START THE GATHERED TXS
!     ======================
      IF(IPRINT.EQ.0) GOTO 5001
      WRITE(OUTD,3001)
3001  FORMAT(36X,29HAND COMPLETES THIS COLLECTION)
5001  CALL UNLOCK(K,IPRINT)
      GATHT(NGATH)=0
      RETURN 1
200   CONTINUE
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3002)
3002  FORMAT(36X,13HAND IS LOCKED)
      RETURN 1
      END SUBROUTINE GATHR2


      SUBROUTINE GAUSS(MEAN,SIGMA,MIN,MAX,RNUM,RANDOM)
!     ***
!     ***          CALL GAUSS(MEAN,SIGMA,MIN,MAX,RNUM,RANDOM)
!     ***
!     *** FUNCTION  : GENERATE A GAUSS-DISTRIBUTED RANDOM SEQUENCE
!     *** PARAMETERS: MEAN   = MEAN VALUE
!     ***             SIGMA  = STANDARD DEVIATION
!     ***             MIN    = INTERVAL'S LOWER BOUND
!     ***             MAX    = INTERVAL'S UPPER BOUND
!     ***             RNUM   = GENERATOR'S IDENTIFIER
!     ***             RANDOM = RESULTING RANDOM NUMBER
!     ***
      INTEGER RNUM
      REAL MAX,MIN,MEAN,SIGMA,RANDOM
	  REAL R,V
!	  
!     GENERATE A GAUSS-DISTRIBUTED RANDOM NUMBER
!     ==========================================
100   R = RN(RNUM)
      IF(R.EQ.0.) GOTO 100
      V = (-2.0 *ALOG(R))**0.5*COS(6.2832*RN(RNUM))
      RANDOM = V * SIGMA + MEAN
!
!     INTERVAL CHECK
!     ==============
      IF(RANDOM.LT.MIN.OR.RANDOM.GT.MAX) GOTO 100
      RETURN  
      END SUBROUTINE GAUSS



      SUBROUTINE GENERA(ET,ZTX,PR,ID,*,IPRINT)
!
!     ***
!     ***          CALL GENERA(ET,LTX,PR,ID,&1006,IPRINT)
!     *** FUNCTION  : CREATE TXS
!     *** PARAMETERS: ET  = EMISSION INTERVAL
!     ***             ZTX = GENERATION LIMIT
!     ***             PR  = PRIORITY
!     ***             ID  = GENERA CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     FIND THE SOURCE
!     ===============
      DO  100 J = 1,SRC1
      IF(SRC(J,1).EQ.LEV) GOTO 150
100   CONTINUE
      WRITE(OUTD,3000) T,LEV
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE GENERA ,30(1H+)/1X,&
      24(1H+),3H T=,I7,5H LINE,I3,19H OF THE EVENT LIST ,&
      19(1H+)/1X,24(1H+),20H CONTAINS NO SOURCE ,36(1H+)/)
      RETURN 1
!
!     ADVANCE THE COUNTERS
!     ====================
150   SRC(J,1) = 0
      SRC(J,2) = SRC(J,2) + 1
      NTXC = NTXC + 1
!
!     GENERATE
!     ========
 !     DO 200 LTX = 1 , TX1
 !     IF(TX(LTX,1).EQ.0) GOTO 250
!200   CONTINUE
      DO LTX = 1 , TX1
         IF(TX(LTX,1).EQ.0) GOTO 250
      END DO
      WRITE(OUTD,3001) T,J
3001  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE GENERA ,30(1H+)/1X,&
      24(1H+),3H T=,I7,20H TX MATRIX OVERFLOW ,26(1H+)/1X,24(1H+),&
      11H FOR SOURCE,I3,1X,41(1H+)/)
      RETURN 1
250   TX(LTX,1) = NTXC
      TX(LTX,3) = T
      TX(LTX,4) = PR
      IF(IPRINT.EQ.0) GOTO 5002
      WRITE(OUTD,3002) T,NTXC,J
3002  FORMAT(8H GENERA:,3X,2HT=,I7,2X,2HTX,I5,6H,  0  ,&
      24H IS GENERATED FOR SOURCE,I3)
!
!     RESET THE EVENT LIST
!     ====================
5002  IF(LTX.GT.LAL) LAL = LTX
!     write(*,*) 'J=',j,' SRC(J,2)=',SRC(J,2),' ZTX=',ZTX
      IF(SRC(J,2).GE.ZTX) RETURN
!     write(*,*) 'passed'
      SRC(J,1) = LEV
      EL(LEV,1) = ID
      EL(LEV,2) = T + ET
      IF(LEV.GT.LEL) LEL = LEV
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3003) EL(LEV,2)
3003  FORMAT(36X,40HNEXT TX GENERATION FOR THIS SOURCE AT T=,I7)
      RETURN
      END SUBROUTINE GENERA



      SUBROUTINE GRAPH(TAB,NTAB,IU,IO,Y,*)
!     ***
!     ***          CALL GRAPH(TAB,NTAB,IU,IO,Y,*1006)
!     ***
!     *** FUNCTION  : PRINT A HISTOGRAM OF A FREQUENCY TABLE
!     *** PARAMETERS: TAB  = FREQUENCY TABLE
!     ***             NTAB = TABLE IDENTIFIER
!     ***             IU   = NUMBER OF LEFTMOST INTERVAL
!     ***             IO   = NUMBER OF RIGHTMOST INTERVAL
!     ***             Y    = CODE
!     ***                    = 0: RELATIVE FREQUENCY
!     ***                    = 1: ABSOLUTE FREQUENCY
!     ***
      INTEGER Y,NTAB,IU,IO
      REAL MAX
      CHARACTER*28 L(15)
      REAL TAB(max_tab,4)
      character*1 LINE(MAX_TAB)
      real M(MAX_TAB)
!	  INTEGER*1 KSTERN,KBLANC,KSTR,KPL
      integer i,j,ib,int,int1,int2,i1,i2,imax,inty,modul1,intb,ipl,mult,is
      real sum,ymax,yvall,ycorr
      CHARACTER*1 KSTERN,KBLANC,KSTR,KPL
      DATA KSTERN/'*'/, KBLANC/' '/,KSTR/'-'/,KPL/'+'/
      DATA L(1)/&
      '(18X,5(F6.0,14X))           '/,&
      L(2)/'(16X,6(F6.0,10X))           '/,&
      L(3)/'(15X,7(F6.0,8X))            '/,&
      L(4)/'(14X,8(F6.0,6X))            '/,&
      L(5)/'(14X,9(F6.0,5X))            '/,&
      L(6)/'(13X,F6.0,34X,F6.0,44X,F6.0)'/,&
      L(7)/'(13X,F6.0,30X,F6.0,39X,F6.0)'/,&
      L(8)/'(12X,F6.0,26X,F6.0,34X,F6.0)'/,&
      L(9)/'(12X,F6.0,22X,F6.0,29X,F6.0)'/,&
      L(10)/'(11X,F6.0,18X,3(F6.0,24X))  '/,&
      L(11)/'(11X,F6.0,14X,4(F6.0,19X))  '/,&
      L(12)/'(10X,F6.0,10X,5(F6.0,14X))  '/,&
      L(13)/'(10X,F6.0,6X,6(F6.0,9X))    '/,&
      L(14)/'(9X,F6.0,2X,10(F6.0,4X))    '/,&
      L(15)/'(9X,F6.0,3X,10(F6.0,4X))    '/
!
!     ERROR CHECK
!C     ===========
      IF(IO-IU.GT.99.OR.IO-IU.LT.4) GOTO 900
!
!     INITIALIZE THE VARIABLES
!     ========================
      DO 1 I=1,MAX_TAB
      M(I)=0.
1     LINE(I)=KBLANC
      INT = 100/(IO-IU+1)
      INT1 = 1
      IF(INT.EQ.1) GOTO 2
      INT1=INT/2
      INT2=INT1
      IF(INT1+INT2.LT.INT) INT1=INT1+1
2     MAX=0.
!
!     FIND THE TOPS OF THE COLUMNS
!     ============================
      I1=0
      IF(Y.GT.0) GOTO 5
      SUM=0.
      DO 3 I=1,TAB1
3     SUM=SUM+TAB(I,2)
5     DO 50 I = IU,IO
      I1=I1+1
      IF(Y.GT.0) M(I1)=TAB(I,2)
      IF(Y.EQ.0) M(I1)=TAB(I,2)/SUM*100.
      IF(M(I1).GT.MAX) MAX=M(I1)
      IF(INT.EQ.1) GOTO 50
      IF(INT1.EQ.1) GOTO 20
      DO 10 I2=2,INT1
      I1=I1+1
10    M(I1)=M(I1-1)
20    I1=I1+INT2
50    CONTINUE
!
!     SCALE THE FREQUENCY AXIS
!     ========================
      IND = 0
      IMAX = MAX
      DO 60 I=1,20
      IF(IMAX.LT.10.OR.IMAX.EQ.10.AND.MAX.EQ.10**(IND+1)) GOTO 70
      IND=IND+1
      IMAX=IMAX/10
60    CONTINUE
70    IF(MAX.GT.IMAX*10**IND) IMAX=IMAX+1
      YMAX=IMAX*10.**IND
      IF(IMAX.LE.5.AND.YMAX.GT.10) IMAX=IMAX*2
!
!     PRINT A TITLE FOR THE GRAPH
!     ===========================
100   WRITE(OUTD,3000) NTAB
3000  FORMAT(1H1,15X,6HTABLE ,I3//)
      IF(Y.EQ.0) WRITE(OUTD,3001)
3001  FORMAT(4X,7HPERCENT/)
!
!     PRINT THE HISTOGRAM COLUMNS
!     ===========================
      INTY = 50
	  if (IMAX.EQ.0) THEN 
	     MODUL1=0
	  ELSE
         MODUL1 = MOD(50,IMAX)
	  END IF
      INTY=INTY-MODUL1
      YVALL=YMAX/INTY
      YCORR=YVALL/1000.
	  IF (IMAX.EQ.0) THEN
        INTB=10
      ELSE 
        INTB=INTY/IMAX
      END IF
      J=0
      DO 150 I=1,INTY
      DO 120 I1=1,MAX_TAB
      IF(LINE(I1).NE.KBLANC) GOTO 120
      IF(M(I1)+YCORR.GE.YMAX) LINE(I1)=KSTERN
120   CONTINUE
      J=J+1
      IF(J.GT.1) GOTO 130
      IMAX=IFIX(YMAX+0.5)
      WRITE(OUTD,3002) IMAX,(LINE(I1),I1=1,100)
3002  FORMAT(1X,I10,2H I,100A1)
      GOTO 140
130   WRITE(OUTD,3003) (LINE(I1),I1=1,100)
3003  FORMAT(12X,1HI,100A1)
140   YMAX=YMAX-YVALL
      IF(J.EQ.INTB) J=0
150   CONTINUE
!
!     PRINT THE X AXIS
!     ================
      IPL = INT1
      MULT=1
      IF(INT.GT.1.AND.INT.LT.11) MULT=4
      IF(INT.EQ.1) MULT=9
      DO 170 I=1,100
      IF(I.EQ.IPL) GOTO 160
      LINE(I)=KSTR
      GOTO 170
160   LINE(I)=KPL
      IPL=IPL+MULT*INT
      IF(MULT.EQ.4.OR.MULT.EQ.9) MULT=MULT+1
170   CONTINUE
      WRITE(OUTD,3004) (LINE(I1),I1=1,100)
3004  FORMAT(10X,3H0 I,100A1)
      IS=16-INT
      IF(INT.EQ.14) IS=3
      IF(INT.EQ.16) IS=2
      IF(INT.EQ.20) IS=1
      IB=IU+1
      IS=1
      IF(INT.GT.10) GOTO 190
      IF(INT.GT.1) IS=5
      IF(INT.EQ.1) IS=10
      IB=IU+IS-1
190   WRITE(OUTD,L(IS)) TAB(IU,1),(TAB(I,1),I=IB,IO,IS)
      RETURN
!
!     ERROR EXIT
!     =========
900   WRITE(OUTD,3020) IU,IO
3020  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE GRAPH ,31(1H+)/1X,&
      24(1H+),4H IO=,I3,25H IO-IU > 99 OR IO-IU < 4 ,17(1H+)/)
      RETURN 1
      END SUBROUTINE GRAPH


      SUBROUTINE INIT1
!     ***
!     ***          CALL INIT1
!     *** FUNCTION  : INITIALIZE VALUES FOR RANDOM NUMBER GENERATORS
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     SET THE MULTIPLIERS
!     ===================
      DFACT( 1) =  10753813.D0
      DFACT( 2) = 228181237.D0
      DFACT( 3) = 348984893.D0
      DFACT( 4) = 590634277.D0
      DFACT( 5) =  10841757.D0
      DFACT( 6) = 107408213.D0
      DFACT( 7) = 469770781.D0
      DFACT( 8) = 832159853.D0
      DFACT( 9) =  10842045.D0
      DFACT(10) = 228244373.D0
      DFACT(11) = 348969437.D0
      DFACT(12) =  10763125.D0
      DFACT(13) = 107452989.D0
      DFACT(14) = 469773661.D0
      DFACT(15) = 590648085.D0
      DFACT(16) =  10784189.D0
      DFACT(17) = 228206429.D0
      DFACT(18) = 469793853.D0
      DFACT(19) = 832230813.D0
      DFACT(20) = 107464661.D0
      DFACT(21) = 469808301.D0
      DFACT(22) = 590585901.D0
      DFACT(23) = 228178005.D0
      DFACT(24) = 349045949.D0
      DFACT(25) =  10845149.D0
      DFACT(26) = 107441485.D0
      DFACT(27) = 469795933.D0
      DFACT(28) = 349062285.D0
      DFACT(29) = 107380645.D0
      DFACT(30) =  10767581.D0
!
!     SPECIFY CONSTANTS AND MODULO
!     ============================
      DMODUL = 2.**30
      DO 10  I = 1 , DRN1
10    DCONST(I) = 227623267.D0
!
!     SET THE STARTING VALUES
!     =======================
      DO 20 I=1,DRN1
20    DRN(I)=1.D0
      RETURN
      END SUBROUTINE INIT1

      SUBROUTINE INIT2(*)
!     ***
!     ***          CALL INIT2(&9999)
!     ***
!     *** FUNCTION  : SET UP DATA AREAS FOR MULTIFACILITIES
!     *** PARAMETER : EXIT1 = LIST-END EXIT
!     ***

      IMPLICIT INTEGER (A-Z)
!
!     SPECIFIY THE MULTIFACILITIES USED
!     =================================
      LSE = 1
      DO  100 I=1,MFAC1
      IF(MFAC(I,2).EQ.0) GOTO 100
!
!     SET UP A SECTION OF THE SE MATRIX
!     =================================
      MBV(I) = LSE
      LSE = LSE + MFAC(I,2)
      IF(LSE-1.GT.SE1) GOTO 200
100   CONTINUE
      RETURN
200   WRITE(OUTD,3000) I
3000  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE INIT2 ,31(1H+)/1X,&
      24(1H+),28H SE MATRIX OVERFLOW FOR MFAC,I3,1X,24(1H+)/)
      RETURN 1
      END SUBROUTINE INIT2

      SUBROUTINE INIT3(*)
!     ***
!     ***          CALL INIT3(&9999)
!     ***
!     *** FUNCTION  : SET UP DATA AREAS FOR ADDRESSIBLE STORAGES
!     *** PARAMETER : EXIT1 = LIST-END EXIT
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     CHEK FOR ADDRESSIBILITY
!     =======================
      LSM = 1
      DO 100 I=1,STO1
      IF(STRAMA(I,1).EQ.0) GOTO 100
!
!     SET UP A SECTION IN THE SEGMENT MATRIX
!     ======================================
      SBV(I) = LSM
      SM(LSM,1) = STO(I,2)
      SM(LSM,2) = -1
      LSM = LSM + STO(I,2)
      IF(LSM-1.GT.SM1) GOTO 200
100   CONTINUE
      RETURN
200   WRITE(OUTD,3000) I
3000  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE INIT3 ,31(1H+),/1X,&
      24(1H+),27H SM MATRIX OVERFLOW FOR STO,I3,1X,25(1H+)/)
      RETURN 1
      END SUBROUTINE INIT3



      subroutine initd
      implicit integer (a-z)
 
      el1    = MAX_EL      !  30
      TX1    = MAX_TX
      TX2    = MAX_TX2
      FAM1   = MAX_FAM     ! 300
      SRC1   = MAX_SRC     !  20
!      ASM1   = MAX_ASM     !  10
      ASM1   = 10          !  10
      FAC1   = MAX_FAC     !  40
      MFAC1  = MAX_MFAC    !  20
      SE1    = MAX_SE      !  40
      STO1   = MAX_STO     !  40
      SM1    = MAX_SM      !1024
      NGATE1 = 25
      NGATE2 = 25
      GATHT1 = 10
      GATHF1 = 10
      UCHT1  = 10
      UCHF1  = 10
      KEND   = 9200
      STATE1 = 9200
      EFAC   = FAC1
      EMFAC  = EFAC+MFAC1
      ESTO   = EMFAC+STO1
      EGATE1 = ESTO+NGATE1
      EGATE2 = EGATE1+NGATE2
      EGATHF = EGATE2+FAM1*GATHF1
      EGATHT = EGATHF+GATHT1
      EUCHF  = EGATHT+UCHF1*FAM1*2
      BMFAC  = EFAC+1
      BSTO   = EMFAC+1
      BGATE2 = EGATE1+1
      POL1   = max_pol      !  20
      BIN1   = max_bin      !  50  
      TAB1   = max_tab      ! 100
      DRN1   = max_drn      ! 
      IND    =  5
      OUTD   =  7
      SAVI   = 25
      SAVO   = 26
      RETURN
      END subroutine initd


      SUBROUTINE KNOCKD(KT,NFA,IDN,*,*,IPRINT)
!     ***
!     ***          CALL KNOCKD(KT,NFA,IDN,&1005,&1006,IPRINT)
!     ***
!     *** FUNCTION  :  KNOCK DOWN A FACILITY
!     *** PARAMETERS:  KT  = KNOCKDOWN TIME
!     ***              NFA = FACILITY NUMBER
!     ***              IDN = TARGET
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     ERROR CHECK
!     ===========
      IF(IABS(FAC(NFA,1)).EQ.LTX) GOTO 100
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NFA
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE KNOCKD ,30(1H+)/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,19H DOESN'T OCCUPY FAC,&
      1X,10(1H+)/)
      RETURN 2
!
!     KNOCK DOWN THE FACILITY
!     =======================
100   FAC(NFA,1)=-LTX
      FAC(NFA,3)=3
      AL(LTX,1)=IDN
      AL(LTX,2)=T+KT
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NFA,AL(LTX,2)
3001  FORMAT(8H KNOCKD:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      4H FAC,I3,30H WILL BE KNOCKED DOWN UNTIL T=,I7)
      RETURN 1
      END SUBROUTINE KNOCKD



      SUBROUTINE LEAVE(NST,NE,*,IPRINT)
!     ***
!     ***          CALL LEAVE(NST,NE,EXIT1,IPRINT)
!     ***
!     *** FUNCTION  : FREE SPACE IN A NON-ADDRESSIBLE STORAGE
!     *** PARAMETERS: NST   = STORAGE NUMBER
!     ***             NE    = NUMBER OF ELEMENTS TO BE FREED
!     ***             EXIT1 = ABORTED-FREEING EXIT
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     FREE
!     ====
      IF(STO(NST,1).LT.NE) GOTO 10
      STO(NST,1) = STO(NST,1) - NE
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NST,NE
3000  FORMAT(8H LEAVE :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      6H FREES,I5,16H ELEMENTS IN STO,I3)
      RETURN
10    WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NST
3001  FORMAT(1H0,24(1H+),27H WARNING: SUBROUTINE LEAVE ,29(1H+),/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,10H TOO MUCH ,&
      23(1H+),/1X,24(1H+),25H SPACE BEING FREED IN STO,&
      I3,1X,27(1H+)/)
      RETURN 1
      END SUBROUTINE LEAVE


      SUBROUTINE LFIRST(MFA,*)
!     ***
!     ***          CALL LFIRST(MFA, EXIT1)
!     ***
!     *** FUNCTION  : FIND THE FIRST FREE ELEMENT IN A MULTIFAC
!     *** PARAMETERS: MFA   = MULTIFACILITY NUMBER
!     ***             EXIT1 = PLAN-ERROR EXIT
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     FIND A FREE SERVICE ELEMENT
!     ===========================
      I1 = MBV(MFA)
      I2 = I1 + MFAC(MFA,2) - 1
      DO 100  LSE = I1 , I2
      IF(SE(LSE,1).EQ.0) RETURN
100   CONTINUE
!
!     ERROR
!     =====
      WRITE(OUTD,3000) MFA
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE LFIRST ,30(1H+)/1X,&
      24(1H+),32H NO FREE SERVICE ELEMENT IN MFAC,I3,1X,20(1H+)/)
      RETURN 1
      END SUBROUTINE LFIRST



      SUBROUTINE LINK1(NUCHN, ID, *, IPRINT)
!     ***
!     ***          CALL LINK1(NUCHN, ID, &1005, IPRINT)
!     ***
!     *** FUNCTION  : BLOCK KINSMEN AT A USER CHAIN
!     *** PARAMETERS: NUCHN = TYPE-1 USER CHAIN NUMBER
!     ***             ID    = LINK1 CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     TEST FOR FAMILY MEMBERSHIP
!     ==========================
      IF(LFAM.EQ.0) RETURN
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EGATHT+2*FAM1*(NUCHN-1)+2*LFAM-1
!
!     BLOCK
!     =====
      IF(TX(LTX,8).NE.0) GOTO 100
      UCHF(LFAM,NUCHN,1)=UCHF(LFAM,NUCHN,1)+1
      IF(UCHF(LFAM,NUCHN,2).GT.0) GOTO 50
      STATE(K)=0
      STATE(K+1)=1
50    AL(LTX,1)=ID
      AL(LTX,2)=-K
      TX(LTX,8)=T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NUCHN
3000  FORMAT(8H LINK1 :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      29H IS BLOCKED AT USER CHAIN1 NO,I3)
      RETURN 1
!
!     REMOVE
!     ======
100   UCHF(LFAM,NUCHN,1)=UCHF(LFAM,NUCHN,1)-1
      UCHF(LFAM,NUCHN,2)=UCHF(LFAM,NUCHN,2)-1
      IF(UCHF(LFAM,NUCHN,1).EQ.0) UCHF(LFAM,NUCHN,2)=0
      IF(UCHF(LFAM,NUCHN,2).GT.0) GOTO 150
      STATE(K)=0
      STATE(K+1)=1
150   TX(LTX,8)=0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NUCHN
3001  FORMAT(8H LINK1 :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      22H LEAVES USER CHAIN1 NO,I3)
      RETURN
      END SUBROUTINE LINK1




      SUBROUTINE LINK2(NUCHN, ID, *, IPRINT)
!     ***
!     ***          CALL LINK2(NUCHN, ID, &1005, IPRINT)
!     ***
!     *** FUNCTION  : BLOCK TXS IN A USER CHAIN
!     *** PARAMETERS: NUCHN = TYPE-2 USER CHAIN NUMBER
!     ***             ID    = LINK2 CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EUCHF+2*NUCHN-1
!
!     BLOCK
!     =====
      IF(TX(LTX,8).NE.0) GOTO 100
      UCHT(NUCHN,1)=UCHT(NUCHN,1)+1
      IF(UCHT(NUCHN,2).GT.0) GOTO 50
      STATE(K)=0
      STATE(K+1)=1
50    AL(LTX,1)=ID
      AL(LTX,2)=-K
      TX(LTX,8)=T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NUCHN
3000  FORMAT(8H LINK2 :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      29H IS BLOCKED AT USER CHAIN2 NO,I3)
      RETURN 1
!
!     REMOVE
!     ======
100   UCHT(NUCHN,1)=UCHT(NUCHN,1)-1
      UCHT(NUCHN,2)=UCHT(NUCHN,2)-1
      IF(UCHT(NUCHN,1).EQ.0) UCHT(NUCHN,2)=0
      IF(UCHT(NUCHN,2).GT.0) GOTO 150
      STATE(K)=0
      STATE(K+1)=1
150   TX(LTX,8)=0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NUCHN
3001  FORMAT(8H LINK2 :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      22H LEAVES USER CHAIN2 NO,I3)
      RETURN
      END SUBROUTINE LINK2


      SUBROUTINE LOGNOR(MEAN,SIGMA,MIN,MAX,RNUM,RANDOM)
!     ***
!     ***          CALL LOGNOR(MEAN,SIGMA,MIN,MAX,RNUM,RANDOM)
!     ***
!     *** FUNCTION  : GENERATE A LOGNORMAL RANDOM SEQUENCE
!     *** PARAMETERS: MEAN   = MEAN VALUE
!     ***             SIGMA  = STANDART DEVIATION
!     ***             MIN    = INTERVAL'S LOWER BOUND
!     ***             MAX    = INTERVAL'S UPPER BOUND
!     ***             RNUM   = GENERATOR'S IDENTIFIER
!     ***             RANDOM = RESULTING RANDOM SEQUENCE
!     ***
      INTEGER RNUM
      REAL MAX,MIN,MEANX,MEAN,RANDOM,SIGMA
	  REAL SIGMX, SIGMX2, R, V
!
!     GENERATE A LOGNORMALLY DISTRIBUTED RANDOM NUMBER
!     ================================================
      SIGMX2 = ALOG(SIGMA**2/MEAN**2+1.)
      MEANX  = ALOG(MEAN)-0.5*SIGMX2
      SIGMX  = SQRT(SIGMX2)
100   R=RN(RNUM)
      IF(R.EQ.0.) GOTO 100
      V=(-2.0*ALOG(R))**0.5*COS(6.2832*RN(RNUM))
      RANDOM=EXP(V*SIGMX+MEANX)
!
!     INTERVAL CHECK
!     ==============
      IF(RANDOM.LT.MIN.OR.RANDOM.GT.MAX) GOTO 100
      RETURN
      END SUBROUTINE LOGNOR


      SUBROUTINE MCLEAR(MFA,*,*,IPRINT)
!     ***
!     ***          CALL MCLEAR(MFA ,EXIT1, &1006, IPRINT)
!     ***
!     *** FUNCTION  : FREE A SERVICE ELEMENTE IN A
!     ***             MULTIFACILITY
!     *** PARAMETERS: MFA   = MULTIFACILITY NUMBER
!     ***             EXIT1 = PREEMPTION EXIT
!     ***
      IMPLICIT INTEGER(A - Z)
!
!     FIND THE TX'S SERVICE ELEMENT
!     =============================
      I1 = MBV(MFA)
      I2 = I1 + MFAC(MFA,2) - 1
      DO 50  LSE = I1 , I2
      IF(IABS(SE(LSE,1)).EQ.LTX) GOTO 100
50    CONTINUE
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),MFA
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE MCLEAR ,30(1H+)/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,20H DOESN'T OCCUPY MFAC,&
      I3,1X,9(1H+)/)
      RETURN 2
!
!     FREE
!     ====
100   DO 150  I = 1 , 3
150   SE(LSE,I) = 0
      MFAC(MFA,1) = MFAC(MFA,1) - 1
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EFAC + MFA
!
!     TERMINATE THE OCCUPATION
!     ========================
      IF(TX(LTX,6).EQ.0) GOTO 200
      AL(LTX,1) = TX(LTX,5)
      AL(LTX,2) = - K
      TX(LTX,8) = T
200   TX(LTX,5) = 0
      IF(TX(LTX,6).EQ.0) TX(LTX,7) = 0
      IF(IPRINT.EQ.0) GOTO 5001
      NSE = LSE - MBV(MFA) + 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NSE,MFA
3001  FORMAT(8H MCLEAR:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      23H LEAVES SERVICE ELEMENT,I3,8H IN MFAC,I3)
!
!     START THE LOCKED  TXS
!     =====================
5001  STATE(K) = 1
      CALL UNLOCK(K,IPRINT)
      IF(TX(LTX,6).EQ.0) RETURN
      RETURN 1
      END SUBROUTINE MCLEAR



      SUBROUTINE MKNOCK(KT,MFA,IDN,*,*,IPRINT)
!     ***
!     ***          CALL MKNOCK(KT,MFA,IDN,&1005,&1006,IPRINT)
!     ***
!     *** FUNCTION  : KNOCK DOWN A SERVICE ELEMENT
!     ***             IN A MULTIFACILITY
!     *** PARAMETERS: KT  = KNOCKDOWN TIME
!     ***             MFA = MULTIFACILITY NUMBER
!     ***             IDN = TARGET
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     FIND THE TX'S SERVICE ELEMENT
!     =============================
      I1=MBV(MFA)
      I2=I1+MFAC(MFA,2)-1
      DO 100 LSE=I1,I2
      IF(IABS(SE(LSE,1)).EQ.LTX) GOTO 200
100   CONTINUE
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),MFA
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE MKNOCK ,30(1H+)/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,20H DOESN'T OCCUPY MFAC,&
      I3,1X,9(1H+)/)
      RETURN 2
!
!     KNOCK DOWN THE SERVICE ELEMENT
!     ==============================
200   SE(LSE,1)=-LTX
      SE(LSE,3)=3
      AL(LTX,1)=IDN
      AL(LTX,2)=T+KT
      IF(IPRINT.EQ.0) RETURN 1
      NSE = LSE - MBV(MFA)+1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NSE,MFA,AL(LTX,2)
3001  FORMAT(8H MKNOCK:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      16H SERVICE ELEMENT,I3,8H IN MFAC,I3,21H WILL BE KNOCKED DOWN,&
      9H UNTIL T=,I7)
      RETURN 1
      END SUBROUTINE MKNOCK


      SUBROUTINE MPREEM(MFA,ID,REP,*,*,IPRINT)
!     ***
!     ***          CALL MPREEM(MFA,ID,REP,&1005,&1006,IPRINT)
!     ***
!     *** FUNCTION  : PREFERRED ACQUISITION OF A SERVICE ELEMENT
!     ***             IN A MULTIFACILITY WITH PRIORITY COMPARISON
!     *** PARAMETERS: MFA  = MULTIFACILITY NUMBER
!     ***             ID   = MPREEM CALL'S STATEMENT NUMBER
!     ***             REP  = REPEAT CODE
!     ***                    = 0: THE TX NEED NOT BE BOUND TO ITS
!     ***                         INITIAL SERVICE  ELEMENT
!     ***                    = 1: TX MUST BE BOUND TO ITS
!     ***                         INITIAL SERVICE ELEMENT
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EFAC+MFA
!
!     BLOCKING DICISION
!     =================
      K1 = 1
      IF(MFAC(MFA,1).EQ.MFAC(MFA,2)) GOTO 300
      K1 = 0
      IF(OK.EQ.0) GOTO 300
      OK = 0
!
!     FIND A FREE SERVICE ELEMENT
!     ===========================
      IF(TX(LTX,7).EQ.0) GOTO 100
      LSE=TX(LTX,7)
      IF(SE(LSE,1).NE.0) GOTO 200
      GOTO 150
100   CALL PLANI(MFA,*400)
      IF(LSE.EQ.0) GOTO 200
!
!     ACQUIRE
!     =======
150   MFAC(MFA,1)=MFAC(MFA,1)+1
      IF(MFAC(MFA,1).EQ.MFAC(MFA,2)) STATE(K)=0
      SE(LSE,1)=LTX
      TX(LTX,5)=ID
      IF(REP.EQ.1) TX(LTX,7)=LSE
      TX(LTX,8) = 0
      IF(IPRINT.EQ.0) RETURN
      NSE=LSE-MBV(MFA)+1
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NSE,MFA
3000  FORMAT(8H MPREEM:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      25H ACQUIRES SERVICE ELEMENT,I3,8H IN MFAC,I3)
      RETURN
!
!     LOCK
!     ====
200   AL(LTX,1) = ID
      AL(LTX,2) = -KEND - K
      IF(TX(LTX,8).EQ.0) TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),MFA
3001  FORMAT(8H MPREEM:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      18H IS LOCKED AT MFAC,I3)
      RETURN 1
!
!     BLOCK
!     =====
300   AL(LTX,1) = ID
      AL(LTX,2) = -K
      TX(LTX,8) = T
      IF(IPRINT.EQ.0) GOTO 5002
      WRITE(OUTD,3002) T,TX(LTX,1),TX(LTX,2),MFA
3002  FORMAT(8H MPREEM:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      19H IS BLOCKED AT MFAC,I3)
5002  IF(K1.EQ.0) RETURN 1
!
!     FIND A PREEMPTION CANDIDATE
!     ===========================
      CALL PLANO(MFA,*400)
      IF(LSE.EQ.0.OR.SE(LSE,1).LT.0.AND.SE(LSE,3).EQ.2) RETURN 1
      IF(SE(LSE,3).EQ.3) RETURN 1
!
!     CALL THE POLICY
!     ===============
      POLVEC(1) = IABS(SE(LSE,1))
      POLC = 1
      DO 350  I = 1 , LAL
      IF(AL(I,2).NE.-K) GOTO 350
      POLC = POLC + 1
      POLVEC(POLC) = I
350   CONTINUE
      CALL POLICY(K)
      IF(LTX.EQ.IABS(SE(LSE,1))) RETURN 1
!
!     PREEMPT
!     =======
      LTX = IABS(SE(LSE,1))
      IF(IPRINT.EQ.0) GOTO 5003
      NSE = LSE - MBV(MFA) + 1
      WRITE(OUTD,3003) T,TX(LTX,1),TX(LTX,2),NSE,MFA
3003  FORMAT(8H MPREEM:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      34H IS PREEMPTED FROM SERVICE ELEMENT,I3,8H IN MFAC,I3)
5003  SE(LSE,2) = 1
      IF(SE(LSE,3).EQ.1) RETURN 1
      TX(LTX,6) = AL(LTX,2) - T
      AL(LTX,2) = T
      RETURN 1
!
!     INCORRECT MULTIFACILITY OCCUPATION
!     ==================================
400   RETURN 2
      END SUBROUTINE MPREEM


      SUBROUTINE MSEIZE(MFA,ID,REP,*,*,IPRINT)
!     ***
!     ***          CALL MSEIZE(MFA,ID,REP,&1005,&1006,IPRINT)
!     ***
!     *** FUNCTION  : ACQUIRE A SERVICE ELEMENT IN A
!     ***             MULTIFACILITY
!     *** PARAMETERS: MFA = MULTIFACILITY
!     ***             ID = MSEIZE CALL'S STATEMENT NUMBER
!     ***             REP = REPEAT CODE
!     ***                 = 0:  THE NEED NOT BE BOUND TO IT
!     ***                       INITIAL SERVICE ELEMENT
!     ***                 = 1:  THE TX MUST BE BOUND TO ITS
!     ***                       INITIAL SERVICE ELEMENT
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EFAC + MFA
!
!     BLOCKING DECISION
!     =================
      IF(OK.EQ.0) GOTO 300
      OK = 0
      IF(MFAC(MFA,1).EQ.MFAC(MFA,2)) GOTO 300
!
!     ASSIGN A SERVICE ELEMENT ACCORDING TO THE PLAN
!     ==============================================
      IF(TX(LTX,7).EQ.0) GOTO 100
      LSE = TX(LTX,7)
      IF(SE(LSE,1).NE.0) GOTO 200
      GOTO 150
100   CALL PLANI(MFA,*400)
      IF(LSE.EQ.0) GOTO 200
!
!     ACQUIRE
!     =======
150   MFAC(MFA,1) = MFAC(MFA,1) + 1
      IF(MFAC(MFA,1).EQ.MFAC(MFA,2)) STATE(K) = 0
      SE(LSE,1) = LTX
      TX(LTX,5) = ID
      IF(REP.EQ.1) TX(LTX,7) = LSE
      TX(LTX,8) = 0
      IF(IPRINT.EQ.0) RETURN
      NSE = LSE - MBV(MFA) + 1
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NSE,MFA
3000  FORMAT(8H MSEIZE:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      25H ACQUIRES SERVICE ELEMENT,I3,8H IN MFAC,I3)
      RETURN
!
!     LOCK
!     ====
200   AL(LTX,1) = ID
      AL(LTX,2) = -KEND-K
      IF(TX(LTX,8).EQ.0) TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),MFA
3001  FORMAT(8H MSEIZE:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      18H IS LOCKED AT MFAC,I3)
      RETURN 1
!
!     BLOCK
!     =====
300   AL(LTX,1) = ID
      AL(LTX,2) = -K
      TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE (OUTD,3002) T,TX(LTX,1),TX(LTX,2),MFA
3002  FORMAT(8H MSEIZE:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      19H IS BLOCKED AT MFAC,I3)
      RETURN 1
!
!     INCORRECT MULTIFACILITY OCCUPATION
!     ==================================
400   RETURN 2
      END SUBROUTINE MSEIZE


      SUBROUTINE MSETUP(ST,MFA,IDN,*,*,IPRINT)
!     ***
!     ***          CALL MSETUP (ST,MFA,IDN,&1005,&1006,IPRINT)
!     ***
!     *** FUNCTION  : SET UP A SERVICE ELEMENT IN A MULTIFACILITI
!     *** PARAMETERS: ST  = SETUP TIME
!     ***             MFA = MULTIFACILITY NUMBER
!     ***             IDN = TARGET
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     FIND THE TX'S SERVICE ELEMENT
!     =============================
      I1 = MBV(MFA)
      I2 = I1 + MFAC(MFA,2) - 1
      DO 100  LSE = I1 , I2
      IF(IABS(SE(LSE,1)).EQ.LTX) GOTO 200
100   CONTINUE
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),MFA
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE MSETUP ,30(1H+)/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,20H DOESN'T OCCUPY MFAC,&
      I3,1X,9(1H+)/)
      RETURN 2
!
!     SET UP THE SERVICE ELEMENT
!     ==========================
200   SE(LSE,1) = - LTX
      SE(LSE,3) = 1
      AL(LTX,1) = IDN
      AL(LTX,2) = T + ST
      IF(IPRINT.EQ.0) RETURN 1
      NSE = LSE - MBV(MFA) + 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NSE,MFA,AL(LTX,2)
3001  FORMAT(8H MSETUP:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      16H SERVICE ELEMENT,I3,8H IN MFAC,I3,21H WILL BE SET UP UNTIL,&
      3H T=,I7)
      RETURN 1
      END SUBROUTINE MSETUP



      subroutine mwork(wt,mfa,id,iex,*,*,iprint)
!     ***
!     ***          call mwork(wt, mfa, id, iex, &1005, &1006, iprint)
!     ***
!     *** function   :  serve a tx at a service element in a mult
!     *** parameters :  wt   = work time
!     ***               mfa  = mulitfacility number
!     ***               id   = mwork call's statement number
!     ***               iex  = preemption code
!     ***                    = 0: the tx may be preempted
!     ***                    = 1: the tx may not be preempted
!     ***
      implicit integer (a - z)
!
!     find the tx's service element
!     =============================
      i1 = mbv(mfa)
      i2 = i1 + mfac(mfa,2) - 1
      do 100 lse = i1 , i2
      if(iabs(se(lse,1)).eq.ltx) goto 200
100   continue
      write(outd,3000) t,tx(ltx,1),tx(ltx,2),mfa
3000  format(1h0,24(1h+),25h error: subroutine mwork ,31(1h+)/1x&
      24(1h+),3h t=,i7,2x,2htx,i5,1h,,i3,20h doesn't occupy mfac,&
      i3,1x,9(1h+)/)
      return 2
!
!     serving decision
!     ================
200   if(se(lse,2).eq.0) goto 250
      if(se(lse,3).eq.0) tx(ltx,6) = wt
      return
250   if(se(lse,3).eq.2) return
!
!     serve
!     =====
      if(iex.eq.0) se(lse,1) = ltx
      if(iex.eq.1) se(lse,1) = - ltx
      se(lse,3) = 2
      if(tx(ltx,6).ne.0) goto 300
      al(ltx,1) = id
      al(ltx,2) = t + wt
      if(iprint.eq.0) return 1
      write(outd,3001) t,tx(ltx,1),tx(ltx,2),al(ltx,2)
3001  format(8h mwork :,3x,2ht=,i7,2x,2htx,i5,1h,,i3,2x,&
      24h will be served until t=,i7)
      return 1
!
!     serve reacquiring txs
!     =====================
300   al(ltx,1) = id
      al(ltx,2) = t + tx(ltx,6)
      tx(ltx,6) = 0
      if(iprint.eq.0) return 1
      write(outd,3002) t,tx(ltx,1),tx(ltx,2),al(ltx,2)
3002  format(8h mwork :,3x,2ht=,i7,2x,2htx,i5,1h,,i3,2x,&
      25h is again served until t=,i7)
      return 1
      end subroutine mwork



      subroutine pfifo
!     ***
!     ***          call pfifo
!     ***
!     *** function :   select the tx with the highest priority
!     ***
      implicit integer (a - z)
!
!     initialize the search
!     =====================
      ipr = 1
      ltx = polvec(1)
      if(polc.eq.1) return
!
!     compare
!     =======
      do 100 i = 2 , polc
      ltx1 = polvec(i)
      if(tx(ltx,4)-tx(ltx1,4)) 110, 120, 100
110   ltx = ltx1
      ipr = 1
      polvec(1) = ltx
      goto 100
120   ipr = ipr + 1
      polvec(ipr) = ltx1
100   continue
!
!     call fifo
!     =========
      polc = ipr
      if(ipr.le.1) return
      call fifo
      return
      end subroutine pfifo




      subroutine plani(mfa,*)
!     ***
!     ***          call plani(mfa, exit1)
!     ***
!     *** function  :  find the plan responsible for service-element
!     ***              assignment at this multifacility
!     *** parameters:  mfa   = multifacility number
!     ***              exit1 = plan-error exit
!     ***
      implicit integer(a-z)
!
!     find the plan
!     =============
      if(plama(mfa,1).ne.0) goto 100
      write(outd,3000) mfa
3000  format(1h0,24(1h+),25h error: subroutine plani ,31(1h+)/1x,&
      24(1h+),9h for mfac,i3,9h no plan ,35(1h+)/)
50    return 1
100   iaddr = plama(mfa,1)
      goto (1,2,3,4,5) , iaddr
!
!     call the plani subroutine
!     =========================
1     call lfirst(mfa,*50)
      return
2     call plani2(mfa)
      return
3     call plani3(mfa)
      return
4     call plani4(mfa)
      return
5     call plani5(mfa)
      return
      end subroutine plani
      subroutine plani2(mfa)
	  implicit integer (a-z)
      return
      end subroutine plani2
      subroutine plani3(mfa)
	  implicit integer (a-z)
      return
      end subroutine plani3
      subroutine plani4(mfa)
	  implicit integer (a-z)
      return
      end subroutine plani4
      subroutine plani5(mfa)
	  implicit integer (a-z)
      return 
      end subroutine plani5



      subroutine plano(mfa,*)
!     ***
!     ***          call plano(mfa, exit1)
!     ***
!     *** function  :  find the plan responsible for deciding where
!     ***              service element to preempt
!     *** parameters:  mfa   = multifacility number
!     ***              exit1 = plan-error exit
!     ***
      implicit integer(a-z)
!
!     find the plan
!     =============
      if(plama(mfa,2).ne.0) goto 100
      write(outd,3000) mfa
3000  format(1h0,24(1h+),25H ERROR: SUBROUTINER PLANO ,31(1h+)/1x,&
      24(1h+),9h FOR MFAC,i3,20H NO PREEMPTION PLAN ,35(1h+)/)
      return 1
100   iaddr = plama(mfa,2)
      goto (1,2,3,4,5) , iaddr
!
!     call the plano subroutine
!     =========================
1     call prior(mfa)
      return
2     call plano2(mfa)
      return
3     call plano3(mfa)
      return
4     call plano4(mfa)
      return
5     call plano5(mfa)
      return 
      end subroutine plano
      subroutine plano2(mfa)
      implicit integer(a-z)
      return
      end subroutine plano2
      subroutine plano3(mfa)
      implicit integer(a-z)
      return
      end subroutine plano3
      subroutine plano4(mfa)
      implicit integer(a-z)
      return
      end subroutine plano4
      subroutine plano5(mfa)
      implicit integer(a-z)
      return
      end subroutine plano5


      SUBROUTINE POLICY(K)
!     ***
!     ***          CALL POLICY(K)
!     ***
!     *** FUNCTION  : DETERMINE WHICH POLICY PROCESSES THE TXS WITH
!     ***             AT THE STATION WHOSE STATION NUMBER IS K
!     *** PARAMETER : K = STATION NUMBER
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     FIND THE POLICY
!     ===============
      DO  100 I = 1 , POL1
      IF(POL(I,1).NE.K) GOTO 100
      IADDR = POL(I,2)
      GOTO 120
100   CONTINUE
      IADDR = 1
120   CONTINUE
      GOTO (1,2,3,4,5),IADDR
!
!     CALL THE POLICY-IMPLEMENTING SUBROUTINE
!     =======================================
1     CALL PFIFO
      RETURN
2     CALL FIFO
      RETURN
3     CALL POLI3
      RETURN
4     CALL POLI4
      RETURN
5     CALL POLI5
      RETURN
      END SUBROUTINE POLICY

      SUBROUTINE POLI3
      IMPLICIT INTEGER (A - Z)
      RETURN
      END SUBROUTINE POLI3
      SUBROUTINE POLI4
      IMPLICIT INTEGER (A - Z)
      RETURN
      END SUBROUTINE POLI4
      SUBROUTINE POLI5
      IMPLICIT INTEGER (A - Z)
      RETURN
      END SUBROUTINE POLI5



      SUBROUTINE PREEMP(NFA,ID,*,IPRINT)
!     ***
!     ***         CALL PREEPM(NFA,ID,&1005,IPRINT)
!     ***
!     *** FUNCTION  : PREFERRED ACQUISITION OF A A FACILITY
!     ***             WITH PRIORITY COMPARISON
!     *** PARAMETERS: NFA = FACILITY NUMBER
!     ***             ID  = PREEMP CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     BLOCKING DECISION
!     =================
      K1 = 1
      IF(FAC(NFA,1).NE.0) GOTO 200
      K1 = 0
      IF(OK.EQ.0) GOTO 200
      OK = 0
!
!     ACQUIRE
!     =======
      FAC(NFA,1) = LTX
      STATE(NFA) = 0
      TX(LTX,8) = 0
      TX(LTX,5) = ID
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NFA
3000  FORMAT(8H PREEMP:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      13H ACQUIRES FAC,I3)
      RETURN
!
!     BLOCK
!     =====
200   AL(LTX,1) = ID
      AL(LTX,2) = - NFA
      TX(LTX,8) = T
      IF(IPRINT.EQ.0) GOTO 5001
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NFA
3001  FORMAT(8H PREEMP:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      18H IS BLOCKED AT FAC,I3)
5001  IF(K1.EQ.0) RETURN 1
!
!     SUPRESS THE PREEMPTION
!     ======================
300   IF(FAC(NFA,1).LT.0.AND.FAC(NFA,3).EQ.2) RETURN 1
      IF(FAC(NFA,3).EQ.3.OR.FAC(NFA,2).EQ.1) RETURN 1
!
!     CALL THE POLICY
!     ===============
      POLVEC(1)=IABS(FAC(NFA,1))
      POLC = 1
      DO 350 I=1,LAL
      IF(AL(I,2).NE.-NFA) GOTO 350
      POLC=POLC+1
      POLVEC(POLC)=I
350   CONTINUE
      CALL POLICY(NFA)
      IF(LTX.EQ.IABS(FAC(NFA,1))) RETURN 1
!
!     PREEMPT
!     =======
      LTX = IABS(FAC(NFA,1))
      IF(IPRINT.EQ.0) GOTO 5002
      WRITE(OUTD,3002) T,TX(LTX,1),TX(LTX,2),NFA
3002  FORMAT(8H PREEMP:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      22H IS PREEMTED FROM FAC,I3)
5002  FAC(NFA,2) = 1
      IF(FAC(NFA,3).EQ.1) RETURN 1
      TX(LTX,6) = AL(LTX,2) - T
      AL(LTX,2) = T
      RETURN 1
      END SUBROUTINE PREEMP

       subroutine prior(mfa)
!      ***
!      ***         CALL PRIOR(MFA)
!      ***
!      *** FUNCTION  : FIND THE TXS OCCUPYING THE MFAC WHOSE
!      ***             PRIORITY IS LOWEST
!      *** PARAMETER : MFA = MULTIFACILITY NUMBER
!      ***
       IMPLICIT INTEGER (A-Z)
!
!     INITIALIZE THE SEARCH
!     =====================
      LSE = 0
      PR = TX(LTX,4)
      I1 = MBV(MFA)
      I2 = I1 + MFAC(MFA,2) - 1
!
!     SEARCH OUT THE SERVICE ELEMENT
!     ==============================
      DO 100  I = I1 , I2
      IF(SE(I,1).LT.0.AND.SE(I,3).EQ.2) GOTO 100
      IF(SE(I,2).EQ.1) GOTO 100
      LTX1 = IABS(SE(I,1))
      IF(TX(LTX1,4).GE.PR) GOTO 100
      PR = TX(LTX1,4)
      LSE = I
100   CONTINUE
      RETURN
      END subroutine prior



      SUBROUTINE REPORT(IBIN,IFAC,IMFAC,ISTO,IEL,IAL,ITX,IFAM)
!     ***
!     ***          CALL REPORT(IBIN,IFAC,IMFAC,ISTO,IEL,IAL,ITX,
!     ***                      IFAM)
!     ***
!     *** FUNCTION  : PRINT THE DATA AREAS FOR BINS, FACILITIES,
!     ***             MULTIFACILITIES, STORAGES, AND FAMILIES,
!     ***             THE EVENT LIST, THE ACTIVATION LIST AND
!     ***             TRANSACTION MATRIX
!     *** PARAMETERS: IBIN  = NUMBER OF BINS
!     ***             IFAC  = NUMBER OF FACILITIES
!     ***             IMFAC = NUMBER OF MULIFACILITIES
!     ***             ISTO  = NUMBER OF STORAGES
!     ***             IEL   = NUMBER OF LINES IN THE EVENT LIST
!     ***             IAL   = NUMBER OF LINES IN TH ACTIVATION LIST
!     ***             ITX   = NUMBER OF LINES IN THE TX MATRIX
!     ***             IFAM  = NUMBER OF LINES IN THE FAMILY MATRIX
      IMPLICIT INTEGER(A-Z)
          WRITE(OUTD,3000) T,RT
3000      FORMAT(10X,4HT = ,I7,2X,5HRT = ,I7)
!
!     PRINT THE BIN MATRIX
!     ====================
      IF(IBIN.EQ.0) GOTO 200
      WRITE(OUTD,3001)
3001  FORMAT(/9X,21H CONTENTS OF THE BINS)
      WRITE(OUTD,3002)
3002  FORMAT(/10X,3HNUM,7X,4HSIZE,8X,3HMAX,6X,5HARRIV,5X,6HDEPART,&
      6X,5HFLEET,6X,5HATIME,6X,5HBTIME,7X,4HLCHG,6X,5HMTIME,6X,&
      5HWSIZE)
      DO 100 I = 1 , IBIN
      IF(BIN(I,3).EQ.0) GOTO 100
      WRITE(OUTD,3003) I,(BIN(I,J1),J1=1,8),BINSTA(I,1),BINSTA(I,2)
3003  FORMAT(2X,9(1X,I10),2(1X,F10.5))
100   CONTINUE
!
!     PRINT THE FAC MATRIX
!     ====================
200   IF(IFAC.EQ.0) GOTO 300
      WRITE(OUTD,3004)
3004  FORMAT(///9X,27H CONTENTS OF THE FACILITIES)
      WRITE(OUTD,3005)
3005  FORMAT(/10X,24HNUM   CONT   PRE   PHASE)
      DO 250  I = 1 , IFAC
      IF(FAC(I,1).EQ.0.AND.FAC(I,2).EQ.0) GOTO 250
      WRITE(OUTD,3006) I,(FAC(I,J),J=1,3)
3006  FORMAT(11X,I2,I7,I6,I8)
250   CONTINUE
!
!     PRINT THE MFAC MATRIX
!     =====================
300   IF(IMFAC.EQ.0) GOTO 400
      WRITE(OUTD,3007)
3007  FORMAT(///9X,32H CONTENTS OF THE MULTIFACILITIES)
      WRITE(OUTD,3008)
3008  FORMAT(/10X,16HNUM   CONT   CAP)
      DO 350  I = 1 , IMFAC
      IF(MFAC(I,2).EQ.0) GOTO 350
      WRITE(OUTD,3006) I,MFAC(I,1),MFAC(I,2)
350   CONTINUE
!
!     PRINT THE STO MATRIX
!     ====================
400   IF(ISTO.EQ.0) GOTO 500
      WRITE(OUTD,3009)
3009  FORMAT(///9X,25H CONTENTS OF THE STORAGES)
      WRITE(OUTD,3010)
3010  FORMAT(/10X,16HNUM   CONT   CAP)
      DO 450 I=1,ISTO
      IF(STO(I,2).EQ.0) GOTO 450
      WRITE(OUTD,3006) I,STO(I,1),STO(I,2)
450   CONTINUE
!
!     PRINT THE EVENT LIST
!     ====================
500   IF(IEL.EQ.0) GOTO 600
      WRITE(OUTD,3011)
3011  FORMAT(///9X,11H EVENT LIST)
      WRITE(OUTD,3012)
3012  FORMAT(/10X,24HLINE     TARGET     TIME)
      DO 550 I=1,IEL
      IF(EL(I,1).EQ.0) GOTO 550
      WRITE(OUTD,3013) I,(EL(I,J1),J1=1,2)
3013  FORMAT(4X,3I10)
550   CONTINUE
!
!     PRINT THE ACTIVATION LIST
!     =========================
600   IF(IAL.EQ.0) GOTO 700
      WRITE(OUTD,3014)
3014  FORMAT(///9X,16H ACTIVATION LIST)
      WRITE(OUTD,3015)
3015  FORMAT(/10X,24HLINE    TARGET     STATE)
      DO 650 I=1,IAL
      IF(AL(I,1).EQ.0) GOTO 650
      WRITE(OUTD,3013) I,(AL(I,J1),J1=1,2)
650   CONTINUE
!
!     PRINT THE TX MATRIX
!     ===================
700   IF(ITX.EQ.0) GOTO 800
      WRITE(OUTD,3016)
3016  FORMAT(///9X,19H TRANSACTION MATRIX)
      WRITE(OUTD,3017)
3017  FORMAT(/10X,3HLTX,6X,5HTX NO,5X,6HDUP NO,8X,3HGEN,7X,&
      4HPRIO,6X,5HPTARG,7X,4HREST,5X,6HRET EL,6X,5HBLOCK)
      DO 710 I=1,ITX
      IF(TX(I,1).EQ.0) GOTO 710
      WRITE(OUTD,3018) I,(TX(I,J1),J1=1,8)
3018  FORMAT(2X,9(1X,I10))
710   CONTINUE
      WRITE(OUTD,3019)
3019  FORMAT(/10X,3HLTX,6X,5HBIN 1,6X,5HBIN 2,6X,5HBIN 3,6X,5HBIN 4,&
      6X,5HBIN 5,4X,7HENTRY 1,4X,7HENTRY 2,4X,7HENTRY 3,4X,&
      7HENTRY 4,4X,7HENTRY 5)
      DO 720 I=1,ITX
      IF(TX(I,1).EQ.0) GOTO 720
      WRITE(OUTD,3020) I,(TX(I,J1),J1=9,18)
3020  FORMAT(2X,11(1X,I10))
720   CONTINUE
!
!     PRINT THE FAM MATRIX
!     ====================
800   IF(IFAM.EQ.0) RETURN
      WRITE(OUTD,3021)
3021  FORMAT(///9X,9H FAMILIES)
      WRITE(OUTD,3022)
3022  FORMAT(/10X,24HLINE      FAMNO     SIZE)
      DO 850 I=1,IFAM
      IF(FAM(I,1).EQ.0) GOTO 850
      WRITE(OUTD,3013) I,FAM(I,1),FAM(I,2)
850   CONTINUE
      RETURN
      END SUBROUTINE REPORT


      SUBROUTINE RESET
!     ***
!     ***          CALL RESET
!     ***
!     *** FUNCTION  : CLEAR GPSS-F DATA AREAS
!     ***
      IMPLICIT INTEGER(A-Z)
 
!     CLEAR THE SIMULATIN CLOCKS
!     ==========================
      T = 0
      RT = 0
!
!     SET THE OK AND IT MECHANISMS
!     ============================
      OK = 0
      IT = 0
!
!     RESET THE NTXC COUNTER
!     ======================
      NTXC = 0
!
!     RESET THE LIST-END POINTERS
!     ===========================
      LEL = 0
      LAL = 0
!
!     CLEAR THE DATA AREAS
!     ====================
      DO 111 I=1,SRC1
      DO 111 J=1,2
111   SRC(I,J)=0
      DO 222 I=1,EL1
      DO 222 J=1,2
222   EL(I,J)=0
      DO 444 I=1,TX1
      DO 333 J=1,2
333   AL(I,J)=0
      DO 444 J=1,TX2
444   TX(I,J)=0
      DO 555 I=1,FAC1
      DO 555 J=1,3
555   FAC(I,J)=0
      DO 666 I=1,MFAC1
      MBV(I)=0
      DO 666 J=1,2
      MFAC(I,J)=0
666   PLAMA(I,J)=0
      DO 777 I=1,SE1
      DO 777 J=1,3
777   SE(I,J)=0
      DO 888 I=1,STO1
      SBV(I)=0
      DO 888 J=1,2
      STO(I,J)=0
888   STRAMA(I,J)=0
      DO 999 I=1,SM1
      DO 999 J=1,2
999   SM(I,J)=0
      DO 2222 I=1,BIN1
      BIN(I,8)=1
      DO 1111 J=1,7
1111  BIN(I,J)=0
      DO 2222 J=1,2
2222  BINSTA(I,J)=0
      DO 3333 I=1,STATE1
3333  STATE(I)=1
      DO 7777 I=1,FAM1
      DO 4444 J=1,3
4444  FAM(I,J)=0
      DO 5555 J=1,ASM1
5555  ASM(I,J)=0
      DO 6666 J=1,GATHF1
6666  GATHF(I,J)=0
      DO 7777 J=1,UCHF1
      DO 7777 K=1,2
7777  UCHF(I,J,K)=0
      DO 8888 I=1,UCHT1
      DO 8888 J=1,2
8888  UCHT(I,J)=0
      DO 9999 I=1,GATHT1
9999  GATHT(I)=0
      DO 11111 I=1,POL1
      DO 11111 J=1,2
11111 POL(I,J)=0
      RETURN
      END SUBROUTINE RESET


      FUNCTION RN(RNUM)
!     ***
!     *** FUNCTION  : GENERATE A UNIFORM RANDOM NUMBER ON
!     ***             (O,1)
!     *** PARAMETER : RNUM = GENERATOR'S IDENTIFIER
!     ***
      IMPLICIT NONE
      INTEGER RNUM
	  REAL RN
    !  REAL*8 DRN,DFACT,DMODUL,DCONST
!
!     GENERATE A UNIFORM RANDOM SEQUENCE
!     ==================================
      DRN(RNUM)=DMOD(DFACT(RNUM)*DRN(RNUM)+DCONST(RNUM),DMODUL)
      RN=DRN(RNUM)/(DMODUL-1.)
      RETURN
      END FUNCTION RN




      subroutine save
!     ***
!     ***          call save
!     ***
!     *** function  : save all data areas on a file with
!     ***             logical device number 'savo'
!     ***
      implicit integer (a-z)
      rewind savo
      write(savo,1) lal,lel
1     format(10i10)
      k=5
      do 2 i=1 , lel , 5
      if(k.gt.lel) k = lel
      write(savo,1) ((el(j,m),m=1,2),j=i,k)
2     k = k +5
      k = 5
      do 3 i = 1 , lal , 5
      if(k.gt.lal) k = lal
      write(savo,1) ((al(j,m),m = 1 , 2 ) , j = i , k )
3     k = k + 5
      do 4 i=1,lal
      k = 10
      do 4 j = 1,tx2,10
      if(k.gt.tx2) k = tx2
      write(savo,1) (tx(i,m),m=j,k)
4     k = k + 10
      k = 10
      do 5 i = 1 , state1 , 10
      if(k.gt.state1) k = state1
      write(savo,1) (state(m),m=i,k)
5     k = k + 10
      write(savo,1) t
      write(savo,1) rt
      do 6 i = 1 , fam1
      k = 10
      do 6 j = 1 , asm1 , 10
      if(k.gt.asm1) k = asm1
      write(savo,1) (asm(i,m),m=j,k)
6     k = k + 10
      k = 3
      do 7 i = 1 , fam1 , 3
      if(k.gt.fam1) k = fam1
      write(savo,1) ((fam(j,m),m=1,3),j=i,k)
7     k = k + 3
      do 8 i = 1 , fam1
      k = 10
      do 8 j = 1 , gathf1 , 10
      if(k.gt.gathf1) k = gathf1
      write(savo,1) (gathf(i,m),m=j,k)
8     k = k + 10
      do 9 i = 1 , fam1
      k = 5
      do 9 j = 1 , uchf1 , 5
      if(k.gt.uchf1) k = uchf1
      write(savo,1) ((uchf(i,m,l),l=1,2),m=j,k)
9     k = k + 5
      k = 3
      do 10 i = 1 , fac1 , 3
      if(k.gt.fac1) k = fac1
      write(savo,1) ((fac(j,m),m=1,3),j=i,k)
10    k = k + 3
      k = 5
      do 11  i = 1 , mfac1 , 5
      if(k.gt.mfac1) k = mfac1
      write(savo,1) ((mfac(j,m),m=1,2),j=i,k)
11    k = k + 5
      k = 10
      do 12  i = 1 , mfac1 , 10
      if(k.gt.mfac1) k = mfac1
      write(savo,1) (mbv(m),m=i,k)
12    k = k + 10
      k = 5
      do 13  i = 1 , mfac1 , 5
      if(k.gt.mfac1) k = mfac1
      write(savo,1) ((plama(j,m),m=1,2),j=i,k)
13    k = k + 5
      k = 3
      do 14  i = 1 , se1 , 3
      if(k.gt.se1) k = se1
      write(savo,1) ((se(j,m),m=1,3),j=i,k)
14    k = k + 3
      k = 5
      do 15  i = 1 , sto1 , 5
      if(k.gt.sto1) k = sto1
      write(savo,1) ((sto(j,m),m=1,2),j=i,k)
15    k = k + 5
      k = 10
      do 16 i = 1 , sto1 , 10
      if(k.gt.sto1) k = sto1
      write(savo,1) (sbv(m),m=i,k)
16    k = k + 10
      k = 5
      do 17  i = 1 , sm1 , 5
      if(k.gt.sm1) k = sm1
      write(savo,1) ((sm(j,m),m=1,2),j=i,k)
17    k = k + 5
      k = 5
      do 18  i = 1 , sto1 , 5
      if(k.gt.sto1) k = sto1
      write(savo,1) ((strama(j,m),m=1,2),j=i,k)
18    k = k + 5
      k = 5
      do 19  i = 1 , pol1 , 5
      if(k.gt.pol1) k = pol1
      write(savo,1) ((pol(j,m),m=1,2),j=i,k)
19    k = k + 5
      k = 10
      do 20  i = 1 , gatht1 , 10
      if(k.gt.gatht1) k = gatht1
      write(savo,1) (gatht(m),m=i,k)
20    k = k + 10
      k = 5
      do 21  i = 1 , ucht1 , 5
      if(k.gt.ucht1) k = ucht1
      write(savo,1) ((ucht(j,m),m=1,2),j=i,k)
21    k = k + 10
      write(savo,1) ntxc
      k = 5
      do 22  i = 1 , src1 , 5
      if(k.gt.src1) k = src1
      write(savo,1) ((src(j,m),m=1,2),j=i,k)
22    k = k + 5
      do 23  i = 1 , bin1
23    write(savo,1) (bin(i,m),m=1,8)
      k = 3
      do 25  i = 1 , bin1 , 3
      if(k.gt.bin1) k = bin1
      write(savo,24) ((binsta(j,m),m=1,2),j=i,k)
24    format(6d15.8)
25    k = k + 3
      k = 6
      do 26  i = 1 , drn1 , 6
      if(k.gt.drn1) k = drn1
      write(savo,24) (drn(m),m=i,k)
26    k = k + 6
      rewind savo
      return
      end subroutine save



      SUBROUTINE SEIZE(NFA,ID,*,IPRINT)
!     ***
!     ***          CALL SEIZE(NFA,ID,&1005,IPRINT)
!     ***
!     *** FUNCTION  : ACUQUIRE A FACILITY
!     *** PARAMETERS: NFA = FACILITY NUMBER
!     ***             ID  = SEIZE CALL'S STATEMENT NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     BLOCKING DECISION
!     =================
      IF(OK.EQ.0) GOTO 100
      OK = 0
      IF(FAC(NFA,1).NE.0) GOTO 100
!
!     ACQUIRE
!     =======
      FAC(NFA,1) = LTX
      STATE(NFA) = 0
      TX(LTX,5) = ID
      TX(LTX,8) = 0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NFA
3000  FORMAT(8H SEIZE :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      13H ACQUIRES FAC,I3)
      RETURN
!
!     BLOCK
!     =====
100   AL(LTX,1) = ID
      AL(LTX,2) = -NFA
      TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NFA
3001  FORMAT(8H SEIZE :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      18H IS BLOCKED AT FAC,I3)
      RETURN 1
      END SUBROUTINE SEIZE


      SUBROUTINE SELIST(MFA)
!     ***
!     ***          CALL SELIST(NFA)
!     ***
!     *** FUNCTION  : PRINT THE CONTENTS OF THE MULTIFACILITY MFA
!     *** PARAMETER : MFA = MULTIFACILITY NUMBER
!     ***
      IMPLICIT INTEGER(A-Z)
!
!     PRINT THE SE METRIX
!     ===================
      WRITE(OUTD,3000) MFA
3000  FORMAT(1H0,10X,24HS E   M A T R I X : MFAC,I3,//10X,&
      32HLINE   ELEM   CONT    PRC  PHASE)
      I1 = MBV(MFA)
      I2 = I1 + MFAC(MFA,2) - 1
      I3 = 1
      DO 100 I=I1,I2
      WRITE(OUTD,3001) I,I3,(SE(I,J),J=1,3)
3001  FORMAT(10X,5(I4,3X))
100   I3=I3+1
      WRITE(OUTD,3002)
3002  FORMAT(1H0)
      RETURN
      END SUBROUTINE SELIST


      SUBROUTINE SETUP(ST,NFA,IDN,*,*,IPRINT)
!     ***
!     ***          CALL SETUP(ST,NFA,IDN,&1005,&1006,IPRINT)
!     ***
!     *** FUNCTION  : SET UP A FACILITY
!     *** PARAMETERS: ST  = SETUP TIME
!     ***             NFA = FACILITY NUMBER
!     ***             IDN = TARGET
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     ERROR CHECK
!     ===========
      IF(IABS(FAC(NFA,1)).EQ.LTX) GOTO 100
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NFA
3000  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE SETUP ,31(1H+)/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,19H DOESN'T OCCUPY FAC,&
      I3,1X,10(1H+)/)
      RETURN 2
!
!     SET UP THE FACILITY
!     ===================
100   FAC(NFA,1) = -LTX
      FAC(NFA,3) = 1
      AL(LTX,1) = IDN
      AL(LTX,2) = T + ST
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NFA,AL(LTX,2)
3001  FORMAT(8H SETUP :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      4H FAC,I3,24H WILL BE SET UP UNTIL T=,I7)
      RETURN 1
      END SUBROUTINE SETUP


      SUBROUTINE SIMEND(NBN,NE,P,*)
!     ***
!     ***          CALL SIMEND(NBN,NE,P,&1006)
!     ***
!     *** FUNCTION  : FIND THE END OF SIMULATION WITH HELP OF
!     ***             THE MEAN WAITING TIME OF BIN NBN
!     *** PARAMETERS: NBN   = BIN NUMBER
!     ***             NE    = NUMBER OF WITHDRAWN TOKENS BETWEEN
!     ***                     CHECK
!     ***             P     = PERCENTAGE DEVIATION
!     ***             EXIT1 = TERMINATION EXIT
!     ***
      IMPLICIT INTEGER (A-Z)
      REAL RMEAN(20),PR,MU,MO
      DATA PM/0/,NA/0/,RMEAN/20*0./
	  INTEGER pm
!
!     CHECK FOR A NEW TEST
!     ====================
      NT = NA + NE
      IF(BIN(NBN,4).LT.NT) RETURN
!
!     COMPUTE AND ENTER THE NEW MEAN
!     ==============================
      NA = NT
      PM = MOD(PM,20) + 1
      RMEAN(PM) = FLOAT(BIN(NBN,6))/FLOAT(BIN(NBN,4))
!
!     COMPUTE THE ADMISSIBLE INTERVAL
!     ===============================
      PR = (RMEAN(PM)*P)/100.
      MU = RMEAN(PM) - PR
      MO = RMEAN(PM) + PR
!
!     TEST THE MEAN VALUES
!     ====================
      IM = PM
10    IM = MOD(IM,20) + 1
      IF(RMEAN(IM).EQ.0.OR.RMEAN(IM).LT.MU.OR.RMEAN(IM).GT.MO) RETURN
      IF(IM.NE.PM) GOTO 10
      WRITE(OUTD,3000) T,NTXC,NBN,RMEAN(PM),BIN(NBN,4)
3000  FORMAT(1H1,///20X,43H****** S I M U L A T I O N   H A L T ******,&
      //27X,7HT    = ,I10/27X,7HNTXC = ,I10//20X,17HMEAN WAITING TIME,&
      /20X,8HFOR BIN ,I2,4H IS ,F10.5,5H FOR ,I10,8H DEPARTS)
      RETURN 1
      END SUBROUTINE SIMEND 


      SUBROUTINE SMLIST(NST)
!     ***
!     ***          CALL SMLIST(NST)
!     ***
!     *** FUNCTION  : PRINT THE CONTENTS OF THE STORAGE NUMBERED
!     *** PARAMETER : NST = STORAGE NUMBER
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     PRINT THE SM MATRIX
!     ===================
      WRITE(OUTD,3000) NST
3000  FORMAT(1H0,10X,27HS M   M A T R I X : STORAGE,I3//10X,&
      25HLINE   ADDR   SEGL   CODE)
      I = SBV(NST)
      IE = I + STO(NST,2) - 1
10    IADDR = I - SBV(NST) + 1
      WRITE(OUTD,3001) I,IADDR,SM(I,1),SM(I,2)
3001  FORMAT(10X,4(I4,3X))
      I = I + SM(I,1)
      IF(I.LE.IE) GOTO 10
      WRITE(OUTD,3002)
3002  FORMAT(1H0)
      RETURN
      END SUBROUTINE SMLIST


      SUBROUTINE SPLIT(NDUP,DUPADD,*,IPRINT)
!     ***
!     ***          CALL SPLIT(NDUP,DUPADD,&1006,IPRINT)
!     ***
!     *** FUNCTION  : REPLICATE A TX
!     *** PARAMETERS: NDUP   = NUMBER OF DUPLICATES
!     ***             DUPADD = TARGET FOR OFFSPRING
!     ***
      IMPLICIT INTEGER (A - Z)
!
!     TEST
!     ====
      IF(NDUP.EQ.0) RETURN
      IF(IPRINT.EQ.0) GOTO 5000
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NDUP
3000  FORMAT(8H SPLIT :,3X,2HTX,I7,2X,2HTX,I5,1H,,I3,2X,&
      14H IS SPLIT INTO,I3,8H KINSMEN)
5000  IF(LFAM.NE.0) GOTO 300
!
!     FOUND A FAMILY
!     ==============
      DO 100 LFAM=1,FAM1
      IF(FAM(LFAM,1).EQ.0) GOTO 200
100   CONTINUE
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2)
3001  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE SPLIT ,31(1H+),1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,9H OVERRUN 24(1H+)/&
      1X,24(1H+),19H IN THE FAM MATRIX,37(1H+)/)
      RETURN 1
200   FAM(LFAM,1) = TX(LTX,1)
      FAM(LFAM,2) = 1
      FAM(LFAM,3) = 1
      TX(LTX,2) = 1
!
!     RESET
!     =====
300   DO 340  I = 1 , NDUP
      FAM(LFAM,2) = FAM(LFAM,2) + 1
      FAM(LFAM,3) = FAM(LFAM,3) + 1
      DO 310  IA = 1 , TX1
      IF(TX(IA,1).EQ.0) GOTO 320
310   CONTINUE
      WRITE(OUTD,3002) T,TX(LTX,1),TX(LTX,2)
3002  FORMAT(1H0,24(1H+),25H ERROR: SUBROUTINE SPLIT ,31(1H+)/1X,&
      24(1H+),3H T=,I7,2X,2HTX,I5,1H,,I3,9H OVERRUN ,22(1H+)/&
      1X,24(1H+),18H IN THE TX MATRIX ,38(1H+)/)
      RETURN 1
320   LTXD = IA
      TX(LTXD,1) = TX(LTX,1)
      TX(LTXD,2) = FAM(LFAM,3)
      TX(LTXD,3) = T
      DO 330  JC = 4 , TX2
330   TX(LTXD,JC) = TX(LTX,JC)
      IF(LTXD.GT.LAL) LAL = LTXD
      AL(LTXD,1) = DUPADD
      AL(LTXD,2) = T
340   CONTINUE
      RETURN
      END SUBROUTINE SPLIT



      SUBROUTINE STRATA(NST,NE,*)
!     ***
!     ***          CALL STRATA(NST,NE,EXIT1)
!     ***
!     *** FUNCTION  : FIND  THE A-STRATEGY RESPOSIBLE FOR GIVEN
!     ***             STORAGE
!     *** PARAMETERS: NST   =  STORAGE NUMBER
!     ***             NE    =  NUMBER OF ELEMENTS REQUESTED
!     ***             EXIT1 =  STRATEGY-ERROR EXIT
!     ***
      IMPLICIT INTEGER(A-Z)
!
!     FIND THE STRATEGY
!     =================
      IF(STRAMA(NST,1).NE.0) GOTO 100
      WRITE(OUTD,3000) NST
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE STRATA ,30(1H+)/1X,&
      24(1H+),8H FOR STO,I3,24H NO ALLOCATION STRATEGY,21(1H+)/)
      RETURN 1
100   IADDR = STRAMA(NST,1)
      GOTO (1,2,3,4,5),IADDR
!
!     CALL THE STRAA SUBROUTINE
!     =========================
1     CALL FFIT(NST,NE)
      RETURN
2     CALL BFIT(NST,NE)
      RETURN
3     CALL STRAA3(NST,NE)
      RETURN
4     CALL STRAA4(NST,NE)
      RETURN
5     CALL STRAA5(NST,NE)
      RETURN
      END SUBROUTINE STRATA
      SUBROUTINE STRAA3(NST,NE)
      IMPLICIT INTEGER(A-Z)
      RETURN
      END SUBROUTINE STRAA3
      SUBROUTINE STRAA4(NST,NE)
      IMPLICIT INTEGER(A-Z)
      RETURN
      END SUBROUTINE STRAA4
      SUBROUTINE STRAA5(NST,NE)
      IMPLICIT INTEGER(A-Z)
      RETURN
      END SUBROUTINE STRAA5




      SUBROUTINE STRATF(NST,KEY,*)
!     ***
!     ***          CALL STRATF(NST,KEY,EXIT1)
!     ***
!     *** FUNCTION  : FIND THE F-STRATEGY RESPOSIBLE FOR A GIVEN
!     ***             STORAGE
!     *** PARAMETERS: NST   = STORAGE NUMBER
!     ***             KEY   = FREEING KEY
!     ***             EXIT1 = STRATEGY-ERROR EXIT
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     FIND THE STRATEGY
!     =================
      IF(STRAMA(NST,2).NE.0) GOTO 100
      WRITE(OUTD,3000) NST
3000  FORMAT(1H0,24(1H+),26H ERROR: SUBROUTINE STRATF ,30(1H+)/1X,&
      24(1H+),8H FOR STO,I3,21H NO FREEING STRATEGY ,24(1H+)/)
      RETURN 1
100   IADDR = STRAMA(NST,2)
      GOTO (1,2,3,4,5),IADDR
!
!     CALL THE STRAF SUBROUTINE
!     =========================
1     CALL STRAF1(NST,KEY)
      RETURN
2     CALL STRAF2(NST,KEY)
      RETURN
3     CALL STRAF3(NST,KEY)
      RETURN
4     CALL STRAF4(NST,KEY)
      RETURN
5     CALL STRAF5(NST,KEY)
      RETURN
      END SUBROUTINE STRATF
      SUBROUTINE STRAF1(NST,KEY)
      IMPLICIT INTEGER (A-Z)
      RETURN
      END SUBROUTINE STRAF1
      SUBROUTINE STRAF2(NST,KEY)
      IMPLICIT INTEGER (A-Z)
      RETURN
      END SUBROUTINE STRAF2
      SUBROUTINE STRAF3(NST,KEY)
      IMPLICIT INTEGER (A-Z)
      RETURN
      END SUBROUTINE STRAF3
      SUBROUTINE STRAF4(NST,KEY)
      IMPLICIT INTEGER (A-Z)
      RETURN
      END SUBROUTINE STRAF4
      SUBROUTINE STRAF5(NST,KEY)
      IMPLICIT INTEGER (A-Z)
      RETURN
      END SUBROUTINE STRAF5


      SUBROUTINE TABULA(X,Y,OG1,GBR,NG,TAB)
!     ***
!     ***          CALL TABULA(X,Y,OG1,GBR,NG,TAB)
!     ***
!     *** FUNCTION  : ACCUMULATE A FREQUENCY TABLE FOR X
!     ***             AND SUM A COORDINATED VALUE Y
!     ***             ACCORDING TO THE INTERVAL ON WHICH X FALLS
!     *** PARAMETERS: X    = PRINSIPAL DATUM
!     ***             Y    = CORELATED DATUM
!     ***             OG1  = UPPER BOUND OF THE FIRST INTERVAL
!     ***             GBR  = INTERVAL LENGTH
!     ***             NG   = NUMBER OF INTERVALS
!     ***             TAB  = FREQUENCY TABLE
!     ***
      REAL X,Y,OG1,GBR
	  INTEGER NG
      REAL TAB(100,4)
	  REAL G
	  INTEGER I,J
!
!     SET UP THE FREQUENCY TABLE
!     ==========================
      IF(NG.GT.TAB1) NG = TAB1
      IF(TAB(1,1).NE.0.OR.TAB(2,1).NE.0) GOTO 50
      G = OG1
      DO 10 J = 1 , NG
      TAB(J,1) = G
10    G = G + GBR
!
!     TABULATE THE DATA
!     =================
50    DO 100 J = 1, NG
      IF(X.LE.TAB(J,1)) GOTO 150
100   CONTINUE
      I = NG
150   TAB(J,2) = TAB(J,2) + 1
      TAB(J,4) = TAB(J,4) + Y
      RETURN
      END SUBROUTINE TABULA








      SUBROUTINE TERMIN(*,IPRINT)
!     ***
!     ***          CALL TERMIN(&1005,IPRINT)
!     ***
!     *** FUNCTION :  ANNIHILATE A TX
!     ***
      IMPLICIT INTEGER (A-Z)
      IF(IPRINT.EQ.0) GOTO 5000
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2)
3000  FORMAT(8H TERMIN:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      15H IS ANNIHILATED)
!
!     FAMILY
!     ======
5000  IF(LFAM.EQ.0) GOTO 200
      FAM(LFAM,2) = FAM(LFAM,2) - 1
      IF(FAM(LFAM,2).GT.0) GOTO 200
      DO 100 IB = 1 , 3
100   FAM(LFAM,IB) = 0
      DO 110 IB = 1,ASM1
110   ASM(LFAM,IB) = 0
      LFAM = 0
!
!     ANNIHILATE
!     ==========
200   DO 201 IA=1,TX2
201   TX(LTX,IA) = 0
      DO 202 IA=1,2
202   AL(LTX,IA) = 0
!
!     RESET THE LIST-END POINTER, LAL
!     ===============================
250   IF(TX(LAL,1).NE.0.OR.LAL.EQ.1) RETURN 1
      LAL = LAL - 1
      GOTO 250
      END SUBROUTINE TERMIN


      SUBROUTINE TRANSF(RATIO,*,RNUM,IPRINT)
!     ***
!     ***          CALL TRANSF(RATIO,EXIT1,RNUM,IPRINT)
!     ***
!     *** FUNCTION  : STOHASTIC SPLITING OF A TX STREAM
!     ***             ACCORDING TO A SPECIFIED RATIO
!     *** PARAMETERS: RATIO = SELECTION PROBABILITY
!     ***             EXIT1 = EXIT FOR THE SELECTED TXS
!     ***             RNUM  = GENERATOR'S IDENTIFIER
!     ***
      IMPLICIT INTEGER (A-Z)
      REAL RATIO,Z
      Z = RN(RNUM)
      IF(Z.LE.RATIO.AND.RATIO.GT.0) GOTO 100
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2)
3000  FORMAT(8H TRANSF:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      30H CONTINUES TO THE NEXT STATION)
      RETURN
100   CONTINUE
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2)
3001  FORMAT(8H TRANSF:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      34H BRANCHES TO THE SPECIFIED STATION)
      RETURN 1
      END SUBROUTINE TRANSF


      SUBROUTINE UNIFRM(A,B,RNUM,RANDOM)
!     ***
!     ***          CALL UNIFRM(A,B,RNUM,RANDOM)
!     ***
!     *** FUNCTION  : GENERATE A UNIFORM RANDOM SEQUENCE ON THE
!     ***             INTERVAL (A,B)
!     *** PARAMETERS: A      = INTERVAL'S LOWER BOUND
!     ***             B      = INTERVAL'S UPPER BOUND
!     ***             RNUM   = GENERATOR'S IDENTIFIER
!     ***             RENDAOM= RESULTING RANDOM NUMBER
!     ***
      implicit none 
      INTEGER RNUM
	  real a,b,random
      RANDOM = A + (B - A) * RN(RNUM)
      RETURN
      END SUBROUTINE UNIFRM



      SUBROUTINE UNLIN1(NUCHN,NUMIN,NUMAX,ID,*,IPRINT)
!     ***
!     ***          CALL UNLIN1(NUCHN,NUMIN,NUMAX,ID,&1005,IPRINT)
!     ***
!     *** FUNCTION  : REMOVE A GROUP OF TXS FROM A
!     ***             USER CHAIN
!     *** PARAMETERS: NUCHN  = TYPE-1 USER CHAIN NUMBER
!     ***             NUMIN  = MINIMUM REQUEST
!     ***             NUMAX  = MAXIMUM REQUEST
!     ***             ID     = UNLIN1 CALL'S STATEMENNT NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     TEST FOR FAMILY MEMBERSHIP
!     ==========================
      IF(LFAM.EQ.0) RETURN
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EGATHT + 2 * FAM1 * (NUCHN -1) + 2 * LFAM - 1
!
!     BLOCKING DECISION
!     =================
      IF(UCHF(LFAM,NUCHN,1).GE.NUMIN.AND.STATE(K+1).EQ.1) GOTO 100
!
!     BLOCK
!     =====
      STATE(K+1) = 0
      AL(LTX,1) = ID
      AL(LTX,2) = - (K+1)
      IF(TX(LTX,8).EQ.0) TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NUCHN
3000  FORMAT(8H UNLIN1:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      34H IS BLOCKED AT TRIGGER STATION1 NO,I3)
      RETURN 1
!
!     START A REMOVAL PROCEDURE
!     =========================
100   UCHF(LFAM,NUCHN,2) = NUMAX
      STATE(K) = 1
      STATE(K+1) = 0
      IF(UCHF(LFAM,NUCHN,1).EQ.0) UCHF(LFAM,NUCHN,2) = 0
      TX(LTX,8) = 0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NUCHN,NUMIN,NUMAX
3001  FORMAT(8H UNLIN1:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      28H REMOVES FROM USER CHAIN1 NO,I3,/36X,8HAT LEAST,&
      I3,12H BUT AT MOST,I3,4H TXS)
      RETURN
      END SUBROUTINE UNLIN1

      SUBROUTINE UNLIN2(NUCHN,NUMIN,NUMAX,ID,*,IPRINT)
!     ***
!     ***          CALL UNLIN2(NUCHN,NUMIN,NUMAX,ID,&1005,IPRINT)
!     ***
!     *** FUNCTION  : REMOVAL OF TXS AT A USER CHAIN
!     *** PARAMETERS: NUCHN  = TYPE-2 USER CHAIN NUMBER
!     ***             NUMIN  = MINIMUM REQUEST
!     ***             NUMAX  = MAXIMUM REQUEST
!     ***             ID     = UNLIN2 CALL'S STATEMENNT NUMBER
!     ***
      IMPLICIT INTEGER (A-Z)
!
!     DETERMINE THE STATION NUMBER
!     ============================
      K = EUCHF + 2 * NUCHN - 1
!
!     BLOCKING DECISION
!     =================
      IF(UCHT(NUCHN,1).GE.NUMIN.AND.STATE(K+1).EQ.1) GOTO 100
!
!     BLOCK
!     =====
      STATE(K+1) = 0
      AL(LTX,1) = ID
      AL(LTX,2) = - (K+1)
      IF(TX(LTX,8).EQ.0) TX(LTX,8) = T
      IF(IPRINT.EQ.0) RETURN 1
      WRITE(OUTD,3000) T,TX(LTX,1),TX(LTX,2),NUCHN
3000  FORMAT(8H UNLIN2:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      34H IS BLOCKED AT TRIGGER STATION2 NO,I3)
      RETURN 1
!
!     START A REMOVAL PROCEDURE
!     =========================
100   UCHT(NUCHN,2) = NUMAX
      STATE(K) = 1
      STATE(K+1) = 0
      IF(UCHT(NUCHN,1).EQ.0) UCHT(NUCHN,2) = 0
      TX(LTX,8) = 0
      IF(IPRINT.EQ.0) RETURN
      WRITE(OUTD,3001) T,TX(LTX,1),TX(LTX,2),NUCHN,NUMIN,NUMAX
3001  FORMAT(8H UNLIN2:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      28H REMOVES FROM USER CHAIN2 NO,I3,/36X,8HAT LEAST,&
      I3,12H BUT AT MOST,I3,4H TXS)
      RETURN
      END SUBROUTINE UNLIN2

      SUBROUTINE UNLOCK(K,IPRINT)
!     ***
!     ***          CALL UNLOCK(K,IPRINT)
!     ***
!     *** FUNCTION  : START THE LOCKED TXS AT A STATION
!     *** PARAMETER : K = STATION NUMBER
!     ***
      IMPLICIT INTEGER(A-Z)
!
!     START
!     =====
      K1 = - KEND - K
      DO 100 I = 1 , LAL
      IF(AL(I,2).NE.K1) GOTO 100
      AL(I,2) = - K
      IF(IPRINT.EQ.0) GOTO 100
      WRITE(OUTD,3000) T,TX(I,1),TX(I,2)
3000  FORMAT(8H UNLOCK:,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      11H IS STARTED)
100   CONTINUE
      RETURN
      END SUBROUTINE UNLOCK


      SUBROUTINE VERSIO
!     ***
!     ***          CALL VERSIO
!     ***
!     *** FUNCTION  : PRINT A GRSS-F - VIRSION HEADING
!     ***
      IMPLICIT INTEGER(A-Z)
      WRITE(OUTD,1)
1     FORMAT(/1X,22HGPSS FORTRAN SIMULATOR,/   &
              1X,22H======================,//  &
              1X,20HVERSION FORM 9.06.89,/     &
              1X,20H====================,/     &
              1X,20HEXTENDED MFAC=20    ,/     &
              1X,20H====================///)
      RETURN
      END SUBROUTINE VERSIO


end module GPSS