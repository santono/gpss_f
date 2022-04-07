PROGRAM V_LEV_12
USE GPSS               
!     GPSS FORTRAN SIMULATION PROGRAM
!     ===============================
!
!     MODEL ROTOR
!
!     1. GENERAL FORTRAN DECLARATIONS
!
      IMPLICIT INTEGER(A-Z)
!
!     ****** USER'S TYPE DECLARATIONS ******
!
      REAL RAN
      REAL MEAN1 , MIN1 , MAX1 , MEAN2 , MIN2 , MAX2
      REAL SMEAN1 , SMEAN2
      REAL ME,MI,MA,PRCNT,PROIZ,A,B
!     Потери от выхода из строя позиций
      INTEGER PP
      REAL    BPP
!     Потери от простоев машины из-за ремонта всех позиций
      INTEGER PN
      REAL    BPN
      REAL STEP1,STEP2
      REAL PROCP
      INTEGER NSTEP1,NSTEP2,UNLIM,STEPLM
      LOGICAL LS1,LS2,LSREM,LSWORK
 !    CHARACTER*1 C
 !    CHARACTER*1 INBYTE
 !    LOGICAL INSTAT
      
!
!     ****** END OF USER'S TYPE DECLARATIONS ******
!
!
!     ****** USER'S DIMENSION STARTEMENTS ******
!
!    ****** WARNING: ALL FREQUENCY TABLES ******
!     ****** MUST HAVE THE DIMENSION (100,4) ******
!
!     ****** END OF USER'S DIMENSION STATEMENTS ******
!
!
!     ****** USER'S COMMON STATEMENT'S ******
!
!     ****** END OF USER'S COMMON STATEMENTS ******
!
!     2. STATEMERNT FUNCTIONS
!
!    ****** USER'S STATEMENT FUNCTIONS ******
!    ****** END OF USER'S STATEMENT FUNCTIONS ******
!
!    3. CLEAR THE DATA AREA
!
3000  CONTINUE
      CALL INITD
      CALL RESET
!
!     ****** CLEAR USER'S VARIABLES AND TABLES ******
!
!     ****** END OF CLEARING USER'S VARIABLES AND TABLES ******
!
!     4. INITIALIZE CONSTANTS AND CONTROL VARIABLES
!
!     Время моделирования
!     ===================
      NLIM  = 0
!
!     Тип стратегии 1 - до указанного числа позиций
!                   2 - остановки на BRK через ТСМ
!     ============================================= 
      ISTR  = 1
      MEAN2 = 0.
!
!     Средняя длительность ремонта 1-й раб. позиций
!     =============================================
      MEAN1 = 0.
!
!     Число позиций в роторной машине
!     ===============================
      NPOZ  = 0.
!
!     Число позиций при поломке которых производится остановка
!     ========================================================
      NLIM  = 0
!
!     Продолжительность смены
!     =======================
      TSM   = 0
!
!     Длительность перерыва
!     =====================
      BRK   = 0
!
!     Длительность кинематического цикла
!     ==================================
      PROIZ = 0
!   
!     Доля нарастания длительности ремонта при поломке каждого новой позиций
!     ======================================================================
      PRCNT = 0
      OPEN(OUTD,FILE='F2.TXT',STATUS='REPLACE')
4000  CONTINUE
!
!     ****** SET POLICY, STRATEGY AND PLAN MATRICES ******
!     ****** END OF SETTING POLICY, STRATEGY AND PLAN MATRICES ******
!
!     ****** SET STORAGE CAPACITIES *******
!     ****** END OF SETTING STORAGE CAPACITIES
!
!     ****** SET MULTIFACILITY CAPACITIES ******
!     ****** END OF SETTING MULTIFACILITY CAPACITIES *******
!

      CALL INITIO(N,NPOZ,PROIZ,MEAN1,STEP1,NSTEP1,&
      ISTR,IPRINT,SNLIM,UNLIM,STEPLM,&
      PRCNT,MEAN2,STEP2,NSTEP2)
      SMEAN1 = MEAN1
      SMEAN2 = MEAN2
      NLIM   = SNLIM
      IPRINT = 0
      OUTR   = 8
      OPEN(OUTR,FILE='RES_12.TXT',STATUS='REPLACE')
      write(OUTR,88882)
88882 FORMAT(' ',4HNPOZ,4HNLIM,14H    MEAN1     ,14H    MEAN2     ,&
      12H     P      ,7H   B   ,8HME1/MEA2,&
      4HNLIM,7H    B   ,7H 100-B ,7H   PN  ,10H     PP   ,&
      7H  BPN  ,7H  BPP   ,10H     T    ,10H    TP    )

      DO 6660 JNLIM=1,UNLIM
      PRINT *, 'Step 1 =',JNLIM,' total amnt of step ',UNLIM
      DO 7770 JMEAN2=1,NSTEP2
!      WRITE(*,*) 'Шаг 2 =',JMEAN2,' всего шагов ',NSTEP2
      DO 8880 JMEAN1=1,NSTEP1
      CALL RESET
!      WRITE(*,*) 'Шаг 3 ',JMEAN1,' всего шагов ',NSTEP1
      
      CALL INIT1
      CALL INIT2(*9999)
      CALL INIT3(*9999)
!
!     5. READ IN CONTROL VALUES, INITUIAL VALUES AND FUNCTION-VAL-
!        TABLES
!
5000  CONTINUE
!      CALL VERSIO
       CONTIN=0
       SECURE=0
      GOTO 5500
5097  WRITE(*,*) '**** BAD INPUT ****'
      GOTO 9999
!
5500  CONTINUE
!
!     ****** READ IN USER'S VALUES ******
!     ****** END OF READING IN USER'S VALUES ******
!
!     ****** SET USER'S VALUES ******
!      CALL SETTMR
      P=0
      TP=0
!      LS1=.TRUE.
      LS1    = .FALSE.
      LSWORK = .FALSE.
      LSREM  = .FALSE.
      MIN1   = 0.17*MEAN1
      MAX1   = 5.*MEAN1
      MIN2   = 0.17*MEAN2
      MAX2   = 5.*MEAN2
      TN     = 0
      NFLR   = 0
      TB     = TSM
      PN     = 0            !К-во не изготовленных деталей
      PP     = 0

!     ****** END OF SETTING USER'S VALUES ******
!
      IF(CONTIN.GT.0) CALL CONT
!
!     6. SCHEDULE THE FIRST EVENTS
!
6000  CONTINUE
      CALL EVENT(1,1,1,*1006,IPRINT)
      CALL EVENT(1,4,2,*1006,IPRINT)
      CALL EVENT(N,7,3,*1006,IPRINT)
!
!     7. FLOW MANAGEMENT 1 : EVENTS AND SCHEDULED ACTIVATIONS
!
1001  CONTINUE
      CALL ACTIV1(*1006)
      LT = T
      NT = RT
!
!     ****** TARGET SELECTOR ******
!
1003  CONTINUE
!     IF(.NOT.INSTAT(DUMMY)) GOTO 2003
!     C=INBYTE(DUMMY)
!     IF(C.EQ.'T') CALL PRTTMR
!     IF(C.EQ.'O') IPRINT=1-IPRINT
!     IF(C.EQ.'Q') GOTO 1006
2003  CONTINUE
      GOTO (1,2,3,4,5,6,7,8,9,10),NADDR
      WRITE(*,*) 'ERROR NADDR=',NADDR
      STOP
!     ****** END OF TARGET SELEKTOR ******
!
!     9.  MODEL
!
!     ****** MODEL ******
!
!     CREATE THE TXS
!     ==============
1     CALL GENERA(0,NPOZ,1,1,*1006,IPRINT)
1101   CALL ERLANG(MEAN1,1,MIN1,MAX1,1,RAN,*1006)
      if (JNLIM.NE.1)  IPRINT=0
      if (JMEAN2.NE.1) IPRINT=0
      if (JMEAN1.NE.1) IPRINT=0
      if (T.GT.5000) IPRINT=0

      IZ=IFIX(RAN+.5)
      CALL ADVANC(IZ,2,*1005,IPRINT)
2     IF (FAC(1,1).GT.0) GOTO 3
!     Количество вышедших из строя раб позиций
      NFLR=NFLR+1
!     Общее количество изготовленных деталей
      P  = FLOAT(P)+FLOAT(T-TN)*FLOAT(NPOZ-NFLR+1)*PROIZ
!     Количество не изготовленных деталей
      PN = FLOAT(PN) + FLOAT(T-TN)*FLOAT(NFLR-1)*PROIZ
      AI=FLOAT(T)*FLOAT(NPOZ)*PROIZ
      IF (P.GT.AI) WRITE(outd,3033) T,P,AI
3033  FORMAT(8H P>AI  :,3X,2HT=,I7,2X,2HP=,F15.0,4H,AI=,F15.2)
!     Отметка времени текущего состояния (количество вышедших из строя позиций)
      TN=T
8     IF(.NOT.LS2(ISTR,TB,NFLR,NLIM)) GOTO 31
!     Выпустить транзакт "ремонтник"
      IF (IPRINT.EQ.1) WRITE(OUTD,3006) T,TX(LTX,1),TX(LTX,2)
3006  FORMAT(8H UNL1  :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      13H UNLOCK REMMM)
      LSREM=.TRUE.
      CALL UNLOCK(ESTO+2,IPRINT)
!     Для первой стратегии обход
31    IF(ISTR.EQ.1) GOTO 3
      IF(NFLR.LT.NPOZ) GO TO 3
      IF(T.GE.TB) GOTO 3
      IZ=TB-T
!     TP=TP+IZ
      CALL ADVANC(IZ,8,*1005,IPRINT)
3     CALL GATE1(LSWORK,0,1,3,1,*1005,IPRINT)
      GOTO 1101
!
!     Секцмя моделирования ремонта
!     ============================
!
4     CALL GENERA(0,1,4,2,*1006,IPRINT)
6     CALL GATE1(LSREM,0,2,6,1,*1005,IPRINT)
10    CALL SEIZE(1,10,*1005,IPRINT)
      LSREM=.FALSE.
      IF (IPRINT.EQ.1) WRITE(OUTD,3008) T,TX(LTX,1),TX(LTX,2)
3008  FORMAT(8H UNL1  :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      13H MOVED-------)
      TB=TB+TSM+BRK
!     P=FLOAT(P)+FLOAT(T-TN)*FLOAT(NPOZ-NFLR+1)*PROIZ
      IF (NFLR.GT.1) THEN
         A=PRCNT*FLOAT(NFLR-1)
      ELSE
         A=0
      END IF
      ME=MEAN2*(1.+A)
      MI=MIN2 *(1.+A)
      MA=MAX2 *(1.+A)
      IF(ISTR.NE.1) GOTO 29
      CALL ERLANG(MEAN2,1,MIN2,MAX2,1,RAN,*1006)
      IZ=IFIX(RAN+.5)
      GOTO 19
29    IZ=BRK
19    CONTINUE
      DO 11 I=1,LAL
!     IF(I.EQ.LTX) GOTO 11
!     IF(AL(I,2).GE.T) AL(I,2)=AL(I,2)+IZ
11    CONTINUE
      CALL ADVANC(IZ,5,*1005,IPRINT)
!     Время простоя 
5     TP=TP+IZ
      I=ESTO+1
      IF (IPRINT.EQ.1) WRITE(OUTD,3007) T,TX(LTX,1),TX(LTX,2)
3007  FORMAT(8H UNL1  :,3X,2HT=,I7,2X,2HTX,I5,1H,,I3,2X,&
      13H UNLOCK QQQQQ)
      LSWORK=.TRUE.
      TN=T
      CALL UNLOCK(I,IPRINT)
      CALL ADVANC(1,9,*1005,IPRINT)
9     LSWORK=.FALSE.
      NFLR=0
      TN=T
      CALL CLEAR(1,*1005,*1006,IPRINT)
      GOTO 6
!
!     10. FLOW MANAGEMENT 2 : CONDITIONED ACTIVATIONS
!
1005  CONTINUE
      CALL ACTIV2(*1001,*1003)
!
!     11. FINAL ANALYSIS
!
7     CONTINUE
1006  CONTINUE
      T=LT
      RT=NT
      CALL ENDBIN
!
!     ****** USER'S FINAL ANALYSIS *******
!
!
!     ******** END OF USER'S  FINAL ANALYSIS ********
!
!     12. PRINT RESULTS
!
8000  CONTINUE
       AI=FLOAT(T)*FLOAT(NPOZ)*PROIZ
       B=FLOAT(P)/FLOAT(AI)*100.
       PP  = TP*FLOAT(NPOZ)*PROIZ
       BPN = FLOAT(PN) / FLOAT(AI) * 100. 
       BPP = FLOAT(PP) / FLOAT(AI) * 100. 
       IF (ABS(BPN+BPP).LT.0.01) PROCP=1
       IF (ABS(BPN+BPP).GE.0.01) PROCP=BPP/(BPN+BPP)
    !   BPN=(100-B)*(1-PROCP)
    !   BPP=(100-B)*PROCP 
       WRITE(OUTR,8881) NPOZ,NLIM,MEAN1,MEAN2,P,B,MEAN1/MEAN2,NLIM,B,&
       100-B,PN,PP,BPN,BPP,T,TP
8881   FORMAT(' ',I4,I4,F14.4,F14.2,I12,F7.2,F8.2,I4,F7.2,&
       F7.2,I7,I10,F7.2,F7.2,I10,I10)
!      WRITE(*,'(A15)') 'Внутренний цикл'
       MEAN1=MEAN1+STEP1
8880  CONTINUE
      MEAN1=SMEAN1
      MEAN2=MEAN2+STEP2
7770  CONTINUE
      MEAN2=SMEAN2
      MEAN1=SMEAN1
      NLIM=NLIM+STEPLM
6660  CONTINUE
      IF(SECURE.EQ.1) CALL SAVE
!
!     ****** PRINT USER'S RESULTS ******
!
!     ****** END OF PRINTING USER'S RESULTS ******
!
9999  CLOSE(OUTR)
      STOP
      END PROGRAM V_LEV_12
      FUNCTION LS2(ISTR,TSM,NFLR,NLIM)
!
!     Функция возвращает TRUE можно выпускать транзакт-ремонтник для обоих стратегий
!
      IMPLICIT INTEGER(A-Z)
      LOGICAL LS2
      LS2=.FALSE.
      IF(ISTR.NE.1) GOTO 1
!     Первая стратегия - по количеству вышедших из строя позиций 
      IF(NLIM.LE.NFLR) LS2=.TRUE.
      RETURN
!     Вторая стратеги по дительности смены 
1     IF(T.EQ.TSM) LS2=.TRUE.
      RETURN
      END FUNCTION LS2
      SUBROUTINE INITIO(N,NPOZ,PROIZ,MEAN1,STEP1,NSTEP1,&
      ISTR,IPRINT,NLIM,UNLIM,STEPLM,&
      PRCNT,MEAN2,STEP2,NSTEP2)
      IMPLICIT INTEGER(A-Z)
      REAL MEAN1,MEAN2,PROIZ,STEP1,STEP2,PRCNT
      OPEN(1,FILE='ROTOR11.TXT')
      READ(1,*) N
      WRITE(*,*) ' TIME OF MODELING (MIN.)                  -->',N
      READ(1,*) NPOZ
      WRITE(*,*) ' NUMBER OF WORK POSITION                  -->',NPOZ
      READ(1,*) PROIZ
      WRITE(*,*) ' TIME OF KINEMATIC CYCLE (MIN.)           -->',PROIZ
      READ(1,*) MEAN1,STEP1,NSTEP1
      WRITE(*,*) ' MEAN NARABOTKA NA OTKAZ SINGLE R. POZ.   -->',MEAN1
      WRITE(*,*) ' STEP AND NMB OF STEPS                    -->',STEP1,NSTEP1
      READ(1,*) ISTR
      WRITE(*,*) ' KIND OF STRATEGY WORK (1 OR 2)',ISTR
      READ(1,*) IPRINT
      WRITE(*,*) ' TARCING FLAG (0 OR 1)                    -->',IPRINT
!      IF(ISTR.EQ.2) GOTO 5096
!     FIRST STRATEGY
      READ(1,*) NLIM,STEPLM,UNLIM
      WRITE(*,*) ' NMB OTKAZ BEGIN ',NLIM,' NMB OF STEPS ',STEPLM,' STEP SIZE ',UNLIM
      READ(1,*) PRCNT
      WRITE(*,*) ' PART OF INCREMENTING TIME OF REPAIR (0..1)-->',PRCNT
      READ(1,*) MEAN2
      READ(1,*) STEP2
      READ(1,*) NSTEP2
      WRITE(*,*) ' MEAN TIME OF REPAIR                       -->',MEAN2
      WRITE(*,*) ' Step and Nmb Of STep                      -->',STEP2,NSTEP2
    !  PAUSE    'Press Any Key For Continue for Starting Modeling Process'
      CLOSE(1)
      RETURN
      END  SUBROUTINE INITIO