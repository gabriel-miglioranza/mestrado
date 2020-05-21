
      SUBROUTINE colocacao(N,INF,SUP,A,B,Raiz,W)
        !
	IMPLICIT NONE
	Integer ND,N,i,j,iter,niter,NDIM,IHLF, Nt
	External Peso
	Integer, Parameter :: Size=20
	REAL(kind(1.d0)) ALFA(Size),BETA(Size),W(Size),Raiz(Size)
	REAL(kind(1.d0)) A(Size,Size),B(Size,Size)
	REAL(kind(1.d0)) INF, SUP
		
      ND=Size

        ! C�lculo das ra�zes 
      CALL ORTSOM(ND,INF,SUP,N,PESO,ALFA,BETA,RAIZ)
      CALL PESSOMA(ND,INF,SUP,N,PESO,RAIZ,W)
	DO i=N,1,-1
	   Raiz(i+1)=Raiz(i)
	   W(i+1)=W(i)
      END DO
      Nt=N+2
      Raiz(1)=0.0d0
	Raiz(Nt)=1.0d0
	W(1)=0.0d0
	W(Nt)=0.0d0

      ! Calc�lo das derivadas dos polinomios
      call DER1(ND,Nt,RAIZ,A)
	call DER2(ND,Nt,RAIZ,B)
      
      end
     
       
   
C***********************************************************************

	Subroutine Peso(X,Val)
	Implicit double precision (a-z)
        Val = 1.0
      RETURN
	END SUBROUTINE Peso

C***********************************************************************

      SUBROUTINE COLOC(PESO,ND,INF,M,N,ALFA,BETA,W,RAIZ)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	Double Precision INF,M
      EXTERNAL PESO
      DIMENSION ALFA(ND),BETA(ND),RAIZ(ND),W(ND)
C
      CALL ORTSOM(ND,INF,M,N,PESO,ALFA,BETA,RAIZ)
C
      CALL PESSOMA(ND,INF,M,N,PESO,RAIZ,W)
      RETURN
      END
C**************************************************************
C     Subrotina que calcula as raizes do polinomio ortogonal.
C     A integracao interna e feita por Simpson.
C**************************************************************
      SUBROUTINE ORTSOM(ND,INF,M,N,PESO,ALFA,BETA,RAIZ)
      IMPLICIT Double Precision (A-H,O-Z)
	Double Precision INF,M
	INTEGER SizeSR,z
	Parameter (SizeSR=40)
      DIMENSION ALFA(ND),BETA(ND),RAIZ(ND)
      DIMENSION S(SizeSR),R(SizeSR)
C
C     O numero de raizes tem que ser positivo.
C
      IF(N.LE.0) call Error(01)
      IF(M.LE.INF) call Error(02)
	IF(SizeSR.lt.N) call Error(03)
	IF(N.GT.ND) call Error(04)
C
C
C ******* CALCULO DE ALFA(1) *******
C
      BETA(1)=0.
      LIM=M-INF
C     
       I = 10
       ZNODIV = 2**I
       ZNOSIM = ZNODIV/2
       SOMAS = 0.
       WIDTH = LIM/ZNODIV
        DO 40 J = 1, ZNOSIM
        XL = INF + 2.*(J-1)*WIDTH
        XR = INF + 2.*J*WIDTH
        XM = (XL + XR)/2.
        CALL PESO(XL,VALL)
        CALL PESO(XR,VALR)  
        CALL PESO(XM,VALM)  
        A = (WIDTH/3.)*(VALL + 4.*VALM + VALR)
        SOMAS = SOMAS + A
 40     continue
C
       I = 10
       ZNODIV = 2**I
       ZNOSIM = ZNODIV/2
       SOMAR = 0.
       WIDTH = LIM/ZNODIV
        DO 60 J = 1, ZNOSIM
        XL = INF + 2.*(J-1)*WIDTH
        XR = INF + 2.*J*WIDTH
        XM = (XL + XR)/2.
        CALL PESO(XL,VALL)
        vall=vall*xl
        CALL PESO(XR,VALR)  
        valr=valr*xr
        CALL PESO(XM,VALM)
        valm=valm*xm
        A = (WIDTH/3.)*(VALL + 4.*VALM + VALR)
        SOMAR = SOMAR + A
 60     continue
C        
      S(1)=SOMAS
      R(1)=SOMAR
      ALFA(1)=R(1)/S(1)
      IF(N.GT.1) GOTO 5010
      RAIZ(1)=ALFA(1)
      RETURN
C
C     Se N > 1, demais raizes sao calculadas aqui.
C
 5010 CONTINUE
      DO 6020 Z=2,N
        SOMAR=0.
        SOMAS=0.
C
       I = 10
       ZNODIV = 2**I
       ZNOSIM = ZNODIV/2
       WIDTH = LIM/ZNODIV
        DO 90 J = 1, ZNOSIM
        XL = INF + 2.*(J-1)*WIDTH
        XR = INF + 2.*J*WIDTH
        XM = (XL + XR)/2.
        CALL PESO(XL,VAL1L)
        CALL PESO(XR,VAL1R)  
        CALL PESO(XM,VAL1M)  
        K = Z - 1
        CALL POLI(nd,XL,K,ALFA,BETA,VAL2L)
        CALL POLI(nd,XR,K,ALFA,BETA,VAL2R)  
        CALL POLI(nd,XM,K,ALFA,BETA,VAL2M)  
        CALL FUN1(XL,Z,VAL1L,VAL2L,F1XL)
        CALL FUN1(XR,Z,VAL1R,VAL2R,F1XR) 
        CALL FUN1(XM,Z,VAL1M,VAL2M,F1XM) 
        AS = (WIDTH/3.)*(F1XL + 4.*F1XM + F1XR)
        CALL FUN2(XL,Z,VAL1L,VAL2L,F2XL)
        CALL FUN2(XR,Z,VAL1R,VAL2R,F2XR) 
        CALL FUN2(XM,Z,VAL1M,VAL2M,F2XM) 
        AR = (WIDTH/3.)*(F2XL + 4.*F2XM + F2XR)  
        SOMAR = SOMAR + AR
        SOMAS = SOMAS + AS
 90     continue
C
        S(Z)=SOMAS
        R(Z)=SOMAR
        BETA(Z)=S(Z)/S(Z-1)
        ALFA(Z)=(R(Z) - BETA(Z)*R(Z-1))/S(Z)
 6020 CONTINUE
C
      CALL ROOT(ND,INF,M,N,ALFA,BETA,RAIZ)
      RETURN
      END
C**************************************************************
C     Subrotina que calcula as raizes distintas de um polinomio
C     ortogonal, conhecendo-se os coeficiente de recursao.
C**************************************************************
      SUBROUTINE ROOT(ND,INF,M,N,ALFA,BETA,RAIZ)
      IMPLICIT Double Precision (A-H,O-Z)
	Double Precision INF,M
      DIMENSION ALFA(ND),BETA(ND),RAIZ(ND)
      IF(N.LE.0) call Error(01)
      X=INF
      DEL=0.01
      DO 6010 I=1,N
        SOMA=0.
        X=X+DEL
        NITER=0
 5010   CONTINUE
        NITER=NITER+1
        IF(NITER.GT.100) call Error(05)
        CALL POLI(ND,X,N,ALFA,BETA,VAL)
        CALL POLIN(ND,X,N,ALFA,BETA,DER)
        IF(I.EQ.1) GOTO 5020
        SOMA=0.
        DO 6020 J=1,I-1
          SOMA=SOMA + 1/(X-RAIZ(J))
 6020   CONTINUE
 5020   CONTINUE
        X1=X - (VAL/DER)/(1-(VAL/DER)*SOMA)
        DESV=DABS(X1-X)
        X=X1
        IF(DSQRT(DESV).GT.1.0D-6) GOTO 5010
        RAIZ(I)=X
 6010 CONTINUE
      RETURN
      END
C**************************************************************
C     Subrotina que calcula o valor do polinomio, dados os
C     coeficientes da recursao.
C**************************************************************
      SUBROUTINE POLI(ND,X,N,ALFA,BETA,VAL)
      IMPLICIT Double Precision (A-H,O-Z)
      DIMENSION ALFA(ND),BETA(ND),PN(3)
      IF(N.LE.0) GOTO 5010
      PN(1)=1.
      PN(2)=1.
      PN(3)=1.
      DO 6010 J=1,N
        PN(1)=PN(2)
        PN(2)=PN(3)
        PN(3)=(X-ALFA(J))*PN(2) - BETA(J)*PN(1)
 6010 CONTINUE
      VAL=PN(3)
      RETURN
 5010 CONTINUE
      VAL=1.
      RETURN
      END
C**************************************************************
C     Subrotina que calcula o valor da derivada polinomio, 
C     dados os coeficientes da recursao.
C**************************************************************
      SUBROUTINE POLIN(ND,X,N,ALFA,BETA,DER)
      IMPLICIT Double Precision(A-H,O-Z)
      DIMENSION ALFA(ND),BETA(ND),PLIN(3)
      IF(N.LE.0) GOTO 5010
      PLIN(1)=0.
      PLIN(2)=0.
      PLIN(3)=0.
      DO 6010 I=1,N
        K=I-1
        CALL POLI(ND,X,K,ALFA,BETA,VAL)
        PLIN(1)=PLIN(2)
        PLIN(2)=PLIN(3)
        PLIN(3)=VAL + (X-ALFA(I))*PLIN(2) - BETA(I)*PLIN(1)
 6010 CONTINUE
      DER=PLIN(3)
      RETURN
 5010 CONTINUE
      DER=0.
      RETURN
      END
C**************************************************************
C     Subrotina que faz a interpolacao com os polinomios 
C     interpoladores de Lagrange.
C     
C     K e o indice da funcao interpoladora, comecando em 1.
C     VAL e o valor da funcao interpoladora K em Y.
C**************************************************************
      SUBROUTINE AGRAN(ND,N,K,RAIZ,Y,VAL)
      IMPLICIT Double Precision (A-H,O-Z)
      DIMENSION RAIZ(ND)
      PROD1=1.
      PROD2=1.
      DO 6010 I=1,N
        IF(I.EQ.K) GOTO 6010
        PROD1=PROD1*(Y-RAIZ(I))
        PROD2=PROD2*(RAIZ(K)-RAIZ(I))
 6010 CONTINUE
      VAL=PROD1/PROD2
      RETURN
      END
C**************************************************************
C     Subrotina que faz o calculo dos pesos da quadratura.
C**************************************************************
      SUBROUTINE PESSOMA(ND,INF,M,N,PESO,RAIZ,W)
      IMPLICIT Double Precision (A-H,O-Z)
	Double Precision INF,M
      DIMENSION RAIZ(ND),W(ND)
C
      IF(N.LE.0) call Error(01)
      IF(M.LE.INF) call Error(02)
	IF(N.GT.ND) call Error(04)
C
      ZK=M-INF
      DO 6010 J=1,N
        LIM=J
        SOMA=0.
        I = 10
        ZNODIV = 2**I
        ZNOSIM = ZNODIV/2
        WIDTH = ZK/ZNODIV
        DO 110 L = 1, ZNOSIM
        YL = INF + 2.*(L-1)*WIDTH
        YR = INF + 2.*L*WIDTH
        YM = (YL + YR)/2.
        CALL AGRAN(ND,N,LIM,RAIZ,YL,VALL)
        call peso(yl,vallp)
        vall=vall*vallp
        CALL AGRAN(ND,N,LIM,RAIZ,YR,VALR) 
        call peso(yr,valrp)
        valr=valr*valrp
        CALL AGRAN(ND,N,LIM,RAIZ,YM,VALM) 
        call peso(ym,valmp)
        valm=valm*valmp
        A = (WIDTH/3.)*(VALL + 4.*VALM + VALR)
        SOMA = SOMA + A
 110    continue
        W(J)=SOMA
 6010 CONTINUE
      RETURN
      END   
C**************************************************************
C     Subrotina que faz o calculo dos pesos da quadratura para
C     integracoes de Raiz(i) a M.
C**************************************************************
      SUBROUTINE PEMAIS(ND,M,N,RAIZ,I,WMAIS)
      IMPLICIT Double Precision (A-H,O-Z)
	DOUBLE PRECISION M
      DIMENSION RAIZ(ND),WMAIS(ND,ND)
c
      IF(N.LE.0) call Error(01)
	IF(N.GT.ND) call Error(04)
c
      DO 6010 J = 1,N
        ZK = M - RAIZ(I)
        LIM=J
        SOMAM=0.
        II = 10
        ZNODIV = 2**II
        ZNOSIM = ZNODIV/2
        WIDTH = ZK/ZNODIV
        DO 110 L = 1, ZNOSIM
        YL = RAIZ(I) + 2.*(L-1)*WIDTH
        YR = RAIZ(I) + 2.*L*WIDTH
        YM = (YL + YR)/2.
        CALL AGRAN(ND,N,LIM,RAIZ,YL,VALL)
        call peso(yl,vallp)
        vall=vall*vallp
        CALL AGRAN(ND,N,LIM,RAIZ,YR,VALR) 
        call peso(yr,valrp)
        valr=valr*valrp
        CALL AGRAN(ND,N,LIM,RAIZ,YM,VALM) 
        call peso(ym,valmp)
        valm=valm*valmp
        A = (WIDTH/3.)*(VALL + 4.*VALM + VALR)
        SOMAM = SOMAM + A
 110     continue
        WMAIS(I,J)=SOMAM
 6010 CONTINUE
      RETURN
      END   
C**************************************************************
C     Subrotina que faz o calculo dos pesos da quadratura para
C     integracoes de INF ate Raiz(i).
C**************************************************************
      SUBROUTINE PEMENOS(ND,N,W,I,WMAIS,WMENOS) 
      IMPLICIT Double Precision (A-H,O-Z)
      DIMENSION W(ND),WMAIS(ND,ND),WMENOS(ND,ND)
c
      IF(N.LE.0) call Error(01)
	IF(N.GT.ND) call Error(04)
c
      DO 6010 J=1,N
         WM = W(J) - WMAIS(I,J)
         WMENOS(I,J)=WM
 6010 CONTINUE  
      RETURN
      END
C**************************************************************
C     Subrotina que faz o calculo auxiliar para regra de 
C     recursao.
C**************************************************************
      SUBROUTINE FUN1(X,Z,VAL1,VAL2,F1X)
      IMPLICIT Double Precision (A-H,O-Z)
      integer z
      F1X = (X**(Z-1))*VAL1*VAL2
      RETURN 
      END
C**************************************************************
C     Subrotina que faz o calculo auxiliar para regra de 
C     recursao.
C**************************************************************
      SUBROUTINE FUN2(X,Z,VAL1,VAL2,F2X)
      IMPLICIT Double Precision (A-H,O-Z)
      integer z
      F2X = (X**(Z))*VAL1*VAL2
      RETURN 
      END
C**************************************************************
C     Subrotina que faz o calculo da matriz de derivadas 
C     primeiras nos pontos de interpolacao.
C**************************************************************
      SUBROUTINE MDGRAN(ND,N,RAIZ,DLK)
      IMPLICIT Double Precision (A-H,O-Z)
      DIMENSION RAIZ(ND)
      DIMENSION DLK(ND,ND)
       DO 1 K=1,N
        DO 2 J=1,N
        IF(J.EQ.K) THEN 
        CALL DGRAN1(ND,N,K,J,RAIZ,DLK)     
        ELSE 
        CALL DGRAN(ND,N,K,J,RAIZ,DLK)
        ENDIF
 2      CONTINUE
 1     CONTINUE
      RETURN 
      END
C**************************************************************
C     Subrotina que faz o calculo da derivada do interpolador
C     Lagrangeano para K diferente de J.
C**************************************************************
      SUBROUTINE DGRAN(ND,N,K,J,RAIZ,DLK)
      IMPLICIT Double Precision (A-H,O-Z)
      DIMENSION RAIZ(ND)
      DIMENSION DLK(ND,ND)
      PROD1=1.
      PROD2=1.
      DO 6010 I=1,N
        IF(I.EQ.K) GOTO 6010
        IF(I.EQ.J) GOTO 6010    
        PROD1=PROD1*(RAIZ(J)-RAIZ(I))
        PROD2=PROD2*(RAIZ(K)-RAIZ(I))
 6010 CONTINUE
      DLK(K,J)=(PROD1/PROD2)*(1/(RAIZ(K)-RAIZ(J)))
      RETURN
      END
C**************************************************************
C     Subrotina que faz o calculo da derivada do interpolador
C     Lagrangeano para K igual a J.
C**************************************************************
      SUBROUTINE DGRAN1(ND,N,K,J,RAIZ,DLK)
      IMPLICIT Double Precision (A-H,O-Z)
      DIMENSION RAIZ(ND)
      DIMENSION DLK(ND,ND)
      SOMA=0.
      DO 6010 I=1,N
      IF(I.EQ.J) GOTO 6010    
      SOMA=SOMA + 1/(RAIZ(J)-RAIZ(I))
 6010 CONTINUE
      DLK(K,J) = SOMA
      RETURN
      END
C**************************************************************
C     Subrotina que emite as mensagens de erro
C**************************************************************
      Subroutine Error(Nerror)
	Integer Nerror
C
      If(Nerror.eq.01) write(*,1010)
 1010 format(1x,' Erro01 - Numero de raizes deve ser positivo.')
C
      If(Nerror.eq.02) write(*,1020)
 1020 format(1x,' Erro02 - O limite superior e menor que o inferior.')
C
      If(Nerror.eq.03) write(*,1030)
 1030 format(1x,' Erro03 - Vetores S e R subdimensionados em ORTOG.')
C
      If(Nerror.eq.04) write(*,1040)
 1040 format(1x,' Erro04 - Vetores de entrada subdimensionados (N>ND).')
C
      If(Nerror.eq.05) write(*,1050)
 1050 format(1x,' Erro05 - Nao convergiu para as raizes.')
C
	stop
	end
C**************************************************************
C     Subrotina que calcula as derivadas dos polinomios
C     interpoladores de Lagrange.
C**************************************************************
      Subroutine DER1(ND,N,ROOT,A)
	Implicit Double Precision(a-z)
	Integer ND,N,I,J,K
	Dimension Root(N),A(ND,ND)
C
C     O numero de raizes tem que ser positivo.
C
      IF(N.LE.0) call Error(01)
	IF(N.GT.ND) call Error(04)
C
C     Para um polin�mio de grau zero.
C
      IF(N.EQ.1) then
	             A(1,1)=0.0d0
	             return
	           endif
C
C     Para polinomios de grau igual ou superior a 1.
C
      do 6010 i=1,N
         do 6020 j=1,N
C
C        Bases para a derivada e para o interpolador sao iguais.
C
	      if (i.eq.j) then
	         Soma=0.0d0
	         do 6030 k=1,N
	            if(k.eq.i) goto 6030
	            soma=soma+1.0d0/(ROOT(i)-ROOT(k))
 6030          continue
               A(J,I)=soma
C
C        Bases para a derivada e para o interpolador sao distintas.
C
	                  else
	         Prod=1.d0
	         do 6040 k=1,N
	            if(k.eq.i) goto 6040
	            if(k.eq.j) goto 6040
	            Prod=Prod*(ROOT(j)-ROOT(k))/(ROOT(i)-ROOT(k))
 6040          continue
               A(J,I)=Prod/(ROOT(i)-ROOT(j))
	                  endif
 6020    continue
 6010 continue
C
	return
	end
C**************************************************************
C     Subrotina que calcula as derivadas segundas dos polinomios
C     interpoladores de Lagrange.
C**************************************************************
      Subroutine DER2(ND,N,ROOT,B)
	Implicit Double Precision(a-z)
	Integer ND,N,I,J,K1,K2
	Dimension Root(N),B(ND,ND)
C
C     O numero de raizes tem que ser positivo.
C
      IF(N.LE.0) call Error(01)
	IF(N.GT.ND) call Error(04)
C
C     Para um polin�mio de grau zero.
C
      IF(N.EQ.1) then
	             B(1,1)=0.0d0
	             return
	           endif
C
C     Para um polin�mio de grau um.
C
      IF(N.EQ.2) then
	             B(1,1)=0.0d0
	             B(1,2)=0.0d0
	             B(2,1)=0.0d0
	             B(2,2)=0.0d0
	             return
	           endif
C
C     Para polinomios de grau igual ou superior a 2.
C
      do 6010 i=1,N
         do 6020 j=1,N
C
C        Bases para a derivada e para o interpolador sao iguais.
C
	      if (i.eq.j) then
	         Soma=0.0d0
	         do 6030 k1=1,N
	            if(k1.eq.i) goto 6030
                  do 6040 k2=1,N 
	               if(k2.eq.i) goto 6040
	               if(k2.eq.k1) goto 6040
	               soma=soma+1.0d0/(ROOT(i)-ROOT(k1))
     1                              /(ROOT(i)-ROOT(k2))
 6040             continue
 6030          continue
               B(J,I)=soma
C
C        Bases para a derivada e para o interpolador sao distintas.
C
	                  else
	         Soma=0.d0
	         do 6050 k1=1,N
	            if(k1.eq.i) goto 6050
	            if(k1.eq.j) goto 6050
	            Prod=1.0d0
	            do 6060 k2=1,N
	               if(k2.eq.i) goto 6060
	               if(k2.eq.j) goto 6060
	               if(k2.eq.k1) goto 6060
	               Prod=Prod*(ROOT(j)-ROOT(k2))/(ROOT(i)-ROOT(k2))
 6060             continue
                  Soma=Soma+Prod/(ROOT(i)-ROOT(k1))
 6050          continue
               B(J,I)=2*Soma/(ROOT(i)-ROOT(j))
	                  endif
 6020    continue
 6010 continue
C
	return
	end
C
