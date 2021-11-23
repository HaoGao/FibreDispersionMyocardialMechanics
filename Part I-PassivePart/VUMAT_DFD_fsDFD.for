	 !DEC$ ATTRIBUTES ALIAS:"vumat"::VUMAT 
 	  subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock),C(ndir,ndir),
     *  F(ndir,ndir),B(ndir,ndir),ftf(ndir,ndir),sts(ndir,ndir),ftfa(ndir,ndir),
     *  pntn(ndir,ndir),ftstf(ndir,ndir),sig(ndir,ndir),ftntf(ndir,ndir),U(ndir,ndir),UI(ndir,ndir),Ro(ndir,ndir),
     *  Fg(ndir,ndir), pIden(ndir,ndir),signew(ndir, ndir),p_value(ndir,ndir),
     *  f0(ndir),fNn(ndir),fn0(ndir),RT(ndir, ndir),fNnNn(ndir,ndir),
     *  fv(ndir),fvf(ndir),fcar(ndir),stnts(ndir,ndir)
C  
      character*80 cmname
	  
      parameter ( half = 0.5d0,zero = 0.0d0, one  = 1.d0, 
     *            two  = 2.d0, three = 3.d0,index_J = 3,
     *            asmall  = 2.d-16,pi=3.1415926d0)	 
      real*8 FibreSheet(10000,7), PDFdata(1500,5), Fb(3,3), sigmaa(3,3)
      integer*4 PDFelement, StartAnalysis
      common /FibreSheetblk/FibreSheet, PDFdata, PDFelement, StartAnalysis
	  
	  
C
C Read material properties
C
      P_a  = props(1)
      P_b  = props(2)
      P_af = props(3)
      P_bf = props(4)
	  P_as = props(5)
      P_bs = props(6)
	  P_an = props(7)
      P_bn = props(8)
	  P_afs = props(9)
      P_bfs = props(10)
	  P_afn = props(11)
      P_bfn = props(12)
	  P_asn = props(13)
      P_bsn = props(14)

C - Compressible case
      P_D  = props(15)
      dinv=1/P_D
C - active 	 
	  t0= props(16)
	  p_m= props(17)
	  b_len= props(18)
	  p_l0= props(19)
	  b_sha= props(20)
	  Ca0_m= props(21)
	  Ca0= props(22)
	  Tmax= props(23)
	  a_thr= props(24)
	  p_nf= props(25)
	  p_ns= props(26)
	  p_nn= props(27)
	  

	  if (totalTime .EQ. zero .AND. StartAnalysis .EQ. 0) then
			call vexter
	  end if
	  
	  if (totalTime .GT. zero .AND. totalTime .LE. dt .AND. StartAnalysis .EQ. 0) then
			call vexter
	  end if

C
C Please note here, initial value is not zero
      do k = 1, nblock

C--Compute deformation gradient tensor F
		F(1,1)=defgradNew(k,1)
		F(2,2)=defgradNew(k,2)
		F(3,3)=defgradNew(k,3)
		F(1,2)=defgradNew(k,4)
		F(2,3)=defgradNew(k,5)
		F(3,1)=defgradNew(k,6)
		F(2,1)=defgradNew(k,7)
		F(3,2)=defgradNew(k,8)
		F(1,3)=defgradNew(k,9)
	

			U(1,1)=stretchNew(k,1)
			U(2,2)=stretchNew(k,2)
			U(3,3)=stretchNew(k,3)
			U(1,2)=stretchNew(k,4)
			U(2,3)=stretchNew(k,5)
			U(3,1)=stretchNew(k,6)
			U(2,1)=stretchNew(k,4)
			U(3,2)=stretchNew(k,5)
			U(1,3)=stretchNew(k,6)
			
		DETU = U(1,1)*U(2,2)*U(3,3) 
     *      - U(1,1)*U(2,3)*U(3,2)  
     *      - U(1,2)*U(2,1)*U(3,3)  
     *      + U(1,2)*U(2,3)*U(3,1)  
     *      + U(1,3)*U(2,1)*U(3,2)  
     *      - U(1,3)*U(2,2)*U(3,1)		
		
		UI(1,1)=1.0D0/DETU*(U(2,2)*U(3,3)-U(2,3)*U(3,2))
		UI(1,2)=1.0D0/DETU*(U(1,3)*U(3,2)-U(1,2)*U(3,3))
		UI(1,3)=1.0D0/DETU*(U(1,2)*U(2,3)-U(1,3)*U(2,2))
		UI(2,1)=1.0D0/DETU*(U(2,3)*U(3,1)-U(2,1)*U(3,3))
		UI(2,2)=1.0D0/DETU*(U(1,1)*U(3,3)-U(1,3)*U(3,1))
		UI(2,3)=1.0D0/DETU*(U(1,3)*U(2,1)-U(1,1)*U(2,3))
		UI(3,1)=1.0D0/DETU*(U(2,1)*U(3,2)-U(2,2)*U(3,1))
		UI(3,2)=1.0D0/DETU*(U(1,2)*U(3,1)-U(1,1)*U(3,2))
		UI(3,3)=1.0D0/DETU*(U(1,1)*U(2,2)-U(1,2)*U(2,1))
		
C--Compute the determinant 	 
		DET = F(1,1)*F(2,2)*F(3,3) 
     *      - F(1,1)*F(2,3)*F(3,2)  
     *      - F(1,2)*F(2,1)*F(3,3)  
     *      + F(1,2)*F(2,3)*F(3,1)  
     *      + F(1,3)*F(2,1)*F(3,2)  
     *      - F(1,3)*F(2,2)*F(3,1)

C- F is decomposed
        do i=1, ndir
			do j=1, ndir			
			Fb(i,j)=F(i,j)*DET**(-1.0/3.0)
			end do
		end do
C-Compute the C, B, f tensor f, s tensor s			
		do i=1, ndir
		  do j=1, ndir
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+Fb(i,ij)*Fb(j,ij)
			end do
			B(i,j)=tmp
					
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+Fb(ij,i)*Fb(ij,j)
			end do
			C(i,j)=tmp
			
			tmp=0.0d0
			do ij=1, ndir
				tmp=tmp+F(i,ij)*UI(ij,j)
			end do
			Ro(i,j)=tmp
			
			tmp=0.0d0
			pIden(i,j)=tmp
			sigmaa(i,j)=tmp 

		  end do
		end do 
C--Identity tensor			
			pIden(1,1)=1.0d0
			pIden(2,2)=1.0d0
			pIden(3,3)=1.0d0
C--Invariant 
		    p_I1=C(1,1)+C(2,2)+C(3,3)

	    call fibre(Fb,sigmaa,P_af,P_bf,P_as,P_bs,P_an,P_bn,P_afs,P_bfs,P_afn,P_bfn,P_asn,P_bsn)

C--Cauchy stress	
		do i=1, ndir
			do j=1, ndir
				sig(i,j)=P_a*exp(P_b*(p_I1-3.0d0))*(B(i,j)-p_I1/3.0d0*pIden(i,j))
     *        		+sigmaa(i,j)+(DET**2-1.0d0)*dinv*pIden(i,j)
			end do
		end do  			
			
c  ***********************************************************
      Time=totalTime
      t=0.0d0
		 p_lr=0.00185d0
         p_lff=p_lr*sqrt(p_I4f)
		 
      if (tempOld(k).GE.a_thr .AND. stateOld(k,4).LT.asmall) then
	     stateNew(k,4)=Time
		 stateNew(k,3)=p_m*p_lff+b_len
      elseif (stateOld(k,4).GE.asmall) then
	     t=Time-stateOld(k,4)
	     stateNew(k,4)=stateOld(k,4)
		 stateNew(k,3)=stateOld(k,3)
      else
	     t=0.0d0
      end if 
	  
c 	  - Ta

         t_r=p_m*p_lff+b_len
         t1=t_r+t0
         w_t=zero
		 
         if(t .GT. 0.0d0 .and. t .LE. t0)then 
            w_t=pi*t/t0
         end if
         if(t .GT. t0 .and. t .LE. t1)then
            w_t=pi*(t-t0+t_r)/t_r
         end if	 
         if(t .GT. t1)then
            w_t=0.0d0
			stateNew(k,4)=zero
			stateNew(k,3)=zero
         end if
		 	

         if(t .GT. zero .AND. p_lff .GT. p_l0)then
            term1=exp(b_sha*(p_lff-p_l0))
            Eca=Ca0_m/sqrt(term1-1.0d0)
			term3=Ca0**2+Eca**2
            term2=Ca0**2/term3

            Ta=half*Tmax*term2*(1.0d0-cos(w_t))

            stateNew(k,2)=Ta
	
c-------------------------------------

			do i=1, ndir
				do j=1, ndir
					sig(i,j)=sig(i,j)/DET+p_nf*Ta/p_I4f*ftf(i,j)+p_ns*Ta/p_I4s*sts(i,j)+p_nn*Ta/p_I4n*pntn(i,j)
				end do
			end do
			 
         end if
		  
c- Do finial rotation back to initial configuration		 
		  do i=1, ndir
			do j=1, ndir
				tmp=0.0d0
				do ia=1, ndir
					do ib=1, ndir
						tmp=tmp+Ro(ia,i)*sig(ia,ib)*Ro(ib,j)
					end do
				end do
				signew(i,j)=tmp
			end do
		  end do
c  ***********************************************************
		
			stressNew(k,1) = signew(1,1)
			stressNew(k,2) = signew(2,2)
			stressNew(k,3) = signew(3,3)
			stressNew(k,4) = signew(1,2)
			stressNew(k,5) = signew(2,3)
			stressNew(k,6) = signew(3,1)


c                   write(*,*) p_I4fn,fNnNn(1,1)		
C
C
C Update the specific internal energy -
C

			stressPower = half * (
     * 		( stressOld(k,1) + stressNew(k,1) ) * strainInc(k,1) +
     * 		( stressOld(k,2) + stressNew(k,2) ) * strainInc(k,2) +
     * 		( stressOld(k,3) + stressNew(k,3) ) * strainInc(k,3) ) +
     * 		( stressOld(k,4) + stressNew(k,4) ) * strainInc(k,4) +
     * 		( stressOld(k,5) + stressNew(k,5) ) * strainInc(k,5) +
     * 		( stressOld(k,6) + stressNew(k,6) ) * strainInc(k,6)
	 
			enerInternNew(k) = enerInternOld(k) + stressPower / density(k)
C
C Update the dissipated inelastic specific energy -
C
C   
	 
	        plasticWorkInc = 0.0d0
			enerInelasNew(k) = enerInelasOld(k)
     * 		+ plasticWorkInc / density(k)
	 
       end do
C
      return
      end
	  
c-------------------------------------
c-------------------------------------
c read fibre direction and used for tranform stress into local material coordinate 	  
      subroutine fibre(pFb,psigmaa,Pp_af,Pp_bf,Pp_as,Pp_bs,Pp_an,Pp_bn,Pp_afs,Pp_bfs,Pp_afn,Pp_bfn,Pp_asn,Pp_bsn)
      
      include 'vaba_param.inc'   
C     Possible values for the lOp argument
      integer*4 i, j, ij
      real*8  pFb(3,3),pIden(3,3), psigmaa(3,3),pMFD(3),ps0(3),pn0(3),pNn(3),pk1(3),pk(3)
      real*8  pks(3),pkn(3),psr(3),pnr(3),pfNn(3),psNn(3),pnNn(3),pftf(3,3),psts(3,3),pntn(3,3),pftstf(3,3)
	  real*8  pftntf(3,3),pstnts(3,3),pC(3,3)
      real*8 FibreSheet(10000,7), PDFdata(1500,5) 
      integer*4 PDFelement, StartAnalysis
      common /FibreSheetblk/FibreSheet, PDFdata, PDFelement, StartAnalysis

c ****************************	
      
C-Compute the C, B, f tensor f, s tensor s			
			do i=1, 3
				do j=1, 3
					tmp=0.0d0
					pIden(i,j)=tmp
					psigmaa(i,j)=tmp 
				end do
			end do 
			pIden(1,1)=1.0d0
			pIden(2,2)=1.0d0
			pIden(3,3)=1.0d0
			
c-			MFD=[1 0 0], s0=[0 1 0], n0=[0 0 1]	
			pMFD(1)=1.0d0
			pMFD(2)=0.0d0
			pMFD(3)=0.0d0
			ps0(1)=0.0d0
			ps0(2)=1.0d0
			ps0(3)=0.0d0
			pn0(1)=0.0d0
			pn0(2)=0.0d0
			pn0(3)=1.0d0			

C-dispersed fibre 
		do kk=1, PDFelement
			do i=1,3
				pNn(i)=PDFdata(kk,i+2)
				psr(i)=PDFdata(kk,i+2)
c-				pnr(i)=PDFdata(kk,i+2)
			end do


			do i=1,3
				pfNn(i)=pFb(i,1)*pNn(1)+pFb(i,2)*pNn(2)+pFb(i,3)*pNn(3)	
				psNn(i)=pFb(i,1)*psr(1)+pFb(i,2)*psr(2)+pFb(i,3)*psr(3)
c-				pnNn(i)=pFb(i,1)*pnr(1)+pFb(i,2)*pnr(2)+pFb(i,3)*pnr(3)
			end do
			
c-			call pdot(pfNn,pfNn,p_I4fa)
c-			call pdot(psNn,psNn,p_I4sa)
c-			call pdot(pnNn,pnNn,p_I4na)
			p_I4fa=pfNn(1)*pfNn(1)+pfNn(2)*pfNn(2)+pfNn(3)*pfNn(3)
			p_I4sa=psNn(1)*psNn(1)+psNn(2)*psNn(2)+psNn(3)*psNn(3)
c-			p_I4na=pnNn(1)*pnNn(1)+pnNn(2)*pnNn(2)+pnNn(3)*pnNn(3)


			
			do i=1, 3
				do j=1, 3
					pftf(i,j)=pfNn(i)*pfNn(j)
					psts(i,j)=psNn(i)*psNn(j)
c-					pntn(i,j)=pnNn(i)*pnNn(j)
				end do
			end do 
					
        
			if (p_I4fa .GE. 1.0d0) then
				psigmaa=psigmaa+PDFdata(kk,1)*2.0d0*Pp_af*(p_I4fa-1)*exp(Pp_bf*((p_I4fa-1)**2))*(pftf-p_I4fa/3.0d0*pIden)
			end if
        
			if (p_I4sa .GE. 1.0d0) then
				psigmaa=psigmaa+PDFdata(kk,2)*2.0d0*Pp_as*(p_I4sa-1)*exp(Pp_bs*((p_I4sa-1)**2))*(psts-p_I4sa/3.0d0*pIden)
			end if
       
c-			if (p_I4na .GE. 1.0d0) then
c-				psigmaa=psigmaa+PDFdata(kk,1)*2.0d0*Pp_an*(p_I4na-1)*exp(Pp_bn*((p_I4na-1)**2))*(pntn-p_I4na/3.0d0*pIden)
c-			end if		

        
	 
c-        		write(*,*) 'pdf' , PDFdata(kk,1),PDFdata(kk,2),PDFdata(kk,3),PDFdata(kk,4)
		end do
		
c- I_8



			do i=1,3
				do j=1,3
					tmp=0.0d0
					do ij=1, ndir
						tmp=tmp+pFb(ij,i)*pFb(ij,j)
					end do
					pC(i,j)=tmp
					
					pftstf(i,j)=pFb(i,1)*pFb(j,2)+pFb(i,2)*pFb(j,1)
c-					pftntf(i,j)=pFb(i,1)*pFb(j,3)+pFb(i,3)*pFb(j,1)
c-					pstnts(i,j)=pFb(i,2)*pFb(j,3)+pFb(i,3)*pFb(j,2)
				end do
			end do 
		
			p_I8fs=pC(1,2)
c-			p_I8fn=pC(1,3)
c-			p_I8sn=pC(2,3)

			psigmaa=psigmaa+(Pp_afs*p_I8fs*exp(Pp_bfs*p_I8fs**2))*(pftstf-2.0d0*p_I8fs/3.0d0*pIden)
c-     *          	+(Pp_afn*p_I8fn*exp(Pp_bfn*p_I8fn**2))*(pftntf-2.0d0*p_I8fn/3.0d0*pIden)
c-     *	        	+(Pp_asn*p_I8sn*exp(Pp_bsn*p_I8sn**2))*(pstnts-2.0d0*p_I8sn/3.0d0*pIden)
	 
c-		write(*,*) 'invarinat' , p_I4fa,p_I4sa,p_I4na,p_I8fs,p_I8fn,p_I8sn
c-		write(*,*) 'parameter' , Pp_af,Pp_bf,Pp_as,Pp_bs,Pp_an,Pp_bn,Pp_afs,Pp_bfs,Pp_afn,Pp_bfn,Pp_asn,Pp_bsn

		
      return 
      end			
	

c-------------------------------------


c read fibre direction and used for tranform stress into local material coordinate 	  
      subroutine vexter
      
      include 'vaba_param.inc'   
C     Possible values for the lOp argument
      integer*4 i, elemNTotal
      real*8 FibreSheet(10000,7), PDFdata(1500,5) 
      integer*4 PDFelement, StartAnalysis
      common /FibreSheetblk/FibreSheet, PDFdata, PDFelement, StartAnalysis


		StartAnalysis=StartAnalysis+1
!        if( lOp .eq. j_int_StartAnalysis) then
!	    write(*,*) 'uexternaldb:  beginning'
			
!            open(401,FILE='/home/pgrad1/2306902g/DBGuan/DB/UMAT/DBG_DiscreteFibreDispersion/fibre.txt',STATUS='OLD')
!            read(401,*) elemNTotal
!            write(*,*) 'total data to read: n = ', elemNTotal

!            do i = 1 , elemNTotal
!                read(401, *) FibreSheet(i,1),FibreSheet(i,2),FibreSheet(i,3),FibreSheet(i,4)
!            end do 
!            write(*,*) 'read fibre value finished!'
!            close(401)		
			
			
			
			open(402,FILE='D:/DBGuan/HaoSpecialIssue/PassiveFit/ABAQUS/PDF_b1b24539_fsDFD.txt',STATUS='OLD')
            read(402,*) PDFelement
            write(*,*) 'total data to read: n = ', PDFelement

            do i = 1 , PDFelement
                read(402, *) PDFdata(i,1),PDFdata(i,2),PDFdata(i,3),PDFdata(i,4),PDFdata(i,5)
				write(*, *) PDFdata(i,1),PDFdata(i,2),PDFdata(i,3),PDFdata(i,4),PDFdata(i,5)
            end do 
            write(*,*) 'read ,PDFdata(i,3) value finished!'
            close(402)	
	   
!        end if 	  
	  
      return 
      end
	  
c-------------------------------------	  
      subroutine pcross(pa, pb,pcross_r)
      include 'vaba_param.inc'   
      real*8  pcross_r(3)
      real*8  pa(3), pb(3)

		pcross_r(1) = pa(2) * pb(3) - pa(3) * pb(2)
		pcross_r(2) = pa(3) * pb(1) - pa(1) * pb(3)
		pcross_r(3) = pa(1) * pb(2) - pa(2) * pb(1)
      return
      end subroutine 
	  
	  
      subroutine pnorm(pa,pnorm_r)
      include 'vaba_param.inc'   
      real*8   pnorm_r
      real*8   pa(3)

		pnorm_r = pa(1)**2+pa(2)**2+pa(3)**2
		pnorm_r = sqrt(pnorm_r)
      return
      end subroutine 

      subroutine pdot(pa, pb, pdot_r)
      include 'vaba_param.inc'   
      real*8   pdot_r
      real*8   pa(3), pb(3)

		pdot_r = pa(1) * pb(1) + pa(2) * pb(2) + pa(3) * pb(3)
      return
      end subroutine
	  
	  
      subroutine pmatvec(pa, pb, pmatvec_r)
      include 'vaba_param.inc'   
      real*8  pmatvec_r(3)
      real*8  pa(3,3)
      real*8  pb(3)
		
		
		pmatvec_r(1)=pa(1,1)*pb(1)+pa(1,2)*pb(2)+pa(1,3)*pb(3)
		pmatvec_r(2)=pa(2,1)*pb(1)+pa(2,2)*pb(2)+pa(2,3)*pb(3)
		pmatvec_r(3)=pa(3,1)*pb(1)+pa(3,2)*pb(2)+pa(3,3)*pb(3)
      return
      end subroutine 
	  
      subroutine pmatmat(pa, pb, pmatmat_r)
      include 'vaba_param.inc'   
      real*8  tmp
      real*8  pmatmat_r(3,3)
      real*8  pa(3,3), pb(3,3)
      integer :: i,j,ij
		
		
		do i=1, 3
		  do j=1, 3
			tmp=0.0d0
			do ij=1, 3
				tmp=tmp+pa(i,ij)*pb(ij,j)
			end do
			pmatmat_r(i,j)=tmp

		  end do
		end do 
      return
      end subroutine 

      subroutine pvecvec(pa, pb, pvecvec_r)
      include 'vaba_param.inc'   
      real*8  tmp
      real*8  pvecvec_r(3,3)
      real*8  pa(3), pb(3)
      integer :: i,j,ij
		
		
		do i=1, 3
		  do j=1, 3
			pvecvec_r(i,j)=pa(i)*pb(j)

		  end do
		end do 
      return
      end subroutine 