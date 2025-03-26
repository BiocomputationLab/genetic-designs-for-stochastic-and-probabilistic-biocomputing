!Genetic designs for stochastic and probabilistic biocomputing
!Lewis Grozinger,Jesús Miró-Bueno,Ángel Goñi-Moreno

program pbits_NOT_gate
   implicit none

   real(8)::t,a0,r1,r2,tau,v,aux1,aux2,aux_write,incremento_t
   real(8),dimension(:),allocatable::a,c,h
   integer,dimension(:),allocatable::y
   integer::n,m,mu,i
   logical::p0_1,p1_1,p0_2,p1_2
   real(8)::time_0_1,time_1_1,time_0_2,time_1_2
   real(8)::time00,time01,time10,time11

   call random_seed()

   open(1,file='time1.dat')
   open(2,file='time2.dat')
   open(3,file='bars_genes.dat')

   n=14
   m=30
   v=1.d0
   allocate(a(1:m),h(1:m),c(1:m),y(1:n))

!GENE1
   time_0_1=0.d0
   time_1_1=0.d0
   p0_1=.false.
   p1_1=.false.

!GENE2
   time_0_2=0.d0
   time_1_2=0.d0
   p0_2=.false.
   p1_2=.true.

   time00=0.d0
   time01=0.d0
   time10=0.d0
   time11=0.d0

   write(1,*)'time',',','P_1',',','Pa_1',',','Pr_1',',','Pra_1',',','M_1',',','TF_1',',','TF2_1'
   write(2,*)'time',',','P_2',',','Pa_2',',','Pr_2',',','Pra_2',',','M_2',',','TF_2',',','TF2_2'
   write(3,*)'time',',','00',',','01',',','10',',','11'

   a0=0.d0
   t=0.d0

   aux1=0.d0
   aux2=0.d0
   aux_write=100.d0
   incremento_t=0.01d0

   y=0
   a=0.d0
   c=0.d0
   h=0.d0

!GENE 1
   c(1)=10.d0   !reaction 1: P_1-->Pr_1
   c(2)=1.d0   !reaction 2: Pr_1-->P_1
   c(3)=0.01d0   !reaction 3: Pa_1-->Pra_1
   c(4)=c(2)   !reaction 4: Pra_1-->Pa_1
   c(5)=1.0d0/v   !reaction 5: TF2_2+P_1-->Pa_1
   c(6)=1.d0   !reaction 6: Pa_1-->TF2_2+P_1
   c(7)=c(5)   !reaction 7: TF2_2+Pr_1-->Pra_1
   c(8)=c(6)   !reaction 8: Pra_1-->TF2_2+Pr_1
   c(9)=10.d0   !reaction 9: Pr_1-->Pr_1+M_1k2 
   c(10)=c(9)   !reaction 10: Pra_1-->Pra_1+M_1
   c(11)=10.d0   !reaction 11: M_1-->...k4
   c(12)=50.d0   !reaction 12: M_1-->M_1+TF_1k5
   c(13)=1.d0   !reaction 13: TF_1-->...k6
   c(14)=1.d0/v*0.1   !reaction 14: TF_1+TF_1--> TF2_1
   c(15)=1.d0   !reaction 15: TF2_1 -->TF_1+TF_1

   if (p0_1.eqv..true.) then
      c(1)=0.01d0   !reaction 1: P_1-->Pr_1
      c(3)=0.01d0   !reaction 3: Pa_1-->Pra_1
   else if (p1_1.eqv..true.) then
      c(1)=10.d0   !reaction 1: P_1-->Pr_1
      c(3)=10.d0   !reaction 3: Pa_1-->Pra_1
   end if

!GENE 2
   c(16)=10.d0   !reaction 16: P_2-->Pr_2
   c(17)=1.d0   !reaction 17: Pr_2-->P_2
   c(18)=0.01d0   !reaction 18: Pa_2-->Pra_2
   c(19)=c(17)   !reaction 19: Pra_2-->Pa_2
   c(20)=1.0d0/v   !reaction 20: TF2_1+P_2-->Pa_2
   c(21)=1.d0   !reaction 21: Pa_2-->TF2_1+P_2
   c(22)=c(20)   !reaction 22: TF2_1+Pr_2-->Pra_2
   c(23)=c(21)   !reaction 23: Pra_2-->TF2_1+Pr_2
   c(24)=10.d0   !reaction 24: Pr_2-->Pr_2+M_2k2 
   c(25)=c(24)   !reaction 25: Pra_2-->Pra_2+M_2
   c(26)=10.d0   !reaction 26: M_2-->...k4
   c(27)=50.d0   !reaction 27: M_2-->M_2+TF_2k5
   c(28)=1.d0   !reaction 28: TF_2-->...k6
   c(29)=1.d0/v*0.1   !reaction 29: TF_2+TF_2--> TF2_2
   c(30)=1.d0   !reaction 30: TF2_2 -->TF_2+TF_2

   if (p0_2.eqv..true.) then
      c(16)=0.01d0
      c(18)=0.01d0
   else if (p1_2.eqv..true.) then
      c(16)=10.d0
      c(18)=10.d0
   end if

   y(1)=1   !P_1
   y(2)=0   !Pa_1
   y(3)=0   !Pr_1
   y(4)=0   !Pra_1
   y(5)=0   !M_1
   y(6)=50   !TF_1
   y(7)=100   !TF2_1

   y(8)=1   !P_2
   y(9)=0   !Pa_2
   y(10)=0   !Pr_2
   y(11)=0   !Pra_2
   y(12)=0   !M_2
   y(13)=50   !TF_2
   y(14)=100   !TF2_2

   do while(t<4000)   !BUCLE PRINCIPAL

      if (t>=aux_write) then 
         write(1,*)t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',',y(6),',',y(7)
         write(2,*)t,',',y(8),',',y(9),',',y(10),',',y(11),',',y(12),',',y(13),',',y(14)
         aux_write=t+incremento_t
      end if

      h(1)=real(y(1),8)   !reaction 1
      h(2)=real(y(3),8)   !reaction 2
      h(3)=real(y(2),8)   !reaction 3
      h(4)=real(y(4),8)   !reaction 4
      h(5)=real(y(1),8)*real(y(14),8)   !reaction 5
      h(6)=real(y(2),8)   !reaction 6
      h(7)=real(y(3),8)*real(y(14),8)   !reaction 7
      h(8)=real(y(4),8)   !reaction 8
      h(9)=real(y(3),8)   !reaction 9
      h(10)=real(y(4),8)   !reaction 10
      h(11)=real(y(5),8)   !reaction 11
      h(12)=real(y(5),8)   !reaction 12
      h(13)=real(y(6),8)   !reaction 13
      h(14)=real(y(6),8)*(real(y(6),8)-1)/2   !reaction 14
      h(15)=real(y(7),8)   !reaction 15
      h(16)=real(y(8),8)   !reaction 16
      h(17)=real(y(10),8)   !reaction 17
      h(18)=real(y(9),8)   !reaction 18
      h(19)=real(y(11),8)   !reaction 19
      h(20)=real(y(8),8)*real(y(7),8)   !reaction 20
      h(21)=real(y(9),8)   !reaction 21
      h(22)=real(y(10),8)*real(y(7),8)   !reaction 22
      h(23)=real(y(11),8)   !reaction 23
      h(24)=real(y(10),8)   !reaction 24
      h(25)=real(y(11),8)   !reaction 25
      h(26)=real(y(12),8)   !reaction 26
      h(27)=real(y(12),8)   !reaction 27
      h(28)=real(y(13),8)   !reaction 28
      h(29)=real(y(13),8)*(real(y(13),8)-1)/2   !reaction 29
      h(30)=real(y(14),8)   !reaction 30

      do i=1,m
         a(i)=h(i)*c(i)
      end do

      a0=0.d0
      do i=1,m
         a0=a0+a(i)
      end do

      call random_number(r1)
      do while (r1==0.d0)
         call random_number(r1)
      end do
      call random_number(r2)

      tau=(1.d0/a0)*dlog(1.d0/r1)

      if ((y(3)+y(4)==1).and.(y(10)+y(11)==1)) then
         time11=time11+tau
      else if ((y(3)+y(4)==1).and.(y(10)+y(11)==0)) then
         time10=time10+tau
      else if ((y(3)+y(4)==0).and.(y(10)+y(11)==1)) then
         time01=time01+tau
      else if ((y(3)+y(4)==0).and.(y(10)+y(11)==0)) then
         time00=time00+tau
      end if

      aux1=r2*a0
      aux2=0.0
      mu=0
      do i=1,m
         aux2=aux2+a(i)
         if (aux1<aux2) then 
            mu=i
            exit
         end if
      end do

      t=t+tau

      if (mu==1) then !reaction 1
         y(1)=y(1)-1
         y(3)=y(3)+1
      else if (mu==2) then !reaction 2
         y(1)=y(1)+1
         y(3)=y(3)-1
      else if (mu==3) then !reaction 3
         y(2)=y(2)-1
         y(4)=y(4)+1
      else if (mu==4) then !reaction 4
         y(2)=y(2)+1
         y(4)=y(4)-1
      else if (mu==5) then !reaction 5
         y(1)=y(1)-1
         y(14)=y(14)-1
         y(2)=y(2)+1
      else if (mu==6) then !reaction 6
         y(1)=y(1)+1
         y(14)=y(14)+1
         y(2)=y(2)-1
      else if (mu==7) then !reaction 7
         y(3)=y(3)-1
         y(14)=y(14)-1
         y(4)=y(4)+1
      else if (mu==8) then !reaction 8
         y(3)=y(3)+1
         y(14)=y(14)+1
         y(4)=y(4)-1
      else if (mu==9) then !reaction 9
         y(5)=y(5)+1
      else if (mu==10) then !reaction 10
         y(5)=y(5)+1
      else if (mu==11) then !reaction 11
         y(5)=y(5)-1
      else if (mu==12) then !reaction 12
         y(6)=y(6)+1
      else if (mu==13) then !reaction 13
         y(6)=y(6)-1
      else if (mu==14) then !reaction 14
         y(6)=y(6)-2
         y(7)=y(7)+1
      else if (mu==15) then !reaction 15
         y(6)=y(6)+2
         y(7)=y(7)-1
      else if (mu==16) then !reaction 16
         y(8)=y(8)-1
         y(10)=y(10)+1
      else if (mu==17) then !reaction 17
         y(8)=y(8)+1
         y(10)=y(10)-1
      else if (mu==18) then !reaction 18
         y(9)=y(9)-1
         y(11)=y(11)+1
      else if (mu==19) then !reaction 19
         y(9)=y(9)+1
         y(11)=y(11)-1
      else if (mu==20) then !reaction 20
         y(8)=y(8)-1
         y(7)=y(7)-1
         y(9)=y(9)+1
      else if (mu==21) then !reaction 21
         y(8)=y(8)+1
         y(7)=y(7)+1
         y(9)=y(9)-1
      else if (mu==22) then !reaction 22
         y(10)=y(10)-1
         y(7)=y(7)-1
         y(11)=y(11)+1
      else if (mu==23) then !reaction 23
         y(10)=y(10)+1
         y(7)=y(7)+1
         y(11)=y(11)-1
      else if (mu==24) then !reaction 24
         y(12)=y(12)+1
      else if (mu==25) then !reaction 25
         y(12)=y(12)+1
      else if (mu==26) then !reaction 26
         y(12)=y(12)-1
      else if (mu==27) then !reaction 27
         y(13)=y(13)+1
      else if (mu==28) then !reaction 28
         y(13)=y(13)-1
      else if (mu==29) then !reaction 29
         y(13)=y(13)-2
         y(14)=y(14)+1
      else if (mu==30) then !reaction 30
         y(13)=y(13)+2
         y(14)=y(14)-1
      end if

   end do 

   write(3,*)t,',',time00,',',time01,',',time10,',',time11

end program
