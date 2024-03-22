!Genetic designs for stochastic and probabilistic biocomputing
!Lewis Grozinger,Jesús Miró-Bueno,Ángel Goñi-Moreno

program pbits_AND_gate
   implicit none

   real(8)::t,a0,r1,r2,tau,v,aux1,aux2,aux_write,incremento_t
   real(8),dimension(:),allocatable::a,c,h 
   integer,dimension(:),allocatable::y
   integer::n,m,mu,i
   logical::p0_1,p1_1,p0_2,p1_2,p0_3,p1_3
   real(8)::time_0_1,time_1_1,time_0_2,time_1_2,time_0_3,time_1_3
   real(8)::time000,time001,time010,time011,time100,time101,time110,time111

   call random_seed()

   open(1,file='time1.dat')
   open(2,file='time2.dat')
   open(3,file='time3.dat')
   open(4,file='bars_genes.dat')

   n=33
   m=99
   v=1.d0
   allocate(a(1:m),h(1:m),c(1:m),y(1:n))


!GENE1
   time_0_1=0.d0
   time_1_1=0.d0
   p0_1=.true.
   p1_1=.false.

!GENE2
   time_0_2=0.d0
   time_1_2=0.d0
   p0_2=.false.
   p1_2=.false.

!GENE3
   time_0_3=0.d0
   time_1_3=0.d0
   p0_3=.false.
   p1_3=.false.

!GENE1
   time_0_TF2_1=0.d0
   time_1_TF2_1=0.d0

!GENE2
   time_0_TF2_2=0.d0
   time_1_TF2_2=0.d0

!GENE3
   time_0_TF2_3=0.d0
   time_1_TF2_3=0.d0

!three genes
   time000=0.d0
   time001=0.d0
   time010=0.d0
   time011=0.d0
   time100=0.d0
   time101=0.d0
   time110=0.d0
   time111=0.d0


write(1,*)'time',',','P_1',',','Pa2_1',',','Pa3_1',',','Pa23_1',',','Pr_1',',','Pra2_1',',','Pra3_1',',','Pra23_1',',','M_1',','&
   ,'TF_1',',','TF2_1'
   write(2,*)'time',',','P_2',',','Pa2_2',',','Pa3_2',',','Pa23_2',',','Pr_2',',','Pra2_2',',','Pra3_2',',','Pra23_2',',','M_2',','&
   ,'TF_2',',','TF2_2'
   write(3,*)'time',',','P_3',',','Pa2_3',',','Pa3_3',',','Pa23_3',',','Pr_3',',','Pra2_3',',','Pra3_3',',','Pra23_3',',','M_3',','&
   ,'TF_3',',','TF2_3'
   write(4,*)'time',',','000',',','001',',','010',',','011',',','100',',','101',',','110',',','111'
   write(5,*)'time',',','000',',','001',',','010',',','011',',','100',',','101',',','110',',','111'

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
   c(3)=0.01d0   !reaction 3: Pa2_1-->Pra2_1 
   c(4)=c(2)   !reaction 4: Pra2_1-->Pa2_1
   c(5)=10.d0   !reaction 5: Pa3_1-->Pra3_1 
   c(6)=c(2)   !reaction 6: Pra3_1-->Pa3_1
   c(7)=10.d0   !reaction 7: Pa23_1-->Pra23_1 
   c(8)=c(2)   !reaction 8: Pra23_1-->Pa23_1
   c(9)=1.0d0/v   !reaction 9: TF2_2+P_1-->Pa2_1
   c(10)=1.d0   !reaction 10: Pa2_1-->TF2_2+P_1
   c(11)=c(9)   !reaction11: TF2_2+Pr_1-->Pra2_1
   c(12)=c(10)   !reaction 12: Pra2_1-->TF2_2+Pr_1
   c(13)=1.0d0/v   !reaction 13: TF2_3+P_1-->Pa3_1
   c(14)=1.d0   !reaction 14: Pa3_1-->TF2_3+P_1
   c(15)=c(9)   !reaction 15: TF2_3+Pr_1-->Pra3_1
   c(16)=c(10)   !reaction 16: Pra2_1-->TF2_3+Pr_1
   c(17)=1.0d0/v   !reaction 17: TF2_2+Pa3_1-->Pa23_1
   c(18)=1.d0   !reaction 18: Pa23_1-->TF2_2+Pa3_1
   c(19)=c(9)   !reaction 19: TF2_2+Pra3_1-->Pra23_1
   c(20)=c(10)   !reaction 20: Pra23_1-->TF2_2+Pra3_1
   c(21)=1.0d0/v   !reaction 21: TF2_3+Pa2_1-->Pa23_1
   c(22)=1.d0   !reaction 22: Pa23_1-->TF2_3+Pa2_1
   c(23)=c(9)   !reaction 23: TF2_3+Pra2_1-->Pra23_1
   c(24)=c(10)   !reaction 24: Pra23_1-->TF2_3+Pra2_1
   c(25)=10.d0   !reaction 25: Pr_1-->Pr_1+M_1
   c(26)=c(25)   !reaction 26: Pra2_1-->Pra2_1+M_1
   c(27)=c(25)   !reaction 27: Pra3_1-->Pra3_1+M_1
   c(28)=c(25)   !reaction 28: Pra23_1-->Pra23_1+M_1
   c(29)=10.d0   !reaction 29: M_1-->...k4
   c(30)=50.d0   !reaction 30: M_1-->M_1+TF_1
   c(31)=1.d0   !reaction 31: TF_1-->...
   
!module TF levels
   c(32)=1.d0/v*0.1   !reaction 32: TF_1+TF_1--> TF2_1
   c(33)=1.d0   !reaction 33: TF2_1 -->TF_1+TF_1


   if (p0_1.eqv..true.) then
      c(1)=0.01d0   !reaction 1: P_1-->Pr_1
      c(3)=c(1)   !reaction 3: Pa2_1-->Pra2_1
      c(5)=c(1)   !reaction 5: Pa3_1-->Pra3_1
      c(7)=c(1)   !reaction 7: Pa23_1-->Pra23_1
   else if (p1_1.eqv..true.) then
      c(1)=10.d0   !reaction 1: P_1-->Pr_1
      c(3)=c(1)   !reaction 3: Pa2_1-->Pra2_1
      c(5)=c(1)   !reaction 5: Pa3_1-->Pra3_1
      c(7)=c(1)   !reaction 7: Pa23_1-->Pra23_1
   end if


!GENE 2
   c(1+33)=10.d0   !reaction 1+33: P_2-->Pr_2 
   c(2+33)=c(2)   !reaction 2+33: Pr_2-->P_2
   c(3+33)=0.01d0   !reaction 3+33: Pa1_2-->Pra1_2 
   c(4+33)=c(2)   !reaction 4+33: Pra1_2-->Pa1_2
   c(5+33)=10.d0   !reaction 5+33: Pa3_2-->Pra3_2 
   c(6+33)=c(2)   !reaction 6+33: Pra3_2-->Pa3_2
   c(7+33)=10.d0   !reaction 7+33: Pa13_2-->Pra13_2 
   c(8+33)=c(2)   !reaction 8+33: Pra13_2-->Pa13_2
   c(9+33)=1.0d0/v   !reaction 9: TF2_1+P_2-->Pa1_2
   c(10+33)=1.d0   !reaction 10: Pa1_2-->TF2_1+P_2
   c(11+33)=c(9)   !reaction11: TF2_1+Pr_2-->Pra1_2
   c(12+33)=c(10)   !reaction 12: Pra1_2-->TF2_1+Pr_2
   c(13+33)=1.0d0/v   !reaction 13: TF2_3+P_2-->Pa3_2
   c(14+33)=1.d0   !reaction 14: Pa3_2-->TF2_3+P_2
   c(15+33)=c(9)   !reaction 15: TF2_3+Pr_2-->Pra3_2
   c(16+33)=c(10)   !reaction 16: Pra1_2-->TF2_3+Pr_2
   c(17+33)=1.0d0/v   !reaction 17: TF2_1+Pa3_2-->Pa13_2
   c(18+33)=1.d0   !reaction 18: Pa13_2-->TF2_1+Pa3_2
   c(19+33)=c(9)   !reaction 19: TF2_1+Pra3_2-->Pra13_2
   c(20+33)=c(10)   !reaction 20: Pra13_2-->TF2_1+Pra3_2
   c(21+33)=1.0d0/v   !reaction 21: TF2_3+Pa1_2-->Pa13_2
   c(22+33)=1.d0   !reaction 22: Pa13_2-->TF2_3+Pa1_2
   c(23+33)=c(9)   !reaction 23: TF2_3+Pra1_2-->Pra13_2
   c(24+33)=c(10)   !reaction 24: Pra13_2-->TF2_3+Pra1_2
   c(25+33)=10.d0   !reaction 25: Pr_2-->Pr_2+M_2
   c(26+33)=c(25)   !reaction 26: Pra1_2-->Pra1_2+M_2
   c(27+33)=c(25)   !reaction 27: Pra3_2-->Pra3_2+M_2
   c(28+33)=c(25)   !reaction 27: Pra13_2-->Pra13_2+M_2
   c(29+33)=10.d0   !reaction 11: M_2-->...k4
   c(30+33)=50.d0   !reaction 12: M_2-->M_2+TF_2
   c(31+33)=1.d0   !reaction 13: TF_2-->...
   c(32+33)=1.d0/v*0.1   !reaction 14: TF_2+TF_2--> TF2_2
   c(33+33)=1.d0   !reaction 15: TF2_2 -->TF_2+TF_2

   if (p0_2.eqv..true.) then
      c(1+33)=0.01d0   !reaction 1: P_2-->Pr_2
      c(3+33)=c(1+33)   !reaction 3: Pa1_2-->Pra1_2
      c(5+33)=c(1+33)   !reaction 5: Pa3_2-->Pra3_2
      c(7+33)=c(1+33)   !reaction 7: Pa13_2-->Pra13_2
   else if (p1_2.eqv..true.) then
      c(1+33)=10.d0   !reaction 1: P_2-->Pr_2
      c(3+33)=c(1+33)   !reaction 3: Pa1_2-->Pra1_2
      c(5+33)=c(1+33)   !reaction 5: Pa3_2-->Pra3_2
      c(7+33)=c(1+33)   !reaction 7: Pa13_2-->Pra13_2
   end if

!GENE 3
   c(1+66)=0.01d0   !reaction 1: P_3-->Pr_3 
   c(2+66)=c(2)   !reaction 2: Pr_3-->P_3
   c(3+66)=0.01d0   !reaction 3: Pa1_3-->Pra1_3 
   c(4+66)=c(2)   !reaction 4: Pra1_3-->Pa1_3
   c(5+66)=0.01d0   !reaction 5: Pa2_3-->Pra2_3  
   c(6+66)=c(2)   !reaction 6: Pra2_3-->Pa2_3
   c(7+66)=10.d0   !reaction 7: Pa12_3-->Pra12_3 
   c(8+66)=c(2)   !reaction 8: Pra12_3-->Pa12_3
   c(9+66)=1.0d0/v   !reaction 9: TF2_1+P_3-->Pa1_3
   c(10+66)=1.d0   !reaction 10: Pa1_3-->TF2_1+P_3
   c(11+66)=c(9)   !reaction11: TF2_1+Pr_3-->Pra1_3
   c(12+66)=c(10)   !reaction 12: Pra1_3-->TF2_1+Pr_3
   c(13+66)=1.0d0/v   !reaction 13: TF2_2+P_3-->Pa2_3
   c(14+66)=1.d0   !reaction 14: Pa2_3-->TF2_2+P_3
   c(15+66)=c(9)   !reaction 15: TF2_2+Pr_3-->Pra2_3
   c(16+66)=c(10)   !reaction 16: Pra1_3-->TF2_2+Pr_3
   c(17+66)=1.0d0/v   !reaction 17: TF2_1+Pa2_3-->Pa12_3
   c(18+66)=1.d0   !reaction 18: Pa12_3-->TF2_1+Pa2_3
   c(19+66)=c(9)   !reaction 19: TF2_1+Pra2_3-->Pra12_3
   c(20+66)=c(10)   !reaction 20: Pra12_3-->TF2_1+Pra2_3
   c(21+66)=1.0d0/v   !reaction 21: TF2_2+Pa1_3-->Pa12_3
   c(22+66)=1.d0   !reaction 22: Pa12_3-->TF2_2+Pa1_3
   c(23+66)=c(9)   !reaction 23: TF2_2+Pra1_3-->Pra12_3
   c(24+66)=c(10)   !reaction 24: Pra12_3-->TF2_2+Pra1_3
   c(25+66)=10.d0   !reaction 25: Pr_3-->Pr_3+M_3
   c(26+66)=c(25)   !reaction 26: Pra1_3-->Pra1_3+M_3
   c(27+66)=c(25)   !reaction 27: Pra2_3-->Pra2_3+M_3
   c(28+66)=c(25)   !reaction 27: Pra12_3-->Pra12_3+M_3
   c(29+66)=10.d0   !reaction 11: M_3-->...k4
   c(30+66)=50.d0   !reaction 12: M_3-->M_3+TF_3
   c(31+66)=1.d0   !reaction 13: TF_3-->...
   c(32+66)=1.d0/v*0.1   !reaction 14: TF_3+TF_3--> TF2_3
   c(33+66)=1.d0   !reaction 15: TF2_3 -->TF_3+TF_3

   if (p0_3.eqv..true.) then
      c(1+66)=0.01d0   !reaction 1: P_3-->Pr_3
      c(3+66)=c(1+66)   !reaction 3: Pa1_3-->Pra1_3
      c(5+66)=c(1+66)   !reaction 5: Pa2_3-->Pra2_3
      c(7+66)=c(1+66)   !reaction 7: Pa12_3-->Pra12_3
   else if (p1_3.eqv..true.) then
      c(1+66)=10.d0   !reaction 1: P_3-->Pr_3
      c(3+66)=c(1+66)   !reaction 3: Pa1_3-->Pra1_3
      c(5+66)=c(1+66)   !reaction 5: Pa2_3-->Pra2_3
      c(7+66)=c(1+66)   !reaction 7: Pa12_3-->Pra12_3
   end if



   y(1)=1   !P_1
   y(2)=0   !Pa2_1
   y(3)=0   !Pa3_1
   y(4)=0   !Pa23_1
   y(5)=0   !Pr_1
   y(6)=0   !Pra2_1
   y(7)=0   !Pra3_1
   y(8)=0   !Pra23_1
   y(9)=0   !M_1
   y(10)=50   !TF_1
   y(11)=100   !TF2_1

   y(12)=1   !P_2
   y(13)=0   !Pa1_2
   y(14)=0   !Pa3_2
   y(15)=0   !Pa13_2
   y(16)=0   !Pr_1
   y(17)=0   !Pra1_2
   y(18)=0   !Pra3_2
   y(19)=0   !Pra13_2
   y(20)=0   !M_2
   y(21)=50   !TF_2
   y(22)=100   !TF2_2

   y(23)=1   !P_3
   y(24)=0   !Pa1_3
   y(25)=0   !Pa2_3
   y(26)=0   !Pa12_3
   y(27)=0   !Pr_3
   y(28)=0   !Pra1_3
   y(29)=0   !Pra2_3
   y(30)=0   !Pra12_3
   y(31)=0   !M_3
   y(32)=50   !TF_3
   y(33)=100   !TF2_3

   do while(t<4000)

         if (t>=aux_write) then 
            write(1,*)t,',',y(1),',',y(2),',',y(3),',',y(4),',',y(5),',',y(6),',',y(7),',',y(8),',',y(9),',',y(10),',',y(11)
            write(2,*)t,',',y(12),',',y(13),',',y(14),',',y(15),',',y(16),',',y(17),',',y(18),',',y(19),',',y(20),',',y(21),','&
            ,y(22)
            write(3,*)t,',',y(23),',',y(24),',',y(25),',',y(26),',',y(27),',',y(28),',',y(29),',',y(30),',',y(31),',',y(32),','&
            ,y(33)

            aux_write=t+incremento_t

         end if

      h(1)=real(y(1),8)   !reaction 1: P_1-->Pr_1
      h(2)=real(y(5),8)   !reaction 2: Pr_1-->P_1
      h(3)=real(y(2),8)   !reaction 3: Pa2_1-->Pra2_1
      h(4)=real(y(6),8)   !reaction 4: Pra2_1-->Pa2_1
      h(5)=real(y(3),8)   !reaction 5: Pa3_1-->Pra3_1
      h(6)=real(y(7),8)   !reaction 6: Pra3_1-->Pa3_1
      h(7)=real(y(4),8)   !reaction 7: Pa23_1-->Pra23_1
      h(8)=real(y(8),8)   !reaction 8: Pra23_1-->Pa23_1

      h(9)=real(y(22),8)*real(y(1),8)   !reaction 9: TF2_2+P_1-->Pa2_1
      h(10)=real(y(2),8)   !reaction 10: Pa2_1-->TF2_2+P_1
      h(11)=real(y(22),8)*real(y(5),8)   !reaction11: TF2_2+Pr_1-->Pra2_1
      h(12)=real(y(6),8)   !reaction 12: Pra2_1-->TF2_2+Pr_1

      h(13)=real(y(33),8)*real(y(1),8)   !reaction 13: TF2_3+P_1-->Pa3_1
      h(14)=real(y(3),8)   !reaction 14: Pa3_1-->TF2_3+P_1
      h(15)=real(y(33),8)*real(y(5),8)   !reaction 15: TF2_3+Pr_1-->Pra3_1
      h(16)=real(y(7),8)   !reaction 16: Pra2_1-->TF2_3+Pr_1

      h(17)=real(y(22),8)*real(y(3),8)   !reaction 17: TF2_2+Pa3_1-->Pa23_1
      h(18)=real(y(4),8)   !reaction 18: Pa23_1-->TF2_2+Pa3_1
      h(19)=real(y(22),8)*real(y(7),8)   !reaction 19: TF2_2+Pra3_1-->Pra23_1
      h(20)=real(y(8),8)   !reaction 20: Pra23_1-->TF2_2+Pra3_1

      h(21)=real(y(33),8)*real(y(2),8)   !reaction 21: TF2_3+Pa2_1-->Pa23_1
      h(22)=real(y(4),8)   !reaction 22: Pa23_1-->TF2_3+Pa2_1
      h(23)=real(y(33),8)*real(y(6),8)   !reaction 23: TF2_3+Pra2_1-->Pra23_1
      h(24)=real(y(8),8)   !reaction 24: Pra23_1-->TF2_3+Pra2_1

      h(25)=real(y(5),8)   !reaction 25: Pr_1-->Pr_1+M_1
      h(26)=real(y(6),8)   !reaction 26: Pra2_1-->Pra2_1+M_1
      h(27)=real(y(7),8)   !reaction 27: Pra3_1-->Pra3_1+M_1
      h(28)=real(y(8),8)   !reaction 28: Pra23_1-->Pra23_1+M_1

      h(29)=real(y(9),8)   !reaction 29: M_1-->...
      h(30)=real(y(9),8)   !reaction 30: M_1-->M_1+TF_1
      h(31)=real(y(10),8)   !reaction 31: TF_1-->...

      h(32)=real(y(10),8)*(real(y(10),8)-1)/2   !reaction 32: TF_1+TF_1--> TF2_1
      h(33)=real(y(11),8)   !reaction 33: TF2_1 -->TF_1+TF_1

      h(1+33)=real(y(1+11),8)   !reaction 1: P_1-->Pr_1
      h(2+33)=real(y(5+11),8)   !reaction 2: Pr_1-->P_1
      h(3+33)=real(y(2+11),8)   !reaction 3: Pa2_1-->Pra2_1
      h(4+33)=real(y(6+11),8)   !reaction 4: Pra2_1-->Pa2_1
      h(5+33)=real(y(3+11),8)   !reaction 5: Pa3_1-->Pra3_1
      h(6+33)=real(y(7+11),8)   !reaction 6: Pra3_1-->Pa3_1
      h(7+33)=real(y(4+11),8)   !reaction 7: Pa23_1-->Pra23_1
      h(8+33)=real(y(8+11),8)   !reaction 8: Pra23_1-->Pa23_1

      h(9+33)=real(y(11),8)*real(y(1+11),8)   !reaction 9: TF2_2+P_1-->Pa2_1
      h(10+33)=real(y(2+11),8)   !reaction 10: Pa2_1-->TF2_2+P_1
      h(11+33)=real(y(11),8)*real(y(5+11),8)   !reaction11: TF2_2+Pr_1-->Pra2_1
      h(12+33)=real(y(6+11),8)   !reaction 12: Pra2_1-->TF2_2+Pr_1

      h(13+33)=real(y(33),8)*real(y(1+11),8)   !reaction 13: TF2_3+P_1-->Pa3_1
      h(14+33)=real(y(3+11),8)   !reaction 14: Pa3_1-->TF2_3+P_1
      h(15+33)=real(y(33),8)*real(y(5+11),8)   !reaction 15: TF2_3+Pr_1-->Pra3_1
      h(16+33)=real(y(7+11),8)   !reaction 16: Pra2_1-->TF2_3+Pr_1

      h(17+33)=real(y(11),8)*real(y(3+11),8)   !reaction 17: TF2_2+Pa3_1-->Pa23_1
      h(18+33)=real(y(4+11),8)   !reaction 18: Pa23_1-->TF2_2+Pa3_1
      h(19+33)=real(y(11),8)*real(y(7+11),8)   !reaction 19: TF2_2+Pra3_1-->Pra23_1
      h(20+33)=real(y(8+11),8)   !reaction 20: Pra23_1-->TF2_2+Pra3_1

      h(21+33)=real(y(33),8)*real(y(2+11),8)   !reaction 21: TF2_3+Pa2_1-->Pa23_1
      h(22+33)=real(y(4+11),8)   !reaction 22: Pa23_1-->TF2_3+Pa2_1
      h(23+33)=real(y(33),8)*real(y(6+11),8)   !reaction 23: TF2_3+Pra2_1-->Pra23_1
      h(24+33)=real(y(8+11),8)   !reaction 24: Pra23_1-->TF2_3+Pra2_1

      h(25+33)=real(y(5+11),8)   !reaction 25: Pr_1-->Pr_1+M_1
      h(26+33)=real(y(6+11),8)   !reaction 26: Pra2_1-->Pra2_1+M_1
      h(27+33)=real(y(7+11),8)   !reaction 27: Pra3_1-->Pra3_1+M_1
      h(28+33)=real(y(8+11),8)   !reaction 28: Pra23_1-->Pra23_1+M_1

      h(29+33)=real(y(9+11),8)   !reaction 29: M_1-->...
      h(30+33)=real(y(9+11),8)   !reaction 30: M_1-->M_1+TF_1
      h(31+33)=real(y(10+11),8)   !reaction 31: TF_1-->...

      h(32+33)=real(y(10+11),8)*(real(y(10+11),8)-1)/2   !reaction 32: TF_1+TF_1--> TF2_1
      h(33+33)=real(y(11+11),8)   !reaction 33: TF2_1 -->TF_1+TF_1

      h(1+66)=real(y(1+22),8)   !reaction 1: P_1-->Pr_1
      h(2+66)=real(y(5+22),8)   !reaction 2: Pr_1-->P_1
      h(3+66)=real(y(2+22),8)   !reaction 3: Pa2_1-->Pra2_1
      h(4+66)=real(y(6+22),8)   !reaction 4: Pra2_1-->Pa2_1
      h(5+66)=real(y(3+22),8)   !reaction 5: Pa3_1-->Pra3_1
      h(6+66)=real(y(7+22),8)   !reaction 6: Pra3_1-->Pa3_1
      h(7+66)=real(y(4+22),8)   !reaction 7: Pa23_1-->Pra23_1
      h(8+66)=real(y(8+22),8)   !reaction 8: Pra23_1-->Pa23_1

      h(9+66)=real(y(11),8)*real(y(1+22),8)   !reaction 9: TF2_2+P_1-->Pa2_1
      h(10+66)=real(y(2+22),8)   !reaction 10: Pa2_1-->TF2_2+P_1
      h(11+66)=real(y(11),8)*real(y(5+22),8)   !reaction11: TF2_2+Pr_1-->Pra2_1
      h(12+66)=real(y(6+22),8)   !reaction 12: Pra2_1-->TF2_2+Pr_1

      h(13+66)=real(y(22),8)*real(y(1+22),8)   !reaction 13: TF2_3+P_1-->Pa3_1
      h(14+66)=real(y(3+22),8)   !reaction 14: Pa3_1-->TF2_3+P_1
      h(15+66)=real(y(22),8)*real(y(5+22),8)   !reaction 15: TF2_3+Pr_1-->Pra3_1
      h(16+66)=real(y(7+22),8)   !reaction 16: Pra2_1-->TF2_3+Pr_1

      h(17+66)=real(y(11),8)*real(y(3+22),8)   !reaction 17: TF2_2+Pa3_1-->Pa23_1
      h(18+66)=real(y(4+22),8)   !reaction 18: Pa23_1-->TF2_2+Pa3_1
      h(19+66)=real(y(11),8)*real(y(7+22),8)   !reaction 19: TF2_2+Pra3_1-->Pra23_1
      h(20+66)=real(y(8+22),8)   !reaction 20: Pra23_1-->TF2_2+Pra3_1

      h(21+66)=real(y(22),8)*real(y(2+22),8)   !reaction 21: TF2_3+Pa2_1-->Pa23_1
      h(22+66)=real(y(4+22),8)   !reaction 22: Pa23_1-->TF2_3+Pa2_1
      h(23+66)=real(y(22),8)*real(y(6+22),8)   !reaction 23: TF2_3+Pra2_1-->Pra23_1
      h(24+66)=real(y(8+22),8)   !reaction 24: Pra23_1-->TF2_3+Pra2_1

      h(25+66)=real(y(5+22),8)   !reaction 25: Pr_1-->Pr_1+M_1
      h(26+66)=real(y(6+22),8)   !reaction 26: Pra2_1-->Pra2_1+M_1
      h(27+66)=real(y(7+22),8)   !reaction 27: Pra3_1-->Pra3_1+M_1
      h(28+66)=real(y(8+22),8)   !reaction 28: Pra23_1-->Pra23_1+M_1

      h(29+66)=real(y(9+22),8)   !reaction 29: M_1-->...
      h(30+66)=real(y(9+22),8)   !reaction 30: M_1-->M_1+TF_1
      h(31+66)=real(y(10+22),8)   !reaction 31: TF_1-->...

      h(32+66)=real(y(10+22),8)*(real(y(10+22),8)-1)/2   !reaction 32: TF_1+TF_1--> TF2_1
      h(33+66)=real(y(11+22),8)   !reaction 33: TF2_1 -->TF_1+TF_1

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

      if ((y(5)+y(6)+y(7)+y(8)==0).and.(y(16)+y(17)+y(18)+y(19)==0).and.(y(27)+y(28)+y(29)+y(30)==0)) then
         time000=time000+tau
      else if ((y(5)+y(6)+y(7)+y(8)==0).and.(y(16)+y(17)+y(18)+y(19)==0).and.(y(27)+y(28)+y(29)+y(30)==1)) then
         time001=time001+tau
      else if ((y(5)+y(6)+y(7)+y(8)==0).and.(y(16)+y(17)+y(18)+y(19)==1).and.(y(27)+y(28)+y(29)+y(30)==0)) then
         time010=time010+tau
      else if ((y(5)+y(6)+y(7)+y(8)==0).and.(y(16)+y(17)+y(18)+y(19)==1).and.(y(27)+y(28)+y(29)+y(30)==1)) then
         time011=time011+tau
      else if ((y(5)+y(6)+y(7)+y(8)==1).and.(y(16)+y(17)+y(18)+y(19)==0).and.(y(27)+y(28)+y(29)+y(30)==0)) then
         time100=time100+tau
      else if ((y(5)+y(6)+y(7)+y(8)==1).and.(y(16)+y(17)+y(18)+y(19)==0).and.(y(27)+y(28)+y(29)+y(30)==1)) then
         time101=time101+tau
      else if ((y(5)+y(6)+y(7)+y(8)==1).and.(y(16)+y(17)+y(18)+y(19)==1).and.(y(27)+y(28)+y(29)+y(30)==0)) then
         time110=time110+tau
      else if ((y(5)+y(6)+y(7)+y(8)==1).and.(y(16)+y(17)+y(18)+y(19)==1).and.(y(27)+y(28)+y(29)+y(30)==1)) then
         time111=time111+tau
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

      if (mu==1) then !reaction 1: P_1-->Pr_1
         y(1)=y(1)-1
         y(5)=y(5)+1
      else if (mu==2) then !reaction 2: Pr_1-->P_1
         y(1)=y(1)+1
         y(5)=y(5)-1
      else if (mu==3) then !reaction 3: Pa2_1-->Pra2_1
         y(2)=y(2)-1
         y(6)=y(6)+1
      else if (mu==4) then !reaction 4: Pra2_1-->Pa2_1
         y(2)=y(2)+1
         y(6)=y(6)-1
      else if (mu==5) then !reaction 5: Pa3_1-->Pra3_1
         y(3)=y(3)-1
         y(7)=y(7)+1
      else if (mu==6) then !reaction 6: Pra3_1-->Pa3_1
         y(3)=y(3)+1
         y(7)=y(7)-1
      else if (mu==7) then !reaction 7: Pa23_1-->Pra23_1
         y(4)=y(4)-1
         y(8)=y(8)+1
      else if (mu==8) then !reaction 8: Pra23_1-->Pa23_1
         y(4)=y(4)+1
         y(8)=y(8)-1

      else if (mu==9) then !reaction 9: TF2_2+P_1-->Pa2_1
         y(1)=y(1)-1
         y(22)=y(22)-1
         y(2)=y(2)+1
      else if (mu==10) then !reaction 10: Pa2_1-->TF2_2+P_1
         y(1)=y(1)+1
         y(22)=y(22)+1
         y(2)=y(2)-1
      else if (mu==11) then !reaction11: TF2_2+Pr_1-->Pra2_1
         y(5)=y(5)-1
         y(22)=y(22)-1
         y(6)=y(6)+1
      else if (mu==12) then !reaction 12: Pra2_1-->TF2_2+Pr_1
         y(5)=y(5)+1
         y(22)=y(22)+1
         y(6)=y(6)-1

      else if (mu==13) then !reaction 13: TF2_3+P_1-->Pa3_1
         y(1)=y(1)-1
         y(33)=y(33)-1
         y(3)=y(3)+1
      else if (mu==14) then !reaction 14: Pa3_1-->TF2_3+P_1
         y(1)=y(1)+1
         y(33)=y(33)+1
         y(3)=y(3)-1
      else if (mu==15) then !reaction 15: TF2_3+Pr_1-->Pra3_1
         y(5)=y(5)-1
         y(33)=y(33)-1
         y(7)=y(7)+1
      else if (mu==16) then !reaction 16: Pra2_1-->TF2_3+Pr_1
         y(5)=y(5)+1
         y(33)=y(33)+1
         y(7)=y(7)-1

      else if (mu==17) then   !reaction 17: TF2_2+Pa3_1-->Pa23_1
         y(3)=y(3)-1
         y(22)=y(22)-1
         y(4)=y(4)+1
      else if (mu==18) then   !reaction 18: Pa23_1-->TF2_2+Pa3_1
         y(3)=y(3)+1
         y(22)=y(22)+1
         y(4)=y(4)-1
      else if (mu==19) then   !reaction 19: TF2_2+Pra3_1-->Pra23_1
         y(7)=y(7)-1
         y(22)=y(22)-1
         y(8)=y(8)+1
      else if (mu==20) then   !reaction 20: Pra23_1-->TF2_2+Pra3_1
         y(7)=y(7)+1
         y(22)=y(22)+1
         y(8)=y(8)-1

      else if (mu==21) then   !reaction 21: TF2_3+Pa2_1-->Pa23_1
         y(2)=y(2)-1
         y(33)=y(33)-1
         y(4)=y(4)+1
      else if (mu==22) then   !reaction 22: Pa23_1-->TF2_3+Pa2_1
         y(2)=y(2)+1
         y(33)=y(33)+1
         y(4)=y(4)-1
      else if (mu==23) then   !reaction 23: TF2_3+Pra2_1-->Pra23_1
         y(6)=y(6)-1
         y(33)=y(33)-1
         y(8)=y(8)+1
      else if (mu==24) then   !reaction 24: Pra23_1-->TF2_3+Pra2_1
         y(6)=y(6)+1
         y(33)=y(33)+1
         y(8)=y(8)-1

      else if (mu==25) then   !reaction 25: Pr_1-->Pr_1+M_1
         y(9)=y(9)+1
      else if (mu==26) then   !reaction 26: Pra2_1-->Pra2_1+M_1
         y(9)=y(9)+1
      else if (mu==27) then   !reaction 27: Pra3_1-->Pra3_1+M_1
         y(9)=y(9)+1
      else if (mu==28) then   !reaction 28: Pra23_1-->Pra23_1+M_1
         y(9)=y(9)+1

      else if (mu==29) then   !reaction 29: M_1-->...
         y(9)=y(9)-1
      else if (mu==30) then   !reaction 30: M_1-->M_1+TF_1
         y(10)=y(10)+1
      else if (mu==31) then   !reaction 31: TF_1-->...
         y(10)=y(10)-1

      else if (mu==32) then   !reaction 32: TF_1+TF_1--> TF2_1
         y(10)=y(10)-2
         y(11)=y(11)+1
      else if (mu==33) then   !reaction 33: TF2_1 -->TF_1+TF_1
         y(10)=y(10)+2
         y(11)=y(11)-1

      else if (mu==1+33)  then !reaction 1: P_1-->Pr_1
         y(1+11)=y(1+11)-1
         y(5+11)=y(5+11)+1
      else if (mu==2+33)  then !reaction 2: Pr_1-->P_1
         y(1+11)=y(1+11)+1
         y(5+11)=y(5+11)-1
      else if (mu==3+33)  then !reaction 3: Pa2_1-->Pra2_1
         y(2+11)=y(2+11)-1
         y(6+11)=y(6+11)+1
      else if (mu==4+33)  then !reaction 4: Pra2_1-->Pa2_1
         y(2+11)=y(2+11)+1
         y(6+11)=y(6+11)-1
      else if (mu==5+33)  then !reaction 5: Pa3_1-->Pra3_1
         y(3+11)=y(3+11)-1
         y(7+11)=y(7+11)+1
      else if (mu==6+33)  then !reaction 6: Pra3_1-->Pa3_1
         y(3+11)=y(3+11)+1
         y(7+11)=y(7+11)-1
      else if (mu==7+33)  then !reaction 7: Pa23_1-->Pra23_1
         y(4+11)=y(4+11)-1
         y(8+11)=y(8+11)+1
      else if (mu==8+33)  then !reaction 8: Pra23_1-->Pa23_1
         y(4+11)=y(4+11)+1
         y(8+11)=y(8+11)-1

      else if (mu==9+33)  then !reaction 9: TF2_2+P_1-->Pa2_1
         y(1+11)=y(1+11)-1
         y(11)=y(11)-1
         y(2+11)=y(2+11)+1
      else if (mu==10+33)  then !reaction 10: Pa2_1-->TF2_2+P_1
         y(1+11)=y(1+11)+1
         y(11)=y(11)+1
         y(2+11)=y(2+11)-1
      else if (mu==11+33)  then !reaction11: TF2_2+Pr_1-->Pra2_1
         y(5+11)=y(5+11)-1
         y(11)=y(11)-1
         y(6+11)=y(6+11)+1
      else if (mu==12+33)  then !reaction 12: Pra2_1-->TF2_2+Pr_1
         y(5+11)=y(5+11)+1
         y(11)=y(11)+1
         y(6+11)=y(6+11)-1

      else if (mu==13+33)  then !reaction 13: TF2_3+P_1-->Pa3_1
         y(1+11)=y(1+11)-1
         y(33)=y(33)-1
         y(3+11)=y(3+11)+1
      else if (mu==14+33)  then !reaction 14: Pa3_1-->TF2_3+P_1
         y(1+11)=y(1+11)+1
         y(33)=y(33)+1
         y(3+11)=y(3+11)-1
      else if (mu==15+33)  then !reaction 15: TF2_3+Pr_1-->Pra3_1
         y(5+11)=y(5+11)-1
         y(33)=y(33)-1
         y(7+11)=y(7+11)+1
      else if (mu==16+33)  then !reaction 16: Pra2_1-->TF2_3+Pr_1
         y(5+11)=y(5+11)+1
         y(33)=y(33)+1
         y(7+11)=y(7+11)-1

      else if (mu==17+33)  then   !reaction 17: TF2_2+Pa3_1-->Pa23_1
         y(3+11)=y(3+11)-1
         y(11)=y(11)-1
         y(4+11)=y(4+11)+1
      else if (mu==18+33)  then   !reaction 18: Pa23_1-->TF2_2+Pa3_1
         y(3+11)=y(3+11)+1
         y(11)=y(11)+1
         y(4+11)=y(4+11)-1
      else if (mu==19+33)  then   !reaction 19: TF2_2+Pra3_1-->Pra23_1
         y(7+11)=y(7+11)-1
         y(11)=y(11)-1
         y(8+11)=y(8+11)+1
      else if (mu==20+33)  then   !reaction 20: Pra23_1-->TF2_2+Pra3_1
         y(7+11)=y(7+11)+1
         y(11)=y(11)+1
         y(8+11)=y(8+11)-1

      else if (mu==21+33)  then   !reaction 21: TF2_3+Pa2_1-->Pa23_1
         y(2+11)=y(2+11)-1
         y(33)=y(33)-1
         y(4+11)=y(4+11)+1
      else if (mu==22+33)  then   !reaction 22: Pa23_1-->TF2_3+Pa2_1
         y(2+11)=y(2+11)+1
         y(33)=y(33)+1
         y(4+11)=y(4+11)-1
      else if (mu==23+33)  then   !reaction 23: TF2_3+Pra2_1-->Pra23_1
         y(6+11)=y(6+11)-1
         y(33)=y(33)-1
         y(8+11)=y(8+11)+1
      else if (mu==24+33)  then   !reaction 24: Pra23_1-->TF2_3+Pra2_1
         y(6+11)=y(6+11)+1
         y(33)=y(33)+1
         y(8+11)=y(8+11)-1

      else if (mu==25+33)  then   !reaction 25: Pr_1-->Pr_1+M_1
         y(9+11)=y(9+11)+1
      else if (mu==26+33)  then   !reaction 26: Pra2_1-->Pra2_1+M_1
         y(9+11)=y(9+11)+1
      else if (mu==27+33)  then   !reaction 27: Pra3_1-->Pra3_1+M_1
         y(9+11)=y(9+11)+1
      else if (mu==28+33)  then   !reaction 28: Pra23_1-->Pra23_1+M_1
         y(9+11)=y(9+11)+1

      else if (mu==29+33)  then   !reaction 29: M_1-->...
         y(9+11)=y(9+11)-1
      else if (mu==30+33)  then   !reaction 30: M_1-->M_1+TF_1
         y(10+11)=y(10+11)+1
      else if (mu==31+33)  then   !reaction 31: TF_1-->...
         y(10+11)=y(10+11)-1

      else if (mu==32+33)  then   !reaction 32: TF_1+TF_1--> TF2_1
         y(10+11)=y(10+11)-2
         y(11+11)=y(11+11)+1
      else if (mu==33+33)  then   !reaction 33: TF2_1 -->TF_1+TF_1
         y(10+11)=y(10+11)+2
         y(11+11)=y(11+11)-1

      else if (mu==1+66)  then !reaction 1: P_1-->Pr_1
         y(1+22)=y(1+22)-1
         y(5+22)=y(5+22)+1
      else if (mu==2+66)  then !reaction 2: Pr_1-->P_1
         y(1+22)=y(1+22)+1
         y(5+22)=y(5+22)-1
      else if (mu==3+66)  then !reaction 3: Pa2_1-->Pra2_1
         y(2+22)=y(2+22)-1
         y(6+22)=y(6+22)+1
      else if (mu==4+66)  then !reaction 4: Pra2_1-->Pa2_1
         y(2+22)=y(2+22)+1
         y(6+22)=y(6+22)-1
      else if (mu==5+66)  then !reaction 5: Pa3_1-->Pra3_1
         y(3+22)=y(3+22)-1
         y(7+22)=y(7+22)+1
      else if (mu==6+66)  then !reaction 6: Pra3_1-->Pa3_1
         y(3+22)=y(3+22)+1
         y(7+22)=y(7+22)-1
      else if (mu==7+66)  then !reaction 7: Pa23_1-->Pra23_1
         y(4+22)=y(4+22)-1
         y(8+22)=y(8+22)+1
      else if (mu==8+66)  then !reaction 8: Pra23_1-->Pa23_1
         y(4+22)=y(4+22)+1
         y(8+22)=y(8+22)-1

      else if (mu==9+66)  then !reaction 9: TF2_2+P_1-->Pa2_1
         y(1+22)=y(1+22)-1
         y(11)=y(11)-1
         y(2+22)=y(2+22)+1
      else if (mu==10+66)  then !reaction 10: Pa2_1-->TF2_2+P_1
         y(1+22)=y(1+22)+1
         y(11)=y(11)+1
         y(2+22)=y(2+22)-1
      else if (mu==11+66)  then !reaction11: TF2_2+Pr_1-->Pra2_1
         y(5+22)=y(5+22)-1
         y(11)=y(11)-1
         y(6+22)=y(6+22)+1
      else if (mu==12+66)  then !reaction 12: Pra2_1-->TF2_2+Pr_1
         y(5+22)=y(5+22)+1
         y(11)=y(11)+1
         y(6+22)=y(6+22)-1

      else if (mu==13+66)  then !reaction 13: TF2_3+P_1-->Pa3_1
         y(1+22)=y(1+22)-1
         y(22)=y(22)-1
         y(3+22)=y(3+22)+1
      else if (mu==14+66)  then !reaction 14: Pa3_1-->TF2_3+P_1
         y(1+22)=y(1+22)+1
         y(22)=y(22)+1
         y(3+22)=y(3+22)-1
      else if (mu==15+66)  then !reaction 15: TF2_3+Pr_1-->Pra3_1
         y(5+22)=y(5+22)-1
         y(22)=y(22)-1
         y(7+22)=y(7+22)+1
      else if (mu==16+66)  then !reaction 16: Pra2_1-->TF2_3+Pr_1
         y(5+22)=y(5+22)+1
         y(22)=y(22)+1
         y(7+22)=y(7+22)-1

      else if (mu==17+66)  then   !reaction 17: TF2_2+Pa3_1-->Pa23_1
         y(3+22)=y(3+22)-1
         y(11)=y(11)-1
         y(4+22)=y(4+22)+1
      else if (mu==18+66)  then   !reaction 18: Pa23_1-->TF2_2+Pa3_1
         y(3+22)=y(3+22)+1
         y(11)=y(11)+1
         y(4+22)=y(4+22)-1
      else if (mu==19+66)  then   !reaction 19: TF2_2+Pra3_1-->Pra23_1
         y(7+22)=y(7+22)-1
         y(11)=y(11)-1
         y(8+22)=y(8+22)+1
      else if (mu==20+66)  then   !reaction 20: Pra23_1-->TF2_2+Pra3_1
         y(7+22)=y(7+22)+1
         y(11)=y(11)+1
         y(8+22)=y(8+22)-1

      else if (mu==21+66)  then   !reaction 21: TF2_3+Pa2_1-->Pa23_1
         y(2+22)=y(2+22)-1
         y(22)=y(22)-1
         y(4+22)=y(4+22)+1
      else if (mu==22+66)  then   !reaction 22: Pa23_1-->TF2_3+Pa2_1
         y(2+22)=y(2+22)+1
         y(22)=y(22)+1
         y(4+22)=y(4+22)-1
      else if (mu==23+66)  then   !reaction 23: TF2_3+Pra2_1-->Pra23_1
         y(6+22)=y(6+22)-1
         y(22)=y(22)-1
         y(8+22)=y(8+22)+1
      else if (mu==24+66)  then   !reaction 24: Pra23_1-->TF2_3+Pra2_1
         y(6+22)=y(6+22)+1
         y(22)=y(22)+1
         y(8+22)=y(8+22)-1

      else if (mu==25+66)  then   !reaction 25: Pr_1-->Pr_1+M_1
         y(9+22)=y(9+22)+1
      else if (mu==26+66)  then   !reaction 26: Pra2_1-->Pra2_1+M_1
         y(9+22)=y(9+22)+1
      else if (mu==27+66)  then   !reaction 27: Pra3_1-->Pra3_1+M_1
         y(9+22)=y(9+22)+1
      else if (mu==28+66)  then   !reaction 28: Pra23_1-->Pra23_1+M_1
         y(9+22)=y(9+22)+1

      else if (mu==29+66)  then   !reaction 29: M_1-->...
         y(9+22)=y(9+22)-1
      else if (mu==30+66)  then   !reaction 30: M_1-->M_1+TF_1
         y(10+22)=y(10+22)+1
      else if (mu==31+66)  then   !reaction 31: TF_1-->...
         y(10+22)=y(10+22)-1

      else if (mu==32+66)  then   !reaction 32: TF_1+TF_1--> TF2_1
         y(10+22)=y(10+22)-2
         y(11+22)=y(11+22)+1
      else if (mu==33+66)  then   !reaction 33: TF2_1 -->TF_1+TF_1
         y(10+22)=y(10+22)+2
         y(11+22)=y(11+22)-1
      end if

   end do

   write(4,*)t,',',time000,',',time001,',',time010,',',time011,',',time100,',',time101,',',time110,',',time111

end program

