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

   n=33 !number of molecular species
   m=99 !number of reactions
   v=1.d0 !volume
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
   c(32)=1.d0/v*0.1   !reaction 32: TF_1+TF_1--> TF2_1
   c(33)=1.d0   !reaction 33: TF2_1 -->TF_1+TF_1

   if (p0_1.eqv..true.) then
      c(1)=0.01d0
      c(3)=c(1)
      c(5)=c(1)
      c(7)=c(1)
   else if (p1_1.eqv..true.) then
      c(1)=10.d0
      c(3)=c(1)
      c(5)=c(1)
      c(7)=c(1)
   end if

!GENE 2
   c(34)=10.d0   !reaction 34: P_2-->Pr_2 
   c(35)=c(2)   !reaction 35: Pr_2-->P_2
   c(36)=0.01d0   !reaction 36: Pa1_2-->Pra1_2 
   c(37)=c(2)   !reaction 37: Pra1_2-->Pa1_2
   c(38)=10.d0   !reaction 38: Pa3_2-->Pra3_2 
   c(39)=c(2)   !reaction 39: Pra3_2-->Pa3_2
   c(40)=10.d0   !reaction 40: Pa13_2-->Pra13_2 
   c(41)=c(2)   !reaction 41: Pra13_2-->Pa13_2
   c(42)=1.0d0/v   !reaction 42: TF2_1+P_2-->Pa1_2
   c(43)=1.d0   !reaction 43: Pa1_2-->TF2_1+P_2
   c(44)=c(9)   !reaction 44: TF2_1+Pr_2-->Pra1_2
   c(45)=c(10)   !reaction 45: Pra1_2-->TF2_1+Pr_2
   c(46)=1.0d0/v   !reaction 46: TF2_3+P_2-->Pa3_2
   c(47)=1.d0   !reaction 47: Pa3_2-->TF2_3+P_2
   c(48)=c(9)   !reaction 48: TF2_3+Pr_2-->Pra3_2
   c(49)=c(10)   !reaction 49: Pra1_2-->TF2_3+Pr_2
   c(50)=1.0d0/v   !reaction 50: TF2_1+Pa3_2-->Pa13_2
   c(51)=1.d0   !reaction 51: Pa13_2-->TF2_1+Pa3_2
   c(52)=c(9)   !reaction 52: TF2_1+Pra3_2-->Pra13_2
   c(53)=c(10)   !reaction 53: Pra13_2-->TF2_1+Pra3_2
   c(54)=1.0d0/v   !reaction 54: TF2_3+Pa1_2-->Pa13_2
   c(55)=1.d0   !reaction 55: Pa13_2-->TF2_3+Pa1_2
   c(56)=c(9)   !reaction 56: TF2_3+Pra1_2-->Pra13_2
   c(57)=c(10)   !reaction 57: Pra13_2-->TF2_3+Pra1_2
   c(58)=10.d0   !reaction 58: Pr_2-->Pr_2+M_2
   c(59)=c(25)   !reaction 59: Pra1_2-->Pra1_2+M_2
   c(60)=c(25)   !reaction 60: Pra3_2-->Pra3_2+M_2
   c(61)=c(25)   !reaction 61: Pra13_2-->Pra13_2+M_2
   c(62)=10.d0   !reaction 62: M_2-->...k4
   c(63)=50.d0   !reaction 63: M_2-->M_2+TF_2
   c(64)=1.d0   !reaction 64: TF_2-->...
   c(65)=1.d0/v*0.1   !reaction 65: TF_2+TF_2--> TF2_2
   c(66)=1.d0   !reaction 66: TF2_2 -->TF_2+TF_2

   if (p0_2.eqv..true.) then
      c(34)=0.01d0
      c(36)=c(34)
      c(38)=c(34)
      c(40)=c(34)
   else if (p1_2.eqv..true.) then
      c(34)=10.d0
      c(36)=c(34)
      c(38)=c(34)
      c(40)=c(34)
   end if

!GENE 3
   c(67)=0.01d0   !reaction 67: P_3-->Pr_3 
   c(68)=c(2)   !reaction 68: Pr_3-->P_3
   c(69)=0.01d0   !reaction 69: Pa1_3-->Pra1_3 
   c(70)=c(2)   !reaction 70: Pra1_3-->Pa1_3
   c(71)=0.01d0   !reaction 71: Pa2_3-->Pra2_3  
   c(72)=c(2)   !reaction 72: Pra2_3-->Pa2_3
   c(73)=10.d0   !reaction 73: Pa12_3-->Pra12_3 
   c(74)=c(2)   !reaction 74: Pra12_3-->Pa12_3
   c(75)=1.0d0/v   !reaction 75: TF2_1+P_3-->Pa1_3
   c(76)=1.d0   !reaction 76: Pa1_3-->TF2_1+P_3
   c(77)=c(9)   !reaction 77: TF2_1+Pr_3-->Pra1_3
   c(78)=c(10)   !reaction 78: Pra1_3-->TF2_1+Pr_3
   c(79)=1.0d0/v   !reaction 79: TF2_2+P_3-->Pa2_3
   c(80)=1.d0   !reaction 80: Pa2_3-->TF2_2+P_3
   c(81)=c(9)   !reaction 81: TF2_2+Pr_3-->Pra2_3
   c(82)=c(10)   !reaction 82: Pra1_3-->TF2_2+Pr_3
   c(83)=1.0d0/v   !reaction 83: TF2_1+Pa2_3-->Pa12_3
   c(84)=1.d0   !reaction 84: Pa12_3-->TF2_1+Pa2_3
   c(85)=c(9)   !reaction 85: TF2_1+Pra2_3-->Pra12_3
   c(86)=c(10)   !reaction 86: Pra12_3-->TF2_1+Pra2_3
   c(87)=1.0d0/v   !reaction 87: TF2_2+Pa1_3-->Pa12_3
   c(88)=1.d0   !reaction 88: Pa12_3-->TF2_2+Pa1_3
   c(89)=c(9)   !reaction 89: TF2_2+Pra1_3-->Pra12_3
   c(90)=c(10)   !reaction 90: Pra12_3-->TF2_2+Pra1_3
   c(91)=10.d0   !reaction 91: Pr_3-->Pr_3+M_3
   c(92)=c(91)   !reaction 92: Pra1_3-->Pra1_3+M_3
   c(93)=c(91)   !reaction 93: Pra2_3-->Pra2_3+M_3
   c(94)=c(91)   !reaction 94: Pra12_3-->Pra12_3+M_3
   c(95)=10.d0   !reaction 95: M_3-->...k4
   c(96)=50.d0   !reaction 96: M_3-->M_3+TF_3
   c(97)=1.d0   !reaction 97: TF_3-->...
   c(98)=1.d0/v*0.1   !reaction 98: TF_3+TF_3--> TF2_3
   c(99)=1.d0   !reaction 99: TF2_3 -->TF_3+TF_3

   if (p0_3.eqv..true.) then
      c(67)=0.01d0   !reaction 67: P_3-->Pr_3
      c(69)=c(67)   !reaction 69: Pa1_3-->Pra1_3
      c(71)=c(67)   !reaction 71: Pa2_3-->Pra2_3
      c(73)=c(67)   !reaction 73: Pa12_3-->Pra12_3
   else if (p1_3.eqv..true.) then
      c(67)=10.d0   !reaction 67: P_3-->Pr_3
      c(69)=c(67)   !reaction 69: Pa1_3-->Pra1_3
      c(71)=c(67)   !reaction 71: Pa2_3-->Pra2_3
      c(73)=c(67)   !reaction 73: Pa12_3-->Pra12_3
   end if

   !molecular species
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

      !molecular reactant combinations
      h(1)=real(y(1),8)   !reaction 1
      h(2)=real(y(5),8)   !reaction 2
      h(3)=real(y(2),8)   !reaction 3
      h(4)=real(y(6),8)   !reaction 4
      h(5)=real(y(3),8)   !reaction 5
      h(6)=real(y(7),8)   !reaction 6
      h(7)=real(y(4),8)   !reaction 7
      h(8)=real(y(8),8)   !reaction 8
      h(9)=real(y(22),8)*real(y(1),8)   !reaction 9
      h(10)=real(y(2),8)   !reaction 10
      h(11)=real(y(22),8)*real(y(5),8)   !reaction 11
      h(12)=real(y(6),8)   !reaction 12
      h(13)=real(y(33),8)*real(y(1),8)   !reaction 13
      h(14)=real(y(3),8)   !reaction 14
      h(15)=real(y(33),8)*real(y(5),8)   !reaction 15
      h(16)=real(y(7),8)   !reaction 16
      h(17)=real(y(22),8)*real(y(3),8)   !reaction 17
      h(18)=real(y(4),8)   !reaction 18
      h(19)=real(y(22),8)*real(y(7),8)   !reaction 19
      h(20)=real(y(8),8)   !reaction 20
      h(21)=real(y(33),8)*real(y(2),8)   !reaction 21
      h(22)=real(y(4),8)   !reaction 22
      h(23)=real(y(33),8)*real(y(6),8)   !reaction 23
      h(24)=real(y(8),8)   !reaction 24
      h(25)=real(y(5),8)   !reaction 25
      h(26)=real(y(6),8)   !reaction 26
      h(27)=real(y(7),8)   !reaction 27
      h(28)=real(y(8),8)   !reaction 28
      h(29)=real(y(9),8)   !reaction 29
      h(30)=real(y(9),8)   !reaction 30
      h(31)=real(y(10),8)   !reaction 31
      h(32)=real(y(10),8)*(real(y(10),8)-1)/2   !reaction 32
      h(33)=real(y(11),8)   !reaction 33
      h(34)=real(y(12),8)   !reaction 34
      h(35)=real(y(16),8)   !reaction 35
      h(36)=real(y(13),8)   !reaction 36
      h(37)=real(y(17),8)   !reaction 37
      h(38)=real(y(14),8)   !reaction 38
      h(39)=real(y(18),8)   !reaction 39
      h(40)=real(y(15),8)   !reaction 40
      h(41)=real(y(19),8)   !reaction 41
      h(42)=real(y(11),8)*real(y(12),8)   !reaction 42
      h(43)=real(y(13),8)   !reaction 43
      h(44)=real(y(11),8)*real(y(16),8)   !reaction 44
      h(45)=real(y(17),8)   !reaction 45
      h(46)=real(y(33),8)*real(y(12),8)   !reaction 46
      h(47)=real(y(14),8)   !reaction 47
      h(48)=real(y(33),8)*real(y(16),8)   !reaction 48
      h(49)=real(y(18),8)   !reaction 49
      h(50)=real(y(11),8)*real(y(14),8)   !reaction 50
      h(51)=real(y(15),8)   !reaction 51
      h(52)=real(y(11),8)*real(y(18),8)   !reaction 52
      h(53)=real(y(19),8)   !reaction 53
      h(54)=real(y(33),8)*real(y(13),8)   !reaction 54
      h(55)=real(y(15),8)   !reaction 55
      h(56)=real(y(33),8)*real(y(17),8)   !reaction 56
      h(57)=real(y(19),8)   !reaction 57
      h(58)=real(y(16),8)   !reaction 58
      h(59)=real(y(17),8)   !reaction 59
      h(60)=real(y(18),8)   !reaction 60
      h(61)=real(y(19),8)   !reaction 61
      h(62)=real(y(20),8)   !reaction 62
      h(63)=real(y(20),8)   !reaction 63
      h(64)=real(y(21),8)   !reaction 64
      h(65)=real(y(21),8)*(real(y(21),8)-1)/2   !reaction 65
      h(66)=real(y(22),8)   !reaction 66
      h(67)=real(y(23),8)   !reaction 67
      h(68)=real(y(27),8)   !reaction 68
      h(69)=real(y(24),8)   !reaction 69
      h(70)=real(y(28),8)   !reaction 70
      h(71)=real(y(25),8)   !reaction 71
      h(72)=real(y(29),8)   !reaction 72
      h(73)=real(y(26),8)   !reaction 73
      h(74)=real(y(30),8)   !reaction 74
      h(75)=real(y(11),8)*real(y(23),8)   !reaction 75
      h(76)=real(y(24),8)   !reaction 76
      h(77)=real(y(11),8)*real(y(27),8)   !reaction 77
      h(78)=real(y(28),8)   !reaction 78
      h(79)=real(y(22),8)*real(y(23),8)   !reaction 79
      h(80)=real(y(25),8)   !reaction 80
      h(81)=real(y(22),8)*real(y(27),8)   !reaction 81
      h(82)=real(y(29),8)   !reaction 82
      h(83)=real(y(11),8)*real(y(25),8)   !reaction 83
      h(84)=real(y(26),8)   !reaction 84
      h(85)=real(y(11),8)*real(y(29),8)   !reaction 85
      h(86)=real(y(30),8)   !reaction 86
      h(87)=real(y(22),8)*real(y(24),8)   !reaction 87
      h(88)=real(y(26),8)   !reaction 88
      h(89)=real(y(22),8)*real(y(28),8)   !reaction 89
      h(90)=real(y(30),8)   !reaction 90
      h(91)=real(y(27),8)   !reaction 91
      h(92)=real(y(28),8)   !reaction 92
      h(93)=real(y(29),8)   !reaction 93
      h(94)=real(y(30),8)   !reaction 94
      h(95)=real(y(31),8)   !reaction 95
      h(96)=real(y(31),8)   !reaction 96
      h(97)=real(y(32),8)   !reaction 97
      h(98)=real(y(32),8)*(real(y(32),8)-1)/2   !reaction 98
      h(99)=real(y(33),8)   !reaction 99

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
      mu=0 !mu represents the chemical reaction that occurs in each iteration
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
         y(5)=y(5)+1
      else if (mu==2) then !reaction 2
         y(1)=y(1)+1
         y(5)=y(5)-1
      else if (mu==3) then !reaction 3
         y(2)=y(2)-1
         y(6)=y(6)+1
      else if (mu==4) then !reaction 4
         y(2)=y(2)+1
         y(6)=y(6)-1
      else if (mu==5) then !reaction 5
         y(3)=y(3)-1
         y(7)=y(7)+1
      else if (mu==6) then !reaction 6
         y(3)=y(3)+1
         y(7)=y(7)-1
      else if (mu==7) then !reaction 7
         y(4)=y(4)-1
         y(8)=y(8)+1
      else if (mu==8) then !reaction 8
         y(4)=y(4)+1
         y(8)=y(8)-1
      else if (mu==9) then !reaction 9
         y(1)=y(1)-1
         y(22)=y(22)-1
         y(2)=y(2)+1
      else if (mu==10) then !reaction 10
         y(1)=y(1)+1
         y(22)=y(22)+1
         y(2)=y(2)-1
      else if (mu==11) then !reaction 11
         y(5)=y(5)-1
         y(22)=y(22)-1
         y(6)=y(6)+1
      else if (mu==12) then !reaction 12
         y(5)=y(5)+1
         y(22)=y(22)+1
         y(6)=y(6)-1
      else if (mu==13) then !reaction 13
         y(1)=y(1)-1
         y(33)=y(33)-1
         y(3)=y(3)+1
      else if (mu==14) then !reaction 14
         y(1)=y(1)+1
         y(33)=y(33)+1
         y(3)=y(3)-1
      else if (mu==15) then !reaction 15
         y(5)=y(5)-1
         y(33)=y(33)-1
         y(7)=y(7)+1
      else if (mu==16) then !reaction 16
         y(5)=y(5)+1
         y(33)=y(33)+1
         y(7)=y(7)-1
      else if (mu==17) then !reaction 17
         y(3)=y(3)-1
         y(22)=y(22)-1
         y(4)=y(4)+1
      else if (mu==18) then !reaction 18
         y(3)=y(3)+1
         y(22)=y(22)+1
         y(4)=y(4)-1
      else if (mu==19) then !reaction 19
         y(7)=y(7)-1
         y(22)=y(22)-1
         y(8)=y(8)+1
      else if (mu==20) then !reaction 20
         y(7)=y(7)+1
         y(22)=y(22)+1
         y(8)=y(8)-1
      else if (mu==21) then !reaction 21
         y(2)=y(2)-1
         y(33)=y(33)-1
         y(4)=y(4)+1
      else if (mu==22) then !reaction 22
         y(2)=y(2)+1
         y(33)=y(33)+1
         y(4)=y(4)-1
      else if (mu==23) then !reaction 23
         y(6)=y(6)-1
         y(33)=y(33)-1
         y(8)=y(8)+1
      else if (mu==24) then !reaction 24
         y(6)=y(6)+1
         y(33)=y(33)+1
         y(8)=y(8)-1
      else if (mu==25) then !reaction 25
         y(9)=y(9)+1
      else if (mu==26) then !reaction 26
         y(9)=y(9)+1
      else if (mu==27) then !reaction 27
         y(9)=y(9)+1
      else if (mu==28) then !reaction 28
         y(9)=y(9)+1
      else if (mu==29) then !reaction 29
         y(9)=y(9)-1
      else if (mu==30) then !reaction 30
         y(10)=y(10)+1
      else if (mu==31) then !reaction 31
         y(10)=y(10)-1
      else if (mu==32) then !reaction 32
         y(10)=y(10)-2
         y(11)=y(11)+1
      else if (mu==33) then !reaction 33
         y(10)=y(10)+2
         y(11)=y(11)-1
      else if (mu==34) then   !reaction 34
         y(12)=y(12)-1
         y(16)=y(16)+1
      else if (mu==35) then   !reaction 35
         y(12)=y(12)+1
         y(16)=y(16)-1
      else if (mu==36) then   !reaction 36
         y(13)=y(13)-1
         y(17)=y(17)+1
      else if (mu==37) then   !reaction 37
         y(13)=y(13)+1
         y(17)=y(17)-1
      else if (mu==38) then   !reaction 38
         y(14)=y(14)-1
         y(18)=y(18)+1
      else if (mu==39) then   !reaction 39
         y(14)=y(14)+1
         y(18)=y(18)-1
      else if (mu==40) then   !reaction 40
         y(15)=y(15)-1
         y(19)=y(19)+1
      else if (mu==41) then   !reaction 41
         y(15)=y(15)+1
         y(19)=y(19)-1
      else if (mu==42) then   !reaction 42
         y(12)=y(12)-1
         y(22)=y(22)-1
         y(13)=y(13)+1
      else if (mu==43) then   !reaction 43
         y(12)=y(12)+1
         y(22)=y(22)+1
         y(13)=y(13)-1
      else if (mu==44) then   !reaction 44
         y(16)=y(16)-1
         y(22)=y(22)-1
         y(17)=y(17)+1
      else if (mu==45) then   !reaction 45
         y(16)=y(16)+1
         y(22)=y(22)+1
         y(17)=y(17)-1
      else if (mu==46) then   !reaction 46
         y(12)=y(12)-1
         y(33)=y(33)-1
         y(14)=y(14)+1
      else if (mu==47) then   !reaction 47
         y(12)=y(12)+1
         y(33)=y(33)+1
         y(14)=y(14)-1
      else if (mu==48) then   !reaction 48
         y(16)=y(16)-1
         y(33)=y(33)-1
         y(18)=y(18)+1
      else if (mu==49) then   !reaction 49
         y(16)=y(16)+1
         y(33)=y(33)+1
         y(18)=y(18)-1
      else if (mu==50) then   !reaction 50
         y(14)=y(14)-1
         y(22)=y(22)-1
         y(15)=y(15)+1
      else if (mu==51) then   !reaction 51
         y(14)=y(14)+1
         y(22)=y(22)+1
         y(15)=y(15)-1
      else if (mu==52) then   !reaction 52
         y(18)=y(18)-1
         y(22)=y(22)-1
         y(19)=y(19)+1
      else if (mu==53) then   !reaction 53
         y(18)=y(18)+1
         y(22)=y(22)+1
         y(19)=y(19)-1
      else if (mu==54) then   !reaction 54
         y(13)=y(13)-1
         y(33)=y(33)-1
         y(15)=y(15)+1
      else if (mu==55) then   !reaction 55
         y(13)=y(13)+1
         y(33)=y(33)+1
         y(15)=y(15)-1
      else if (mu==56) then   !reaction 56
         y(17)=y(17)-1
         y(33)=y(33)-1
         y(19)=y(19)+1
      else if (mu==57) then   !reaction 57
         y(17)=y(17)+1
         y(33)=y(33)+1
         y(19)=y(19)-1
      else if (mu==58) then   !reaction 58
         y(16)=y(16)-1
         y(20)=y(20)+1
      else if (mu==59) then   !reaction 59
         y(17)=y(17)-1
         y(20)=y(20)+1
      else if (mu==60) then   !reaction 60
         y(18)=y(18)-1
         y(20)=y(20)+1
      else if (mu==61) then   !reaction 61
         y(19)=y(19)-1
         y(20)=y(20)+1
      else if (mu==62) then   !reaction 62
         y(20)=y(20)-1
      else if (mu==63) then   !reaction 63
         y(21)=y(21)+1
      else if (mu==64) then   !reaction 64
         y(21)=y(21)-1
      else if (mu==65) then   !reaction 65
         y(21)=y(21)-2
         y(22)=y(22)+1
      else if (mu==66) then   !reaction 66
         y(21)=y(21)+2
         y(22)=y(22)-1
      else if (mu==67) then !reaction 67
         y(23)=y(23)-1
         y(27)=y(27)+1
      else if (mu==68) then !reaction 68
         y(23)=y(23)+1
         y(27)=y(27)-1
      else if (mu==69) then !reaction 69
         y(24)=y(24)-1
         y(28)=y(28)+1
      else if (mu==70) then !reaction 70
         y(24)=y(24)+1
         y(28)=y(28)-1
      else if (mu==71) then !reaction 71
         y(25)=y(25)-1
         y(29)=y(29)+1
      else if (mu==72) then !reaction 72
         y(25)=y(25)+1
         y(29)=y(29)-1
      else if (mu==73) then !reaction 73
         y(26)=y(26)-1
         y(30)=y(30)+1
      else if (mu==74) then !reaction 74
         y(26)=y(26)+1
         y(30)=y(30)-1
      else if (mu==75) then !reaction 75
         y(23)=y(23)-1
         y(33)=y(33)-1
         y(24)=y(24)+1
      else if (mu==76) then !reaction 76
         y(23)=y(23)+1
         y(33)=y(33)+1
         y(24)=y(24)-1
      else if (mu==77) then !reaction 77
         y(27)=y(27)-1
         y(33)=y(33)-1
         y(28)=y(28)+1
      else if (mu==78) then !reaction 78
         y(27)=y(27)+1
         y(33)=y(33)+1
         y(28)=y(28)-1
      else if (mu==79) then !reaction 79
         y(23)=y(23)-1
         y(45)=y(45)-1
         y(25)=y(25)+1
      else if (mu==80) then !reaction 80
         y(23)=y(23)+1
         y(45)=y(45)+1
         y(25)=y(25)-1
      else if (mu==81) then !reaction 81
         y(27)=y(27)-1
         y(45)=y(45)-1
         y(29)=y(29)+1
      else if (mu==82) then !reaction 82
         y(27)=y(27)+1
         y(45)=y(45)+1
         y(29)=y(29)-1
      else if (mu==83) then !reaction 83
         y(25)=y(25)-1
         y(33)=y(33)-1
         y(26)=y(26)+1
      else if (mu==84) then !reaction 84
         y(25)=y(25)+1
         y(33)=y(33)+1
         y(26)=y(26)-1
      else if (mu==85) then !reaction 85
         y(29)=y(29)-1
         y(33)=y(33)-1
         y(30)=y(30)+1
      else if (mu==86) then !reaction 86
         y(29)=y(29)+1
         y(33)=y(33)+1
         y(30)=y(30)-1
      else if (mu==87) then !reaction 87
         y(24)=y(24)-1
         y(45)=y(45)-1
         y(26)=y(26)+1
      else if (mu==88) then !reaction 88
         y(24)=y(24)+1
         y(45)=y(45)+1
         y(26)=y(26)-1
      else if (mu==89) then !reaction 89
         y(28)=y(28)-1
         y(45)=y(45)-1
         y(30)=y(30)+1
      else if (mu==90) then !reaction 90
         y(28)=y(28)+1
         y(45)=y(45)+1
         y(30)=y(30)-1
      else if (mu==91) then !reaction 91
         y(31)=y(31)+1
      else if (mu==92) then !reaction 92
         y(31)=y(31)+1
      else if (mu==93) then !reaction 93
         y(31)=y(31)+1
      else if (mu==94) then !reaction 94
         y(31)=y(31)+1
      else if (mu==95) then !reaction 95
         y(31)=y(31)-1
      else if (mu==96) then !reaction 96
         y(32)=y(32)+1
      else if (mu==97) then !reaction 97
         y(32)=y(32)-1
      else if (mu==98) then !reaction 98
         y(32)=y(32)-2
         y(33)=y(33)+1
      else if (mu==99) then !reaction 99
         y(32)=y(32)+2
         y(33)=y(33)-1
      end if

   end do

   write(4,*)t,',',time000,',',time001,',',time010,',',time011,',',time100,',',time101,',',time110,',',time111

end program
