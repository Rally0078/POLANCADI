            PROGRAM  POLRUN
c------                  MAINLINE FOR N(H) ANALYSIS USING 'POLAN'.    Feb 1985.
c 28.1.93  Basic heads to 80 chr; Data,Peak.. 98/100; debug 132.  Del init '0'.
c 28.1.93  Read input file name;
c     Read (1):   field  (and station heading);
c     Read (2):   data,  data,  data,   . . .
c
c     Use 1 blank line to reread a station,field line (1);   2 blanks to exit.
c
c     An initial data line with FH = 9.  gives quick-check without data lists.
c                          then FH =-9.  reverts to normal output.
C-------------------------------------------------------------------------------
      type :: input_data
            real :: freq
            real :: height
      end type input_data

      DIMENSION FV(399), HT(399), QQ(50)
!     Initialize FV array to have ones(any nonzero value)
!     Not initializing causes platform-dependent branching 
!      due to uninitialized values
      data fv /399*1/
      data ndim / 399 /
      character head*25, dat*8, datin*15
      common /test/ test1, test2, test3              !put in POLAN... for debugs
    8 format (/, a,/ ( f8.2, 9f10.2) )
    9 format (/, a,/ (f6.2,f6.1, 6(f8.2, f6.1)) )
      
                                                       !set array size
      integer, parameter :: N = 399
      type(input_data) :: inputs(N)
      integer :: size = 0
      integer :: endoffile = 0
      call date_and_time (dat)

!     Fixed input name to receive input from IDL
!     See in.dat for the compatible format

      datin = 'in.dat'
      OPEN (UNIT=1, FILE=datin, STATUS='OLD')
      OPEN (UNIT=2, FILE='out.dat')
      write (2,40) datin, dat
40    format (' P O L A N    of  FEBRUARY 1993.     ',
     &        'Data file: ',a15, '   Run: ',A8,'.')
cc      QQ(1) = -1.          ! no coefficients back
10       quik = 0.         ! quick-look off
         GO TO 100
20          quik = 1.      ! quick-look on
c                                      read field and mode
100   READ (1,120) HEAD, FH,DIP, AMODE, VALLEY, LIST                    read (1)
120   format (a25, 4f5.0 , i5, f5.0)
              if (fh.eq. 0.)  stop
              if (fh.eq. 9.)  go to 20
              if (fh.eq.-9.)  go to 10
c
      WRITE (2,160)
      WRITE (2,160) HEAD, FH,DIP, AMODE, VALLEY, LIST
      write (*,160) head, fh,dip, amode,valley,list                         show
160   format (a25,'   FH',f5.2, '  Dip',f5.1,
     1       '    Amode',f5.1,'  Valy',f6.2,'  List',i2)
c
c-------------------------------------- read data - loop for new ionogram
200   READ (1,220) HEAD, START
220   format (a25, f5.3)
!     Below line is not needed for single data source
         if (fv(1).eq.0.) go to 100                    !read new field constants
      WRITE (2,240)
      WRITE (2,240) HEAD, START
      write (*,240) HEAD, START                                             show
240   format (a25,'   Start =',F7.3)
!     Read from file into the type input_data
!     Stop reading and jump to 330 if end of file is reached
      do i = 1, N
            read(1, *, end=330) inputs(i)
            size = size + 1
      end do
330   endoffile=1
!     Initialize the FV and HT arrays with the received input
      FV = inputs%freq
      HT = inputs%height
!     Old code for old format, not needed anymore
!         READ (1,320) (FV(I),HT(I), I= NH+1, NH+8)                      -read(2)
!320      format (8(f5.3,f5.2))
!         NH = NH+8
!     Write output to out.dat
      if (quik<1.) write (2,9) 'Input data', (fv(i),ht(i), i= 1,size)         list
c-------------------------------------- analysis
c
            if (quik.ge.1.)  list = 0
!Verify output with real results, dont touch this subroutine and ruin the blood sacrifice
      CALL POLAN (size,FV,HT, QQ, FH,DIP, START, AMODE, VALLEY, LIST)
c
c-------------------------------------- output
               if (quik.ge.1.) go to 200
      write (2,9) 'Real Heights', (fv(i),ht(i), i= 1,size+2)
c
      if (qq(1).gt.2.) write(2,8) 'Coefficients QQ',
     &                            (qq(i), i=1, ifix(qq(1))+1)
      write (2,*)'===================================================='
      END
