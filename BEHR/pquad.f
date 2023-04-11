      subroutine pquad (f, a, b, absacc, relacc, minrul, maxrul,ruleno, 
     &quads, nevals, ercode)
      double precision a, absacc, b, relacc
      double precision f, quads(9)
      external f
      integer ercode, minrul, maxrul, nevals, ruleno
      double precision abserr, length, midpt, p(0:255), w(0:510), zero
      double precision sum, value, values(0:127)
      integer ptnum, ptstrt, pstops(8), ptstop, wtnum
      common /ppw511/ p, w
      data pstops /0, 1, 3, 7, 15, 31, 63, 127/
      data zero /0d0/
      if(.not.( (absacc .le. zero .and. relacc .le. zero).or. minrul 
     &.lt. 2 .or. minrul .gt. maxrul .or. maxrul .gt. 9 ))goto 23000
      do 23002 ruleno = 1,9
      quads(ruleno) = zero
23002 continue
      ruleno = min(9, max(2, minrul))
      nevals = 0
      ercode = 2
      return
23000 continue
      length = b - a
      midpt = (a + b) / 2
      ercode = 1
      ptstop = pstops(minrul - 1)
      wtnum = ptstop
      values(0) = f(midpt)
      sum = values(0) * w(wtnum)
      do 23004 ptnum = 1, ptstop 
      wtnum = wtnum + 1
      value = f(a + length * p(ptnum)) + f(b - length * p(ptnum))
      sum = sum + value * w(wtnum)
      values(ptnum) = value
23004 continue
      quads(minrul - 1) = sum * length
      do 23006 ruleno = minrul, maxrul 
      sum = 0
      do 23008 ptnum = 0, ptstop 
      wtnum = wtnum + 1
      sum = sum + values(ptnum) * w(wtnum)
23008 continue
      ptstrt = ptstop + 1
      ptstop = ptstrt + ptstop
      do 23010 ptnum = ptstrt, ptstop 
      wtnum = wtnum + 1
      value = f(a + length * p(ptnum)) + f(b - length * p(ptnum))
      sum = sum + value * w(wtnum)
      if(.not.(ruleno .lt. 9))goto 23012
      values(ptnum) = value
23012 continue
23010 continue
      quads(ruleno) = sum * length
      abserr = abs( quads(ruleno) - quads(ruleno - 1) )
      if(.not.(abserr .le. max(absacc,relacc * abs(quads(ruleno)))))
     &goto 23014
      ercode = 0
      goto 23007
23014 continue
23006 continue
23007 continue
      ruleno = min (ruleno, maxrul)
      nevals = 2 * ptstop + 1
      return
      end
