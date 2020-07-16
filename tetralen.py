from numpy import *
import time

# Created by Mohammed Abdellateef

def tetraLen(x1, x2, x3, ph1, ph2, ph3, precision):
   """
   Tetrahedreon lengths:
   	finding missing lengths in tetrahedreon from known tetrahedreon
      base lengths (3 lengths) and head-base angles (3 angles).

   In: 3 lengths, 3 angles
   Out: 3 lengths
   """
   T = time.time()
   # Interpolation function
   def intrPol(rng, rngd, wantd):
        topd = rngd[0]; downd = rngd[1]; top = rng[0]; down = rng[1]
        difr = top - down
        difrd = topd - downd
        difrw = wantd - downd
	# print(topd, downd, difr, difrd, difrw)
        rangeSize = ((difrw*difr)/difrd)*0.3
        newDown = down + rangeSize
        newTop = down + 3*rangeSize
        if i == 1 or i ==2:  # for the mid you've to flip upside-down
          if newTop < top:
             top = newTop
        else:
          if newTop > top:
             top = newTop
        rngs = [top, newTop, newDown, down]
	# print rngs
        return rngs  # left lengths

     # we need this for rng 1, 2 first
   # Middle lengths calculated from outside lengths
   def midFromOut(i, j, rng, ph):
      if i==0:
         lmid = midUp(rng[j], x1, ph1)
         rmid = midUp(rng[j], x3, ph3)
      elif i==1:
         lmid = midUp(rng[j], x1, ph1)
         rmid = midDown(rng[j], x3, ph3)
      elif i==2:
         lmid = midDown(rng[j], x1, ph1)
         rmid = midUp(rng[j], x3, ph3)
      else:
         lmid = midDown(rng[j], x1, ph1)
         rmid = midDown(rng[j], x3, ph3)
      xbar = thirdLen(lmid, rmid, ph)
      return xbar, lmid, rmid
   # Comparing ranges to mid-mid known length
   def compRanges(rngs, rngd, x, ph):
      xbar, lmid, rmid = midFromOut(i, 1, rngs, ph)
      xbar2, lmid2, rmid2 = midFromOut(i, 2, rngs, ph)

      if xbar > x:
         if xbar2 > x:
           top =  rngs[2]
           down = rngs[3]
           topd = xbar2
           downd = rngd[1]
           topMid = [lmid2, rmid2]
           downMid = []
         else:
           top = rngs[1]
           down = rngs[2]
           topd = xbar
           downd = xbar2
           topMid = [lmid, rmid]
           downMid = [lmid2, rmid2]
      else:
        top = rngs[0]
        down = rngs[1]
        topd = rngd[0]
        downd = xbar
        topMid = []
        downMid = [lmid, rmid]
      rng = [top, down]
      rngd = [topd, downd]
      mids = [topMid, downMid]
      # print(rngd, rng, x2)
      return rng, rngd, mids

   # Farthest Points (longest lengths):
   # Equally sides Tetrahedreon
   def isoTri(x1, ph1):
      th1 = (pi-ph1)/2
      lr = x1*sin(th1)/sin(ph1)
      return lr

   # longest length Projection
   def riTri(x1, ph1):
      inside = x1/sin(ph1)     # hypotenuse
      thp1 = pi/2-ph1
      outside = inside*sin(thp1)
      return outside, inside

   # min as applied for a one from bottom
   def twoMid(l, x, ph):
      # find lower and upper lengths (Isosceles Tri)
      # given 2 lengths and angle
      ath = l*sin(ph)/x
      if ath > 1:
         th = nan          # uncontained case
      else:
         th = arcsin(ath)
         if pi-th > th:
            th = pi - th
      # th = arcsin(l*sin(ph)/x)
      thu = th - ph
      thl = th + ph
      up = l*sin(thu)/sin(th)
      lo = l*sin(thl)/sin(th)
      return up, lo 

   def midUp(l, x, ph):
      ath = l*sin(ph)/x
      if ath > 1:
         th = nan          # uncontained case
      else:
         th = arcsin(ath)
         if pi-th > th:
            th = pi - th
      thu = th - ph
      up = l*sin(thu)/sin(th)
      return up

   def midDown(l, x, ph):
      ath = l*sin(ph)/x
      if ath > 1:
         th = nan          # uncontained case
      else:
         th = arcsin(ath)
         if pi-th > th:
            th = pi - th
      thl = th + ph
      lo = l*sin(thl)/sin(th)
      return lo 

   # cos rule to find missing length
   def thirdLen(l1, l2, th):
      l3 = sqrt(l1**2 + l2**2 - 2*l1*l2*cos(th))
      return l3

   # Mid-mid possible lengths for a candidate
   def fourLen(leftu, leftl, rightu, rightl, ph):
      uu = thirdLen(leftu, rightu, ph)
      ul = thirdLen(leftu, rightl, ph)
      lu = thirdLen(leftl, rightu, ph)
      ll = thirdLen(leftl, rightl, ph)
      minimax = [uu, ul, lu, ll]
      return minimax

   # Acceptance Criterea
   def accepCrit(pericision, rngd, rng,  mids, x):
      topd = rngd[0]
      downd = rngd[1]
      top = rng[0]
      down = rng[1]
      # when should I accept
      dud = topd - downd
      solval = 0
      if round(topd, 4) == x:
         solval = [top, mids[0][0], mids[0][1]]                          # left
         return True, solval
      elif round(downd, 4) == x:
         solval = [down, mids[1][0], mids[1][1]]
         return True, solval
      else:
         return False, solval

   def compLen(rng, rngd, wantd):
   # decimal digits equality delmma round(leng/maxi/mini, 4)
      if rngd[0] < wantd or rngd[1] > wantd:
         return False
      else:
         return True

   # Maximum (left/right) length (upper bound)
   outl, inl = riTri(x1, ph1)
   outr, inr = riTri(x3, ph3)
   lm = isoTri(x1, ph1)
   rm = isoTri(x3, ph3)
   # taking smaller length from left/right lengths
   if lm > rm:
      isoOut = rm  # all 4 diagonal lines are equal
   else:
      isoOut = lm

   # upper region
   if outl > outr:             # right is smaller in length
      out = outr   # no change in inr
      inl = midUp(out, x1, ph1)
   else:
      out = outl   # no change in inl
      inr = midUp(out, x3, ph3)

   inll = midDown(out, x1, ph1) # max left lower
   inrl = midDown(out, x3, ph3) # max right lower
   min_ul = thirdLen(inl, inrl, ph2)
   min_lu = thirdLen(inll, inr, ph2)

   if x1 < x3:     # take the shortes side
      lol = x1    # lower left
      ul = sin(pi-2*ph1)*x1/sin(ph1) # min: upper left
      lr = midDown(x1, x3, ph3)  # min: lower right
      mx_ul = thirdLen(ul, lr, ph2)
      mx_lu = midUp(x1, x3, ph3)
   else:
      lol = x3
      lr = sin(pi-2*ph3)*x3/sin(ph3)
      ul = midDown(x3, x1, ph1)
      mx_lu = thirdLen(lr, ul, ph2)
      mx_ul = midUp(x3, x1, ph1) # directly equal to mid-mid length

   if inl > inr:
      ull = inr
      umr = outr
      uml = midDown(inr, x1, ph1)
   else:
      ull = inl
      uml = outl
      umr = midDown(inl, x3, ph3)

   # out, mx_ul, mx_lu
   outside = [out, lol, lol, ull, 0,inl, inll, 0]  # flipped for ul/lu
   # outside = [out, lol, lol, ull, 0, inl, inll, 0]

   mx_uu = thirdLen(inl,    inr,    ph2)
   mx_ll = thirdLen(uml, umr, ph2)
   # mx_lu/ul is close to the unknown point and min_lu/ul is far
   # work upside-down
   minimax = [mx_uu, mx_ul, mx_lu, mx_ll, 0, min_ul, min_lu, 0] # upper-lower should apparently be more than
   # minimax = [mx_uu, min_ul, min_lu, mx_ll, 0, mx_ul, mx_lu, 0] # upper-lower should apparently be more than

   mids = [0, 0]
   # Special Case solution
   sol = []
   th = (pi - ph2)/2
   l = sin(th)*x2/sin(ph2)
   if round(l, 4) == round(isoOut, 4):
      sol.append([isoOut, isoOut, isoOut])

   all_checked = False
   t = 1
   k = 0
   q = 0
   while all_checked is False:
      for i in range(4):
         rng = [outside[i], outside[i+4]]
         rngd = [minimax[i], minimax[i+4]]
         inRange = compLen(rng, rngd, x2)
         if inRange is False:
            continue
         ss, val = accepCrit(precision, rngd, rng, mids, x2)
         while ss is False:
	    # if q == 5:
	    #    break
            # q = q+1
            rngs = intrPol(rng, rngd, x2)  # new search bound
            rng, rngd, mids = compRanges(rngs, rngd, x2, ph2)
            ss, val = accepCrit(precision, rngd, rng, mids, x2)
         if min(val) > 0:
            sol.append(val)
            #print ("Solution No. " + str(t) + ": " + str(val))
            t = t+1
         k = 0
         all_checked = True
   elapsed = time.time() - T
   nSol = t-1 # number of solutions
   return sol, nSol, elapsed

