import numpy as np
from numpy import sin, cos, sqrt, arcsin, arccos, pi, isfinite, array, dot, clip, sort
import time
import random


# Created by Mohammed Abdellateef
# Root-finding driver rewritten for correctness; the geometric primitives
# (midUp, midDown, thirdLen, riTri, isoTri) are unchanged/validated against
# known geometry. The coarse bracket scan is vectorized with NumPy for
# speed; final roots are refined with scalar bisection to `precision`.

def tetraLen(x1, x2, x3, ph1, ph2, ph3, precision, samples=250):
   """
   Tetrahedron missing-edge solver.

   Same geometric algorithm as the original version, with instrumentation
   for profiling the root-finding cost.
   """
   T = time.perf_counter()

   tol = 10 ** (-precision) if precision >= 1 else precision
   dedup_radius = max(1e-6, tol)

   stats = {
       "scan_points": 0,
       "brackets": 0,
       "exact_roots": 0,
       "bisect_roots": 0,
       "bisect_iterations": 0,
       "g_calls": 0,
   }

   def isoTri(x1, ph1):
      th1 = (pi-ph1)/2
      lr = x1*sin(th1)/sin(ph1)
      return lr

   def riTri(x1, ph1):
      inside = x1/sin(ph1)
      thp1 = pi/2-ph1
      outside = inside*sin(thp1)
      return outside, inside

   # l*sin(ph)/x should be exactly 1.0 right at the feasibility
   # boundary (h == hm_max), but hm_max itself was back-computed as
   # x/sin(ph), so the round trip can land a hair past 1.0
   # (e.g. 1.0000000000000002). Without tolerance, that single ULP
   # trips "ath > 1" and silently turns a legitimate boundary sample
   # into NaN -- which can erase the one sample that would have
   # closed a bracket right at the edge of the domain. Clamp instead
   # of rejecting for any overshoot this small; anything larger is a
   # genuinely infeasible triangle and still returns NaN.
   boundary_eps = 1e-9

   def midUp(l, x, ph):
      ath = l*sin(ph)/x
      if ath > 1:
         if ath - 1 > boundary_eps:
            return np.nan
         ath = 1.0
      th = arcsin(ath)
      if pi-th > th:
         th = pi-th
      thu = th-ph
      return l*sin(thu)/sin(th)

   def midDown(l, x, ph):
      ath = l*sin(ph)/x
      if ath > 1:
         if ath - 1 > boundary_eps:
            return np.nan
         ath = 1.0
      th = arcsin(ath)
      if pi-th > th:
         th = pi-th
      thl = th+ph
      return l*sin(thl)/sin(th)

   def midUp_vec(l, x, ph):
      with np.errstate(divide='ignore', invalid='ignore'):
         ath = l*sin(ph)/x
         valid = np.isfinite(ath) & (ath <= 1 + boundary_eps) & (ath >= -1 - boundary_eps)

         ath_c = np.clip(ath, -1.0, 1.0)
         th = arcsin(ath_c)
         th = np.where(pi-th > th, pi-th, th)

         thu = th-ph
         den = sin(th)

         out = np.divide(
             l*sin(thu),
             den,
             out=np.full_like(den, np.nan, dtype=float),
             where=np.abs(den) > 1e-14
         )

         return np.where(valid, out, np.nan)

   def midDown_vec(l, x, ph):
      with np.errstate(divide='ignore', invalid='ignore'):
         ath = l*sin(ph)/x
         valid = np.isfinite(ath) & (ath <= 1 + boundary_eps) & (ath >= -1 - boundary_eps)

         ath_c = np.clip(ath, -1.0, 1.0)
         th = arcsin(ath_c)
         th = np.where(pi-th > th, pi-th, th)

         thl = th+ph
         den = sin(th)

         out = np.divide(
             l*sin(thl),
             den,
             out=np.full_like(den, np.nan, dtype=float),
             where=np.abs(den) > 1e-14
         )

         return np.where(valid, out, np.nan)

   def thirdLen(l1, l2, th):
      return sqrt(l1**2 + l2**2 - 2*l1*l2*cos(th))

   _, inl = riTri(x1, ph1)
   _, inr = riTri(x3, ph3)
   hm_max = min(inl, inr)

   sol = []

   if not (isfinite(hm_max) and hm_max > 0):
      elapsed = time.perf_counter() - T
      stats["time_ms"] = elapsed * 1000
      return sol, 0, elapsed, stats

   branches = [
       (midUp,   midUp,   midUp_vec,   midUp_vec),
       (midUp,   midDown, midUp_vec,   midDown_vec),
       (midDown, midUp,   midDown_vec, midUp_vec),
       (midDown, midDown, midDown_vec, midDown_vec)
   ]

   eps = hm_max * 1e-9
   hs = np.linspace(eps, hm_max, samples)

   # Add extra resolution near the upper geometric boundary.
   upper_width = max(0.01 * hm_max, 10.0 * tol)
   upper_start = max(eps, hm_max - upper_width)

   upper_samples = max(32, samples // 8)
   hs_upper = np.linspace(upper_start, hm_max, upper_samples)

   hs = np.unique(np.concatenate((hs, hs_upper)))
   
   
   
   stats["scan_points"] = samples * len(branches)

   def g_scalar(fL, fR, h):
      stats["g_calls"] += 1

      hl = fL(h, x1, ph1)
      hr = fR(h, x3, ph3)

      if not (isfinite(hl) and isfinite(hr)):
         return np.nan, None, None

      return thirdLen(hl, hr, ph2) - x2, hl, hr

   stats["tangent_roots"] = 0
   stats["refine_windows"] = 0
   stats["refine_points"] = 0

   def bisect_root(fL, fR, a, b, fa):
      # Same scalar bisection used for the coarse brackets, factored
      # out so the adaptive refinement below can reuse it.
      #
      # Near the geometric boundary (h -> hm_max, i.e. ath -> 1) the
      # map from h to the output lengths is nearly singular: arcsin's
      # derivative blows up as its argument approaches 1, so a tiny
      # residual width in h can translate into a large error in the
      # returned lengths. Bisecting only to `0.01*tol` in h is not
      # tight enough there, so converge h itself to (near) machine
      # precision -- 200 halvings is enormously more than needed and
      # simply stops making progress once a/b/m collide in float64.
      h_scale = max(abs(a), abs(b), 1.0)
      width_floor = h_scale * 1e-14

      for _ in range(200):
         stats["bisect_iterations"] += 1
         m = (a + b) / 2.0
         fm, _, _ = g_scalar(fL, fR, m)
         if not isfinite(fm):
            break
         if (fa < 0) == (fm < 0):
            a, fa = m, fm
         else:
            b = m
         if (b - a) < width_floor:
            break
      return (a + b) / 2.0

   def golden_min(fL, fR, a, b, iters=100):
      # Golden-section search for the h in [a, b] that minimizes
      # |g(h)|. Used to pin down genuine tangential (double) roots:
      # points where g grazes zero without ever changing sign, even
      # under fine resampling.
      gr = (sqrt(5.0) - 1.0) / 2.0

      def absg(h):
         val, hl, hr = g_scalar(fL, fR, h)
         if not isfinite(val):
            return np.inf, None, None
         return abs(val), hl, hr

      c = b - gr * (b - a)
      d = a + gr * (b - a)
      fc, _, _ = absg(c)
      fd, _, _ = absg(d)

      for _ in range(iters):
         if (b - a) < 1e-14:
            break
         if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr * (b - a)
            fc, _, _ = absg(c)
         else:
            a, c, fc = c, d, fd
            d = a + gr * (b - a)
            fd, _, _ = absg(d)

      hm = (a + b) / 2.0
      fm, hl_m, hr_m = absg(hm)
      return hm, fm, hl_m, hr_m

   def refine_window(fL, fR, fL_vec, fR_vec, a, b, sub=64):
      # A coarse sample can hide TWO close roots when the function
      # dips below (or rises above) zero and back between two
      # same-signed samples -- the endpoint-only crossing test
      # (v0*v1 < 0) can't see that. Re-sample the suspect window
      # much more finely and re-run the ordinary crossing/exact test
      # on it, recursing one level if the finer grid still shows a
      # turning point that could itself be hiding a pair of roots.
      roots = []

      hs_fine = np.linspace(a, b, sub)
      hl_f = fL_vec(hs_fine, x1, ph1)
      hr_f = fR_vec(hs_fine, x3, ph3)
      valid_f = isfinite(hl_f) & isfinite(hr_f)
      v_f = np.where(valid_f, thirdLen(hl_f, hr_f, ph2) - x2, np.nan)

      stats["refine_points"] += sub

      v0f, v1f = v_f[:-1], v_f[1:]
      both_f = isfinite(v0f) & isfinite(v1f)
      exact_f = both_f & (v0f == 0)
      cross_f = both_f & (v0f * v1f < 0)

      for i in np.nonzero(exact_f | cross_f)[0]:
         if exact_f[i]:
            roots.append(hs_fine[i])
         else:
            roots.append(bisect_root(fL, fR, hs_fine[i], hs_fine[i + 1], v0f[i]))

      if roots:
         return roots

      # Still no crossing at fine resolution: check whether it's a
      # genuine tangency (extremum sitting essentially on zero).
      hm, fm, hl_m, hr_m = golden_min(fL, fR, a, b)
      if isfinite(fm) and fm < tol:
         stats["tangent_roots"] += 1
         roots.append(hm)

      return roots

   found = []

   for fL, fR, fL_vec, fR_vec in branches:

      hl_arr = fL_vec(hs, x1, ph1)
      hr_arr = fR_vec(hs, x3, ph3)

      valid = isfinite(hl_arr) & isfinite(hr_arr)
      v_arr = np.where(
          valid,
          thirdLen(hl_arr, hr_arr, ph2) - x2,
          np.nan
      )

      v0s = v_arr[:-1]
      v1s = v_arr[1:]

      both_finite = isfinite(v0s) & isfinite(v1s)
      exact = both_finite & (v0s == 0)
      crossing = both_finite & (v0s * v1s < 0)

      bracket_idx = np.nonzero(exact | crossing)[0]

      stats["brackets"] += len(bracket_idx)

      for i in bracket_idx:

         if exact[i]:

            root = hs[i]
            stats["exact_roots"] += 1

         else:

            stats["bisect_roots"] += 1
            root = bisect_root(fL, fR, hs[i], hs[i+1], v0s[i])

         _, hl_r, hr_r = g_scalar(fL, fR, root)

         if (
             hl_r is not None
             and hr_r is not None
             and min(root, hl_r, hr_r) > 0
         ):
            found.append([root, hl_r, hr_r])

      # The one place a hidden dip can occur with no sampled point on
      # both sides to reveal it as a "turning point" is the very last
      # coarse interval: h = hm_max is the feasibility edge for
      # whichever branch defines it, and repeated cases show the true
      # root often sits just inside that edge, dipping to zero and
      # partially recovering before hs[-1]. There is no hs[len(hs)]
      # to complete a 3-point turning check there, so always refine
      # that final interval directly (cheap: one more 64-point scan).
      if len(hs) >= 2:
         a0, b0 = v_arr[-2], v_arr[-1]
         if isfinite(a0) and isfinite(b0):
            stats["refine_windows"] += 1
            for root in refine_window(
                fL, fR, fL_vec, fR_vec, hs[-2], hs[-1]
            ):
               _, hl_r, hr_r = g_scalar(fL, fR, root)
               if (
                   hl_r is not None
                   and hr_r is not None
                   and min(root, hl_r, hr_r) > 0
               ):
                  found.append([root, hl_r, hr_r])

      # Same-signed samples can still straddle a pair of close roots
      # (the function dips to/past zero and back within one coarse
      # interval) or a genuine tangency. Any interior point that is a
      # local extremum of the coarse sample is a candidate for this;
      # re-examine that window at much higher resolution.
      for i in range(1, len(v_arr) - 1):

         a0, b0, c0 = v_arr[i - 1], v_arr[i], v_arr[i + 1]

         if not (isfinite(a0) and isfinite(b0) and isfinite(c0)):
            continue

         turning = (b0 - a0) * (c0 - b0) < 0

         if not turning:
            continue

         stats["refine_windows"] += 1

         for root in refine_window(
             fL, fR, fL_vec, fR_vec, hs[i - 1], hs[i + 1]
         ):
            _, hl_r, hr_r = g_scalar(fL, fR, root)
            if (
                hl_r is not None
                and hr_r is not None
                and min(root, hl_r, hr_r) > 0
            ):
               found.append([root, hl_r, hr_r])

   for cand in found:
      if not any(
          all(abs(a-b) < dedup_radius for a,b in zip(cand, existing))
          for existing in sol
      ):
         sol.append(cand)

   elapsed = time.perf_counter() - T
   stats["time_ms"] = elapsed * 1000

   return sol, len(sol), elapsed, stats

# ============================================================
# TETRALEN TEST DATA GENERATOR
# ============================================================

def generate_tetra_test():
    """
    Generate a random tetrahedron and convert it to the
    input format required by tetraLen().

    Returns:
        x1, x2, x3          known base lengths
        ph1, ph2, ph3       head-base angles
        expected            three missing lengths
        points              original 3-D coordinates
    """

    # Random tetrahedron vertices
    A = array([0.0, 0.0, 0.0])

    B = array([
        random.uniform(2.0, 10.0),
        0.0,
        0.0
    ])

    C = array([
        random.uniform(-5.0, 10.0),
        random.uniform(2.0, 10.0),
        0.0
    ])

    D = array([
        random.uniform(-5.0, 10.0),
        random.uniform(-5.0, 10.0),
        random.uniform(2.0, 10.0)
    ])

    # Distance helper
    def dist(P, Q):
        return sqrt(sum((P - Q) ** 2))

    # --------------------------------------------------------
    # Six actual tetrahedron edges
    #
    # Base:
    #   AB = x1
    #   AC = x2
    #   BC = x3
    #
    # Unknown:
    #   AD
    #   BD
    #   CD
    # --------------------------------------------------------

    AB = dist(A, B)
    AC = dist(A, C)
    BC = dist(B, C)

    AD = dist(A, D)
    BD = dist(B, D)
    CD = dist(C, D)

    # --------------------------------------------------------
    # Angles at the head D
    #
    # ph1 = angle A-D-B
    # ph2 = angle A-D-C
    # ph3 = angle B-D-C
    # --------------------------------------------------------

    def angle(P, Q, R):
        """
        Angle P-Q-R.
        """
        v1 = P - Q
        v2 = R - Q

        c = dot(v1, v2) / (dist(P, Q) * dist(R, Q))
        

        # protect against floating-point errors
        c = clip(c, -1.0, 1.0)

        return arccos(c)

    ph1 = angle(A, D, B)
    ph2 = angle(A, D, C)
    ph3 = angle(B, D, C)

    return (
        AB, AC, BC,
        ph1, ph2, ph3,
        array([AD, BD, CD]),
        (A, B, C, D)
    )


def test_tetraLen(n=10, precision=4):
    """
    Generate and test n random tetrahedra against tetraLen().
    """

    print("=" * 70)
    print("TETRALEN RANDOM TEST")
    print("=" * 70)

    passed = 0
    failed = 0

    for test_no in range(1, n + 1):

        (
            x1, x2, x3,
            ph1, ph2, ph3,
            expected,
            points
        ) = generate_tetra_test()

        print("\nTest", test_no)
        print("-" * 70)

        print("Known base:")
        print("x1 =", x1)
        print("x2 =", x2)
        print("x3 =", x3)

        print("\nAngles:")
        print("ph1 =", ph1)
        print("ph2 =", ph2)
        print("ph3 =", ph3)

        print("\nExpected missing lengths:")
        print("AD =", expected[0])
        print("BD =", expected[1])
        print("CD =", expected[2])

        try:
            result, nsol, elapsed, states = tetraLen(
                x1, x2, x3,
                ph1, ph2, ph3,
                precision
            )

            print("\nTetralen result:")
            print(result)
            print("Number of solutions:", nsol)
            print("Time:", elapsed)
            print("States:", states)

            # Check whether any returned solution matches
            # the actual tetrahedron.
            found = False

            for sol in result:

                if len(sol) != 3:
                    continue

                # Because tetraLen may return a different
                # ordering depending on its geometric branch,
                # compare sorted lengths.
                a = sort(array(sol, dtype=float))
                b = sort(expected)

                if all(abs(a - b) < 10 ** (-precision)):
                    found = True
                    break

            if found:
                print("PASS")
                passed += 1
            else:
                print("FAIL")
                print("Expected:", expected)
                print("Returned:", result)
                failed += 1

        except Exception as e:
            print("\nERROR:", type(e).__name__, e)
            failed += 1

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print("Passed:", passed)
    print("Failed:", failed)
    print("Total :", n)

if __name__ == "__main__":
    test_tetraLen(10)
    
